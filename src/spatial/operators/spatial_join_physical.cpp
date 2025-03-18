#include "spatial/operators/spatial_join_physical.hpp"
#include "spatial_join_logical.hpp"
#include "spatial/geometry/sgl.hpp"
#include "spatial/geometry/geometry_type.hpp"
#include "spatial/spatial_types.hpp"

#include "duckdb/common/types/row/tuple_data_iterator.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"
#include "duckdb/common/types/row/tuple_data_collection.hpp"
#include "duckdb/storage/buffer_manager.hpp"

namespace duckdb {

//======================================================================================================================
// Flat RTree
//======================================================================================================================

class FlatRTree {
public:
	using Box = Box2D<float>;

	FlatRTree(uint32_t item_count_p, uint32_t node_size_p) : item_count(item_count_p), node_size(node_size_p) {

		indices.resize(item_count);
		boxes.resize(item_count);
		rows.resize(item_count);

		uint32_t count = item_count;
		uint32_t nodes = item_count;

		layer_bounds.push_back(nodes);

		do {
			count = (count + node_size - 1) / node_size;
			nodes += count;
			layer_bounds.push_back(nodes);
		} while (count > 1);

		boxes.resize(nodes);
		indices.resize(nodes);
		rows.resize(nodes);

		// auto index_size = nodes * sizeof(uint32_t);
		// auto nodes_size = nodes * sizeof(Box);
	}

	uint32_t Count() const {
		return item_count;
	}

	// Return insertion index
	uint32_t Push(const Box &box, data_ptr_t row) {
		// Push the index and the box
		indices[current_position] = current_position;
		boxes[current_position] = box;

		// Update the bounds
		tree_box.Union(box);

		// Store the row pointer
		rows[current_position] = row;

		return current_position++;
	}

	void Sort(vector<uint32_t> &curve) {
		Sort(curve, 0, curve.size() - 1);
	}

	void Sort(vector<uint32_t> &curve, size_t l_idx, size_t r_idx) {
		if (l_idx < r_idx) {
			const auto pivot = curve[(l_idx + r_idx) >> 1];
			auto pivot_l = l_idx - 1;
			auto pivot_r = r_idx + 1;

			while (true) {
				do {
					++pivot_l;
				} while (curve[pivot_l] < pivot);
				do {
					--pivot_r;
				} while (curve[pivot_r] > pivot);

				if (pivot_l >= pivot_r) {
					break;
				}

				// Reorder the curve, boxes and indices
				// TODO: Pass callback here and make static
				std::swap(curve[pivot_l], curve[pivot_r]);
				std::swap(boxes[pivot_l], boxes[pivot_r]);
				std::swap(indices[pivot_l], indices[pivot_r]);
			}

			Sort(curve, l_idx, pivot_r);
			Sort(curve, pivot_r + 1, r_idx);
		}
	}

	void Build() {
		D_ASSERT(item_count == current_position);

		if (item_count <= node_size) {
			boxes[current_position++] = tree_box;
			return;
		}

		// Generate hilbert curve values
		// TODO: Parallelize this with tasks when the number of items is large?
		constexpr auto max_hilbert = std::numeric_limits<uint16_t>::max();
		const auto hw = max_hilbert / (tree_box.max.x - tree_box.min.x);
		const auto hh = max_hilbert / (tree_box.max.y - tree_box.min.y);

		vector<uint32_t> curve(item_count);
		for (idx_t i = 0; i < item_count; i++) {
			const auto &node_box = boxes[i];

			const auto hx = static_cast<uint32_t>(hw * ((node_box.min.x + node_box.max.x) / 2 - tree_box.min.x));
			const auto hy = static_cast<uint32_t>(hh * ((node_box.min.y + node_box.max.y) / 2 - tree_box.min.y));

			curve[i] = sgl::util::hilbert_encode(16, hx, hy);
		}

		// Now, sort the indices based on their curve value
		Sort(curve);

		size_t layer_idx = 0;
		size_t entry_idx = 0;

		while (layer_idx < layer_bounds.size() - 1) {
			const auto entry_end = layer_bounds[layer_idx];

			while (entry_idx < entry_end) {
				auto node_idx = entry_idx;
				auto node_box = boxes[entry_idx];

				size_t child_idx = 0;
				while (child_idx < node_size && entry_idx < entry_end) {

					node_box.Union(boxes[entry_idx]);

					child_idx++;
					entry_idx++;
				}

				// Add a new parent node
				indices[current_position] = node_idx;
				boxes[current_position] = node_box;
				current_position++;
			}

			// Go to the next layer
			layer_idx++;
		}
	}

	size_t UpperBound(size_t node_idx) const {
		const auto it = std::upper_bound(layer_bounds.begin(), layer_bounds.end(), node_idx);
		if (it == layer_bounds.end()) {
			return layer_bounds.back();
		}
		return *it;
	}

	struct LookupState {
		queue<size_t> search_queue;
		Box search_box;
		size_t entry_beg;
		size_t entry_pos;
	};

	void InitLookup(LookupState &state, const Box &box) const {
		state.search_box = box;
		state.entry_beg = boxes.size() - 1;
		state.entry_pos = state.entry_beg;
	}

	idx_t Lookup(LookupState &state, Vector &result) const {
		idx_t count = 0;
		const auto ptr = FlatVector::GetData<data_ptr_t>(result);
		Lookup(state, [&](const data_ptr_t &row) {
			ptr[count++] = row;
			return count == STANDARD_VECTOR_SIZE;
		});
		return count;
	}

	template <class CALLBACK>
	bool Lookup(LookupState &state, CALLBACK &&callback) const {

		while (true) {

			const auto entry_end = std::min(state.entry_beg + node_size, UpperBound(state.entry_beg));

			while (state.entry_pos < entry_end) {
				if (!state.search_box.Intersects(boxes[state.entry_pos])) {
					continue;
				}

				auto yield = false;

				if (state.entry_beg >= item_count) {
					// Internal node
					state.search_queue.push(indices[state.entry_pos]);
				} else {
					// Leaf node
					yield = callback(rows[indices[state.entry_pos]]);
				}

				state.entry_pos++;

				if (yield) {
					// Yield!, return true to signal that there might be more rows!
					return true;
				}
			}

			if (state.search_queue.empty()) {
				// There is no more nodes to search, return false!
				return false;
			}

			state.entry_beg = state.search_queue.front();
			state.entry_pos = state.entry_beg;
			state.search_queue.pop();
		}
	}

	// Return true if any rows were found
	bool Lookup(const Box &box, vector<data_ptr_t> &result) const {
		queue<size_t> search_queue;

		auto node_index = boxes.size() - 1;

		while (true) {
			// find the end index of the node
			const size_t entry_end = std::min(node_index + node_size, UpperBound(node_index));

			// search through child nodes
			for (size_t entry_pos = node_index; entry_pos < entry_end; entry_pos++) {
				// check if node bbox intersects with query bbox
				if (!box.Intersects(boxes[entry_pos])) {
					continue;
				}

				const size_t index = indices[entry_pos];
				if (node_index >= item_count) {
					// Internal node
					search_queue.push(index);
				} else {
					// Leaf node
					result.push_back(rows[index]);
				}
			}

			if (search_queue.empty()) {
				break;
			}

			node_index = search_queue.front();
			search_queue.pop();
		}

		return !result.empty();
	}

private:
	vector<uint32_t> indices;
	vector<Box> boxes;
	vector<uint32_t> layer_bounds;

	vector<data_ptr_t> rows;

	Box tree_box;

	uint32_t item_count = 0;
	uint32_t node_size = 0;
	uint32_t current_position = 0;
};

//======================================================================================================================
// Physical Spatial Join Operator
//======================================================================================================================

PhysicalSpatialJoin::PhysicalSpatialJoin(LogicalOperator &op, unique_ptr<PhysicalOperator> left,
                                         unique_ptr<PhysicalOperator> right, vector<SpatialJoinCondition> conditions_p,
                                         JoinType join_type, idx_t estimated_cardinality)
    : PhysicalJoin(op, PhysicalOperatorType::EXTENSION, join_type, estimated_cardinality),
      conditions(std::move(conditions_p)) {

	children.push_back(std::move(left));
	children.push_back(std::move(right));

	// Only inner joins for now!
	D_ASSERT(join_type == JoinType::INNER);

	const auto &lop = op.Cast<LogicalJoin>();

	// Make sure we have a consistent list of output columns, ids and types, always

	unordered_map<idx_t, idx_t> build_columns_in_conditions;

	// Collect condition types
	for (idx_t i = 0; i < conditions.size(); i++) {
		auto &cond = conditions[i];
		probe_side_condition_types.push_back(cond.left->return_type);
		build_side_condition_types.push_back(cond.right->return_type);

		// If this is just a reference, we dont need to include it as part of the payload
		if (cond.right->GetExpressionClass() == ExpressionClass::BOUND_REF) {
			build_columns_in_conditions.emplace(cond.right->Cast<BoundReferenceExpression>().index, i);
		}
	}

	// Probe-side
	auto &probe_side_input_types = children[0]->GetTypes();

	probe_side_output_columns = lop.left_projection_map;
	if (probe_side_output_columns.empty()) {
		probe_side_output_columns.reserve(probe_side_input_types.size());
		for (idx_t i = 0; i < probe_side_input_types.size(); i++) {
			probe_side_output_columns.emplace_back(i);
		}
	}

	for (const auto &probe_col_idx : probe_side_output_columns) {
		const auto &probe_type = probe_side_input_types[probe_col_idx];
		probe_side_output_types.push_back(probe_type);
	}

	// For ANTI, SEMI and MARK join, we only need to store the keys, so for these the payload/RHS types are empty
	if (join_type == JoinType::ANTI || join_type == JoinType::SEMI || join_type == JoinType::MARK) {
		return;
	}

	// Build-side
	auto &build_side_input_types = children[1]->GetTypes();

	// This is slightly more advanced, as we need to handle the payload columns.
	auto build_side_projection_map = lop.right_projection_map;
	if (build_side_projection_map.empty()) {
		build_side_projection_map.reserve(build_side_input_types.size());
		for (idx_t i = 0; i < build_side_input_types.size(); i++) {
			build_side_projection_map.emplace_back(i);
		}
	}

	for (auto &build_col_idx : build_side_projection_map) {
		const auto &build_type = build_side_input_types[build_col_idx];

		auto it = build_columns_in_conditions.find(build_col_idx);
		if (it == build_columns_in_conditions.end()) {
			// This build side column is not a join key
			payload_columns.push_back(build_col_idx);
			payload_types.push_back(build_type);
			build_side_output_columns.push_back(build_side_condition_types.size() + payload_types.size() - 1);
		} else {
			// This build side column is a join key
			build_side_output_columns.push_back(it->second);
		}
		build_side_output_types.push_back(build_type);
	}
}

InsertionOrderPreservingMap<string> PhysicalSpatialJoin::ParamsToString() const {
	// TODO: Add condition to the result (GetName is wrong)
	auto result = PhysicalOperator::ParamsToString();
	// result["condition"] = condition->GetName();
	return result;
}

string PhysicalSpatialJoin::GetName() const {
	return "SPATIAL_JOIN";
}

//----------------------------------------------------------------------------------------------------------------------
// Sink Interface
//----------------------------------------------------------------------------------------------------------------------
class SpatialJoinGlobalState final : public GlobalSinkState {
public:
	// TODO: Move the layout to the operator instead of the global state.
	TupleDataLayout layout;
	unique_ptr<TupleDataCollection> collection;

	// This is initialized in the finalize state
	unique_ptr<FlatRTree> rtree = nullptr;
};

unique_ptr<GlobalSinkState> PhysicalSpatialJoin::GetGlobalSinkState(ClientContext &context) const {
	vector<LogicalType> layout_types;

	// Joinkey types
	for (const auto &type : build_side_condition_types) {
		layout_types.push_back(type);
	}

	// Any extra payload types
	for (const auto &type : payload_types) {
		layout_types.push_back(type);
	}

	auto gstate = make_uniq<SpatialJoinGlobalState>();
	gstate->layout.Initialize(layout_types);
	gstate->collection = make_uniq<TupleDataCollection>(BufferManager::GetBufferManager(context), gstate->layout);

	return std::move(gstate);
}

class SpatialJoinLocalState final : public LocalSinkState {
public:
	SpatialJoinLocalState(const PhysicalSpatialJoin &op, ClientContext &context, const TupleDataLayout &layout)
	    : join_key_executor(context) {

		// Dont keep the tuples in memory after appending.
		collection = make_uniq<TupleDataCollection>(BufferManager::GetBufferManager(context), layout);
		collection->InitializeAppend(append_state, TupleDataPinProperties::UNPIN_AFTER_DONE);

		// Initialize join key executor/chunk
		// TODO: Use buffer allocator for chunk
		D_ASSERT(op.conditions.size() == 1);
		join_key_executor.AddExpression(*op.conditions[0].right);
		join_key_chunk.Initialize(context, op.build_side_condition_types);

		// If we have a payload, initialize it
		if (!op.payload_types.empty()) {
			payload_chunk.Initialize(context, op.payload_types);
		}
	}

	TupleDataAppendState append_state;
	unique_ptr<TupleDataCollection> collection;
	ExpressionExecutor join_key_executor;
	DataChunk join_key_chunk;
	DataChunk payload_chunk;
};

unique_ptr<LocalSinkState> PhysicalSpatialJoin::GetLocalSinkState(ExecutionContext &context) const {
	auto &gstate = sink_state->Cast<SpatialJoinGlobalState>();
	auto lstate = make_uniq<SpatialJoinLocalState>(*this, context.client, gstate.layout);
	return std::move(lstate);
}

SinkResultType PhysicalSpatialJoin::Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const {

	auto &gstate = input.global_state.Cast<SpatialJoinGlobalState>();
	auto &lstate = input.local_state.Cast<SpatialJoinLocalState>();

	lstate.join_key_chunk.Reset();
	lstate.join_key_executor.Execute(chunk, lstate.join_key_chunk);

	if (payload_types.empty()) {
		// Put an empty chunk into the collection
		lstate.payload_chunk.SetCardinality(chunk.size());
	} else {
		// Reference the input
		lstate.payload_chunk.ReferenceColumns(chunk, payload_columns);
	}

	// Build a chunk to concatenate the join key with the payload columns
	DataChunk source_chunk;
	source_chunk.InitializeEmpty(gstate.layout.GetTypes());

	const auto joinkey_col_count = lstate.join_key_chunk.ColumnCount();
	const auto payload_col_count = lstate.payload_chunk.ColumnCount();

	for (idx_t i = 0; i < joinkey_col_count; i++) {
		source_chunk.data[i].Reference(lstate.join_key_chunk.data[i]);
	}
	for (idx_t i = 0; i < payload_col_count; i++) {
		source_chunk.data[joinkey_col_count + i].Reference(lstate.payload_chunk.data[i]);
	}
	source_chunk.SetCardinality(chunk);

	// Finally, append the combined chunk to the collection
	lstate.collection->Append(lstate.append_state, source_chunk);

	return SinkResultType::NEED_MORE_INPUT;
}

SinkCombineResultType PhysicalSpatialJoin::Combine(ExecutionContext &context, OperatorSinkCombineInput &input) const {
	auto &gstate = input.global_state.Cast<SpatialJoinGlobalState>();
	auto &lstate = input.local_state.Cast<SpatialJoinLocalState>();

	// Flush the append state
	lstate.collection->FinalizePinState(lstate.append_state.pin_state);

	// Append the local collection to the global collection
	gstate.collection->Combine(*lstate.collection);

	return SinkCombineResultType::FINISHED;
}

// This is where we would build the rtree, by iterating through the tupledata collection we've created
SinkFinalizeType PhysicalSpatialJoin::Finalize(Pipeline &pipeline, Event &event, ClientContext &context,
                                               OperatorSinkFinalizeInput &input) const {
	auto &gstate = input.global_state.Cast<SpatialJoinGlobalState>();

	if (gstate.collection->Count() == 0) {
		return SinkFinalizeType::NO_OUTPUT_POSSIBLE;
	}

	// Initialize the flat R-Tree
	static constexpr auto RTREE_NODE_SIZE = 64;
	gstate.rtree = make_uniq<FlatRTree>(gstate.collection->Count(), RTREE_NODE_SIZE);

	// Now, this is where we build the rtree, by iterating over the tuples in the collection.
	// We need to keep everything pinned so that we can probe the pointers later
	TupleDataChunkIterator iterator(*gstate.collection, TupleDataPinProperties::KEEP_EVERYTHING_PINNED, true);

	const auto rows_ptr = iterator.GetRowLocations();
	Vector row_pointer_vector(LogicalType::POINTER, reinterpret_cast<data_ptr_t>(rows_ptr));

	auto &sel = *FlatVector::IncrementalSelectionVector();

	Vector geom_vec(GeoTypes::GEOMETRY());
	auto geom_ptr = FlatVector::GetData<geometry_t>(geom_vec);

	do {
		const auto row_count = iterator.GetCurrentChunkCount();

		// We only need to fetch the LHS GEOMETRY column to build the rtree.
		// The geometry column is always the first column in the layout.
		gstate.collection->Gather(row_pointer_vector, sel, row_count, 0, geom_vec, sel, nullptr);

		// Push the bounding boxes into the R-Tree
		for (idx_t row_idx = 0; row_idx < row_count; row_idx++) {
			const auto &geom = geom_ptr[row_idx];
			Box2D<double> bbox;
			if (!geom.TryGetCachedBounds(bbox)) {
				// Skip empty geometries
				continue;
			}

			Box2D<float> bbox_f;
			bbox_f.min.x = static_cast<float>(bbox.min.x);
			bbox_f.min.y = static_cast<float>(bbox.min.y);
			bbox_f.max.x = static_cast<float>(bbox.max.x);
			bbox_f.max.y = static_cast<float>(bbox.max.y);

			// Push the bounding box into the R-Tree
			gstate.rtree->Push(bbox_f, rows_ptr[row_idx]);
		}
	} while (iterator.Next());

	// Build the R-Tree
	gstate.rtree->Build();

	return SinkFinalizeType::READY;
}

//----------------------------------------------------------------------------------------------------------------------
// Operator Interface
//----------------------------------------------------------------------------------------------------------------------
// This is where we do the probing of the rtree.
class SpatialJoinOperatorState final : public CachingOperatorState {
public:
	idx_t input_idx = 0;
	idx_t match_idx = 0;

	vector<data_ptr_t> matched_rows;
	FlatRTree::LookupState state = {};

	// Expressions Executor
	ExpressionExecutor join_key_executor;
	DataChunk join_key_chunk;
	DataChunk lhs_output;

	ExpressionExecutor predicate_executor;
	DataChunk predicate_chunk;
	unique_ptr<Expression> predicate_expr;

	DataChunk build_side_joinkey_chunk;

	explicit SpatialJoinOperatorState(ClientContext &context)
	    : join_key_executor(context), predicate_executor(context) {
	}
};

class SpatialJoinGlobalOperatorState final : public GlobalOperatorState {
public:
	unique_ptr<FlatRTree> rtree;
	unique_ptr<TupleDataCollection> collection;
};

unique_ptr<OperatorState> PhysicalSpatialJoin::GetOperatorState(ExecutionContext &context) const {
	// auto &gstate = sink_state->Cast<SpatialJoinGlobalState>();
	// Steal the built rtree from the sink state.
	auto lstate = make_uniq<SpatialJoinOperatorState>(context.client);

	// Add the expression executor for the join key
	D_ASSERT(conditions.size() == 1);
	lstate->join_key_executor.AddExpression(*conditions[0].left);
	lstate->join_key_chunk.Initialize(context.client, probe_side_condition_types);

	if (!probe_side_output_types.empty()) {
		lstate->lhs_output.Initialize(context.client, probe_side_output_types);
	}

	// We also need to add the predicate executor
	lstate->predicate_expr = conditions[0].ToExpr(context.client);

	// HACK: Make this reference the join key chunks
	auto &func_expr = lstate->predicate_expr->Cast<BoundFunctionExpression>();
	func_expr.children[0] = make_uniq<BoundReferenceExpression>(conditions[0].left->return_type, 0);
	func_expr.children[1] = make_uniq<BoundReferenceExpression>(conditions[0].right->return_type, 1);

	lstate->predicate_executor.AddExpression(*lstate->predicate_expr);

	// And initalize the chunk to hold the rhs buildside join key
	lstate->build_side_joinkey_chunk.Initialize(context.client, build_side_condition_types);

	return std::move(lstate);
}

unique_ptr<GlobalOperatorState> PhysicalSpatialJoin::GetGlobalOperatorState(ClientContext &context) const {
	auto &gstate = sink_state->Cast<SpatialJoinGlobalState>();

	auto result = make_uniq<SpatialJoinGlobalOperatorState>();
	// Steal the built rtree from the sink state.
	result->rtree = std::move(gstate.rtree);
	// Steal the tuple data collection
	result->collection = std::move(gstate.collection);

	return std::move(result);
}

OperatorResultType PhysicalSpatialJoin::ExecuteInternal(ExecutionContext &context, DataChunk &input, DataChunk &chunk,
                                                        GlobalOperatorState &gstate_p, OperatorState &state_p) const {

	// Probe the rtree
	auto &gstate = gstate_p.Cast<SpatialJoinGlobalOperatorState>();
	auto &lstate = state_p.Cast<SpatialJoinOperatorState>();

	// Return FINISHED if the rtree is empty
	if (gstate.rtree == nullptr || gstate.rtree->Count() == 0) {
		D_ASSERT(join_type == JoinType::INNER);
		return OperatorResultType::FINISHED;
	}

	const auto &rtree = *gstate.rtree;

	// Reference he probe side output columns
	lstate.lhs_output.ReferenceColumns(input, probe_side_output_columns);

	// Execute the probe side join key expression
	// TODO: This is innefficient to do on each cycle
	if (lstate.input_idx == 0) {
		// Reset and execute the join key expression on each new input chunk
		lstate.join_key_chunk.Reset();
		lstate.join_key_executor.Execute(input, lstate.join_key_chunk);
	}

	// Create unified format for the geometry column
	UnifiedVectorFormat probe_geom_format;
	lstate.join_key_chunk.data[0].ToUnifiedFormat(input.size(), probe_geom_format);

	SelectionVector lhs_sel(chunk.GetCapacity());
	SelectionVector target_sel(chunk.GetCapacity());

	idx_t output_count = chunk.GetCapacity();
	idx_t output_idx = 0;

	while (true) {
		const auto remaining_output = output_count - output_idx;
		const auto remaining_matches = lstate.matched_rows.size() - lstate.match_idx;
		const auto remaining_input = lstate.lhs_output.size() - lstate.input_idx;

		//--------------------------------------------------------------------------------------------------------------
		// Case 1: The output chunk is full
		//--------------------------------------------------------------------------------------------------------------
		if (remaining_output == 0) {
			break;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 2: We have matches to output
		//--------------------------------------------------------------------------------------------------------------
		if (remaining_matches != 0) {

			// Output as many as we can in one batch
			const auto batch_size = MinValue<idx_t>(remaining_output, remaining_matches);

			for (idx_t i = 0; i < batch_size; i++) {
				// Offset the selection vector, so that we gather into the target columns at the correct index
				target_sel.set_index(i, output_idx + i);

				// This output idx corresponds to the current input idx
				// TODO: Dont subtract one here
				lhs_sel.set_index(output_idx + i, lstate.input_idx - 1);
			}

			// Fetch each column from the build side
			const auto row_pointers = reinterpret_cast<data_ptr_t>(lstate.matched_rows.data() + lstate.match_idx);
			Vector pointers(LogicalType::POINTER, row_pointers);
			for (idx_t i = 0; i < build_side_output_columns.size(); i++) {

				// TODO: Is this correct?
				auto &build_side_col_idx = build_side_output_columns[i];
				auto &target = chunk.data[lstate.lhs_output.ColumnCount() + i];

				D_ASSERT(gstate.collection->GetLayout().GetTypes()[build_side_col_idx] == target.GetType());

				// Dont forget to offset by 1 here!
				gstate.collection->Gather(pointers, *FlatVector::IncrementalSelectionVector(), batch_size,
				                          build_side_col_idx, target, target_sel, nullptr);
			}

			// Also collect the build-side join key (at index 0) while were here
			gstate.collection->Gather(pointers, *FlatVector::IncrementalSelectionVector(), batch_size, 0,
			                          lstate.build_side_joinkey_chunk.data[0], target_sel, nullptr);

			// Increment the counters
			output_idx += batch_size;
			lstate.match_idx += batch_size;

			continue; // keep going
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 3: We dont have any more input
		//--------------------------------------------------------------------------------------------------------------
		if (remaining_input == 0) {
			// Also reset the input idx
			lstate.input_idx = 0;
			break;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 4: We need to probe the rtree for more matches
		//--------------------------------------------------------------------------------------------------------------
		const auto geom_idx = probe_geom_format.sel->get_index(lstate.input_idx);
		if (!probe_geom_format.validity.RowIsValid(geom_idx)) {
			// If the geometry is invalid, skip it
			lstate.input_idx++;
			continue;
		}

		const auto &geom = UnifiedVectorFormat::GetData<geometry_t>(probe_geom_format)[geom_idx];
		Box2D<double> bbox;
		if (!geom.TryGetCachedBounds(bbox)) {
			// If the geometry is empty, skip it
			lstate.input_idx++;
			continue;
		}

		Box2D<float> bbox_f;
		bbox_f.min.x = static_cast<float>(bbox.min.x);
		bbox_f.min.y = static_cast<float>(bbox.min.y);
		bbox_f.max.x = static_cast<float>(bbox.max.x);
		bbox_f.max.y = static_cast<float>(bbox.max.y);

		// Probe the R-Tree for new matches and increment the input idx
		lstate.matched_rows.clear();
		rtree.Lookup(bbox_f, lstate.matched_rows);
		lstate.match_idx = 0;

		// Increment the input idx
		lstate.input_idx++;
	}

	//------------------------------------------------------------------------------------------------------------------
	// Finalize the output chunk
	//------------------------------------------------------------------------------------------------------------------
	// When we break out of the loop, we're ready to return.

	// Slice the LHS output columns
	chunk.Slice(lstate.lhs_output, lhs_sel, output_idx);

	// All this is kind of a mess...
	// The idea: Combine the join key we gathered from the build side (first column in the layout)
	// and the join key we just computed for the probe side, and then apply the spatial predicate
	// on both keys.

	DataChunk pred_chunk;
	pred_chunk.InitializeEmpty({probe_side_condition_types[0], build_side_condition_types[0]});
	// Note: We have to slice here again so that we get the right join keys matching the columns
	pred_chunk.data[0].Slice(lstate.join_key_chunk.data[0], lhs_sel, output_idx);
	pred_chunk.data[1].Reference(lstate.build_side_joinkey_chunk.data[0]);
	pred_chunk.SetCardinality(output_idx);

	// Actually apply the spatial predicate, and slice again
	SelectionVector valid_sel(output_idx);
	idx_t filtered = lstate.predicate_executor.SelectExpression(pred_chunk, valid_sel);
	chunk.Slice(valid_sel, filtered);

	// TODO: Because we end up slicing quite a bit here, it might make sense to combine multiple batches into one

	// Verify
	chunk.Verify();

	// Return based on if there are more outputs or if we need more input
	return output_idx == output_count ? OperatorResultType::HAVE_MORE_OUTPUT : OperatorResultType::NEED_MORE_INPUT;
}

} // namespace duckdb