#include "spatial/operators/spatial_join_physical.hpp"
#include "spatial/geometry/geometry_type.hpp"
#include "spatial/geometry/sgl.hpp"
#include "spatial/spatial_types.hpp"
#include "spatial_join_logical.hpp"

#include "duckdb/common/types/row/tuple_data_collection.hpp"
#include "duckdb/common/types/row/tuple_data_iterator.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/expression/bound_reference_expression.hpp"
#include "duckdb/execution/operator/join/physical_comparison_join.hpp"
#include "duckdb/storage/buffer_manager.hpp"

namespace duckdb {

//======================================================================================================================
// Flat RTree
//======================================================================================================================

namespace {

class FlatRTree {
public:
	using Box = Box2D<float>;

	FlatRTree(uint32_t item_count_p, uint32_t node_size_p) : item_count(item_count_p), node_size(node_size_p) {

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
		rows.resize(item_count);
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

	class LookupState {
		friend class FlatRTree;

	public:
		explicit LookupState() : matches(LogicalType::POINTER) {
		}

	public:
		Vector matches;
		idx_t matches_count = 0;
		idx_t matches_idx = 0;

	private:
		queue<size_t> search_queue;
		Box search_box;
		size_t entry_beg = 0;
		size_t entry_pos = 0;
		bool exhausted = true;
	};

	void InitScan(LookupState &state, const Box &box) const {
		while (!state.search_queue.empty()) {
			state.search_queue.pop();
		}
		state.search_box = box;
		state.entry_beg = boxes.size() - 1;
		state.entry_pos = state.entry_beg;

		state.exhausted = false;
		state.matches_idx = 0;
		state.matches_count = 0;
	}

	bool Scan(LookupState &state) const {
		if (state.exhausted) {
			return false;
		}

		idx_t count = 0;
		const auto ptr = FlatVector::GetData<data_ptr_t>(state.matches);
		Lookup(state, [&](const data_ptr_t &row) {
			ptr[count++] = row;
			return count == STANDARD_VECTOR_SIZE;
		});
		// Set the count of the result vector
		state.matches_count = count;
		state.matches_idx = 0;

		return count > 0;
	}

	template <class CALLBACK>
	void Lookup(LookupState &state, CALLBACK &&callback) const {

		while (true) {

			const auto entry_end = std::min(state.entry_beg + node_size, UpperBound(state.entry_beg));

			while (state.entry_pos < entry_end) {
				if (!state.search_box.Intersects(boxes[state.entry_pos])) {
					state.entry_pos++;
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
					return;
				}
			}

			if (state.search_queue.empty()) {
				// There is no more nodes to search, return false!
				state.exhausted = true;
				return;
			}

			state.entry_beg = state.search_queue.front();
			state.entry_pos = state.entry_beg;
			state.search_queue.pop();
		}
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

} // namespace

//======================================================================================================================
// Physical Spatial Join Operator
//======================================================================================================================

PhysicalSpatialJoin::PhysicalSpatialJoin(LogicalOperator &op, unique_ptr<PhysicalOperator> left,
                                         unique_ptr<PhysicalOperator> right, unique_ptr<Expression> condition_p,
                                         JoinType join_type, idx_t estimated_cardinality)
    : PhysicalJoin(op, PhysicalOperatorType::EXTENSION, join_type, estimated_cardinality),
      condition(std::move(condition_p)) {

	children.push_back(std::move(left));
	children.push_back(std::move(right));

	auto &func = condition->Cast<BoundFunctionExpression>();

	// Extract the probe side and build side join keys
	probe_side_key = func.children[0].get();
	build_side_key = func.children[1].get();

	// Only simple join types are supported
	D_ASSERT(join_type == JoinType::INNER || join_type == JoinType::LEFT || join_type == JoinType::OUTER ||
	         join_type == JoinType::RIGHT);

	// Always make sure we have a consistent order of the output columns, regardless if we have projection maps or not

	const auto &lop = op.Cast<LogicalJoin>();

	// Probe-side
	const auto &probe_side_input_types = children[0]->types;
	probe_side_output_columns = lop.left_projection_map;
	if (probe_side_output_columns.empty()) {
		probe_side_output_columns.reserve(probe_side_input_types.size());
		for (idx_t i = 0; i < probe_side_input_types.size(); i++) {
			probe_side_output_columns.emplace_back(i);
		}
	}

	for (const auto &probe_col_idx : probe_side_output_columns) {
		const auto type = probe_side_input_types[probe_col_idx];
		probe_side_output_types.push_back(type);
	}

	// Build-side
	const auto &build_side_input_types = children[1]->types;
	build_side_output_columns = lop.right_projection_map;
	if (build_side_output_columns.empty()) {
		build_side_output_columns.reserve(build_side_input_types.size());
		for (idx_t i = 0; i < build_side_input_types.size(); i++) {
			build_side_output_columns.emplace_back(i);
		}
	}
	for (const auto &build_col_idx : build_side_output_columns) {
		const auto type = build_side_input_types[build_col_idx];
		build_side_output_types.push_back(type);
	}
}
// Layout helpers
idx_t PhysicalSpatialJoin::BuildSideKeyColumn() const {
	return 0;
}

idx_t PhysicalSpatialJoin::BuildSideMatchColumn() const {
	D_ASSERT(PropagatesBuildSide(join_type));
	return 1;
}

idx_t PhysicalSpatialJoin::BuildSidePayloadBeg() const {
	if (PropagatesBuildSide(join_type)) {
		return 2;
	} else {
		return 1;
	}
}

idx_t PhysicalSpatialJoin::BuildSidePayloadEnd() const {
	if (PropagatesBuildSide(join_type)) {
		return 2 + build_side_output_columns.size();
	} else {
		return 1 + build_side_output_columns.size();
	}
}

InsertionOrderPreservingMap<string> PhysicalSpatialJoin::ParamsToString() const {
	// TODO: Add condition to the result (GetName is wrong)
	auto result = PhysicalOperator::ParamsToString();
	result["condition"] = condition->GetName();
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

	auto gstate = make_uniq<SpatialJoinGlobalState>();

	// Add the build side key to the layout, which is always the first column
	vector<LogicalType> layout_types;
	layout_types.push_back(build_side_key->return_type);

	// If the join type is outer, we need to add a marker column
	if (PropagatesBuildSide(join_type)) {
		layout_types.push_back(LogicalType::BOOLEAN);
	}

	for (const auto &type : build_side_output_types) {
		layout_types.push_back(type);
	}

	gstate->layout.Initialize(layout_types);
	gstate->collection = make_uniq<TupleDataCollection>(BufferManager::GetBufferManager(context), gstate->layout);

	return std::move(gstate);
}

class SpatialJoinLocalState final : public LocalSinkState {
public:
	SpatialJoinLocalState(const PhysicalSpatialJoin &op, ClientContext &context, const TupleDataLayout &layout)
	    : build_side_key_executor(context) {
		// Dont keep the tuples in memory after appending.
		collection = make_uniq<TupleDataCollection>(BufferManager::GetBufferManager(context), layout);
		collection->InitializeAppend(append_state, TupleDataPinProperties::UNPIN_AFTER_DONE);

		build_side_key_executor.AddExpression(*op.build_side_key);

		build_side_key_chunk.Initialize(context, {op.build_side_key->return_type});
		build_side_row_chunk.InitializeEmpty(layout.GetTypes());
	}

	TupleDataAppendState append_state;
	unique_ptr<TupleDataCollection> collection;

	// Owning DataChunk containing the build side join-key
	DataChunk build_side_key_chunk;
	// Referencing DataChunk chunk, having the same layout as the collection
	DataChunk build_side_row_chunk;
	// Used to execute the build side join key expression
	ExpressionExecutor build_side_key_executor;
};

unique_ptr<LocalSinkState> PhysicalSpatialJoin::GetLocalSinkState(ExecutionContext &context) const {
	auto &gstate = sink_state->Cast<SpatialJoinGlobalState>();
	auto lstate = make_uniq<SpatialJoinLocalState>(*this, context.client, gstate.layout);
	return std::move(lstate);
}

SinkResultType PhysicalSpatialJoin::Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const {

	auto &lstate = input.local_state.Cast<SpatialJoinLocalState>();

	// Reset chunks
	lstate.build_side_key_chunk.Reset();
	lstate.build_side_row_chunk.Reset();

	// Execute the build side join key expression
	lstate.build_side_key_executor.Execute(chunk, lstate.build_side_key_chunk);

	// Make the first column of the build side row chunk a reference to the build side join key
	const auto key_offset = BuildSideKeyColumn();
	lstate.build_side_row_chunk.data[key_offset].Reference(lstate.build_side_key_chunk.data[0]);

	// If we have a match column, we need to add it to the build side row chunk
	if (PropagatesBuildSide(join_type)) {
		const auto match_offset = BuildSideMatchColumn();
		lstate.build_side_row_chunk.data[match_offset].Reference(Value::BOOLEAN(false));
	}

	const auto payload_beg = BuildSidePayloadBeg();
	const auto payload_end = BuildSidePayloadEnd();

	for (idx_t i = payload_beg; i < payload_end; i++) {
		const auto col = build_side_output_columns[i - payload_beg];
		lstate.build_side_row_chunk.data[i].Reference(chunk.data[col]);
	}

	// Set the cardinality to match the input
	lstate.build_side_row_chunk.SetCardinality(chunk.size());

	// Sink the build side chunk
	lstate.collection->Append(lstate.append_state, lstate.build_side_row_chunk);

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

	if (gstate.collection->Count() == 0 && EmptyResultIfRHSIsEmpty()) {
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

	do {
		const auto row_count = iterator.GetCurrentChunkCount();

		// We only need to fetch the build-side key column to build the rtree.
		// The key column is always the first column in the layout.
		gstate.collection->Gather(row_pointer_vector, sel, row_count, BuildSideKeyColumn(), geom_vec, sel, nullptr);

		// Get a pointer to what we just gathered
		const auto geom_ptr = FlatVector::GetData<geometry_t>(geom_vec);
		// Push the bounding boxes into the R-Tree
		for (idx_t row_idx = 0; row_idx < row_count; row_idx++) {
			if (FlatVector::IsNull(geom_vec, row_idx)) {
				// Skip null geometries
				continue;
			}

			const auto &geom = geom_ptr[row_idx];
			Box2D<float> bbox;
			if (!geom.TryGetCachedBounds(bbox)) {
				// Skip empty geometries
				continue;
			}

			// Push the bounding box into the R-Tree
			gstate.rtree->Push(bbox, rows_ptr[row_idx]);
		}
	} while (iterator.Next());

	// Build the R-Tree once we've gathered everything
	gstate.rtree->Build();

	return SinkFinalizeType::READY;
}

//----------------------------------------------------------------------------------------------------------------------
// Operator Interface
//----------------------------------------------------------------------------------------------------------------------
// This is where we do the probing of the rtree.

enum class JoinState { START = 0, INIT, PROBE, SCAN, EMIT, EMIT_LHS };

class SpatialJoinLocalOperatorState final : public CachingOperatorState {
public:
	bool is_initialized = false;

	idx_t input_index = 0;
	JoinState state = JoinState::START;

	DataChunk overflow_matches;

	FlatRTree::LookupState scan;

	DataChunk probe_side_row_chunk; // holds the projected lhs columns
	DataChunk probe_side_key_chunk; // holds the lhs probe key
	DataChunk build_side_key_chunk; // holds the rhs build key
	DataChunk match_pred_arg_chunk; // references the lhs probe key and the rhs build key, used to compute the predicate

	ExpressionExecutor join_probe_executor; // used to compute the probe key
	ExpressionExecutor join_match_executor; // used to compute the predicate

	UnifiedVectorFormat probe_side_key_vformat; // used to access the probe side key, after its been computed

	unique_ptr<Expression> match_expr;

	SelectionVector probe_side_source_sel; // maps what output rows correspond to which input lhs rows
	SelectionVector build_side_source_sel; // used when gathering the build side
	SelectionVector build_side_target_sel; // used when gathering the build side
	SelectionVector match_sel;             // used to select the output rows, by executing the predicate

	SelectionVector lhs_match_sel; // this is only used when emitting the lhs side of the join

	uint8_t left_outer_marker[STANDARD_VECTOR_SIZE] = {};

	idx_t build_side_match_offset = 0;
	unsafe_unique_array<data_ptr_t> build_side_pointers = nullptr;

	explicit SpatialJoinLocalOperatorState(ClientContext &context)
	    : join_probe_executor(context), join_match_executor(context), probe_side_source_sel(STANDARD_VECTOR_SIZE),
	      build_side_source_sel(STANDARD_VECTOR_SIZE), build_side_target_sel(STANDARD_VECTOR_SIZE),
	      match_sel(STANDARD_VECTOR_SIZE), lhs_match_sel(STANDARD_VECTOR_SIZE) {

		build_side_pointers = make_unsafe_uniq_array<data_ptr_t>(STANDARD_VECTOR_SIZE);
	}
};

class SpatialJoinGlobalOperatorState final : public GlobalOperatorState {
public:
	unique_ptr<FlatRTree> rtree;
	unique_ptr<TupleDataCollection> collection;
};

unique_ptr<OperatorState> PhysicalSpatialJoin::GetOperatorState(ExecutionContext &context) const {
	auto lstate = make_uniq<SpatialJoinLocalOperatorState>(context.client);

	// Create a match expression using the condition, that will be used to filter the results
	lstate->match_expr = condition->Copy();
	auto &func_expr = lstate->match_expr->Cast<BoundFunctionExpression>();
	func_expr.children[0] = make_uniq<BoundReferenceExpression>(probe_side_key->return_type, 0);
	func_expr.children[1] = make_uniq<BoundReferenceExpression>(build_side_key->return_type, 1);

	lstate->join_match_executor.AddExpression(*lstate->match_expr);

	// Add the probe side join key expression
	lstate->join_probe_executor.AddExpression(*probe_side_key);

	// The chunks we need for the join
	lstate->probe_side_row_chunk.Initialize(context.client, probe_side_output_types);
	lstate->probe_side_key_chunk.Initialize(context.client, {probe_side_key->return_type});
	lstate->build_side_key_chunk.Initialize(context.client, {build_side_key->return_type});
	lstate->match_pred_arg_chunk.Initialize(context.client, {probe_side_key->return_type, build_side_key->return_type});

	if (PropagatesBuildSide(join_type)) {
		lstate->build_side_match_offset =
		    sink_state->Cast<SpatialJoinGlobalState>().layout.GetOffsets()[BuildSideMatchColumn()];
	}

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
                                                        GlobalOperatorState &gstate_p, OperatorState &lstate_p) const {
	auto &gstate = gstate_p.Cast<SpatialJoinGlobalOperatorState>();
	auto &lstate = lstate_p.Cast<SpatialJoinLocalOperatorState>();

	idx_t output_index = 0;
	idx_t output_count = chunk.GetCapacity();

	while (true) {
		switch (lstate.state) {
		//--------------------------------------------------------------------------------------------------------------
		// START
		//--------------------------------------------------------------------------------------------------------------
		case JoinState::START: {
			// Check if the build side is empty
			if (gstate.rtree == nullptr || gstate.rtree->Count() == 0) {
				if (EmptyResultIfRHSIsEmpty()) {
					return OperatorResultType::FINISHED;
				}

				// TODO: Add test for this (empty rhs with lhs projection map)

				// Slice the input chunk to what we need
				lstate.probe_side_row_chunk.ReferenceColumns(input, probe_side_output_columns);

				// Construct empty join result
				PhysicalComparisonJoin::ConstructEmptyJoinResult(join_type, false, lstate.probe_side_row_chunk, chunk);

				return OperatorResultType::NEED_MORE_INPUT;
			}

			lstate.state = JoinState::INIT;
		} // fall through
		//--------------------------------------------------------------------------------------------------------------
		// INIT
		//--------------------------------------------------------------------------------------------------------------
		case JoinState::INIT: {
			// We have a new fresh input chunk
			lstate.join_probe_executor.Execute(input, lstate.probe_side_key_chunk);
			lstate.probe_side_key_chunk.data[0].ToUnifiedFormat(input.size(), lstate.probe_side_key_vformat);

			// Reference the columns that we actually care about
			lstate.probe_side_row_chunk.ReferenceColumns(input, probe_side_output_columns);

			// zero miss vector
			memset(lstate.left_outer_marker, 0, sizeof(lstate.left_outer_marker));

			// Reset the input index and move on to the next state
			lstate.input_index = 0;
			lstate.state = JoinState::PROBE;

		} // fall through
		//--------------------------------------------------------------------------------------------------------------
		// PROBE
		//--------------------------------------------------------------------------------------------------------------
		case JoinState::PROBE: {
			if (lstate.input_index == input.size()) {
				lstate.state = JoinState::EMIT;
				continue;
			}

			const auto geom_idx = lstate.probe_side_key_vformat.sel->get_index(lstate.input_index);
			if (!lstate.probe_side_key_vformat.validity.RowIsValid(geom_idx)) {
				lstate.input_index++;
				continue;
			}

			const auto geom_ptr = FlatVector::GetData<geometry_t>(lstate.probe_side_key_chunk.data[0]);
			const auto &geom = geom_ptr[lstate.input_index];

			Box2D<float> bbox;
			if (!geom.TryGetCachedBounds(bbox)) {
				lstate.input_index++;
				continue;
			}

			gstate.rtree->InitScan(lstate.scan, bbox);

			if (!gstate.rtree->Scan(lstate.scan)) {
				lstate.input_index++;
				continue;
			}

			lstate.state = JoinState::SCAN;
		} // fall through
		//--------------------------------------------------------------------------------------------------------------
		// SCAN
		//--------------------------------------------------------------------------------------------------------------
		case JoinState::SCAN: {
			const auto matches_remaining = lstate.scan.matches_count - lstate.scan.matches_idx;
			if (matches_remaining == 0) {
				// We are out of matches. Try to get the next probe
				if (gstate.rtree->Scan(lstate.scan)) {
					continue;
				}
				// Otherwise, we are done with this probe
				lstate.input_index++;
				lstate.state = JoinState::PROBE;
				continue;
			}

			const auto output_remaining = output_count - output_index;
			const auto scan_count = MinValue(output_remaining, matches_remaining);

			D_ASSERT(scan_count != 0);

			for (idx_t i = 0; i < scan_count; i++) {
				// These control what we gather from the build side
				lstate.build_side_target_sel.set_index(i, output_index + i);
				lstate.build_side_source_sel.set_index(i, lstate.scan.matches_idx + i);

				// This stores which probe side row we used to gather the build side
				lstate.probe_side_source_sel.set_index(output_index + i, lstate.input_index);
			}

			// Fetch each column from the build side
			const auto probe_side_col_count = probe_side_output_columns.size();
			const auto build_side_key_col = BuildSideKeyColumn();
			const auto payload_offset = BuildSidePayloadBeg();

			auto &row_pointers = lstate.scan.matches;

			// Collect the build side join key (at index 0)
			gstate.collection->Gather(row_pointers, lstate.build_side_source_sel, scan_count, build_side_key_col,
			                          lstate.build_side_key_chunk.data[0], lstate.build_side_target_sel, nullptr);

			// Now, lets collect the rest of the build side columns
			for (idx_t i = 0; i < build_side_output_columns.size(); i++) {

				const auto build_side_col_idx = payload_offset + i;
				auto &target = chunk.data[probe_side_col_count + i];

				D_ASSERT(gstate.collection->GetLayout().GetTypes()[build_side_col_idx] == target.GetType());

				gstate.collection->Gather(row_pointers, lstate.build_side_source_sel, scan_count, build_side_col_idx,
				                          target, lstate.build_side_target_sel, nullptr);
			}

			// Also collect the build side row pointers (if we have a match column)
			if (IsRightOuterJoin(join_type)) {
				const auto ptrs = FlatVector::GetData<data_ptr_t>(row_pointers);
				for (idx_t i = 0; i < scan_count; i++) {
					lstate.build_side_pointers[output_index + i] = ptrs[i];
				}
			}

			// Increment the output and match index
			output_index += scan_count;
			lstate.scan.matches_idx += scan_count;

			if (output_index != output_count) {
				// We still have space left. Scan more!
				continue;
			}

			lstate.state = JoinState::EMIT;
		} // fall through
		//--------------------------------------------------------------------------------------------------------------
		// EMIT
		//--------------------------------------------------------------------------------------------------------------
		case JoinState::EMIT: {

			// Start by adding the lhs columns
			chunk.Slice(lstate.probe_side_row_chunk, lstate.probe_side_source_sel, output_index);

			// Now, lets actually evaluate the predicate
			lstate.match_pred_arg_chunk.data[0].Slice(lstate.probe_side_key_chunk.data[0], lstate.probe_side_source_sel,
			                                          output_index);
			lstate.match_pred_arg_chunk.data[1].Reference(lstate.build_side_key_chunk.data[0]);
			lstate.match_pred_arg_chunk.SetCardinality(output_index);

			const auto filtered =
			    lstate.join_match_executor.SelectExpression(lstate.match_pred_arg_chunk, lstate.match_sel);

			if (IsLeftOuterJoin(join_type)) {
				for (idx_t i = 0; i < filtered; i++) {
					// This is kinda crazy
					// We're first selecting all the output rows that matched.
					// And then figuring out which lhs row that belongs to, by using the probe_side_source_sel vector
					const auto match_source_idx = lstate.match_sel.get_index(i);
					const auto probe_source_idx = lstate.probe_side_source_sel.get_index(match_source_idx);
					lstate.left_outer_marker[probe_source_idx] = 1;
				}
			}

			if (IsRightOuterJoin(join_type)) {
				// Mark each row in the build side as being matched.
				// We need to do this so we dont emit them again in the right-outer join phase
				for (idx_t i = 0; i < filtered; i++) {
					const auto match_source_idx = lstate.match_sel.get_index(i);
					const auto data_ptr =
					    (lstate.build_side_pointers[match_source_idx] + lstate.build_side_match_offset);
					const auto bool_ptr = reinterpret_cast<bool *>(data_ptr);

					*bool_ptr = true;
				}
			}

			chunk.Slice(lstate.match_sel, filtered);

			if (lstate.input_index != input.size()) {
				// We still have more input rows to process
				lstate.state = JoinState::SCAN;
				return OperatorResultType::HAVE_MORE_OUTPUT;
			}

			if (IsLeftOuterJoin(join_type)) {
				// Before we can ask for more input, we need to emit the outer left side
				// But we need a clean output chunk, so we need to return here (cant just fall through)

				lstate.state = JoinState::EMIT_LHS;
				return OperatorResultType::HAVE_MORE_OUTPUT;
			}

			// Otherwise, we need to wait for more input
			lstate.state = JoinState::INIT;
			return OperatorResultType::NEED_MORE_INPUT;
		}
		//--------------------------------------------------------------------------------------------------------------
		// EMIT LEFT OUTER
		//--------------------------------------------------------------------------------------------------------------
		case JoinState::EMIT_LHS: {
			// Emit the outer left side

			idx_t remaining_count = 0;
			for (idx_t i = 0; i < input.size(); i++) {
				if (!lstate.left_outer_marker[i]) {
					lstate.lhs_match_sel.set_index(remaining_count++, i);
				}
			}

			if (remaining_count > 0) {
				chunk.Slice(lstate.probe_side_row_chunk, lstate.lhs_match_sel, remaining_count);

				// Null the RHS columns
				for (idx_t i = 0; i < build_side_output_columns.size(); i++) {
					auto &target = chunk.data[probe_side_output_columns.size() + i];
					target.SetVectorType(VectorType::CONSTANT_VECTOR);
					ConstantVector::SetNull(target, true);
				}
			}

			lstate.state = JoinState::INIT;
			return OperatorResultType::NEED_MORE_INPUT;
		}
		default:
			D_ASSERT(false);
			break;
		}
	}
}

//----------------------------------------------------------------------------------------------------------------------
// Source Interface
//----------------------------------------------------------------------------------------------------------------------
class SpatialJoinGlobalSourceState final : public GlobalSourceState {
public:
	explicit SpatialJoinGlobalSourceState(const PhysicalSpatialJoin &op) : op(op) {
		D_ASSERT(op.sink_state);
		const auto &state = op.op_state->Cast<SpatialJoinGlobalOperatorState>();

		// Initialize a parallel scan
		vector<column_t> column_ids;
		// We dont want to scan the first column (which contains the build-side arg)
		// But we do want to scan the match column
		column_ids.push_back(op.BuildSideMatchColumn());
		for (idx_t i = op.BuildSidePayloadBeg(); i < op.BuildSidePayloadEnd(); i++) {
			column_ids.push_back(i);
		}

		// We dont need to keep the tuples aroun after scanning
		state.collection->InitializeScan(scan_state, std::move(column_ids), TupleDataPinProperties::DESTROY_AFTER_DONE);

		tuples_maximum = state.collection->Count();
	}

	const PhysicalSpatialJoin &op;
	TupleDataParallelScanState scan_state;

	// How many tuples we have scanned so far
	idx_t tuples_maximum = 0;
	atomic<idx_t> tuples_scanned = {0};

public:
	idx_t MaxThreads() override {
		const auto &state = op.op_state->Cast<SpatialJoinGlobalOperatorState>();
		const auto count = state.collection->Count();

		// Rough approximation of the number of threads to use
		return count / (STANDARD_VECTOR_SIZE * 10ULL);
	}
};

class SpatialJoinLocalSourceState final : public LocalSourceState {
public:
	explicit SpatialJoinLocalSourceState(const PhysicalSpatialJoin &op) {
		D_ASSERT(op.sink_state);
		const auto &state = op.op_state->Cast<SpatialJoinGlobalOperatorState>();

		// Initialize a parallel scan
		vector<column_t> column_ids;
		// We dont want to scan the first column (which contains the build-side arg)
		// But we do want to scan the match column
		column_ids.push_back(op.BuildSideMatchColumn());
		for (idx_t i = op.BuildSidePayloadBeg(); i < op.BuildSidePayloadEnd(); i++) {
			column_ids.push_back(i);
		}

		// We dont need to keep the tuples aroun after scanning
		state.collection->InitializeScan(scan_state, std::move(column_ids), TupleDataPinProperties::DESTROY_AFTER_DONE);
		state.collection->InitializeScanChunk(scan_state, scan_chunk);
	}

	TupleDataLocalScanState scan_state;
	DataChunk scan_chunk;
};

unique_ptr<GlobalSourceState> PhysicalSpatialJoin::GetGlobalSourceState(ClientContext &context) const {
	auto gstate = make_uniq<SpatialJoinGlobalSourceState>(*this);
	return std::move(gstate);
}

unique_ptr<LocalSourceState> PhysicalSpatialJoin::GetLocalSourceState(ExecutionContext &context,
                                                                      GlobalSourceState &gstate_p) const {
	auto lstate = make_uniq<SpatialJoinLocalSourceState>(*this);
	return std::move(lstate);
}

SourceResultType PhysicalSpatialJoin::GetData(ExecutionContext &context, DataChunk &chunk,
                                              OperatorSourceInput &input) const {
	D_ASSERT(PropagatesBuildSide(join_type));

	auto &gstate = input.global_state.Cast<SpatialJoinGlobalSourceState>();
	auto &lstate = input.local_state.Cast<SpatialJoinLocalSourceState>();

	const auto &tuples = gstate.op.op_state->Cast<SpatialJoinGlobalOperatorState>().collection;

	SelectionVector sel(STANDARD_VECTOR_SIZE);

	while (tuples->Scan(gstate.scan_state, lstate.scan_state, lstate.scan_chunk)) {
		gstate.tuples_scanned += lstate.scan_chunk.size();

		const auto matches = FlatVector::GetData<bool>(lstate.scan_chunk.data[0]);

		idx_t result_count = 0;
		for (idx_t i = 0; i < lstate.scan_chunk.size(); i++) {
			if (!matches[i]) {
				sel.set_index(result_count++, i);
			}
		}

		if (result_count > 0) {

			const auto lhs_col_count = probe_side_output_columns.size();
			const auto rhs_col_count = build_side_output_columns.size();

			// Null the LHS columns
			for (idx_t i = 0; i < lhs_col_count; i++) {
				auto &target = chunk.data[i];
				target.SetVectorType(VectorType::CONSTANT_VECTOR);
				ConstantVector::SetNull(target, true);
			}

			// Set the RHS columns
			for (idx_t i = 0; i < rhs_col_count; i++) {
				auto &target = chunk.data[lhs_col_count + i];
				// Offset by one here to skip the match column
				target.Slice(lstate.scan_chunk.data[i + 1], sel, result_count);
			}

			chunk.SetCardinality(result_count);
			return SourceResultType::HAVE_MORE_OUTPUT;
		}
	}

	return SourceResultType::FINISHED;
}

//----------------------------------------------------------------------------------------------------------------------
// Misc
//----------------------------------------------------------------------------------------------------------------------

ProgressData PhysicalSpatialJoin::GetProgress(ClientContext &context, GlobalSourceState &gstate) const {
	const auto &state = gstate.Cast<SpatialJoinGlobalSourceState>();
	ProgressData res;
	if (state.tuples_maximum) {
		res.done = state.tuples_scanned.load();
		res.total = state.tuples_maximum;
	} else {
		res.SetInvalid();
	}
	return res;
}

} // namespace duckdb