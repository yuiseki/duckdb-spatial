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
		explicit LookupState() : matches(LogicalType::POINTER) { }
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
		while(!state.search_queue.empty()) {
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
		if(state.exhausted) {
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

	// Only inner joins for now!
	D_ASSERT(join_type == JoinType::INNER);

	// const auto &lop = op.Cast<LogicalJoin>();
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
	vector<LogicalType> layout_types;

	// Add the type of the build-side join key
	layout_types.push_back(build_side_key->return_type);

	// And then add the rest of the types in the RHS
	for (const auto &type : children[1]->types) {
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
	    : join_build_executor(context) {

		// Dont keep the tuples in memory after appending.
		collection = make_uniq<TupleDataCollection>(BufferManager::GetBufferManager(context), layout);
		collection->InitializeAppend(append_state, TupleDataPinProperties::UNPIN_AFTER_DONE);

		// Initialize join key executor/chunk
		// TODO: Use buffer allocator for chunk
		join_build_executor.AddExpression(*op.build_side_key);
		join_build_chunk.Initialize(context, {op.build_side_key->return_type});

		append_chunk.InitializeEmpty(layout.GetTypes());
	}

	TupleDataAppendState append_state;
	unique_ptr<TupleDataCollection> collection;
	ExpressionExecutor join_build_executor;
	DataChunk join_build_chunk;

	// The append chunk just references the vectors from the join key chunk and the input chunk
	DataChunk append_chunk;
};

unique_ptr<LocalSinkState> PhysicalSpatialJoin::GetLocalSinkState(ExecutionContext &context) const {
	auto &gstate = sink_state->Cast<SpatialJoinGlobalState>();
	auto lstate = make_uniq<SpatialJoinLocalState>(*this, context.client, gstate.layout);
	return std::move(lstate);
}

SinkResultType PhysicalSpatialJoin::Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const {

	// auto &gstate = input.global_state.Cast<SpatialJoinGlobalState>();
	auto &lstate = input.local_state.Cast<SpatialJoinLocalState>();

	// Execute the join key expression
	lstate.join_build_chunk.Reset();
	lstate.join_build_executor.Execute(chunk, lstate.join_build_chunk);

	// Now combine the join key with the payload columns
	lstate.append_chunk.Reset();

	// Reference the join key column
	lstate.append_chunk.data[0].Reference(lstate.join_build_chunk.data[0]);

	// Reference the rest of the input columns
	for(idx_t i = 0; i < chunk.ColumnCount(); i++) {
		lstate.append_chunk.data[i + 1].Reference(chunk.data[i]);
	}

	// Append the combined chunk to the collection
	lstate.append_chunk.SetCardinality(chunk);
	lstate.collection->Append(lstate.append_state, lstate.append_chunk);

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

	do {
		const auto row_count = iterator.GetCurrentChunkCount();

		// We only need to fetch the build-side key column to build the rtree.
		// The key column is always the first column in the layout.
		gstate.collection->Gather(row_pointer_vector, sel, row_count, 0, geom_vec, sel, nullptr);

		// Get a pointer to what we just gathered
		const auto geom_ptr = FlatVector::GetData<geometry_t>(geom_vec);
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

	// Build the R-Tree once we've gathered everything
	gstate.rtree->Build();

	return SinkFinalizeType::READY;
}

//----------------------------------------------------------------------------------------------------------------------
// Operator Interface
//----------------------------------------------------------------------------------------------------------------------
// This is where we do the probing of the rtree.
class SpatialJoinLocalOperatorState final : public CachingOperatorState {
public:
	idx_t input_idx = 0;

	FlatRTree::LookupState state;

	DataChunk join_probe_chunk; // holds the lhs probe key
	DataChunk join_build_chunk; // holds the rhs build key
	DataChunk join_match_chunk; // references the lhs probe key and the rhs build key, used to compute the predicate

	ExpressionExecutor join_probe_executor;
	ExpressionExecutor join_match_executor;

	UnifiedVectorFormat join_probe_vformat;

	unique_ptr<Expression> match_expr;

	SelectionVector probe_side_source_sel;
	SelectionVector build_side_source_sel;
	SelectionVector build_side_target_sel;
	SelectionVector match_sel;

	explicit SpatialJoinLocalOperatorState(ClientContext &context)
	    : join_probe_executor(context), join_match_executor(context),
		probe_side_source_sel(STANDARD_VECTOR_SIZE),
		build_side_source_sel(STANDARD_VECTOR_SIZE),
		build_side_target_sel(STANDARD_VECTOR_SIZE),
		match_sel(STANDARD_VECTOR_SIZE) {
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

	// Also add the probe side joinkey expression
	lstate->join_probe_executor.AddExpression(*probe_side_key);

	// The chunks we need for the join
	lstate->join_probe_chunk.Initialize(context.client, {probe_side_key->return_type});
	lstate->join_build_chunk.Initialize(context.client, {build_side_key->return_type});
	lstate->join_match_chunk.Initialize(context.client, {probe_side_key->return_type, build_side_key->return_type});

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

	if (gstate.rtree == nullptr || gstate.rtree->Count() == 0) {
		return OperatorResultType::FINISHED;
	}

	const auto &rtree = *gstate.rtree;

	// Reset and execute join key expression on each new input chunk
	if(lstate.input_idx == 0) {
		const auto remaining_matches = lstate.state.matches_count - lstate.state.matches_idx;
		// Only reset if this chunk is _actually_ new
		if(remaining_matches == 0) {
			lstate.join_probe_chunk.Reset();
			lstate.join_probe_executor.Execute(input, lstate.join_probe_chunk);
			lstate.join_build_chunk.data[0].ToUnifiedFormat(input.size(), lstate.join_probe_vformat);
		}
	}

	idx_t output_count = chunk.GetCapacity();
	idx_t output_index = 0;

	while (true) {
		//--------------------------------------------------------------------------------------------------------------
		// Case 1: The output chunk is full
		//--------------------------------------------------------------------------------------------------------------
		const auto remaining_output = output_count - output_index;
		if(remaining_output == 0) {
			break;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 2: We have matches to output
		//--------------------------------------------------------------------------------------------------------------
		const auto remaining_matches = lstate.state.matches_count - lstate.state.matches_idx;
		if(remaining_matches != 0) {

			const auto batch_size = MinValue<idx_t>(remaining_output, remaining_matches);

			for(idx_t i = 0; i < batch_size; i++) {
				lstate.build_side_target_sel.set_index(i, output_index + 1);
				lstate.build_side_source_sel.set_index(i, lstate.state.matches_idx + i);

				lstate.probe_side_source_sel.set_index(output_index + i, lstate.input_idx - 1);
			}

			auto &pointers = lstate.state.matches;

			// Fetch each column from the build side
			//

			// Also collect the build side join key (at index 0) while were here
			gstate.collection->Gather(pointers,
				lstate.build_side_source_sel,
				batch_size,
				0,
				lstate.join_build_chunk.data[0],
				lstate.build_side_target_sel,
				nullptr
			);

			output_index +=	batch_size;
			lstate.state.matches_idx += batch_size;

			// Keep going
			continue;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 3:
		//--------------------------------------------------------------------------------------------------------------
		if(rtree.Scan(lstate.state)) {
			// We got a next set of matches, so continue to output them
			continue;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 4: We dont have any more input
		//--------------------------------------------------------------------------------------------------------------
		const auto remaining_input = input.size() - lstate.input_idx;
		if(remaining_input == 0) {
			// We are out of input
			lstate.input_idx = 0;
			break;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 5: We need to initialize the next probe batch
		//--------------------------------------------------------------------------------------------------------------
		const auto &probe_geom = lstate.join_probe_vformat;
		const auto geom_idx = probe_geom.sel->get_index(lstate.input_idx);
		// TODO: Make an inner loop here so we dont have to reset the whole state machine
		if (!probe_geom.validity.RowIsValid(geom_idx)) {
			// If the geometry is invalid, skip it
			lstate.input_idx++;
			continue;
		}

		const auto &geom = UnifiedVectorFormat::GetData<geometry_t>(probe_geom)[geom_idx];
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
		rtree.InitScan(lstate.state, bbox_f);
		rtree.Scan(lstate.state);

		// Increment the input idx
		lstate.input_idx++;
	}

	chunk.Slice(lstate.probe_side_source_sel, output_index);

	lstate.join_match_chunk.data[0].Slice(lstate.join_probe_chunk.data[0], lstate.probe_side_source_sel, output_index);
	lstate.join_match_chunk.data[1].Reference(lstate.join_build_chunk.data[0]);
	lstate.join_match_chunk.SetCardinality(output_index);

	const auto filtered = lstate.join_match_executor.SelectExpression(lstate.join_match_chunk, lstate.match_sel);
	chunk.Slice(lstate.match_sel, filtered);

	chunk.Verify();

	return output_index == output_count ? OperatorResultType::HAVE_MORE_OUTPUT : OperatorResultType::NEED_MORE_INPUT;
}

/*
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

	idx_t incoming_index = 0;
	idx_t incoming_count = input.size();

	idx_t outgoing_index = 0;
	idx_t outgoing_count = chunk.GetCapacity();

	idx_t matching_index = 0;
	idx_t matching_count = 0;

	Vector build_side_pointers(LogicalType::POINTER);

	// Used to select the rows from the build side
	SelectionVector build_side_target_sel(STANDARD_VECTOR_SIZE);
	SelectionVector build_side_source_sel(STANDARD_VECTOR_SIZE);
	SelectionVector probe_side_probed_sel(STANDARD_VECTOR_SIZE);

	// Used to apply the exact join predicate on the final chunk
	SelectionVector exact_join_filter_sel(STANDARD_VECTOR_SIZE);
	DataChunk join_chunk;

	while(true) {
		if(outgoing_index == outgoing_count) {
			// There is no more space to output
			break;
		}

		if(matching_index == matching_count) {
			// We are out of matches, try to get another batch
			if(!rtree.Scan(lstate.state)) {

				// We got a next set of matches, so continue to output them
				while(incoming_index != incoming_count) {

				}
				if(incoming_index == incoming_count) {
					// We are out of input
					break;
				}
			}
		}

		// Emit as many matches as we can
		const auto remaining_matching = matching_count - matching_index;
		const auto remaining_outgoing = outgoing_count - outgoing_index;

		const auto remaining = MinValue(remaining_matching, remaining_outgoing);

		for(idx_t i = 0; i < remaining; i++) {
			// We select from our matching index, and output to the outgoing index
			build_side_source_sel.set_index(i, matching_index + i);
			build_side_target_sel.set_index(i, outgoing_index + i);

			// We also need to select the corresponding row from the probe side
			probe_side_probed_sel.set_index(i, incoming_index - 1);
		}

		// Gather the build side join key
		gstate.collection->Gather(build_side_pointers, build_side_source_sel, remaining, 0, join_chunk.data[0], build_side_target_sel, nullptr);

		// Gather the rest of the build side columns
		for (idx_t i = 0; i < build_side_output_columns.size(); i++) {
			// TODO: Is this correct?
			auto &build_side_col_idx = build_side_output_columns[i];
			auto &target = chunk.data[lstate.lhs_output.ColumnCount() + i];

			D_ASSERT(gstate.collection->GetLayout().GetTypes()[build_side_col_idx] == target.GetType());

			gstate.collection->Gather(build_side_pointers, build_side_source_sel, remaining,
									  build_side_col_idx, target, build_side_target_sel, nullptr);
		}

		outgoing_index += remaining;
		matching_index += remaining;
	}

	// Filter the output using the exact predicate
	auto exact_join_count = lstate.predicate_executor.SelectExpression(join_chunk, exact_join_filter_sel);
	chunk.Slice(exact_join_filter_sel, exact_join_count);

	// Return HAVE_MORE_OUTPUT if we have more output
	const auto is_full = outgoing_index == outgoing_count;
	return is_full ? OperatorResultType::HAVE_MORE_OUTPUT : OperatorResultType::NEED_MORE_INPUT;
}
*/
static OperatorResultType EExecuteInternal(ExecutionContext &context, DataChunk &input, DataChunk &chunk,
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
	if (lstate.input_idx == 0 && (lstate.state.matches_count - lstate.state.matches_idx) == 0) {
		// Reset and execute the join key expression on each new input chunk
		lstate.join_key_chunk.Reset();
		lstate.join_key_executor.Execute(input, lstate.join_key_chunk);
	}


	// Create unified format for the geometry column
	UnifiedVectorFormat probe_geom_format;
	lstate.join_key_chunk.data[0].ToUnifiedFormat(input.size(), probe_geom_format);

	SelectionVector lhs_sel(chunk.GetCapacity());
	SelectionVector target_sel(chunk.GetCapacity());
	SelectionVector source_sel(chunk.GetCapacity());

	idx_t output_count = chunk.GetCapacity();
	idx_t output_idx = 0;

	while (true) {
		//--------------------------------------------------------------------------------------------------------------
		// Case 1: The output chunk is full
		//--------------------------------------------------------------------------------------------------------------
		const auto remaining_output = output_count - output_idx;
		if (remaining_output == 0) {
			break;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 2: We have matches to output
		//--------------------------------------------------------------------------------------------------------------
		const auto remaining_matches = lstate.state.matches_count - lstate.state.matches_idx;
		if (remaining_matches != 0) {

			// Output as many as we can in one batch
			const auto batch_size = MinValue<idx_t>(remaining_output, remaining_matches);

			for (idx_t i = 0; i < batch_size; i++) {
				// Offset the selection vector, so that we gather into the target columns at the correct index
				target_sel.set_index(i, output_idx + i);

				// This output idx corresponds to the current input idx
				// TODO: Dont subtract one here
				lhs_sel.set_index(output_idx + i, lstate.input_idx - 1);

				source_sel.set_index(i, lstate.state.matches_idx + i);
			}

			// Fetch each column from the build side
			auto &pointers = lstate.state.matches;
			for (idx_t i = 0; i < build_side_output_columns.size(); i++) {

				// TODO: Is this correct?
				auto &build_side_col_idx = build_side_output_columns[i];
				auto &target = chunk.data[lstate.lhs_output.ColumnCount() + i];

				D_ASSERT(gstate.collection->GetLayout().GetTypes()[build_side_col_idx] == target.GetType());

				gstate.collection->Gather(pointers, source_sel, batch_size,
				                          build_side_col_idx, target, target_sel, nullptr);
			}

			// Also collect the build-side join key (at index 0) while were here
			gstate.collection->Gather(pointers, source_sel, batch_size, 0,
			                          lstate.build_side_joinkey_chunk.data[0], target_sel, nullptr);

			// Increment the counters
			output_idx += batch_size;
			lstate.state.matches_idx += batch_size;

			continue; // keep going
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 3: We need to get the next lookup batch
		//--------------------------------------------------------------------------------------------------------------
		// Get the next set of matches from the R-Tree
		if(rtree.Scan(lstate.state)) {
			// We got a next set of matches, so continue to output them
			continue;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 4: We dont have any more input
		//--------------------------------------------------------------------------------------------------------------
		const auto remaining_input = lstate.lhs_output.size() - lstate.input_idx;
		if (remaining_input == 0) {
			// Also reset the input idx
			lstate.input_idx = 0;
			break;
		}

		//--------------------------------------------------------------------------------------------------------------
		// Case 5: We need to initialize the next probe
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
		rtree.InitScan(lstate.state, bbox_f);
		rtree.Scan(lstate.state);

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