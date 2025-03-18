#pragma once
#include "duckdb/execution/operator/join/physical_join.hpp"
#include "duckdb/planner/operator/logical_join.hpp"

namespace duckdb {

struct SpatialJoinCondition;

class PhysicalSpatialJoin final : public PhysicalJoin {
public:
	static constexpr auto TYPE = PhysicalOperatorType::EXTENSION;

public:
	PhysicalSpatialJoin(LogicalOperator &op, unique_ptr<PhysicalOperator> left, unique_ptr<PhysicalOperator> right,
	                    vector<SpatialJoinCondition> conditions, JoinType join_type, idx_t estimated_cardinality);

	//! The conditions of the join
	vector<SpatialJoinCondition> conditions;

	// The types at each side of the condition, for each condition
	vector<LogicalType> probe_side_condition_types;
	vector<LogicalType> build_side_condition_types;

	//! The indices/types of the left-hand side (probe side) columns that need to be output
	vector<idx_t> probe_side_output_columns;
	vector<LogicalType> probe_side_output_types;

	//! The indices/types of the right-hand side (build side) columns that need to be output
	vector<idx_t> build_side_output_columns;
	vector<LogicalType> build_side_output_types;

	//! The indices/types of the payload columns
	vector<idx_t> payload_columns;
	vector<LogicalType> payload_types;

public:
	// Operator Interface
	unique_ptr<OperatorState> GetOperatorState(ExecutionContext &context) const override;
	unique_ptr<GlobalOperatorState> GetGlobalOperatorState(ClientContext &context) const override;

	bool ParallelOperator() const override {
		return true;
	}

protected:
	// CachingOperatorState Interface
	OperatorResultType ExecuteInternal(ExecutionContext &context, DataChunk &input, DataChunk &chunk,
	                                   GlobalOperatorState &gstate, OperatorState &state) const override;

public:
	// Sink interface
	unique_ptr<GlobalSinkState> GetGlobalSinkState(ClientContext &context) const override;
	unique_ptr<LocalSinkState> GetLocalSinkState(ExecutionContext &context) const override;
	SinkResultType Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const override;
	SinkCombineResultType Combine(ExecutionContext &context, OperatorSinkCombineInput &input) const override;
	SinkFinalizeType Finalize(Pipeline &pipeline, Event &event, ClientContext &context,
	                          OperatorSinkFinalizeInput &input) const override;

	bool IsSink() const override {
		return true;
	}
	bool ParallelSink() const override {
		return true;
	}

public:
	InsertionOrderPreservingMap<string> ParamsToString() const override;
	string GetName() const override;
};

} // namespace duckdb