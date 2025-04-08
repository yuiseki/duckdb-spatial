#pragma once
#include "duckdb/execution/operator/join/physical_join.hpp"
#include "duckdb/planner/operator/logical_join.hpp"

namespace duckdb {

class PhysicalSpatialJoin final : public PhysicalJoin {
public:
	static constexpr auto TYPE = PhysicalOperatorType::EXTENSION;

public:
	PhysicalSpatialJoin(LogicalOperator &op, unique_ptr<PhysicalOperator> left, unique_ptr<PhysicalOperator> right,
	                    unique_ptr<Expression> spatial_predicate, JoinType join_type, idx_t estimated_cardinality);

	//! The condition of the join
	unique_ptr<Expression> condition;
	optional_ptr<Expression> build_side_key;
	optional_ptr<Expression> probe_side_key;

	vector<column_t> build_side_output_columns;
	vector<column_t> probe_side_output_columns;

	vector<LogicalType> probe_side_output_types;
	vector<LogicalType> build_side_output_types;

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


	// CachingOperatorState Interface
	OperatorResultType ExecuteProbeJoin(ExecutionContext &context, DataChunk &input, DataChunk &chunk,
					   GlobalOperatorState &gstate, OperatorState &state) const;
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
	/*
	// Source interface
	unique_ptr<GlobalSourceState> GetGlobalSourceState(ClientContext &context) const override;
	unique_ptr<LocalSourceState> GetLocalSourceState(ExecutionContext &context, GlobalSourceState &gstate) const override;
	SourceResultType GetData(ExecutionContext &context, DataChunk &chunk, OperatorSourceInput &input) const override;

	bool IsSource() const override {
		return PropagatesBuildSide(join_type);
	}

	bool ParallelSource() const override {
		return true;
	}
	*/
public:
	InsertionOrderPreservingMap<string> ParamsToString() const override;
	string GetName() const override;
};

} // namespace duckdb