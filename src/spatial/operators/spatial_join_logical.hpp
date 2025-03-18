#pragma once

#include "duckdb/planner/operator/logical_extension_operator.hpp"

namespace duckdb {

struct SpatialJoinCondition {
	unique_ptr<Expression> left;
	unique_ptr<Expression> right;
	string predicate;

	// TODO: Make non-const, destroy children
	unique_ptr<Expression> ToExpr(ClientContext &context) const;
};

class LogicalSpatialJoin final : public LogicalExtensionOperator {
public:
	static constexpr auto TYPE = LogicalOperatorType::LOGICAL_EXTENSION_OPERATOR;

public:
	explicit LogicalSpatialJoin(JoinType join_type_p);

	vector<ColumnBinding> GetColumnBindings() override;
	void ResolveColumnBindings(ColumnBindingResolver &res, vector<ColumnBinding> &bindings) override;
	bool HasProjectionMap() const override {
		return !left_projection_map.empty() || !right_projection_map.empty();
	}

	string GetName() const override;
	string GetExtensionName() const override;
	unique_ptr<PhysicalOperator> CreatePlan(ClientContext &context, PhysicalPlanGenerator &generator) override;

	//! The type of the join (INNER, OUTER, etc...)
	JoinType join_type;
	//! Table index used to refer to the MARK column (in case of a MARK join)
	idx_t mark_index {};
	//! The conditions of the join
	vector<SpatialJoinCondition> conditions;
	//! The columns of the LHS that are output by the join
	vector<idx_t> left_projection_map;
	//! The columns of the RHS that are output by the join
	vector<idx_t> right_projection_map;
	//! Join Keys statistics (optional)
	vector<unique_ptr<BaseStatistics>> join_stats;

protected:
	void ResolveTypes() override;
};

} // namespace duckdb