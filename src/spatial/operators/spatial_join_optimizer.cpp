#include "spatial_join_optimizer.hpp"

#include "duckdb/main/database.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/operator/logical_any_join.hpp"
#include "spatial_join_logical.hpp"

#include <duckdb/planner/operator/logical_filter.hpp>

namespace duckdb {

static void InsertSpatialJoin(OptimizerExtensionInput &input, unique_ptr<LogicalOperator> &plan) {
	auto &op = *plan;

	// We only care about ANY_JOIN operators
	if (op.type != LogicalOperatorType::LOGICAL_ANY_JOIN) {
		return;
	}

	auto &any_join = op.Cast<LogicalAnyJoin>();

	// We also only support INNER joins for now
	if (any_join.join_type != JoinType::INNER) {
		return;
	}

	// Inspect the join condition
	vector<unique_ptr<Expression>> expressions;
	expressions.push_back(any_join.condition->Copy()); // TODO: Maybe move instead of copy

	// Split by AND
	LogicalFilter::SplitPredicates(expressions);

	// Get the table indexes that are reachable from the left and right children
	auto &left_child = any_join.children[0];
	auto &right_child = any_join.children[1];
	unordered_set<idx_t> left_bindings;
	unordered_set<idx_t> right_bindings;
	LogicalJoin::GetTableReferences(*left_child, left_bindings);
	LogicalJoin::GetTableReferences(*right_child, right_bindings);

	// TODO: Only support a single predicate for now
	if (expressions.size() != 1) {
		return;
	}

	vector<SpatialJoinCondition> join_conditions;

	// Now, check each expression to see if it contains a spatial predicate
	for (auto &expr : expressions) {
		auto total_side = JoinSide::GetJoinSide(*expr, left_bindings, right_bindings);

		if (total_side != JoinSide::BOTH) {
			// Throw?. No, push down the filter
			return;
		}

		// Check if the expression is a spatial predicate
		if (expr->type != ExpressionType::BOUND_FUNCTION) {
			return;
		}

		auto &func = expr->Cast<BoundFunctionExpression>();
		if (func.function.name != "ST_Intersects") {
			return;
		}

		auto left_side = JoinSide::GetJoinSide(*func.children[0], left_bindings, right_bindings);
		auto right_side = JoinSide::GetJoinSide(*func.children[1], left_bindings, right_bindings);

		// The condition can be cleanly split into two sides
		if (left_side != JoinSide::BOTH && right_side != JoinSide::BOTH) {
			SpatialJoinCondition condition;

			// TODO: Support more predicates, and flip/invert them if neccessary
			condition.predicate = "ST_Intersects";

			auto left = std::move(func.children[0]);
			auto right = std::move(func.children[1]);
			if (left_side == JoinSide::RIGHT) {
				std::swap(left, right);
				// TODO: flip predicate here if neccessary
			}

			condition.left = std::move(left);
			condition.right = std::move(right);

			join_conditions.push_back(std::move(condition));
		}
	}

	// Nope!
	if (join_conditions.size() != 1) {
		return;
	}

	// Cool, now we have spatial join conditions. Proceed to create a new LogicalSpatialJoin operator
	auto spatial_join = make_uniq<LogicalSpatialJoin>(JoinType::INNER);

	// Steal the properties from the any join
	spatial_join->conditions = std::move(join_conditions);
	spatial_join->children = std::move(any_join.children);
	spatial_join->expressions = std::move(any_join.expressions);
	spatial_join->types = std::move(any_join.types);
	spatial_join->left_projection_map = std::move(any_join.left_projection_map);
	spatial_join->right_projection_map = std::move(any_join.right_projection_map);
	spatial_join->join_stats = std::move(any_join.join_stats);
	spatial_join->mark_index = any_join.mark_index;
	spatial_join->has_estimated_cardinality = any_join.has_estimated_cardinality;
	spatial_join->estimated_cardinality = any_join.estimated_cardinality;

	// Replace the operator
	plan = std::move(spatial_join);
}

static void TryInsertSpatialJoin(OptimizerExtensionInput &input, unique_ptr<LogicalOperator> &plan) {

	InsertSpatialJoin(input, plan);

	// Recursively call this function on all children
	for (auto &child : plan->children) {
		TryInsertSpatialJoin(input, child);
	}
}

void SpatialJoinRule::Register(DatabaseInstance &db) {

	OptimizerExtension optimizer;
	optimizer.optimize_function = TryInsertSpatialJoin;

	db.config.optimizer_extensions.push_back(optimizer);
}

} // namespace duckdb