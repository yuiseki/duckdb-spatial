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


	// The spatial join condition
	unique_ptr<Expression> spatial_pred_expr = nullptr;

	// Extra predicates that are not spatial predicates
	vector<unique_ptr<Expression>> extra_predicates;

	// Now, check each expression to see if it contains a spatial predicate
	for (auto &expr : expressions) {
		auto total_side = JoinSide::GetJoinSide(*expr, left_bindings, right_bindings);

		if (total_side != JoinSide::BOTH) {
			// Throw?. No, push down the filter
			extra_predicates.push_back(std::move(expr));
			continue;;
		}

		// Check if the expression is a spatial predicate
		if (expr->type != ExpressionType::BOUND_FUNCTION) {
			extra_predicates.push_back(std::move(expr));
			continue;
		}

		auto &func = expr->Cast<BoundFunctionExpression>();
		if (func.function.name != "ST_Intersects") {
			extra_predicates.push_back(std::move(expr));
			continue;
		}

		auto left_side = JoinSide::GetJoinSide(*func.children[0], left_bindings, right_bindings);
		auto right_side = JoinSide::GetJoinSide(*func.children[1], left_bindings, right_bindings);

		// Can the condition can be cleanly split into two sides?
		if(left_side == JoinSide::BOTH || right_side == JoinSide::BOTH) {
			extra_predicates.push_back(std::move(expr));
			continue;
		}

		// TODO: Support more predicates, and flip/invert them if neccessary
		if (left_side == JoinSide::RIGHT) {
			// TODO: Flip function here if needed (if not symmetric)
			std::swap(func.children[0], func.children[1]);
		}

		spatial_pred_expr = std::move(expr);
	}

	// Nope! No spatial predicate found
	if (!spatial_pred_expr) {
		return;
	}

	// TODO: Push a filter for the extra conditions?

	// Cool, now we have spatial join conditions. Proceed to create a new LogicalSpatialJoin operator
	auto spatial_join = make_uniq<LogicalSpatialJoin>(JoinType::INNER);

	// Steal the properties from the any join
	spatial_join->spatial_predicate = std::move(spatial_pred_expr);
	spatial_join->extra_conditions = std::move(extra_predicates);
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