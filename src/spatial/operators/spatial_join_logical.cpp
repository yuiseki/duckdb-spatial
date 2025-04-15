#include "spatial_join_logical.hpp"

#include "duckdb/catalog/catalog_entry/scalar_function_catalog_entry.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "spatial_join_physical.hpp"

#include <duckdb/execution/column_binding_resolver.hpp>

namespace duckdb {

unique_ptr<Expression> SpatialJoinCondition::ToExpr(ClientContext &context) const {
	auto &catalog = Catalog::GetSystemCatalog(context);

	auto &func_set = catalog.GetEntry(context, CatalogType::SCALAR_FUNCTION_ENTRY, DEFAULT_SCHEMA, predicate)
	                     .Cast<ScalarFunctionCatalogEntry>();

	auto func = func_set.functions.GetFunctionByArguments(context, {left->return_type, right->return_type});
	vector<unique_ptr<Expression>> args;
	args.push_back(left->Copy());
	args.push_back(right->Copy());
	return make_uniq_base<Expression, BoundFunctionExpression>(func.return_type, func, std::move(args), nullptr);
}

LogicalSpatialJoin::LogicalSpatialJoin(JoinType join_type_p) : join_type(join_type_p) {
}

vector<ColumnBinding> LogicalSpatialJoin::GetColumnBindings() {
	auto left_bindings = MapBindings(children[0]->GetColumnBindings(), left_projection_map);
	if (join_type == JoinType::SEMI || join_type == JoinType::ANTI) {
		// for SEMI and ANTI join we only project the left hand side
		return left_bindings;
	}

	if (join_type == JoinType::MARK) {
		// for MARK join we project the left hand side plus the MARK column
		left_bindings.emplace_back(mark_index, 0);
		return left_bindings;
	}
	// for other join types we project both the LHS and the RHS
	auto right_bindings = MapBindings(children[1]->GetColumnBindings(), right_projection_map);
	if (join_type == JoinType::RIGHT_SEMI || join_type == JoinType::RIGHT_ANTI) {
		return right_bindings;
	}
	left_bindings.insert(left_bindings.end(), right_bindings.begin(), right_bindings.end());
	return left_bindings;
}

void LogicalSpatialJoin::ResolveColumnBindings(ColumnBindingResolver &res, vector<ColumnBinding> &bindings) {

	auto &cond = spatial_predicate->Cast<BoundFunctionExpression>();

	res.VisitOperator(*children[0]);
	res.VisitExpression(&cond.children[0]);

	// TODO: Duplicate eliminated joins?

	res.VisitOperator(*children[1]);
	res.VisitExpression(&cond.children[1]);

	// Finally, update the bindings
	bindings = GetColumnBindings();
}

void LogicalSpatialJoin::ResolveTypes() {
	types = MapTypes(children[0]->types, left_projection_map);
	if (join_type == JoinType::SEMI || join_type == JoinType::ANTI) {
		// for SEMI and ANTI join we only project the left hand side
		return;
	}
	if (join_type == JoinType::MARK) {
		// for MARK join we project the left hand side, plus a BOOLEAN column indicating the MARK
		types.emplace_back(LogicalType::BOOLEAN);
		return;
	}
	// for any other join we project both sides
	auto right_types = MapTypes(children[1]->types, right_projection_map);
	if (join_type == JoinType::RIGHT_SEMI || join_type == JoinType::RIGHT_ANTI) {
		types = right_types;
		return;
	}
	types.insert(types.end(), right_types.begin(), right_types.end());
}

string LogicalSpatialJoin::GetExtensionName() const {
	return "duckdb_spatial";
}

string LogicalSpatialJoin::GetName() const {
	return "SPATIAL_JOIN";
}

unique_ptr<PhysicalOperator> LogicalSpatialJoin::CreatePlan(ClientContext &context, PhysicalPlanGenerator &generator) {

	// Return a new PhysicalSpatialJoin operator
	auto left = generator.CreatePlan(std::move(children[0]));
	auto right = generator.CreatePlan(std::move(children[1]));

	return make_uniq_base<PhysicalOperator, PhysicalSpatialJoin>(
	    *this, std::move(left), std::move(right), std::move(spatial_predicate), join_type, estimated_cardinality);
}

} // namespace duckdb