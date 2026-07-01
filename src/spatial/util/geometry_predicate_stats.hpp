#pragma once

#include "duckdb/common/types/geometry.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/planner/expression.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"
#include "duckdb/storage/statistics/base_statistics.hpp"
#include "duckdb/storage/statistics/geometry_stats.hpp"

namespace duckdb {

//------------------------------------------------------------------------------
// Geometry predicate statistics pruning
//------------------------------------------------------------------------------
// Shared helpers that let the spatial "ST_" predicates participate in zonemap pruning via the
// scalar-function filter_prune callback. When a predicate is evaluated against a constant geometry, the
// callback reads the constant operand's bounding box from its derived statistics and prunes row groups whose
// column zonemap can never satisfy the predicate.

//! How the bounding boxes of the two operands must relate for a row to possibly match the predicate.
enum class GeometryPredicateBBox : uint8_t {
	//! A match implies bbox(arg0) and bbox(arg1) intersect (symmetric, e.g. ST_Intersects).
	INTERSECTS,
	//! A match implies bbox(arg0) contains bbox(arg1) (e.g. ST_Contains, ST_Covers).
	ARG0_COVERS_ARG1,
	//! A match implies bbox(arg1) contains bbox(arg0) (e.g. ST_Within, ST_CoveredBy).
	ARG1_COVERS_ARG0,
	//! A match implies bbox(arg0) equals bbox(arg1) (e.g. ST_Equals).
	EQUALS,
};

//! How the column zonemap must relate to the constant bounding box for a row to possibly match.
enum class GeometryPrunePredicate : uint8_t {
	//! No constant operand was found; pruning is not possible.
	NONE,
	//! The predicate can never be satisfied (e.g. the constant operand is an empty geometry).
	UNSATISFIABLE,
	//! The column zonemap must intersect the constant bounding box.
	COLUMN_INTERSECTS,
	//! The column zonemap must contain the constant bounding box.
	COLUMN_CONTAINS,
};

//! Identify the operands of a binary geometry predicate: exactly one of the first two arguments must be a
//! (folded) constant. On success sets `column_idx` to the other (column) operand's index. Returns false if the
//! operands don't match that shape (e.g. a spatial join, or both operands constant). The caller prunes against
//! the per-argument statistics, which are derived for us: the constant's geometry stats already carry its
//! bounding box, and the column's look through any CRS-only cast.
inline bool TryGetGeometryPredicateColumn(const BoundFunctionExpression &func, idx_t &column_idx) {
	auto &children = func.GetChildren();
	if (children.size() < 2) {
		return false;
	}
	// Constants are folded by the time we get here, so the constant operand is a plain BoundConstantExpression.
	const bool lhs_const = children[0]->GetExpressionType() == ExpressionType::VALUE_CONSTANT;
	const bool rhs_const = children[1]->GetExpressionType() == ExpressionType::VALUE_CONSTANT;
	if (lhs_const == rhs_const) {
		// Need exactly one constant operand and one column operand.
		return false;
	}
	column_idx = lhs_const ? 1 : 0;
	return true;
}

//! Resolve the runtime pruning predicate given the predicate's bbox semantics and the column operand's index.
inline GeometryPrunePredicate ResolveGeometryPrunePredicate(GeometryPredicateBBox bbox, idx_t column_idx) {
	switch (bbox) {
	case GeometryPredicateBBox::INTERSECTS:
		return GeometryPrunePredicate::COLUMN_INTERSECTS;
	case GeometryPredicateBBox::EQUALS:
		// A match implies identical bounding boxes, so the zonemap must contain the constant.
		return GeometryPrunePredicate::COLUMN_CONTAINS;
	case GeometryPredicateBBox::ARG0_COVERS_ARG1:
		// Containment pruning is only possible when the column is the covering operand (arg0).
		return column_idx == 0 ? GeometryPrunePredicate::COLUMN_CONTAINS : GeometryPrunePredicate::COLUMN_INTERSECTS;
	case GeometryPredicateBBox::ARG1_COVERS_ARG0:
		// Containment pruning is only possible when the column is the covering operand (arg1).
		return column_idx == 1 ? GeometryPrunePredicate::COLUMN_CONTAINS : GeometryPrunePredicate::COLUMN_INTERSECTS;
	default:
		return GeometryPrunePredicate::COLUMN_INTERSECTS;
	}
}

//! Execute the actual pruning given the predicate, the constant bounding box and the column zonemap.
inline FilterPropagateResult ExecuteGeometryPredicatePrune(GeometryPrunePredicate predicate,
                                                           const GeometryExtent &const_extent,
                                                           const BaseStatistics &stats) {
	if (predicate == GeometryPrunePredicate::NONE) {
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}
	if (stats.GetStatsType() != StatisticsType::GEOMETRY_STATS) {
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}
	if (!stats.CanHaveNoNull()) {
		// No non-null values are possible: the predicate is always false.
		return FilterPropagateResult::FILTER_ALWAYS_FALSE;
	}
	if (predicate == GeometryPrunePredicate::UNSATISFIABLE) {
		return FilterPropagateResult::FILTER_ALWAYS_FALSE;
	}

	const auto &col_extent = GeometryStats::GetExtent(stats);
	if (!col_extent.CanPruneXY()) {
		// If neither axis is set (the extent is empty or fully unknown), we cannot prune.
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}

	switch (predicate) {
	case GeometryPrunePredicate::COLUMN_INTERSECTS:
		// Every matching row's bbox intersects the constant, so the union (zonemap) must too.
		return col_extent.IntersectsXY(const_extent) ? FilterPropagateResult::NO_PRUNING_POSSIBLE
		                                             : FilterPropagateResult::FILTER_ALWAYS_FALSE;
	case GeometryPrunePredicate::COLUMN_CONTAINS:
		// Every matching row's bbox contains the constant, so the union (zonemap) must contain it too.
		return col_extent.ContainsXY(const_extent) ? FilterPropagateResult::NO_PRUNING_POSSIBLE
		                                            : FilterPropagateResult::FILTER_ALWAYS_FALSE;
	default:
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}
}

//! filter_prune callback for the symmetric/asymmetric binary geometry predicates.
template <GeometryPredicateBBox BBOX>
FilterPropagateResult GeometryPredicatePruneCallback(const FunctionStatisticsPruneInput &input) {
	idx_t column_idx;
	if (!TryGetGeometryPredicateColumn(input.function, column_idx)) {
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}
	auto column_stats = input.ChildStats(column_idx);
	auto constant_stats = input.ChildStats(1 - column_idx);
	if (!column_stats || !constant_stats || constant_stats->GetStatsType() != StatisticsType::GEOMETRY_STATS) {
		return FilterPropagateResult::NO_PRUNING_POSSIBLE;
	}
	// The constant's derived geometry stats already carry its bounding box; an empty constant has an empty
	// (non-prunable) extent and can never be satisfied.
	const auto &const_extent = GeometryStats::GetExtent(*constant_stats);
	const auto predicate =
	    const_extent.CanPruneXY() ? ResolveGeometryPrunePredicate(BBOX, column_idx) : GeometryPrunePredicate::UNSATISFIABLE;
	return ExecuteGeometryPredicatePrune(predicate, const_extent, *column_stats);
}

} // namespace duckdb
