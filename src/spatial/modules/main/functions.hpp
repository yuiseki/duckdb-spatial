#pragma once

#include "duckdb/common/typedefs.hpp"

namespace duckdb {

class DatabaseInstance;

void RegisterSpatialScalarFunctions(DatabaseInstance &db);
void RegisterSpatialAggregateFunctions(DatabaseInstance &db);
void RegisterSpatialCastFunctions(DatabaseInstance &db);

// TODO: Move these
class Vector;
struct CoreVectorOperations {
public:
	static void Point2DToVarchar(Vector &source, Vector &result, idx_t count);
	static void LineString2DToVarchar(Vector &source, Vector &result, idx_t count);
	static void Polygon2DToVarchar(Vector &source, Vector &result, idx_t count);
	static void Box2DToVarchar(Vector &source, Vector &result, idx_t count);
	static void GeometryToVarchar(Vector &source, Vector &result, idx_t count);
};

} // namespace duckdb