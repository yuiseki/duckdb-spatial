#pragma once

namespace duckdb {

class DatabaseInstance;

void RegisterSpatialOptimizers(DatabaseInstance &db);

} // namespace duckdb