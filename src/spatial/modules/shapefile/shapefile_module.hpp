#pragma once

namespace duckdb {

class DatabaseInstance;

void RegisterShapefileModule(DatabaseInstance &db);

} // namespace duckdb