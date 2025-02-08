#pragma once

namespace duckdb {

class DatabaseInstance;

void RegisterGDALModule(DatabaseInstance &db);

} // namespace duckdb