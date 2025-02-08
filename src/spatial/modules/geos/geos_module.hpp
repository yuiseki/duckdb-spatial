#pragma once

namespace duckdb {

class DatabaseInstance;

void RegisterGEOSModule(DatabaseInstance &db);

}; // namespace duckdb