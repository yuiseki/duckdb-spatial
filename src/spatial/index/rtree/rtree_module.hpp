#pragma once

namespace duckdb {

class DatabaseInstance;

struct RTreeModule {
	static void RegisterIndex(DatabaseInstance &db);
	static void RegisterIndexScan(DatabaseInstance &db);
	static void RegisterIndexPlanScan(DatabaseInstance &db);
	static void RegisterIndexPragmas(DatabaseInstance &db);
};

} // namespace duckdb