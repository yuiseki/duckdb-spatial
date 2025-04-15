#pragma once

namespace duckdb {

class DatabaseInstance;

struct SpatialJoinOptimizer {
	static void Register(DatabaseInstance &db);
};

} // namespace duckdb
