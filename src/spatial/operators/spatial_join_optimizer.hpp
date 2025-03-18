#pragma once

namespace duckdb {

class DatabaseInstance;

struct SpatialJoinRule {
	static void Register(DatabaseInstance &db);
};

} // namespace duckdb
