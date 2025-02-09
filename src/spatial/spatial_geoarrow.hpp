#pragma once

namespace duckdb {

class DatabaseInstance;

struct GeoArrow {
	static void Register(DatabaseInstance &db);
};

} // namespace duckdb