#pragma once

#include "spatial/geometry/geometry_type.hpp"
#include "duckdb/common/types/string_type.hpp"

namespace duckdb {

class ArenaAllocator;

struct WKBWriter {
	// Write a geometry to a WKB blob attached to a vector
	static string_t Write(const geometry_t &geometry, Vector &result);
	static string_t Write(const string_t &geometry, Vector &result);

	// Write a geometry to a WKB blob into a buffer
	static void Write(const geometry_t &geometry, vector<data_t> &buffer);
	static void Write(const string_t &geometry, vector<data_t> &buffer);

	// Write a geometry to a WKB blob into an arena allocator
	static const_data_ptr_t Write(const geometry_t &geometry, uint32_t *size, ArenaAllocator &allocator);
	static const_data_ptr_t Write(const string_t &geometry, uint32_t *size, ArenaAllocator &allocator);
};

} // namespace duckdb