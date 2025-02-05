#pragma once

#include <cstddef>

namespace sgl {
class geometry;
}

namespace duckdb {

class ArenaAllocator;

// todo:
struct Serde {
	static size_t GetRequiredSize(const sgl::geometry &geom);
	static void Serialize(const sgl::geometry &geom, char *buffer, size_t buffer_size);
	static void Deserialize(sgl::geometry &result, ArenaAllocator &arena, const char *buffer, size_t buffer_size);
};

}