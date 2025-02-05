#pragma once

// Wrapper around SGL that injects the DuckDB assert macro
#include "duckdb/common/assert.hpp"
#define SGL_ASSERT(x) D_ASSERT(x)
#include "sgl/sgl.hpp"

#include "duckdb/storage/arena_allocator.hpp"

namespace duckdb {

// sgl::allocator that uses a DuckDB ArenaAllocator to allocate memory
class GeometryAllocator final : public sgl::allocator {
public:
  explicit GeometryAllocator(ArenaAllocator &arena_p) : arena(arena_p) {
  }

  void *alloc(size_t size) override {
    return arena.AllocateAligned(size);
  }
  void dealloc(void *ptr, size_t size) override {
    arena.ReallocateAligned(data_ptr_cast(ptr), size, 0);
  }
  void *realloc(void *ptr, size_t old_size, size_t new_size) override {
    return arena.ReallocateAligned(data_ptr_cast(ptr), old_size, new_size);
  }
private:
  ArenaAllocator &arena;
};

} // namespace duckdb