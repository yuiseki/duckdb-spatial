// TODO: Dont depend on duckdb
#include "duckdb/common/assert.hpp"
#define SGL_ASSERT(condition) D_ASSERT(condition)

#include "sgl/sgl.hpp"

namespace sgl {

namespace ops {

static uint8_t *resize_vertices(allocator &alloc, geometry *geom, bool set_z, bool set_m, double default_z,
                                double default_m) {

	const auto has_z = geom->has_z();
	const auto has_m = geom->has_m();

	const auto source_type = static_cast<vertex_type>(has_z + 2 * has_m);
	const auto target_type = static_cast<vertex_type>(set_z + 2 * set_m);

	const auto source_data = geom->get_vertex_data();
	const auto count = geom->get_count();

	if (source_type == target_type) {
		return source_data;
	}

	switch (source_type) {
	case vertex_type::XY: {
		constexpr auto source_size = sizeof(double) * 2;
		switch (target_type) {
		case vertex_type::XY: {
			// Do nothing
			return source_data;
		}
		case vertex_type::XYZ: {
			constexpr auto target_size = sizeof(double) * 3;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, source_size);
				memcpy(target_data + target_offset + source_size, &default_z, sizeof(double));
			}

			return target_data;
		}
		case vertex_type::XYM: {
			constexpr auto target_size = sizeof(double) * 3;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, source_size);
				memcpy(target_data + target_offset + source_size, &default_m, sizeof(double));
			}

			return target_data;
		}
		case vertex_type::XYZM: {
			constexpr auto target_size = sizeof(double) * 4;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, source_size);
				memcpy(target_data + target_offset + source_size, &default_z, sizeof(double));
				memcpy(target_data + target_offset + source_size + sizeof(double), &default_m, sizeof(double));
			}

			return target_data;
		}
		default:
			SGL_ASSERT(false);
			return nullptr;
		}
	}
	case vertex_type::XYZ: {
		constexpr auto source_size = sizeof(double) * 3;
		switch (target_type) {
		case vertex_type::XY: {
			constexpr auto target_size = sizeof(double) * 2;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, target_size);
			}

			return target_data;
		}
		case vertex_type::XYZ: {
			// Do nothing
			return source_data;
		}
		case vertex_type::XYM: {
			constexpr auto target_size = sizeof(double) * 3;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, target_size);
				memcpy(target_data + target_offset + sizeof(double) * 2, &default_m, sizeof(double));
			}

			return target_data;
		}
		case vertex_type::XYZM: {
			constexpr auto target_size = sizeof(double) * 4;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, target_size);
				memcpy(target_data + target_offset + sizeof(double) * 3, &default_m, sizeof(double));
			}

			return target_data;
		}
		default:
			SGL_ASSERT(false);
			return nullptr;
		}
	}
	case vertex_type::XYM: {
		constexpr auto source_size = sizeof(double) * 3;
		switch (target_type) {
		case vertex_type::XY: {
			constexpr auto target_size = sizeof(double) * 2;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, target_size);
			}

			return target_data;
		}
		case vertex_type::XYZ: {
			constexpr auto target_size = sizeof(double) * 3;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, sizeof(double) * 2);
				memcpy(target_data + target_offset + sizeof(double) * 2, &default_z, sizeof(double));
			}

			return target_data;
		}
		case vertex_type::XYM: {
			// Do nothing
			return source_data;
		}
		case vertex_type::XYZM: {
			constexpr auto target_size = sizeof(double) * 4;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, sizeof(double) * 2);
				memcpy(target_data + target_offset + sizeof(double) * 2, &default_z, sizeof(double));
				memcpy(target_data + target_offset + sizeof(double) * 3,
				       source_data + source_offset + sizeof(double) * 2, sizeof(double));
			}

			return target_data;
		}
		default:
			SGL_ASSERT(false);
			return nullptr;
		}
	}
	case vertex_type::XYZM: {
		constexpr auto source_size = sizeof(double) * 4;
		switch (target_type) {
		case vertex_type::XY: {
			constexpr auto target_size = sizeof(double) * 2;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, sizeof(double) * 2);
			}

			return target_data;
		}
		case vertex_type::XYZ: {
			constexpr auto target_size = sizeof(double) * 3;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, sizeof(double) * 3);
			}

			return target_data;
		}
		case vertex_type::XYM: {
			constexpr auto target_size = sizeof(double) * 3;
			const auto target_data = static_cast<uint8_t *>(alloc.alloc(count * target_size));

			for (size_t i = 0; i < count; i++) {
				const auto source_offset = i * source_size;
				const auto target_offset = i * target_size;
				memcpy(target_data + target_offset, source_data + source_offset, sizeof(double) * 2);
				memcpy(target_data + target_offset + sizeof(double) * 2,
				       source_data + source_offset + sizeof(double) * 3, sizeof(double));
			}

			return target_data;
		}
		case vertex_type::XYZM: {
			// Do nothing
			return source_data;
		}
		default:
			SGL_ASSERT(false);
			return nullptr;
		}
	}
	default:
		SGL_ASSERT(false);
		return nullptr;
	}
}

void force_zm(allocator &alloc, geometry *geom, bool set_z, bool set_m, double default_z, double default_m) {

	auto part = geom;
	if (part == nullptr) {
		return;
	}
	const auto root = part->get_parent();

	while (part != root) {

		switch (part->get_type()) {
		case geometry_type::POINT:
		case geometry_type::LINESTRING: {
			// Convert the vertices
			const auto target_data = resize_vertices(alloc, part, set_z, set_m, default_z, default_m);
			part->set_vertex_data(target_data, part->get_count());
			part->set_z(set_z);
			part->set_m(set_m);
		} break;
		case geometry_type::POLYGON:
		case geometry_type::MULTI_POINT:
		case geometry_type::MULTI_LINESTRING:
		case geometry_type::MULTI_POLYGON:
		case geometry_type::MULTI_GEOMETRY: {
			part->set_z(set_z);
			part->set_m(set_m);
			if (!part->is_empty()) {
				part = part->get_first_part();
			}
		} break;
		default:
			SGL_ASSERT(false);
			break;
		}

		// Now go up/sideways
		while (part != nullptr) {
			const auto parent = part->get_parent();
			if (parent == root) {
				return;
			}

			if (part != parent->get_last_part()) {
				// Go sideways
				part = part->get_next();
				break;
			}

			// Go up
			part = parent;
		}
	}
}

} // namespace ops

} // namespace sgl