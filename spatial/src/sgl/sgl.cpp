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

// TODO: Make non-recursive
size_t to_wkb_size(const geometry *geom) {
	switch (geom->get_type()) {
	case geometry_type::POINT: {
		// order (1)
		// type (4)
		// point (1 * vsize)
		return 1 + 4 + geom->get_vertex_size();
	}
	case geometry_type::LINESTRING: {
		// order (1)
		// type (4)
		// count (4)
		// points (count * vsize)
		return 1 + 4 + 4 + geom->get_vertex_size() * geom->get_count();
	}
	case geometry_type::POLYGON: {
		// order (1)
		// type (4)
		// ring count (4)
		size_t size = 1 + 4 + 4;
		const auto tail = geom->get_last_part();
		auto head = tail;
		if(head) {
			do {
				head = head->get_next();
				// count (4)
				// points (count * vsize)
				size += 4 + head->get_count() * head->get_vertex_size();
			} while (head != tail);
		}
		return size;
	}
	case geometry_type::MULTI_POINT:
	case geometry_type::MULTI_LINESTRING:
	case geometry_type::MULTI_POLYGON:
	case geometry_type::MULTI_GEOMETRY: {
		// order (1)
		// type (4)
		// geometry count (4)
		size_t size = 1 + 4 + 4;
		const auto tail = geom->get_last_part();
		auto head = tail;
		if(head) {
			do {
				head = head->get_next();
				size += to_wkb_size(head);
			} while (head != tail);
		}
		return size;
	}
	default:
		D_ASSERT(false);
		return 0;
	}
}

// TODO: Make non-recursive
size_t to_wkb(const geometry *geom, uint8_t *buffer, size_t size) {

#define WKB_WRITE_U8(PTR, VAL) do { uint8_t v = VAL; memcpy(PTR, &v, sizeof(uint8_t)); PTR += sizeof(uint8_t); } while(0)
#define WKB_WRITE_U32(PTR, VAL) do { uint32_t v = VAL; memcpy(PTR, &v, sizeof(uint32_t)); PTR += sizeof(uint32_t); } while(0)
#define WKB_WRITE_DOUBLE(PTR, VAL) do { double v = VAL; memcpy(PTR, &v, sizeof(double)); PTR += sizeof(double); } while(0)
#define WKB_WRITE_DATA(PTR, SRC, SIZE) do { memcpy(PTR, SRC, SIZE); PTR += SIZE; } while(0)

	auto ptr = buffer;

	// Write header
	const auto type_id = static_cast<uint32_t>(geom->get_type()) + geom->has_z() * 1000 + geom->has_m() * 2000;
	WKB_WRITE_U8(ptr, 1);
	WKB_WRITE_U32(ptr, type_id);

	// Write the body
	switch (geom->get_type()) {
	case geometry_type::POINT: {
		if(geom->is_empty()) {
			// WKB does not support empty points, so we write NaNs instead
			WKB_WRITE_DOUBLE(ptr, std::numeric_limits<double>::quiet_NaN());
			WKB_WRITE_DOUBLE(ptr, std::numeric_limits<double>::quiet_NaN());
		} else {
			WKB_WRITE_DATA(ptr, geom->get_vertex_data(), geom->get_vertex_size());
		}
	} break;
	case geometry_type::LINESTRING: {
		WKB_WRITE_U32(ptr, geom->get_count());
		WKB_WRITE_DATA(ptr, geom->get_vertex_data(), geom->get_vertex_size() * geom->get_count());
	} break;
	case geometry_type::POLYGON: {
		WKB_WRITE_U32(ptr, geom->get_count());
		const auto tail = geom->get_last_part();
		auto head = tail;
		if(head) {
			do {
				head = head->get_next();
				WKB_WRITE_U32(ptr, head->get_count());
				WKB_WRITE_DATA(ptr, head->get_vertex_data(), head->get_vertex_size() * head->get_count());
			} while (head != tail);
		}
	} break;
	case geometry_type::MULTI_POINT:
	case geometry_type::MULTI_LINESTRING:
	case geometry_type::MULTI_POLYGON:
	case geometry_type::MULTI_GEOMETRY: {
		WKB_WRITE_U32(ptr, geom->get_count());
		const auto tail = geom->get_last_part();
		auto head = tail;
		if(head) {
			do {
				head = head->get_next();
				ptr += to_wkb(head, ptr, size - (ptr - buffer));
			} while (head != tail);
		}
	} break;
	default:
		D_ASSERT(false);
		break;
	}

	return ptr - buffer;

#undef WKB_WRITE_U8
#undef WKB_WRITE_U32
#undef WKB_WRITE_DOUBLE
#undef WKB_WRITE_DATA
}

struct wkb_reader {
	const uint8_t* beg = nullptr;
	const uint8_t* ptr = nullptr;
	const uint8_t* end = nullptr;
	bool le = false;
	bool error = false;
};

uint8_t wkb_reader_read_u8(wkb_reader* reader) {
	if(reader->ptr + sizeof(uint8_t) > reader->end) {
		reader->error = true;
		return 0;
	}

	const auto val = *reader->ptr;
	reader->ptr += sizeof(uint8_t);
	return val;
}

uint32_t wkb_reader_read_u32(wkb_reader* reader) {
	if(reader->ptr + sizeof(uint32_t) > reader->end) {
		reader->error = true;
		return 0;
	}

	uint32_t val;
	memcpy(&val, reader->ptr, sizeof(uint32_t));
	reader->ptr += sizeof(uint32_t);
	return val;
}

double wkb_reader_read_f64(wkb_reader* reader) {
	if(reader->ptr + sizeof(double) > reader->end) {
		reader->error = true;
		return 0;
	}

	double val;
	memcpy(&val, reader->ptr, sizeof(double));
	reader->ptr += sizeof(double);
	return val;
}

const uint8_t* wkb_reader_read_data(wkb_reader* reader, size_t size) {
	if(reader->ptr + size > reader->end) {
		reader->error = true;
		return nullptr;
	}

	const auto val = reader->ptr;
	reader->ptr += size;
	return val;
}

geometry from_wkb(allocator *alloc, const uint8_t *buffer, size_t size) {

	// Setup state
	wkb_reader state = {0};
	state.beg = buffer;
	state.ptr = buffer;
	state.end = buffer + size;

	uint32_t stack[256];
	uint32_t depth = 0;

	geometry root;
	geometry *geom = &root;

#define WKB_READ_U8 wkb_reader_read_u8(&state); if (state.error) { goto error; }
#define WKB_READ_U32 wkb_reader_read_u32(&state); if (state.error) { goto error; }
#define WKB_READ_F64 wkb_reader_read_f64(&state); if (state.error) { goto error; }
#define WKB_READ_DATA(SIZE) wkb_reader_read_data(&state, SIZE); if (state.error) { goto error; }

begin:
	const auto le = WKB_READ_U8;
	const auto type_id = WKB_READ_U32;

	const auto type = static_cast<sgl::geometry_type>((type_id & 0xffff) % 1000);
	const auto flags = (type_id & 0xffff) / 1000;
	const auto has_z = (flags == 1) || (flags == 3) || ((type_id & 0x80000000) != 0);
	const auto has_m = (flags == 2) || (flags == 3) || ((type_id & 0x40000000) != 0);

	geom->set_type(type);
	geom->set_z(has_z);
	geom->set_m(has_m);

	switch(geom) {
	case geometry_type::POINT: {
		// Read the point data;
		const auto data = WKB_READ_DATA(geom->get_vertex_size());
		geom->set_vertex_data(data, 1);

		goto next;
	}
	case geometry_type::LINESTRING: {
		// Read the point count
		const auto count = WKB_READ_U32;
		// Read the point data;
		const auto data = WKB_READ_DATA(geom->get_vertex_size() * count);
		geom->set_vertex_data(data, count);

		goto next;
	}
	case geometry_type::POLYGON: {
		// Read the ring count
		const auto count = WKB_READ_U32;

		// Read the point data;
		for(size_t i = 0; i < count; i++) {
			const auto ring_count = WKB_READ_U32;
			const auto data = WKB_READ_DATA(geom->get_vertex_size() * ring_count);

			// create a new ring
			const auto ring = static_cast<geometry*>(alloc->alloc(sizeof(geometry)));
			ring->set_type(geometry_type::LINESTRING);
			ring->set_z(has_z);
			ring->set_m(has_m);
			ring->set_vertex_data(data, ring_count);

			geom->append_part(ring);
		}

		goto next;
	}
	case geometry_type::MULTI_POINT:
	case geometry_type::MULTI_LINESTRING:
	case geometry_type::MULTI_POLYGON:
	case geometry_type::MULTI_GEOMETRY: {
		// Check stack depth

		if(depth == 256) {
			// TODO: Better error handling
			goto error;
		}

		// read the count
		const auto count = WKB_READ_U32;
		stack[depth++] = count;

		// make a new child
		auto new_geom = static_cast<geometry*>(alloc->alloc(sizeof(geometry)));
		geom->append_part(new_geom);
		geom = new_geom;

		goto begin;
	}
	default:
		goto error;
	}

next:
	if(depth == 0) {
		goto done;
	}

	const auto remaining = stack[depth-1];

	if(remaining == 0) {
		depth--;
		geom = geom->get_parent();
		goto next;
	}

	auto new_geom = static_cast<geometry*>(alloc->alloc(sizeof(geometry)));
	auto parent = geom->get_parent();
	parent->append_part(new_geom);
	geom = new_geom;

	stack[depth-1]--;
	goto begin;

done:
	return root;

error:
	return geometry(geometry_type::INVALID);

#undef WKB_READ_U8
#undef WKB_READ_U32
#undef WKB_READ_F64
#undef WKB_READ_DATA

}

} // namespace ops

} // namespace sgl