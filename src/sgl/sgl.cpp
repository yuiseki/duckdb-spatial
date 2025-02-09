// TODO: Dont depend on duckdb
#include "duckdb/common/assert.hpp"
#define SGL_ASSERT(condition) D_ASSERT(condition)

#include "sgl/sgl.hpp"

#include <vector>

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

//------------------------------------------------------------------------------
// WKB Parsing
//------------------------------------------------------------------------------

static uint8_t wkb_reader_read_u8(wkb_reader *state) {
	if (state->pos + sizeof(uint8_t) > state->end) {
		state->error = SGL_WKB_READER_OUT_OF_BOUNDS;
		return 0;
	}

	const auto val = *state->pos;
	state->pos += sizeof(uint8_t);
	return val;
}

static uint32_t wkb_reader_read_u32(wkb_reader *state) {
	if (state->pos + sizeof(uint32_t) > state->end) {
		state->error = SGL_WKB_READER_OUT_OF_BOUNDS;
		return 0;
	}

	uint32_t val;

	if (state->le) {
		memcpy(&val, state->pos, sizeof(uint32_t));
	} else {
		char ibuf[sizeof(uint32_t)];
		char obuf[sizeof(uint32_t)];
		memcpy(ibuf, state->pos, sizeof(uint32_t));
		for (size_t i = 0; i < sizeof(uint32_t); i++) {
			obuf[i] = ibuf[sizeof(uint32_t) - i - 1];
		}
		memcpy(&val, obuf, sizeof(uint32_t));
	}

	state->pos += sizeof(uint32_t);

	return val;
}

static double wkb_reader_read_f64(wkb_reader *state) {
	if (state->pos + sizeof(double) > state->end) {
		state->error = SGL_WKB_READER_OUT_OF_BOUNDS;
		return 0;
	}

	double val;

	if (state->le) {
		memcpy(&val, state->pos, sizeof(double));
	} else {
		char ibuf[sizeof(double)];
		char obuf[sizeof(double)];
		memcpy(ibuf, state->pos, sizeof(double));
		for (size_t i = 0; i < sizeof(double); i++) {
			obuf[i] = ibuf[sizeof(double) - i - 1];
		}
		memcpy(&val, obuf, sizeof(double));
	}

	state->pos += sizeof(double);

	return val;
}

static bool wkb_reader_read_point(wkb_reader *state, geometry *geom) {
	const auto dims = 2 + geom->has_z() + geom->has_m();

	bool all_nan = true;
	double coords[4];

	const auto ptr = state->pos;
	for (size_t i = 0; i < dims; i++) {
		coords[i] = wkb_reader_read_f64(state);
		if (state->error) {
			return false;
		}
		if (!std::isnan(coords[i])) {
			all_nan = false;
		}
	}

	if (state->nan_as_empty && all_nan) {
		geom->set_vertex_data(static_cast<char *>(nullptr), 0);
		return true;
	}
	if (state->le && !state->copy_vertices) {
		geom->set_vertex_data(ptr, 1);
		return true;
	}

	const auto data = static_cast<char *>(state->alloc->alloc(sizeof(double) * dims));
	memcpy(data, coords, sizeof(double) * dims);
	geom->set_vertex_data(data, 1);
	return true;
}

static bool wkb_reader_read_line(wkb_reader *state, geometry *geom) {
	const auto vertex_count = wkb_reader_read_u32(state);
	if (state->error) {
		return false;
	}

	const auto vertex_size = geom->get_vertex_size();
	const auto byte_size = vertex_count * vertex_size;

	if (state->pos + byte_size > state->end) {
		state->error = SGL_WKB_READER_OUT_OF_BOUNDS;
		return false;
	}

	const auto ptr = state->pos;
	state->pos += byte_size;

	// If this is LE encoded, and we dont want to copy the vertices, we can just return the pointer
	if (state->le) {
		if (state->copy_vertices) {
			const auto mem = static_cast<char *>(state->alloc->alloc(byte_size));
			memcpy(mem, ptr, byte_size);
			geom->set_vertex_data(mem, vertex_count);
		} else {
			geom->set_vertex_data(ptr, vertex_count);
		}
	} else {
		// Otherwise, we need to allocate and swap the bytes
		const auto mem = static_cast<char *>(state->alloc->alloc(byte_size));
		for (size_t i = 0; i < vertex_count; i++) {
			const auto src = ptr + i * vertex_size;
			const auto dst = mem + i * vertex_size;

			// Swap doubles within the vertex
			for (size_t j = 0; j < vertex_size; j += sizeof(double)) {
				for (size_t k = 0; k < sizeof(double); k++) {
					dst[j + k] = src[j + sizeof(double) - k - 1];
				}
			}
		}

		geom->set_vertex_data(mem, vertex_count);
	}
	return true;
}

// TODO: Also collect stats?
bool wkb_reader_try_parse(wkb_reader *state, geometry *out) {

// clang-format off
#define read_u8(state) wkb_reader_read_u8(state); if (state->error) { return false; }
#define read_u32(state) wkb_reader_read_u32(state); if (state->error) { return false; }
#define read_u64(state) wkb_reader_read_u64(state); if (state->error) { return false; }
#define read_verts(state, vcount, vsize) wkb_reader_read_vertices(state, vcount, vsize); if (state->error) { return false; }
	// clang-format on

	SGL_ASSERT(state);
	SGL_ASSERT(out);
	SGL_ASSERT(state->buf);
	SGL_ASSERT(state->end);
	SGL_ASSERT(state->alloc);
	SGL_ASSERT(state->stack_buf);
	SGL_ASSERT(state->stack_cap > 0);

	// Setup state
	state->pos = state->buf;
	state->error = SGL_WKB_READER_OK;
	state->depth = 0;
	state->le = false;
	state->type_id = 0;
	state->has_any_m = false;
	state->has_any_z = false;

	geometry *geom = out;

	while (true) {
		state->le = read_u8(state);
		state->type_id = read_u32(state);

		const auto type = static_cast<sgl::geometry_type>((state->type_id & 0xffff) % 1000);
		const auto flags = (state->type_id & 0xffff) / 1000;
		const auto has_z = (flags == 1) || (flags == 3) || ((state->type_id & 0x80000000) != 0);
		const auto has_m = (flags == 2) || (flags == 3) || ((state->type_id & 0x40000000) != 0);
		const auto has_srid = (state->type_id & 0x20000000) != 0;

		if (has_srid) {
			// skip the SRID
			const auto srid = read_u32(state);
			(void)srid;
		}

		geom->set_type(type);
		geom->set_z(has_z);
		geom->set_m(has_m);

		// Compare with root
		if (!state->has_mixed_zm && (out->has_m() != has_m || out->has_z() != has_z)) {
			state->has_any_z |= has_z;
			state->has_any_m |= has_m;
			state->has_mixed_zm = true;
			if (!state->allow_mixed_zm) {
				// Error out!
				state->error = SGL_WKB_READER_MIXED_ZM;
				return false;
			}
		}

		switch (geom->get_type()) {
		case geometry_type::POINT: {
			// Read the point data
			if (!wkb_reader_read_point(state, geom)) {
				return false;
			}
		} break;
		case geometry_type::LINESTRING: {
			if (!wkb_reader_read_line(state, geom)) {
				return false;
			}
		} break;
		case geometry_type::POLYGON: {
			// Read the ring count
			const auto ring_count = read_u32(state);

			// Read the point data;
			for (size_t i = 0; i < ring_count; i++) {
				const auto ring = static_cast<geometry *>(state->alloc->alloc(sizeof(geometry)));
				new (ring) geometry(geometry_type::LINESTRING, has_z, has_m);
				if (!wkb_reader_read_line(state, ring)) {
					return false;
				}
				geom->append_part(ring);
			}
		} break;
		case geometry_type::MULTI_POINT:
		case geometry_type::MULTI_LINESTRING:
		case geometry_type::MULTI_POLYGON:
		case geometry_type::MULTI_GEOMETRY: {

			// Check stack depth
			if (state->depth >= state->stack_cap) {
				state->error = SGL_WKB_READER_RECURSION_LIMIT;
				return false;
			}

			// Read the count
			const auto count = read_u32(state);
			if (count == 0) {
				break;
			}

			state->stack_buf[state->depth++] = count;

			// Make a new child
			auto part = static_cast<geometry *>(state->alloc->alloc(sizeof(geometry)));
			new (part) geometry(geometry_type::INVALID, has_z, has_m);
			geom->append_part(part);

			// Set the new child as the current geometry
			geom = part;

			// Continue to the next iteration in the outer loop
			continue;
		}
		default:
			state->error = SGL_WKB_READER_UNSUPPORTED_TYPE;
			return false;
		}

		// Inner loop
		while (true) {
			const auto parent = geom->get_parent();

			if (state->depth == 0) {
				SGL_ASSERT(parent == nullptr);
				// Done!
				return true;
			}

			SGL_ASSERT(parent != nullptr);

			// Check that we are of the right type
			const auto ptype = parent->get_type();
			const auto ctype = geom->get_type();

			if (ptype == geometry_type::MULTI_POINT && ctype != geometry_type::POINT) {
				state->error = SGL_WKB_INVALID_CHILD_TYPE;
				return false;
			}
			if (ptype == geometry_type::MULTI_LINESTRING && ctype != geometry_type::LINESTRING) {
				state->error = SGL_WKB_INVALID_CHILD_TYPE;
				return false;
			}
			if (ptype == geometry_type::MULTI_POLYGON && ctype != geometry_type::POLYGON) {
				state->error = SGL_WKB_INVALID_CHILD_TYPE;
				return false;
			}

			// Check if we are done with the current part
			state->stack_buf[state->depth - 1]--;

			if (state->stack_buf[state->depth - 1] > 0) {
				// There are still more parts to read
				// Create a new part and append it to the parent
				auto part = static_cast<geometry *>(state->alloc->alloc(sizeof(geometry)));
				new (part) geometry(geometry_type::INVALID, has_z, has_m);
				parent->append_part(part);

				// Go "sideways" to the new part
				geom = part;
				break;
			}

			// Go upwards
			geom = parent;
			state->depth--;
		}
	}

#undef read_u8
#undef read_u32
#undef read_u64
#undef read_verts
}

std::string wkb_reader_get_error_message(const wkb_reader *state) {
	if (!state || state->error == SGL_WKB_READER_OK) {
		return "";
	}

	switch (state->error) {
	case SGL_WKB_READER_OUT_OF_BOUNDS: {
		return "Out of bounds read (is the WKB corrupt?)";
	}
	case SGL_WKB_READER_MIXED_ZM: {
		return "Mixed Z and M values are not allowed";
	}
	case SGL_WKB_READER_RECURSION_LIMIT: {
		return "Recursion limit '" + std::to_string(state->stack_cap) + "' reached";
	}
	case SGL_WKB_READER_UNSUPPORTED_TYPE: {
		// Try to fish out the type anyway
		return "Unsupported geometry type";
	}
	case SGL_WKB_INVALID_CHILD_TYPE: {
		return "Invalid child type";
	}
	default: {
		return "Unknown error";
	}
	}
}

size_t to_wkb_size(const geometry *geom) {
	if (!geom) {
		return 0;
	}

	const auto root = geom->get_parent();

	size_t size = 0;
	const geometry *curr = geom;

	// Main loop
	while (true) {
		switch (curr->get_type()) {
		case geometry_type::POINT: {
			size += 1 + 4 + curr->get_vertex_size();
		} break;
		case geometry_type::LINESTRING: {
			size += 1 + 4 + 4 + curr->get_count() * curr->get_vertex_size();
		} break;
		case geometry_type::POLYGON: {
			size += 1 + 4 + 4;
			const auto tail = curr->get_last_part();
			auto head = tail;
			if (head) {
				do {
					head = head->get_next();
					size += 4 + head->get_count() * head->get_vertex_size();
				} while (head != tail);
			}
		} break;
		case geometry_type::MULTI_POINT: {
			size += 1 + 4 + 4 + curr->get_count() * (1 + 4 + curr->get_vertex_size());
		} break;
		case geometry_type::MULTI_LINESTRING: {
			size += 1 + 4 + 4;
			const auto tail = curr->get_last_part();
			auto head = tail;
			if (head) {
				do {
					head = head->get_next();
					size += 1 + 4 + 4 + head->get_count() * head->get_vertex_size();
				} while (head != tail);
			}
		} break;
		case geometry_type::MULTI_POLYGON: {
			size += 1 + 4 + 4;
			const auto tail = curr->get_last_part();
			auto head = tail;
			if (head) {
				do {
					head = head->get_next();
					size += 1 + 4 + 4;
					const auto rtail = head->get_last_part();
					auto rhead = rtail;
					if (rhead) {
						do {
							rhead = rhead->get_next();
							size += 4 + rhead->get_count() * rhead->get_vertex_size();
						} while (rhead != rtail);
					}
				} while (head != tail);
			}
		} break;
		case geometry_type::MULTI_GEOMETRY: {
			size += 1 + 4 + 4;
			if (!curr->is_empty()) {
				curr = curr->get_first_part();
				continue;
			}
		} break;
		default: {
			SGL_ASSERT(false);
			return 0;
		}
		}

		// Inner loop
		while (true) {
			const auto parent = curr->get_parent();
			if (parent == root) {
				// Done!
				return size;
			}

			if (curr != parent->get_last_part()) {
				// Go sideways
				curr = curr->get_next();
				break;
			}

			// Go upwards
			curr = parent;
		}
	}
}

// TODO: Make non-recursive
size_t to_wkb(const geometry *geom, uint8_t *buffer, size_t size) {

#define WKB_WRITE_U8(PTR, VAL)                                                                                         \
	do {                                                                                                               \
		uint8_t v = VAL;                                                                                               \
		memcpy(PTR, &v, sizeof(uint8_t));                                                                              \
		PTR += sizeof(uint8_t);                                                                                        \
	} while (0)
#define WKB_WRITE_U32(PTR, VAL)                                                                                        \
	do {                                                                                                               \
		uint32_t v = VAL;                                                                                              \
		memcpy(PTR, &v, sizeof(uint32_t));                                                                             \
		PTR += sizeof(uint32_t);                                                                                       \
	} while (0)
#define WKB_WRITE_DOUBLE(PTR, VAL)                                                                                     \
	do {                                                                                                               \
		double v = VAL;                                                                                                \
		memcpy(PTR, &v, sizeof(double));                                                                               \
		PTR += sizeof(double);                                                                                         \
	} while (0)
#define WKB_WRITE_DATA(PTR, SRC, SIZE)                                                                                 \
	do {                                                                                                               \
		memcpy(PTR, SRC, SIZE);                                                                                        \
		PTR += SIZE;                                                                                                   \
	} while (0)

	auto ptr = buffer;

	// Write header
	const auto type_id = static_cast<uint32_t>(geom->get_type()) + geom->has_z() * 1000 + geom->has_m() * 2000;
	WKB_WRITE_U8(ptr, 1);
	WKB_WRITE_U32(ptr, type_id);

	// Write the body
	switch (geom->get_type()) {
	case geometry_type::POINT: {
		if (geom->is_empty()) {
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
		if (head) {
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
		if (head) {
			do {
				head = head->get_next();
				ptr += to_wkb(head, ptr, size - (ptr - buffer));
			} while (head != tail);
		}
	} break;
	default:
		SGL_ASSERT(false);
		break;
	}

	return ptr - buffer;

#undef WKB_WRITE_U8
#undef WKB_WRITE_U32
#undef WKB_WRITE_DOUBLE
#undef WKB_WRITE_DATA
}

//------------------------------------------------------------------------------
// WKT Parsing
//------------------------------------------------------------------------------

static void parse_ws(wkt_reader *state) {
	while (state->pos < state->end && std::isspace(*state->pos)) {
		state->pos++;
	}
}

static bool match_token(wkt_reader *state, const char *token) {
	// case insensitive match
	auto ptr = state->pos;
	while (ptr < state->end && *token != '\0' && std::tolower(*token) == std::tolower(*ptr)) {
		token++;
		ptr++;
	}

	if (*token != '\0') {
		return false;
	}

	state->pos = ptr;
	parse_ws(state);
	return true;
}

static bool match_char(wkt_reader *state, char c) {
	if (state->pos < state->end && std::tolower(*state->pos) == std::tolower(c)) {
		state->pos++;
		parse_ws(state);
		return true;
	}
	return false;
}

static bool match_double(wkt_reader *state, double *result) {
	// Because we care about the length, we cant just use std::strtod straight away without risking
	// out-of-bounds reads. Instead, we will manually parse the number and then use std::strtod to
	// convert the value.

	auto ptr = state->pos;

	// Match sign
	if (ptr < state->end && (*ptr == '+' || *ptr == '-')) {
		ptr++;
	}

	// Match number part
	while (ptr < state->end && std::isdigit(*ptr)) {
		ptr++;
	}

	// Match decimal part
	if (ptr < state->end && *ptr == '.') {
		ptr++;
		while (ptr < state->end && std::isdigit(*ptr)) {
			ptr++;
		}
	}

	// Match exponent part
	if (ptr < state->end && (*ptr == 'e' || *ptr == 'E')) {
		ptr++;
		if (ptr < state->end && (*ptr == '+' || *ptr == '-')) {
			ptr++;
		}

		while (ptr < state->end && std::isdigit(*ptr)) {
			ptr++;
		}
	}

	// Did we manage to parse anything?
	if (ptr == state->pos) {
		return false;
	}

	// If we got here, we know there is something resembling a  number within the bounds of the buffer
	// We can now use std::strtod to actually parse the number
	char *end;
	*result = std::strtod(state->pos, &end);
	if (state->pos == end) {
		return false;
	}
	state->pos = end;
	parse_ws(state);
	return true;
}

struct vertex_buffer {
	allocator *alloc;
	const uint32_t stride;
	double *ptr;
	uint32_t len;
	uint32_t cap;

	vertex_buffer(allocator *alloc, uint32_t stride) : alloc(alloc), stride(stride), len(0), cap(1) {
		ptr = static_cast<double *>(this->alloc->alloc(sizeof(double) * stride * cap));
	}

	void push_back(const double *data) {
		if (len == cap) {
			const auto new_cap = cap * 2;
			const auto old_size = sizeof(double) * stride * cap;
			const auto new_size = sizeof(double) * stride * new_cap;

			ptr = static_cast<double *>(alloc->realloc(ptr, old_size, new_size));
			cap = new_cap;
		}

		memcpy(ptr + len * stride, data, sizeof(double) * stride);
		len++;
	}

	void assign(geometry *geom) {

		// Shrink to fit
		if (cap > len) {
			const auto old_size = sizeof(double) * stride * cap;
			const auto new_size = sizeof(double) * stride * len;
			ptr = static_cast<double *>(alloc->realloc(ptr, old_size, new_size));
		}

		geom->set_vertex_data(reinterpret_cast<const char *>(ptr), len);
	}
};

// TODO: break this up into smaller functions, unify result/state
bool wkt_reader_try_parse(wkt_reader *state, geometry *out) {

	SGL_ASSERT(state != nullptr);
	SGL_ASSERT(out != nullptr);

	// These need to be set by the caller
	SGL_ASSERT(state->alloc != nullptr);
	SGL_ASSERT(state->buf != nullptr);
	SGL_ASSERT(state->end != nullptr);

	// Setup state
	state->pos = state->buf;
	state->error = nullptr;

	allocator *alloc = state->alloc;

	geometry *root = out;
	geometry *geom = root;

	// clang-format off
#define expect_char(STATE, C) do { if(!match_char(STATE, C)) { (STATE)->error = "Expected character: '" #C "'"; return false; } } while(0)
#define expect_number(STATE, RESULT) do { if(!match_double(STATE, RESULT)) { (STATE)->error = "Expected number"; return false; } } while(0)
	// clang-format on

	// Skip whitespace
	parse_ws(state);

	// Skip leading SRID, we dont support it
	// TODO: Parse this and stuff it into the result
	if (match_token(state, "SRID")) {

		while (state->pos < state->end && *state->pos != ';') {
			state->pos++;
		}
		expect_char(state, ';');
	}

	// Main loop
	while (true) {
		// Now we should have a geometry type
		if (match_token(state, "POINT")) {
			geom->set_type(geometry_type::POINT);
		} else if (match_token(state, "LINESTRING")) {
			geom->set_type(geometry_type::LINESTRING);
		} else if (match_token(state, "POLYGON")) {
			geom->set_type(geometry_type::POLYGON);
		} else if (match_token(state, "MULTIPOINT")) {
			geom->set_type(geometry_type::MULTI_POINT);
		} else if (match_token(state, "MULTILINESTRING")) {
			geom->set_type(geometry_type::MULTI_LINESTRING);
		} else if (match_token(state, "MULTIPOLYGON")) {
			geom->set_type(geometry_type::MULTI_POLYGON);
		} else if (match_token(state, "GEOMETRYCOLLECTION")) {
			geom->set_type(geometry_type::MULTI_GEOMETRY);
		} else {
			state->error = "Expected geometry type";
			return false;
		}

		// Match Z and M
		if (match_char(state, 'z')) {
			geom->set_z(true);
		}
		if (match_char(state, 'm')) {
			geom->set_m(true);
		}

		// TODO: make this check configurable
		if ((geom->has_m() != root->has_m()) || (geom->has_z() != root->has_z())) {
			state->error = "Mixed Z and M values are not supported";
			return false;
		}

		const auto vertex_stride = 2 + geom->has_z() + geom->has_m();

		// Parse EMPTY
		if (!match_token(state, "EMPTY")) {
			switch (geom->get_type()) {
			case geometry_type::POINT: {
				expect_char(state, '(');

				vertex_buffer verts(alloc, vertex_stride);
				double vert[4] = {0, 0, 0, 0};
				for (size_t i = 0; i < vertex_stride; i++) {
					expect_number(state, &vert[i]);
				}
				verts.push_back(vert);
				verts.assign(geom);

				expect_char(state, ')');
			} break;
			case geometry_type::LINESTRING: {
				expect_char(state, '(');

				vertex_buffer verts(alloc, vertex_stride);
				do {
					double vert[4] = {0, 0, 0, 0};
					for (size_t i = 0; i < vertex_stride; i++) {
						expect_number(state, &vert[i]);
					}
					verts.push_back(vert);
				} while (match_char(state, ','));

				verts.assign(geom);

				expect_char(state, ')');
			} break;
			case geometry_type::POLYGON: {
				expect_char(state, '(');
				do {
					auto ring = static_cast<geometry *>(alloc->alloc(sizeof(geometry)));
					new (ring) geometry(geometry_type::LINESTRING, geom->has_z(), geom->has_m());
					if (!match_token(state, "EMPTY")) {
						expect_char(state, '(');

						vertex_buffer verts(alloc, vertex_stride);
						do {
							double vert[4] = {0, 0, 0, 0};
							for (size_t i = 0; i < vertex_stride; i++) {
								expect_number(state, &vert[i]);
							}
							verts.push_back(vert);
						} while (match_char(state, ','));

						verts.assign(ring);

						expect_char(state, ')');
					}
					geom->append_part(ring);
				} while (match_char(state, ','));
				expect_char(state, ')');
			} break;
			case geometry_type::MULTI_POINT: {
				expect_char(state, '(');
				// Multipoints are special in that parens around each point is optional.
				do {
					bool has_paren = false;
					if (match_char(state, '(')) {
						has_paren = true;
					}
					auto point = static_cast<geometry *>(alloc->alloc(sizeof(geometry)));
					new (point) geometry(geometry_type::POINT, geom->has_z(), geom->has_m());
					if (!match_token(state, "EMPTY")) {
						// TODO: Do we need to have optional parens to accept EMPTY?

						vertex_buffer verts(alloc, vertex_stride);
						double vert[4] = {0, 0, 0, 0};
						for (size_t i = 0; i < vertex_stride; i++) {
							expect_number(state, &vert[i]);
						}
						verts.push_back(vert);
						verts.assign(point);
					}
					if (has_paren) {
						expect_char(state, ')');
					}
					geom->append_part(point);
				} while (match_char(state, ','));
				expect_char(state, ')');
			} break;
			case geometry_type::MULTI_LINESTRING: {
				expect_char(state, '(');
				do {
					auto line = static_cast<geometry *>(alloc->alloc(sizeof(geometry)));
					new (line) geometry(geometry_type::LINESTRING, geom->has_z(), geom->has_m());
					if (!match_token(state, "EMPTY")) {
						expect_char(state, '(');

						vertex_buffer verts(alloc, vertex_stride);
						do {
							double vert[4] = {0, 0, 0, 0};
							for (size_t i = 0; i < vertex_stride; i++) {
								expect_number(state, &vert[i]);
							}
							verts.push_back(vert);
						} while (match_char(state, ','));

						verts.assign(line);

						expect_char(state, ')');
					}
					geom->append_part(line);
				} while (match_char(state, ','));
				expect_char(state, ')');
			} break;
			case geometry_type::MULTI_POLYGON: {
				expect_char(state, '(');
				do {
					auto poly = static_cast<geometry *>(alloc->alloc(sizeof(geometry)));
					new (poly) geometry(geometry_type::POLYGON, geom->has_z(), geom->has_m());
					if (!match_token(state, "EMPTY")) {
						expect_char(state, '(');
						do {
							auto ring = static_cast<geometry *>(alloc->alloc(sizeof(geometry)));
							new (ring) geometry(geometry_type::LINESTRING, geom->has_z(), geom->has_m());
							if (!match_token(state, "EMPTY")) {
								expect_char(state, '(');

								vertex_buffer verts(alloc, vertex_stride);
								do {
									double vert[4] = {0, 0, 0, 0};
									for (size_t i = 0; i < vertex_stride; i++) {
										expect_number(state, &vert[i]);
									}
									verts.push_back(vert);
								} while (match_char(state, ','));

								verts.assign(ring);

								expect_char(state, ')');
							}
							poly->append_part(ring);
						} while (match_char(state, ','));
						expect_char(state, ')');
					}
					geom->append_part(poly);
				} while (match_char(state, ','));
				expect_char(state, ')');
			} break;
			case geometry_type::MULTI_GEOMETRY: {
				expect_char(state, '(');

				// add another child
				auto new_geom = static_cast<geometry *>(alloc->alloc(sizeof(geometry)));
				new (new_geom) geometry(geometry_type::INVALID);

				geom->append_part(new_geom);
				geom = new_geom;
			}
				continue; // This continue moves us to the next iteration
			default:
				SGL_ASSERT(false);
				state->error = "Invalid geometry type";
				return false;
			}
		}

		while (true) {
			const auto parent = geom->get_parent();
			if (!parent) {
				// Done!
				return true;
			}

			SGL_ASSERT(parent->get_type() == geometry_type::MULTI_GEOMETRY);

			if (match_char(state, ',')) {
				// The geometry collection is not done yet, add another sibling
				auto new_geom = static_cast<geometry *>(alloc->alloc(sizeof(geometry)));
				new (new_geom) geometry(geometry_type::INVALID);

				parent->append_part(new_geom);
				geom = new_geom;

				// goto begin;
				break;
			}

			expect_char(state, ')');
			// The geometry collection is done, go up
			geom = parent;
		}
	}

#undef expect_char
#undef expect_number
}

std::string wkt_reader_get_error_message(const wkt_reader *state) {
	if (!state || !state->error) {
		return "";
	}

	// Return a string of the current position in the input string
	const auto len = 32;
	const auto range_beg = std::max(state->pos - len, state->buf);
	const auto range_End = std::min(state->pos + 1, state->end);
	auto range = std::string(range_beg, range_End);
	if (range_beg != state->buf) {
		range = "..." + range;
	}

	// Add an arrow to indicate the position
	const auto err = std::string(state->error);
	const auto pos = std::to_string(state->pos - state->buf);
	const auto msg = err + " at position '" + pos + "' near: '" + range + "'|<---";

	return msg;
}

//------------------------------------------------------------------------------
// Extract
//------------------------------------------------------------------------------
// TODO: Make these non-recursive

static bool select_points(void *, const sgl::geometry *geom) {
	switch (geom->get_type()) {
	case sgl::geometry_type::POINT:
	case sgl::geometry_type::MULTI_POINT:
	case sgl::geometry_type::MULTI_GEOMETRY:
		return true;
	default:
		return false;
	}
}

static void handle_points(void *state, sgl::geometry *geom) {
	auto &points = *static_cast<sgl::geometry *>(state);

	switch (geom->get_type()) {
	case sgl::geometry_type::POINT:
		points.append_part(geom);
		break;
	case sgl::geometry_type::MULTI_POINT:
	case sgl::geometry_type::MULTI_GEOMETRY:
		geom->filter_parts(state, select_points, handle_points);
		break;
	default:
		SGL_ASSERT(false);
		break;
	}
}

static bool select_lines(void *state, const sgl::geometry *geom) {
	switch (geom->get_type()) {
	case sgl::geometry_type::LINESTRING:
	case sgl::geometry_type::MULTI_LINESTRING:
	case sgl::geometry_type::MULTI_GEOMETRY:
		return true;
	default:
		return false;
	}
}

static void handle_lines(void *state, sgl::geometry *geom) {
	auto &lines = *static_cast<sgl::geometry *>(state);

	switch (geom->get_type()) {
	case sgl::geometry_type::LINESTRING:
		lines.append_part(geom);
		break;
	case sgl::geometry_type::MULTI_LINESTRING:
	case sgl::geometry_type::MULTI_GEOMETRY:
		geom->filter_parts(state, select_lines, handle_lines);
		break;
	default:
		SGL_ASSERT(false);
		break;
	}
}

static bool select_polygons(void *state, const sgl::geometry *geom) {
	switch (geom->get_type()) {
	case sgl::geometry_type::POLYGON:
	case sgl::geometry_type::MULTI_POLYGON:
	case sgl::geometry_type::MULTI_GEOMETRY:
		return true;
	default:
		return false;
	}
}

static void handle_polygons(void *state, sgl::geometry *geom) {
	auto &polygons = *static_cast<sgl::geometry *>(state);

	switch (geom->get_type()) {
	case sgl::geometry_type::POLYGON:
		polygons.append_part(geom);
		break;
	case sgl::geometry_type::MULTI_POLYGON:
	case sgl::geometry_type::MULTI_GEOMETRY:
		geom->filter_parts(state, select_polygons, handle_polygons);
		break;
	default:
		SGL_ASSERT(false);
		break;
	}
}

geometry extract_points(sgl::geometry *geom) {
	auto points = sgl::geometry(sgl::geometry_type::MULTI_POINT, geom->has_z(), geom->has_m());
	geom->filter_parts(&points, select_points, handle_points);
	return points;
}

geometry extract_linestrings(sgl::geometry *geom) {
	auto lines = sgl::geometry(sgl::geometry_type::MULTI_LINESTRING, geom->has_z(), geom->has_m());
	geom->filter_parts(&lines, select_lines, handle_lines);
	return lines;
}

geometry extract_polygons(sgl::geometry *geom) {
	auto polygons = sgl::geometry(sgl::geometry_type::MULTI_POLYGON, geom->has_z(), geom->has_m());
	geom->filter_parts(&polygons, select_polygons, handle_polygons);
	return polygons;
}


//------------------------------------------------------------------------------
// Distance
//------------------------------------------------------------------------------
static double point_point_distance(const sgl::geometry *lhs, const sgl::geometry *rhs) {
	SGL_ASSERT(lhs->get_type() == sgl::geometry_type::POINT);
	SGL_ASSERT(rhs->get_type() == sgl::geometry_type::POINT);

	if(lhs->is_empty() || rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	const auto lhs_vertex = lhs->get_vertex_xy(0);
	const auto rhs_vertex = rhs->get_vertex_xy(0);

	return std::hypot(lhs_vertex.x - rhs_vertex.x, lhs_vertex.y - rhs_vertex.y);
}

/*
function sqr(x) { return x * x }
function dist2(v, w) { return sqr(v.x - w.x) + sqr(v.y - w.y) }
function distToSegmentSquared(p, v, w) {
var l2 = dist2(v, w);
if (l2 == 0) return dist2(p, v);
var t = ((p.x - v.x) * (w.x - v.x) + (p.y - v.y) * (w.y - v.y)) / l2;
t = Math.max(0, Math.min(1, t));
return dist2(p, { x: v.x + t * (w.x - v.x),
y: v.y + t * (w.y - v.y) });
}
function distToSegment(p, v, w) { return Math.sqrt(distToSegmentSquared(p, v, w)); }
 */

static double vertex_distance_squared(const vertex_xy *lhs, const vertex_xy *rhs) {
	return std::pow(lhs->x - rhs->x, 2) + std::pow(lhs->y - rhs->y, 2);
}

static double vertex_distance(const vertex_xy *lhs, const vertex_xy *rhs) {
	return std::hypot(lhs->x - rhs->x, lhs->y - rhs->y);
}

static double point_line_distance(const vertex_xy *p, const vertex_xy *v, const vertex_xy *w) {
	const auto l2 = vertex_distance_squared(v, w);
	if (l2 == 0) {
		// is not better to just compare if w == v?
		return vertex_distance(p, v);
	}

	const auto t = ((p->x - v->x) * (w->x - v->x) + (p->y - v->y) * (w->y - v->y)) / l2;
	const auto t_clamped = std::max(0.0, std::min(1.0, t));
	const auto x = v->x + t_clamped * (w->x - v->x);
	const auto y = v->y + t_clamped * (w->y - v->y);

	const vertex_xy intersection {x, y};

	return vertex_distance(p, &intersection);
}

static double point_linestring_distance(const sgl::geometry *lhs, const sgl::geometry *rhs) {
	SGL_ASSERT(lhs->get_type() == sgl::geometry_type::POINT);
	SGL_ASSERT(rhs->get_type() == sgl::geometry_type::LINESTRING);

	if(lhs->is_empty() || rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	const auto lhs_vertex = lhs->get_vertex_xy(0);
	double min_dist = std::numeric_limits<double>::infinity();
	const auto count = rhs->get_count();

	auto v1 = rhs->get_vertex_xy(0);
	if(count == 1) {
		// Degenerate case, should not happen
		return vertex_distance(&lhs_vertex, &v1);
	}

	for(size_t i = 1; i < count; i++) {
		const auto v2 = rhs->get_vertex_xy(i);
		const auto dist = point_line_distance(&lhs_vertex, &v1, &v2);
		min_dist = std::min(min_dist, dist);
		v1 = v2;
	}

	return min_dist;
}

static double point_polygon_distance(const sgl::geometry *lhs, const sgl::geometry *rhs) {
	SGL_ASSERT(lhs->get_type() == sgl::geometry_type::POINT);
	SGL_ASSERT(rhs->get_type() == sgl::geometry_type::POLYGON);

	if(rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	const auto shell = rhs->get_first_part();
	SGL_ASSERT(shell != nullptr);
	return point_linestring_distance(lhs, shell);
}

static double linestring_linestring_distance(const geometry *lhs, const geometry *rhs) {
	SGL_ASSERT(lhs->get_type() == geometry_type::LINESTRING);
	SGL_ASSERT(rhs->get_type() == geometry_type::LINESTRING);

	if(lhs->is_empty() || rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	double min_dist = std::numeric_limits<double>::infinity();
	// TODO:
	return min_dist;
}

static double linestring_polygon_distance(const geometry *lhs, const geometry *rhs) {
	SGL_ASSERT(lhs->get_type() == geometry_type::LINESTRING);
	SGL_ASSERT(rhs->get_type() == geometry_type::POLYGON);

	if(lhs->is_empty() || rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	if(rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	const auto shell = rhs->get_first_part();
	if(shell->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	return linestring_linestring_distance(lhs, shell);
}

static double polygon_polygon_distance(const geometry *lhs, const geometry *rhs) {
	SGL_ASSERT(lhs->get_type() == geometry_type::POLYGON);
	SGL_ASSERT(rhs->get_type() == geometry_type::POLYGON);

	if(lhs->is_empty() || rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	if(lhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	const auto lhs_shell = lhs->get_first_part();
	if(lhs_shell->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	if(rhs->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	const auto rhs_shell = rhs->get_first_part();
	if(rhs_shell->is_empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	return linestring_linestring_distance(lhs_shell, rhs_shell);
}

static double distance_dispatch(const geometry *lhs_p, const geometry *rhs_p) {
	SGL_ASSERT(!lhs_p->is_collection());
	SGL_ASSERT(!rhs_p->is_collection());

	switch (lhs_p->get_type()) {
	case geometry_type::POINT:
	switch (rhs_p->get_type()) {
		case geometry_type::POINT:
			return point_point_distance(lhs_p, rhs_p);
		case geometry_type::LINESTRING:
			return point_linestring_distance(lhs_p, rhs_p);
		case geometry_type::POLYGON:
			return point_polygon_distance(lhs_p, rhs_p);
		default:
			SGL_ASSERT(false);
			return std::numeric_limits<double>::quiet_NaN();
	}
	case geometry_type::LINESTRING:
		switch (rhs_p->get_type()) {
		case geometry_type::POINT:
			return point_linestring_distance(rhs_p, lhs_p);
		case geometry_type::LINESTRING:
			return linestring_linestring_distance(lhs_p, rhs_p);
		case geometry_type::POLYGON:
			return linestring_polygon_distance(lhs_p, rhs_p);
		default:
			SGL_ASSERT(false);
			return std::numeric_limits<double>::quiet_NaN();
		}
	case geometry_type::POLYGON:
		switch (rhs_p->get_type()) {
		case geometry_type::POINT:
			return point_polygon_distance(rhs_p, lhs_p);
		case geometry_type::LINESTRING:
			return linestring_polygon_distance(rhs_p, lhs_p);
		case geometry_type::POLYGON:
			return polygon_polygon_distance(lhs_p, rhs_p);
		default:
			SGL_ASSERT(false);
			return std::numeric_limits<double>::quiet_NaN();
		}
	default:
		SGL_ASSERT(false);
		return std::numeric_limits<double>::quiet_NaN();
	}
}

double distance(const geometry* lhs_p, const geometry* rhs_p) {
	SGL_ASSERT(lhs_p != nullptr);
	SGL_ASSERT(rhs_p != nullptr);

	auto lhs = lhs_p;
	auto rhs = rhs_p;

	const auto lhs_root = lhs->get_parent();
	const auto rhs_root = rhs->get_parent();

	double min_dist = std::numeric_limits<double>::infinity();

	while(lhs != lhs_root) {

		if(lhs->is_collection() && !lhs->is_empty()) {
			lhs = lhs->get_first_part();
			continue;
		}

		// Otherwise, we have a leaf on the LHS
		// I guess this is where we create an LHS index?
		// I guess it makes sense to re-order lhs and rhs depending on number of parts/verts?
		// Maybe calculate a part/vertex ratio. Although dont count interior polygon rings for that.
		// Alt just cache every calculation.

		while(rhs != rhs_root) {
			if(rhs->is_collection() && !rhs->is_empty()) {
				rhs = rhs->get_first_part();
				continue;
			}

			// If we get here, we have a leaf on both sides!
			min_dist = std::min(min_dist, distance_dispatch(lhs, rhs));

			// Now move the rhs up
			while(rhs != rhs_root) {
				const auto parent = rhs->get_parent();
				if(parent == rhs_root) {
					rhs = parent;
					break;
				}

				if(rhs != parent->get_last_part()) {
					rhs = rhs->get_next();
					break;
				}

				rhs = parent;
			}
		}

		while (lhs != lhs_root) {
			const auto parent = lhs->get_parent();
			if (parent == lhs_root) {
				lhs = parent;
				break;
			}

			if (lhs != parent->get_last_part()) {
				lhs = lhs->get_next();
				break;
			}

			lhs = parent;
		}
	}

	return min_dist;
}


//----------------------------------------------------------------------------------------------------------------------
// Validity
//----------------------------------------------------------------------------------------------------------------------

bool is_valid(const sgl::geometry *geom) {
	if(!geom) {
		return false;
	}

	const auto root = geom->get_parent();
	auto curr = geom;

	while (true) {
		switch(curr->get_type()) {
			case sgl::geometry_type::POINT: {
				// Points cant have more than one vertex
				if(curr->get_count() > 1) {
					return false;
				}
			} break;
			case sgl::geometry_type::LINESTRING: {
				// Linestrings must have zero or at least two vertices
				if(curr->get_count() == 1) {
					return false;
				}
			} break;
			case sgl::geometry_type::POLYGON: {
				const auto tail = curr->get_last_part();
				auto head = tail;
				if (!head) {
					break;
				}
				do {
					head = head->get_next();
					// Polygon rings must have at least four vertices
					if(head->get_count() < 4) {
						return false;
					}
				} while (head != tail);
			} break;
			case sgl::geometry_type::MULTI_POINT:
			case sgl::geometry_type::MULTI_LINESTRING:
			case sgl::geometry_type::MULTI_POLYGON:
			case sgl::geometry_type::MULTI_GEOMETRY: {
				if(!curr->is_empty()) {
					// Go downwards
					curr = curr->get_first_part();
					continue;
				}
			} break;
			default:
				// Just return false!
				return false;
		}

		// Inner loop
		while(true) {
			const auto parent = curr->get_parent();
			if(parent == root) {
				// Done!
				return true;
			}
			if(curr != parent->get_last_part()) {
				// Go sideways
				curr = curr->get_next();
				break;
			}
			// Go upwards
			curr = parent;
		}
	}
}



} // namespace ops

} // namespace sgl