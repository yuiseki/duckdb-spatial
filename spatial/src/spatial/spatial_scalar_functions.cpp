#include "spatial/core/function_builder.hpp"
#include "spatial/core/module.hpp"
#include "spatial/core/types.hpp"
#include "spatial/core/util/math.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"

#define SGL_ASSERT(x) D_ASSERT(x)
#include "sgl/sgl.hpp"
#include "spatial/core/functions/cast.hpp"
#include "spatial/core/geometry/wkb_writer.hpp"
#include "yyjson.h"

#include <duckdb/common/types/blob.hpp>

namespace spatial {
namespace core {

//######################################################################################################################
// Util
//######################################################################################################################

class BinaryReader {
public:
	BinaryReader(const char *ptr, const char *end) : beg(ptr), end(end), ptr(ptr) {
	}
	BinaryReader(const char *buffer, const size_t size) : BinaryReader(buffer, buffer + size) {
	}

	template <class T>
	T Read() {
		static_assert(std::is_trivially_copyable<T>::value, "Type must be trivially copyable");
		CheckSize(sizeof(T));
		T value;
		memcpy(&value, ptr, sizeof(T));
		ptr += sizeof(T);
		return value;
	}

	template <class T>
	T ReadBE() {
		static_assert(std::is_trivially_copyable<T>::value, "Type must be trivially copyable");
		CheckSize(sizeof(T));

		uint8_t in[sizeof(T)];
		uint8_t out[sizeof(T)];
		memcpy(in, ptr, sizeof(T));
		ptr += sizeof(T);

		for (size_t i = 0; i < sizeof(T); i++) {
			out[i] = in[sizeof(T) - i - 1];
		}

		T swapped = 0;
		memcpy(&swapped, out, sizeof(T));
		return swapped;
	}

	const char *Reserve(const size_t size) {
		CheckSize(size);
		const char *result = ptr;
		ptr += size;
		return result;
	}

	void Skip(const size_t size) {
		CheckSize(size);
		ptr += size;
	}

private:
	void CheckSize(const size_t size) const {
		if (ptr + size > end) {
			throw InternalException("Buffer overflow");
		}
	}

	const char *beg;
	const char *end;
	const char *ptr;
};

class BinaryWriter {
public:
	BinaryWriter(char *ptr, char *end) : beg(ptr), end(end), ptr(ptr) {
	}
	BinaryWriter(char *buffer, const size_t size) : BinaryWriter(buffer, buffer + size) {
	}

	template <class T>
	void Write(const T &value) {
		static_assert(std::is_trivially_copyable<T>::value, "Type must be trivially copyable");
		CheckSize(sizeof(T));
		memcpy(ptr, &value, sizeof(T));
		ptr += sizeof(T);
	}

	char *Reserve(const size_t size) {
		CheckSize(size);
		char *result = ptr;
		ptr += size;
		return result;
	}

	void Skip(const size_t size, const bool zero = false) {
		CheckSize(size);
		if (zero) {
			memset(ptr, 0, size);
		}
		ptr += size;
	}

	void Copy(const char *buffer, const size_t size) {
		CheckSize(size);
		memcpy(ptr, buffer, size);
		ptr += size;
	}

private:
	void CheckSize(const size_t size) const {
		if (ptr + size > end) {
			throw InternalException("Buffer overflow");
		}
	}

	char *beg;
	char *end;
	char *ptr;
};

// todo:
struct Serde {
	static size_t GetRequiredSize(const sgl::geometry &geom);
	static void Serialize(const sgl::geometry &geom, char *buffer, size_t buffer_size);
	static void Deserialize(sgl::geometry &result, ArenaAllocator &arena, const char *buffer, size_t buffer_size);
};

static size_t GetRequiredSizeInternal(const sgl::geometry *geom) {
	const auto vertex_size = geom->get_vertex_size();
	const auto part_count = geom->get_count();

	switch (geom->get_type()) {
	case sgl::geometry_type::POINT:
	case sgl::geometry_type::LINESTRING:
		// 4 bytes for the type
		// 4 bytes for the length
		// sizeof(vertex) * count;
		return 4 + 4 + part_count * vertex_size;
	case sgl::geometry_type::POLYGON: {
		// Polygons are special because they may pad between the rings and the ring data
		// 4 bytes for the type
		// 4 bytes for the length
		// sizeof(vertex) * count;
		size_t size = 4 + 4;

		const auto tail = geom->get_last_part();
		if (!tail) {
			return size;
		}
		auto part = tail;
		do {
			part = part->get_next();
			size += 4 + part->get_count() * vertex_size;
		} while (part != tail);

		if (part_count % 2 == 1) {
			size += 4;
		}
		return size;
	}
	case sgl::geometry_type::MULTI_POINT:
	case sgl::geometry_type::MULTI_LINESTRING:
	case sgl::geometry_type::MULTI_POLYGON:
	case sgl::geometry_type::MULTI_GEOMETRY: {
		// 4 bytes for the type
		// 4 bytes for the length
		// recursive call for each part
		size_t size = 4 + 4;
		const auto tail = geom->get_last_part();
		if (!tail) {
			return size;
		}
		auto part = tail;
		do {
			part = part->get_next();
			size += GetRequiredSizeInternal(part);
		} while (part != tail);
		return size;
	}
	default:
		D_ASSERT(false);
		return 0;
	}
}

size_t Serde::GetRequiredSize(const sgl::geometry &geom) {
	const auto type = geom.get_type();

	const auto has_bbox = type != sgl::geometry_type::POINT && !geom.is_empty();
	const auto has_z = geom.has_z();
	const auto has_m = geom.has_m();

	const auto dims = 2 + (has_z ? 1 : 0) + (has_m ? 1 : 0);

	const auto head_size = 4 + 4; // type + props + padding
	const auto geom_size = GetRequiredSizeInternal(&geom);
	const auto bbox_size = has_bbox ? dims * sizeof(float) * 2 : 0;

	const auto full_size = head_size + geom_size + bbox_size;

	// Check that the size is a multiple of 8
	D_ASSERT(full_size % 8 == 0);

	return full_size;
}

static void SerializeVertices(BinaryWriter &cursor, const sgl::geometry *geom, const uint32_t count, const bool has_z,
                              const bool has_m, const bool has_bbox, const uint32_t vsize, sgl::box_xyzm &bbox) {

	const auto verts = geom->get_vertex_data();

	// Copy the vertices to the cursor
	const auto dst = cursor.Reserve(count * vsize);

	if (!has_bbox) {
		// Fast path, issue on memcpy to the cursor
		memcpy(dst, verts, count * vsize);
		return;
	}

	sgl::vertex_xyzm vertex = {0};
	for (uint32_t i = 0; i < count; i++) {

		// Load the vertex from the geometry
		memcpy(&vertex, verts + i * vsize, vsize);

		// Copy the vertex to the cursor
		memcpy(dst + i * vsize, &vertex, vsize);

		bbox.min.x = std::min(bbox.min.x, vertex.x);
		bbox.min.y = std::min(bbox.min.y, vertex.y);
		bbox.max.x = std::max(bbox.max.x, vertex.x);
		bbox.max.y = std::max(bbox.max.y, vertex.y);

		if (has_z) {
			bbox.min.zm = std::min(bbox.min.zm, vertex.zm);
			bbox.max.zm = std::max(bbox.max.zm, vertex.zm);
		}
		if (has_m) {
			bbox.min.m = std::min(bbox.min.m, vertex.m);
			bbox.max.m = std::max(bbox.max.m, vertex.m);
		}
	}
}

static void SerializeRecursive(BinaryWriter &cursor, const sgl::geometry *geom, const bool has_z, const bool has_m,
                               const bool has_bbox, const uint32_t vsize, sgl::box_xyzm &bbox) {
	const auto type = geom->get_type();
	const auto count = geom->get_count();

	if (type < sgl::geometry_type::POINT || type > sgl::geometry_type::MULTI_GEOMETRY) {
		throw InvalidInputException("Cannot serialize geometry of type %d", static_cast<int>(type));
	}

	// The GeometryType enum used to start with POINT = 0
	// but now it starts with INVALID = 0, so we need to subtract 1
	cursor.Write<uint32_t>(static_cast<uint32_t>(type) - 1);
	cursor.Write<uint32_t>(count);

	switch (type) {
	case sgl::geometry_type::POINT:
	case sgl::geometry_type::LINESTRING:
		SerializeVertices(cursor, geom, count, has_z, has_m, has_bbox, vsize, bbox);
		break;
	case sgl::geometry_type::POLYGON: {
		auto ring_cursor = cursor;
		cursor.Skip((count * 4) + (count % 2 == 1 ? 4 : 0), true);

		const auto tail = geom->get_last_part();
		if (!tail) {
			break;
		}

		auto ring = tail;
		do {
			ring = ring->get_next();
			ring_cursor.Write<uint32_t>(ring->get_count());
			SerializeVertices(ring_cursor, ring, ring->get_count(), has_z, has_m, has_bbox, vsize, bbox);
		} while (ring != tail);

	} break;
	case sgl::geometry_type::MULTI_POINT:
	case sgl::geometry_type::MULTI_LINESTRING:
	case sgl::geometry_type::MULTI_POLYGON:
	case sgl::geometry_type::MULTI_GEOMETRY: {
		const auto tail = geom->get_last_part();
		if (!tail) {
			break;
		}

		auto part = tail;
		do {
			part = part->get_next();
			SerializeRecursive(cursor, part, has_z, has_m, has_bbox, vsize, bbox);
		} while (part != tail);
	} break;
	default:
		D_ASSERT(false);
	}
}

void Serde::Serialize(const sgl::geometry &geom, char *buffer, size_t buffer_size) {
	const auto type = geom.get_type();

	const auto has_bbox = type != sgl::geometry_type::POINT && !geom.is_empty();
	const auto has_z = geom.has_z();
	const auto has_m = geom.has_m();

	// Set flags
	uint8_t flags = 0;
	flags |= has_z ? 0x01 : 0;
	flags |= has_m ? 0x02 : 0;
	flags |= has_bbox ? 0x04 : 0;

	BinaryWriter cursor(buffer, buffer_size);

	cursor.Write<uint8_t>(static_cast<uint8_t>(type));
	cursor.Write<uint8_t>(flags);
	cursor.Write<uint16_t>(0); // unused for now
	cursor.Write<uint32_t>(0); // padding

	const auto dims = 2 + (has_z ? 1 : 0) + (has_m ? 1 : 0);
	const auto vert_size = dims * sizeof(double);
	const auto bbox_size = has_bbox ? dims * sizeof(float) * 2 : 0;

	// Setup a bbox to store the min/max values
	sgl::box_xyzm bbox = sgl::box_xyzm::smallest();

	auto bbox_cursor = cursor;
	cursor.Skip(bbox_size, true);

	SerializeRecursive(cursor, &geom, has_z, has_m, has_bbox, vert_size, bbox);

	if (has_bbox) {
		bbox_cursor.Write<float>(MathUtil::DoubleToFloatDown(bbox.min.x)); // xmin
		bbox_cursor.Write<float>(MathUtil::DoubleToFloatDown(bbox.min.y)); // ymin
		bbox_cursor.Write<float>(MathUtil::DoubleToFloatUp(bbox.max.x));   // xmax
		bbox_cursor.Write<float>(MathUtil::DoubleToFloatUp(bbox.max.y));   // ymax

		if (has_z) {
			bbox_cursor.Write<float>(MathUtil::DoubleToFloatDown(bbox.min.zm)); // zmin
			bbox_cursor.Write<float>(MathUtil::DoubleToFloatUp(bbox.max.zm));   // zmax
		}

		if (has_m) {
			bbox_cursor.Write<float>(MathUtil::DoubleToFloatDown(bbox.min.m)); // mmin
			bbox_cursor.Write<float>(MathUtil::DoubleToFloatUp(bbox.max.m));   // mmax
		}
	}
}

static void DeserializeRecursive(BinaryReader &cursor, sgl::geometry &geom, const bool has_z, const bool has_m,
                                 ArenaAllocator &arena) {
	const auto count = cursor.Read<uint32_t>();
	switch (geom.get_type()) {
	case sgl::geometry_type::POINT:
	case sgl::geometry_type::LINESTRING: {
		const auto verts = cursor.Reserve(count * geom.get_vertex_size());
		geom.set_vertex_data(verts, count);
	} break;
	case sgl::geometry_type::POLYGON: {
		auto ring_cursor = cursor;
		cursor.Skip((count * 4) + (count % 2 == 1 ? 4 : 0));
		for (uint32_t i = 0; i < count; i++) {
			const auto ring_count = ring_cursor.Read<uint32_t>();
			const auto verts = cursor.Reserve(ring_count * geom.get_vertex_size());

			auto ring_mem = arena.AllocateAligned(sizeof(sgl::geometry));
			const auto ring = new (ring_mem) sgl::geometry(sgl::geometry_type::LINESTRING);

			ring->set_z(has_z);
			ring->set_m(has_m);
			ring->set_vertex_data(verts, ring_count);

			geom.append_part(ring);
		}
	} break;
	case sgl::geometry_type::MULTI_POINT:
	case sgl::geometry_type::MULTI_LINESTRING:
	case sgl::geometry_type::MULTI_POLYGON:
	case sgl::geometry_type::MULTI_GEOMETRY: {
		for (uint32_t i = 0; i < count; i++) {
			const auto part_type = static_cast<sgl::geometry_type>(cursor.Read<uint32_t>() + 1);
			auto part_mem = arena.AllocateAligned(sizeof(sgl::geometry));
			const auto part = new (part_mem) sgl::geometry(part_type);
			part->set_z(has_z);
			part->set_m(has_m);
			DeserializeRecursive(cursor, *part, has_z, has_m, arena);

			geom.append_part(part);
		}
	} break;
	default:
		break;
	}
}

void Serde::Deserialize(sgl::geometry &result, ArenaAllocator &arena, const char *buffer, size_t buffer_size) {

	BinaryReader cursor(buffer, buffer_size);

	const auto type = static_cast<sgl::geometry_type>(cursor.Read<uint8_t>() + 1);
	const auto flags = cursor.Read<uint8_t>();
	cursor.Skip(sizeof(uint16_t));
	cursor.Skip(sizeof(uint32_t)); // padding

	// Parse flags
	const auto has_z = (flags & 0x01) != 0;
	const auto has_m = (flags & 0x02) != 0;
	const auto has_bbox = (flags & 0x04) != 0;

	const auto format_v1 = (flags & 0x40) != 0;
	const auto format_v0 = (flags & 0x80) != 0;

	if (format_v1 || format_v0) {
		// Unsupported version, throw an error
		throw NotImplementedException(
		    "This geometry seems to be written with a newer version of the DuckDB spatial library that is not "
		    "compatible with this version. Please upgrade your DuckDB installation.");
	}

	if (has_bbox) {
		// Skip past bbox if present
		cursor.Skip(sizeof(float) * 2 * (2 + has_z + has_m));
	}

	// Create root geometry
	result.set_type(type);
	result.set_z(has_z);
	result.set_m(has_m);

	// Read the first type
	cursor.Read<uint32_t>();

	// Deserialize the geometry
	DeserializeRecursive(cursor, result, has_z, has_m, arena);
}

//----------------------------------------------------------------------------------------------------------------------
// LocalState
//----------------------------------------------------------------------------------------------------------------------

namespace {

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

class LocalState final : public FunctionLocalState {
public:
	explicit LocalState(ClientContext &context) : arena(BufferAllocator::Get(context)), allocator(arena) {
	}

	static unique_ptr<FunctionLocalState> Init(ExpressionState &state, const BoundFunctionExpression &expr,
	                                           FunctionData *bind_data);
	static LocalState &ResetAndGet(ExpressionState &state);

	// De/Serialize geometries
	sgl::geometry Deserialize(const string_t &blob);
	string_t Serialize(Vector &vector, const sgl::geometry &geom);

	ArenaAllocator &GetArena() {
		return arena;
	}
	GeometryAllocator &GetAllocator() {
		return allocator;
	}

private:
	ArenaAllocator arena;
	GeometryAllocator allocator;
};

unique_ptr<FunctionLocalState> LocalState::Init(ExpressionState &state, const BoundFunctionExpression &expr,
                                                FunctionData *bind_data) {
	return make_uniq_base<FunctionLocalState, LocalState>(state.GetContext());
}

LocalState &LocalState::ResetAndGet(ExpressionState &state) {
	auto &local_state = ExecuteFunctionState::GetFunctionState(state)->Cast<LocalState>();
	local_state.arena.Reset();
	return local_state;
}

sgl::geometry LocalState::Deserialize(const string_t &blob) {
	sgl::geometry geom;
	Serde::Deserialize(geom, arena, blob.GetDataUnsafe(), blob.GetSize());
	return geom;
}

string_t LocalState::Serialize(Vector &vector, const sgl::geometry &geom) {
	const auto size = Serde::GetRequiredSize(geom);
	auto blob = StringVector::EmptyString(vector, size);
	Serde::Serialize(geom, blob.GetDataWriteable(), size);
	blob.Finalize();
	return blob;
}
} // namespace

namespace {

//######################################################################################################################
// Functions
//######################################################################################################################

//======================================================================================================================
// ST_Area
//======================================================================================================================

struct ST_Area {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			return sgl::ops::area(&geom);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	static void PolygonAreaFunction(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);

		auto &input = args.data[0];
		auto count = args.size();

		auto &ring_vec = ListVector::GetEntry(input);
		auto ring_entries = ListVector::GetData(ring_vec);
		auto &coord_vec = ListVector::GetEntry(ring_vec);
		auto &coord_vec_children = StructVector::GetEntries(coord_vec);
		auto x_data = FlatVector::GetData<double>(*coord_vec_children[0]);
		auto y_data = FlatVector::GetData<double>(*coord_vec_children[1]);

		UnaryExecutor::Execute<list_entry_t, double>(input, result, count, [&](list_entry_t polygon) {
			auto polygon_offset = polygon.offset;
			auto polygon_length = polygon.length;

			bool first = true;
			double area = 0;
			for (idx_t ring_idx = polygon_offset; ring_idx < polygon_offset + polygon_length; ring_idx++) {
				auto ring = ring_entries[ring_idx];
				auto ring_offset = ring.offset;
				auto ring_length = ring.length;

				double sum = 0;
				for (idx_t coord_idx = ring_offset; coord_idx < ring_offset + ring_length - 1; coord_idx++) {
					sum += (x_data[coord_idx] * y_data[coord_idx + 1]) - (x_data[coord_idx + 1] * y_data[coord_idx]);
				}
				sum = std::abs(sum);
				if (first) {
					// Add outer ring
					area = sum * 0.5;
					first = false;
				} else {
					// Subtract holes
					area -= sum * 0.5;
				}
			}
			return area;
		});

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void LineStringAreaFunction(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		UnaryExecutor::Execute<list_entry_t, double>(input, result, args.size(), [](list_entry_t) { return 0; });
	}

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D
	//------------------------------------------------------------------------------------------------------------------
	static void PointAreaFunction(DataChunk &args, ExpressionState &state, Vector &result) {
		using POINT_TYPE = StructTypeBinary<double, double>;
		using AREA_TYPE = PrimitiveType<double>;
		GenericExecutor::ExecuteUnary<POINT_TYPE, AREA_TYPE>(args.data[0], result, args.size(),
		                                                     [](POINT_TYPE) { return 0; });
	}

	//------------------------------------------------------------------------------------------------------------------
	// BOX_2D
	//------------------------------------------------------------------------------------------------------------------
	static void BoxAreaFunction(DataChunk &args, ExpressionState &state, Vector &result) {

		using BOX_TYPE = StructTypeQuaternary<double, double, double, double>;
		using AREA_TYPE = PrimitiveType<double>;

		GenericExecutor::ExecuteUnary<BOX_TYPE, AREA_TYPE>(args.data[0], result, args.size(), [&](BOX_TYPE &box) {
			auto minx = box.a_val;
			auto miny = box.b_val;
			auto maxx = box.c_val;
			auto maxy = box.d_val;
			return AREA_TYPE {(maxx - minx) * (maxy - miny)};
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr const char *DESCRIPTION = R"(
    Compute the area of a geometry.

    Returns `0.0` for any geometry that is not a `POLYGON`, `MULTIPOLYGON` or `GEOMETRYCOLLECTION` containing polygon
	geometries.

	The area is in the same units as the spatial reference system of the geometry.

    The `POINT_2D` and `LINESTRING_2D` overloads of this function always return `0.0` but are included for completeness.
	)";

	static constexpr const char *EXAMPLE = R"(
    select ST_Area('POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))'::geometry);
	-- 1.0
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {

		FunctionBuilder::RegisterScalar(db, "ST_Area", [](ScalarFunctionBuilder &func) {
			// GEOMETRY
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			// POLYGON_2D
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("polygon", GeoTypes::POLYGON_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(PolygonAreaFunction);
			});

			// LINESTRING_2D
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(LineStringAreaFunction);
			});

			// POINT_2D
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(PointAreaFunction);
			});

			// BOX_2D
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(BoxAreaFunction);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			// TODO: Set Example
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_AsGeoJSON
//======================================================================================================================

using namespace duckdb_yyjson_spatial;

class JSONAllocator {
	// Stolen from the JSON extension :)
public:
	explicit JSONAllocator(ArenaAllocator &allocator)
		: allocator(allocator), yyjson_allocator({Allocate, Reallocate, Free, &allocator}) {
	}
	yyjson_alc *GetYYJSONAllocator() {
		return &yyjson_allocator;
	}
	void Reset() {
		allocator.Reset();
	}
private:
	static void *Allocate(void *ctx, size_t size) {
		const auto alloc = static_cast<ArenaAllocator *>(ctx);
		return alloc->AllocateAligned(size);
	}
	static void *Reallocate(void *ctx, void *ptr, size_t old_size, size_t size) {
		const auto alloc = static_cast<ArenaAllocator *>(ctx);
		return alloc->ReallocateAligned(data_ptr_cast(ptr), old_size, size);
	}
	static void Free(void *ctx, void *ptr) {
		// NOP because ArenaAllocator can't free
	}
	ArenaAllocator &allocator;
	yyjson_alc yyjson_allocator;
};

struct ST_AsGeoJSON {

	//------------------------------------------------------------------------------------------------------------------
	// JSON Formatting Functions
	//------------------------------------------------------------------------------------------------------------------
	// TODO: Move these into SGL at some point, make non-recursive
	static void FormatCoords(const sgl::geometry *geom, yyjson_mut_doc *doc, yyjson_mut_val *obj) {
		const auto vertex_type = static_cast<sgl::vertex_type>(geom->has_z() + geom->has_m() * 2);
		const auto vertex_count = geom->get_count();

		// GeoJSON does not support M values, so we ignore them
		switch (vertex_type) {
		case sgl::vertex_type::XY:
		case sgl::vertex_type::XYM: {
			for(uint32_t i = 0; i < vertex_count; i++) {
				const auto coord = yyjson_mut_arr(doc);
				const auto vert = geom->get_vertex_xy(i);
				yyjson_mut_arr_add_real(doc, coord, vert.x);
				yyjson_mut_arr_add_real(doc, coord, vert.y);
				yyjson_mut_arr_append(obj, coord);
			}
		} break;
		case sgl::vertex_type::XYZ:
		case sgl::vertex_type::XYZM: {
			for(uint32_t i = 0; i < vertex_count; i++) {
				const auto coord = yyjson_mut_arr(doc);
				const auto vert = geom->get_vertex_xyzm(i);

				yyjson_mut_arr_add_real(doc, coord, vert.x);
				yyjson_mut_arr_add_real(doc, coord, vert.y);
				yyjson_mut_arr_add_real(doc, coord, vert.zm);
				yyjson_mut_arr_append(obj, coord);

			}
		} break;
		default:
			D_ASSERT(false);
			break;
		}
	}

	static void FormatRecursive(const sgl::geometry *geom, yyjson_mut_doc *doc, yyjson_mut_val *obj) {
		switch(geom->get_type()) {
			case sgl::geometry_type::POINT: {
				yyjson_mut_obj_add_str(doc, obj, "type", "Point");
				const auto coords = yyjson_mut_arr(doc);
				yyjson_mut_obj_add_val(doc, obj, "coordinates", coords);
				FormatCoords(geom, doc, coords);
			} break;
			case sgl::geometry_type::LINESTRING: {
				yyjson_mut_obj_add_str(doc, obj, "type", "LineString");
				const auto coords = yyjson_mut_arr(doc);
				yyjson_mut_obj_add_val(doc, obj, "coordinates", coords);
				FormatCoords(geom, doc, coords);
			} break;
			case sgl::geometry_type::POLYGON: {
				yyjson_mut_obj_add_str(doc, obj, "type", "Polygon");
				const auto coords = yyjson_mut_arr(doc);
				yyjson_mut_obj_add_val(doc, obj, "coordinates", coords);

				const auto tail = geom->get_last_part();
				auto head = tail;
				if(head) {
					do {
						head = head->get_next();
						const auto ring = yyjson_mut_arr(doc);
						FormatCoords(head, doc, ring);
						yyjson_mut_arr_append(coords, ring);
					} while(head != tail);
				}
			} break;
			case sgl::geometry_type::MULTI_POINT: {
				yyjson_mut_obj_add_str(doc, obj, "type", "MultiPolygon");

				const auto coords = yyjson_mut_arr(doc);
				yyjson_mut_obj_add_val(doc, obj, "coordinates", coords);

				const auto tail = geom->get_last_part();
				auto head = tail;

				if(head) {
					do {
						head = head->get_next();
						FormatCoords(head, doc, coords);
					} while(head != tail);
				}
			} break;
			case sgl::geometry_type::MULTI_LINESTRING: {
				yyjson_mut_obj_add_str(doc, obj, "type", "MultiLineString");

				const auto coords = yyjson_mut_arr(doc);
				yyjson_mut_obj_add_val(doc, obj, "coordinates", coords);

				const auto tail = geom->get_last_part();
				auto head = tail;

				if(head) {
					do {
						head = head->get_next();
						const auto line = yyjson_mut_arr(doc);
						FormatCoords(head, doc, line);
						yyjson_mut_arr_append(coords, line);
					} while(head != tail);
				}
			} break;
			case sgl::geometry_type::MULTI_POLYGON: {
				yyjson_mut_obj_add_str(doc, obj, "type", "MultiPolygon");

				const auto coords = yyjson_mut_arr(doc);
				yyjson_mut_obj_add_val(doc, obj, "coordinates", coords);

				const auto tail = geom->get_last_part();
				auto head = tail;

				if(head) {
					do {
						head = head->get_next();
						const auto poly = yyjson_mut_arr(doc);

						const auto ring_tail = head->get_last_part();
						auto ring_head = ring_tail;
						if(ring_head) {
							do {
								ring_head = ring_head->get_next();
								const auto ring = yyjson_mut_arr(doc);
								FormatCoords(ring_head, doc, ring);
								yyjson_mut_arr_append(poly, ring);
							} while(ring_head != ring_tail);
						}
						yyjson_mut_arr_append(coords, poly);
					} while(head != tail);
				}
			} break;
			case sgl::geometry_type::MULTI_GEOMETRY: {
				yyjson_mut_obj_add_str(doc, obj, "type", "GeometryCollection");

				const auto geoms = yyjson_mut_arr(doc);
				yyjson_mut_obj_add_val(doc, obj, "geometries", geoms);

				const auto tail = geom->get_last_part();
				auto head = tail;

				if(head) {
					do {
						head = head->get_next();
						const auto sub_geom = yyjson_mut_obj(doc);
						FormatRecursive(head, doc, sub_geom);
						yyjson_mut_arr_append(geoms, sub_geom);
					} while(head != tail);
				}
			} break;
			default:
				D_ASSERT(false);
				break;
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		JSONAllocator allocator(lstate.GetArena());

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](string_t &blob) {
			const auto geom = lstate.Deserialize(blob);

			const auto doc = yyjson_mut_doc_new(allocator.GetYYJSONAllocator());
			const auto obj = yyjson_mut_obj(doc);
			yyjson_mut_doc_set_root(doc, obj);

			FormatRecursive(&geom, doc, obj);

			size_t json_size = 0;
			char *json_data = yyjson_mut_write_opts(doc, 0, allocator.GetYYJSONAllocator(), &json_size, nullptr);
			// Because the arena allocator only resets after each pipeline invocation, we can safely just point into the
			// arena here without needing to copy the data to the string heap with StringVector::AddString
			return { json_data, json_size };
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
	    Returns the geometry as a GeoJSON fragment

	    This does not return a complete GeoJSON document, only the geometry fragment.
		To construct a complete GeoJSON document or feature, look into using the DuckDB JSON extension in conjunction with this function.
		This function supports geometries with Z values, but not M values. M values are ignored.
	)";

	static constexpr auto EXAMPLE = R"(
		select ST_AsGeoJSON('POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))'::geometry);
		----
		{"type":"Polygon","coordinates":[[[0.0,0.0],[0.0,1.0],[1.0,1.0],[1.0,0.0],[0.0,0.0]]]}

		-- Convert a geometry into a full GeoJSON feature (requires the JSON extension to be loaded)
		SELECT CAST({
			type: 'Feature',
			geometry: ST_AsGeoJSON(ST_Point(1,2)),
			properties: {
				name: 'my_point'
			}
		} AS JSON);
		----
		{"type":"Feature","geometry":{"type":"Point","coordinates":[1.0,2.0]},"properties":{"name":"my_point"}}
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_AsGeoJSON", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::JSON());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "conversion");
		});
	}
};

//======================================================================================================================
// ST_AsText
//======================================================================================================================

struct ST_AsText {

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePoint(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);
		auto &input = args.data[0];
		auto count = args.size();
		CoreVectorOperations::Point2DToVarchar(input, result, count);
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	// TODO: We want to format these to trim trailing zeros
	static void ExecuteLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);
		auto &input = args.data[0];
		auto count = args.size();
		CoreVectorOperations::LineString2DToVarchar(input, result, count);
	}

	//------------------------------------------------------------------------------------------------------------------
	// POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	// TODO: We want to format these to trim trailing zeros
	static void ExecutePolygon(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);
		auto count = args.size();
		auto &input = args.data[0];
		CoreVectorOperations::Polygon2DToVarchar(input, result, count);
	}

	//------------------------------------------------------------------------------------------------------------------
	// BOX_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteBox(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);
		auto count = args.size();
		auto &input = args.data[0];
		CoreVectorOperations::Box2DToVarchar(input, result, count);
	}

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	// TODO: Move this to SGL once we have proper double formatting
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);
		auto count = args.size();
		auto &input = args.data[0];
		CoreVectorOperations::GeometryToVarchar(input, result, count);
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr const char *DESCRIPTION = R"(
		Returns the geometry as a WKT string
	)";

	static constexpr const char *EXAMPLE = R"(
		SELECT ST_AsText(ST_MakeEnvelope(0,0,1,1));
		----
		POLYGON ((0 0, 0 1, 1 1, 1 0, 0 0))
	)";


	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_AsText", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::VARCHAR);

				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::VARCHAR);

				variant.SetFunction(ExecutePoint);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(LogicalType::VARCHAR);

				variant.SetFunction(ExecuteLineString);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("polygon", GeoTypes::POLYGON_2D());
				variant.SetReturnType(LogicalType::VARCHAR);

				variant.SetFunction(ExecutePolygon);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2D());
				variant.SetReturnType(LogicalType::VARCHAR);

				variant.SetFunction(ExecuteBox);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "conversion");
		});
	}
};

//======================================================================================================================
// ST_AsWKB
//======================================================================================================================

struct ST_AsWKB {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(),
			[&](const string_t &input) {
				return WKBWriter::Write(input, result);
			});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = "Returns the geometry as a WKB (Well-Known-Binary) blob";
	static constexpr auto EXAMPLE = R"(
		SELECT ST_AsWKB('POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))'::GEOMETRY)::BLOB;
		----
		\x01\x03\x00\x00\x00\x01\x00\x00\x00\x05...
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_AsWKB", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::WKB_BLOB());

				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "conversion");
		});
	}
};

//======================================================================================================================
// ST_AsHEXWKB
//======================================================================================================================

struct ST_AsHEXWKB {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		vector<data_t> buffer;
		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			buffer.clear();

			WKBWriter::Write(blob, buffer);

			auto blob_size = buffer.size() * 2; // every byte is rendered as two characters
			auto blob_str = StringVector::EmptyString(result, blob_size);
			auto blob_ptr = blob_str.GetDataWriteable();

			idx_t str_idx = 0;
			for (auto byte : buffer) {
				auto byte_a = byte >> 4;
				auto byte_b = byte & 0x0F;
				blob_ptr[str_idx++] = Blob::HEX_TABLE[byte_a];
				blob_ptr[str_idx++] = Blob::HEX_TABLE[byte_b];
			}

			blob_str.Finalize();
			return blob_str;
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr const char *DESCRIPTION = R"(
		Returns the geometry as a HEXWKB string
	)";

	static constexpr const char *EXAMPLE = R"(
		SELECT ST_AsHexWKB('POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))'::geometry);
		----
		01030000000100000005000000000000000000000000000...
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_AsHEXWKB", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::VARCHAR);

				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "conversion");
		});
	}
};

//======================================================================================================================
// ST_AsSVG
//======================================================================================================================

struct ST_AsSVG {

	//------------------------------------------------------------------------------------------------------------------
	// SVG Formatting Functions
	//------------------------------------------------------------------------------------------------------------------
	// TODO: Move this to sgl once we have proper double formatting. And make non-recursive please.

	static void FormatPoint(const sgl::geometry *geom, vector<char> &buffer, int32_t max_digits, bool rel) {
		D_ASSERT(geom->get_type() == sgl::geometry_type::POINT);
		if(geom->is_empty()) {
			return;
		}
		const auto vert = geom->get_vertex_xy(0);
		if(rel) {
			constexpr auto x = "x=\"";
			constexpr auto y = "y=\"";
			buffer.insert(buffer.end(), x, x + 3);
			MathUtil::format_coord(vert.x, buffer, max_digits);
			buffer.push_back('"');
			buffer.push_back(' ');
			buffer.insert(buffer.end(), y, y + 3);
			MathUtil::format_coord(-vert.y, buffer, max_digits);
			buffer.push_back('"');
		} else {
			constexpr auto cx = "cx=\"";
			constexpr auto cy = "cy=\"";
			buffer.insert(buffer.end(), cx, cx + 4);
			MathUtil::format_coord(vert.x, buffer, max_digits);
			buffer.push_back('"');
			buffer.push_back(' ');
			buffer.insert(buffer.end(), cy, cy + 4);
			MathUtil::format_coord(-vert.y, buffer, max_digits);
			buffer.push_back('"');
		}
	}

	static void FormatLineString(const sgl::geometry *geom, vector<char> &buffer, int32_t max_digits, bool rel, bool close) {
		D_ASSERT(geom->get_type() == sgl::geometry_type::LINESTRING);

		const auto vertex_count = geom->get_count();
		if (vertex_count == 0) {
			return;
		}

		sgl::vertex_xy last_vert = geom->get_vertex_xy(0);
		buffer.push_back('M');
		buffer.push_back(' ');
		MathUtil::format_coord(last_vert.x, -last_vert.y, buffer, max_digits);

		if(vertex_count == 1) {
			return;
		}

		buffer.push_back(' ');
		buffer.push_back(rel ? 'l' : 'L');

		if(rel) {
			for(uint32_t i = 1; i < vertex_count; i++) {
				if(i == vertex_count - 1 && close) {
					buffer.push_back(' ');
					buffer.push_back('z');
				} else {
					const auto vert = geom->get_vertex_xy(i);
					const auto delta = vert - last_vert;
					MathUtil::format_coord(delta.x, -delta.y, buffer, max_digits);
				}
			}
		} else {
			for (uint32_t i = 1; i < vertex_count; i++) {
				if (i == vertex_count - 1 && close) {
					buffer.push_back(' ');
					buffer.push_back('Z');
				} else {
					const auto vert = geom->get_vertex_xy(i);
					buffer.push_back(' ');
					MathUtil::format_coord(vert.x, -vert.y, buffer, max_digits);
				}
			}
		}
	}

	static void FormatPolygon(const sgl::geometry *geom, vector<char> &buffer, int32_t max_digits, bool rel) {
		const auto tail = geom->get_last_part();
		auto head = tail;
		if(head) {
			do {
				head = head->get_next();
				FormatLineString(head, buffer, max_digits, rel, true);
			} while(head != tail);
		}
	}

	static void FormatRecursive(const sgl::geometry *geom, vector<char> &buffer, int32_t max_digits, bool rel) {
		switch(geom->get_type()) {
		case sgl::geometry_type::POINT: FormatPoint(geom, buffer, max_digits, rel); break;
		case sgl::geometry_type::LINESTRING: FormatLineString(geom, buffer, max_digits, rel, false); break;
		case sgl::geometry_type::POLYGON: FormatPolygon(geom, buffer, max_digits, rel); break;
		case sgl::geometry_type::MULTI_POINT: {
			const auto tail = geom->get_last_part();
			auto head = tail;
			if(head) {
				do {
					head = head->get_next();
					FormatPoint(head, buffer, max_digits, rel);
					if(head != tail) {
						buffer.push_back(',');
					}
				} while(head != tail);
			}
		} break;
		case sgl::geometry_type::MULTI_LINESTRING: {
			const auto tail = geom->get_last_part();
			auto head = tail;
			if(head) {
				do {
					head = head->get_next();
					FormatLineString(head, buffer, max_digits, rel, false);
					if(head != tail) {
						buffer.push_back(' ');
					}
				} while(head != tail);
			}
		} break;
		case sgl::geometry_type::MULTI_POLYGON: {
			const auto tail = geom->get_last_part();
			auto head = tail;
			if(head) {
				do {
					head = head->get_next();
					FormatPolygon(head, buffer, max_digits, rel);
					if(head != tail) {
						buffer.push_back(' ');
					}
				} while(head != tail);
			}
		} break;
		case sgl::geometry_type::MULTI_GEOMETRY: {
			const auto tail = geom->get_last_part();
			auto head = tail;
			if(head) {
				do {
					head = head->get_next();
					FormatRecursive(head, buffer, max_digits, rel);
					if(head != tail) {
						buffer.push_back(';');
					}
				} while(head != tail);
			}
		} break;
		default:
			D_ASSERT(false);
			break;
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		vector<char> buffer;

		TernaryExecutor::Execute<string_t, bool, int32_t, string_t>(
			args.data[0], args.data[1], args.data[2], result, args.size(),
			[&](const string_t &blob, const bool rel, const int32_t max_digits) {

				// Clear buffer
				buffer.clear();

				// Deserialize geometry
				const auto geom = lstate.Deserialize(blob);

				FormatRecursive(&geom, buffer, max_digits, rel);

				return StringVector::AddString(result, buffer.data(), buffer.size());
			});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Convert the geometry into a SVG fragment or path

	    Convert the geometry into a SVG fragment or path
		The SVG fragment is returned as a string. The fragment is a path element that can be used in an SVG document.
		The second boolean argument specifies whether the path should be relative or absolute.
		The third argument specifies the maximum number of digits to use for the coordinates.

		Points are formatted as cx/cy using absolute coordinates or x/y using relative coordinates.
	)";

	static constexpr auto EXAMPLE = R"(
		SELECT ST_AsSVG('POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))'::GEOMETRY, false, 15);
		----
		M 0 0 L 0 -1 1 -1 1 0 Z
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_AsSVG", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("relative", LogicalType::BOOLEAN);
				variant.AddParameter("precision", LogicalType::INTEGER);

				variant.SetReturnType(LogicalType::VARCHAR);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "conversion");
		});
	}
};

//======================================================================================================================
// ST_Centroid
//======================================================================================================================
// The GEOMETRY version is currently implemented in the GEOS module

struct ST_Centroid {

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D
	//------------------------------------------------------------------------------------------------------------------
	// Provided for completeness sake
	static void ExecutePoint(DataChunk &args, ExpressionState &state, Vector &result) {
		result.Reference(args.data[0]);
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto count = args.size();
		UnifiedVectorFormat format;
		input.ToUnifiedFormat(count, format);

		auto line_vertex_entries = ListVector::GetData(input);
		auto &line_vertex_vec = ListVector::GetEntry(input);
		auto &line_vertex_vec_children = StructVector::GetEntries(line_vertex_vec);
		auto line_x_data = FlatVector::GetData<double>(*line_vertex_vec_children[0]);
		auto line_y_vec = FlatVector::GetData<double>(*line_vertex_vec_children[1]);

		auto &point_vertex_children = StructVector::GetEntries(result);
		auto point_x_data = FlatVector::GetData<double>(*point_vertex_children[0]);
		auto point_y_data = FlatVector::GetData<double>(*point_vertex_children[1]);
		for (idx_t out_row_idx = 0; out_row_idx < count; out_row_idx++) {

			auto in_row_idx = format.sel->get_index(out_row_idx);
			if (format.validity.RowIsValid(in_row_idx)) {
				auto line = line_vertex_entries[in_row_idx];
				auto line_offset = line.offset;
				auto line_length = line.length;

				double total_x = 0;
				double total_y = 0;
				double total_length = 0;

				// To calculate the centroid of a line, we calculate the centroid of each segment
				// and then weight the segment centroids by the length of the segment.
				// The final centroid is the sum of the weighted segment centroids divided by the total length.
				for (idx_t coord_idx = line_offset; coord_idx < line_offset + line_length - 1; coord_idx++) {
					auto x1 = line_x_data[coord_idx];
					auto y1 = line_y_vec[coord_idx];
					auto x2 = line_x_data[coord_idx + 1];
					auto y2 = line_y_vec[coord_idx + 1];

					auto segment_length = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
					total_length += segment_length;
					total_x += (x1 + x2) * 0.5 * segment_length;
					total_y += (y1 + y2) * 0.5 * segment_length;
				}

				point_x_data[out_row_idx] = total_x / total_length;
				point_y_data[out_row_idx] = total_y / total_length;

			} else {
				FlatVector::SetNull(result, out_row_idx, true);
			}
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePolygon(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto count = args.size();
		UnifiedVectorFormat format;
		input.ToUnifiedFormat(count, format);

		auto poly_entries = ListVector::GetData(input);
		auto &ring_vec = ListVector::GetEntry(input);
		auto ring_entries = ListVector::GetData(ring_vec);
		auto &vertex_vec = ListVector::GetEntry(ring_vec);
		auto &vertex_vec_children = StructVector::GetEntries(vertex_vec);
		auto x_data = FlatVector::GetData<double>(*vertex_vec_children[0]);
		auto y_data = FlatVector::GetData<double>(*vertex_vec_children[1]);

		auto &centroid_children = StructVector::GetEntries(result);
		auto centroid_x_data = FlatVector::GetData<double>(*centroid_children[0]);
		auto centroid_y_data = FlatVector::GetData<double>(*centroid_children[1]);

		for (idx_t in_row_idx = 0; in_row_idx < count; in_row_idx++) {
			if (format.validity.RowIsValid(in_row_idx)) {
				auto poly = poly_entries[in_row_idx];
				auto poly_offset = poly.offset;
				auto poly_length = poly.length;

				double poly_centroid_x = 0;
				double poly_centroid_y = 0;
				double poly_area = 0;

				// To calculate the centroid of a polygon, we calculate the centroid of each ring
				// and then weight the ring centroids by the area of the ring.
				// The final centroid is the sum of the weighted ring centroids divided by the total area.
				for (idx_t ring_idx = poly_offset; ring_idx < poly_offset + poly_length; ring_idx++) {
					auto ring = ring_entries[ring_idx];
					auto ring_offset = ring.offset;
					auto ring_length = ring.length;

					double ring_centroid_x = 0;
					double ring_centroid_y = 0;
					double ring_area = 0;

					// To calculate the centroid of a ring, we calculate the centroid of each triangle
					// and then weight the triangle centroids by the area of the triangle.
					// The final centroid is the sum of the weighted triangle centroids divided by the ring area.
					for (idx_t coord_idx = ring_offset; coord_idx < ring_offset + ring_length - 1; coord_idx++) {
						auto x1 = x_data[coord_idx];
						auto y1 = y_data[coord_idx];
						auto x2 = x_data[coord_idx + 1];
						auto y2 = y_data[coord_idx + 1];

						auto tri_area = (x1 * y2) - (x2 * y1);
						ring_centroid_x += (x1 + x2) * tri_area;
						ring_centroid_y += (y1 + y2) * tri_area;
						ring_area += tri_area;
					}
					ring_area *= 0.5;

					ring_centroid_x /= (ring_area * 6);
					ring_centroid_y /= (ring_area * 6);

					if (ring_idx == poly_offset) {
						// The first ring is the outer ring, and the remaining rings are holes.
						// For the outer ring, we add the area and centroid to the total area and centroid.
						poly_area += ring_area;
						poly_centroid_x += ring_centroid_x * ring_area;
						poly_centroid_y += ring_centroid_y * ring_area;
					} else {
						// For holes, we subtract the area and centroid from the total area and centroid.
						poly_area -= ring_area;
						poly_centroid_x -= ring_centroid_x * ring_area;
						poly_centroid_y -= ring_centroid_y * ring_area;
					}
				}
				centroid_x_data[in_row_idx] = poly_centroid_x / poly_area;
				centroid_y_data[in_row_idx] = poly_centroid_y / poly_area;
			} else {
				FlatVector::SetNull(result, in_row_idx, true);
			}
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// BOX_2D/F
	//------------------------------------------------------------------------------------------------------------------
	template<class T>
	static void ExecuteBox(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto count = args.size();
		UnifiedVectorFormat format;
		input.ToUnifiedFormat(count, format);
		auto &box_children = StructVector::GetEntries(input);
		auto minx_data = FlatVector::GetData<T>(*box_children[0]);
		auto miny_data = FlatVector::GetData<T>(*box_children[1]);
		auto maxx_data = FlatVector::GetData<T>(*box_children[2]);
		auto maxy_data = FlatVector::GetData<T>(*box_children[3]);

		auto &centroid_children = StructVector::GetEntries(result);
		auto centroid_x_data = FlatVector::GetData<double>(*centroid_children[0]);
		auto centroid_y_data = FlatVector::GetData<double>(*centroid_children[1]);

		for (idx_t out_row_idx = 0; out_row_idx < count; out_row_idx++) {
			auto in_row_idx = format.sel->get_index(out_row_idx);
			if (format.validity.RowIsValid(in_row_idx)) {
				centroid_x_data[out_row_idx] = (minx_data[in_row_idx] + maxx_data[in_row_idx]) * 0.5;
				centroid_y_data[out_row_idx] = (miny_data[in_row_idx] + maxy_data[in_row_idx]) * 0.5;
			} else {
				FlatVector::SetNull(result, out_row_idx, true);
			}
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	// TODO: add example & desc
	static constexpr auto DESCRIPTION = "";
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Centroid", [&](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.SetReturnType(GeoTypes::POINT_2D());
				variant.SetFunction(ExecutePoint);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(GeoTypes::POINT_2D());
				variant.SetFunction(ExecuteLineString);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("polygon", GeoTypes::POLYGON_2D());
				variant.SetReturnType(GeoTypes::POINT_2D());
				variant.SetFunction(ExecutePolygon);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2D());
				variant.SetReturnType(GeoTypes::POINT_2D());
				variant.SetFunction(ExecuteBox<double>);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2DF());
				variant.SetReturnType(GeoTypes::POINT_2D());
				variant.SetFunction(ExecuteBox<float>);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Collect
//======================================================================================================================

struct ST_Collect {

	//------------------------------------------------------------------------------------------------------------------
	// Execution
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		auto &child_vec = ListVector::GetEntry(args.data[0]);
		auto child_count = ListVector::GetListSize(args.data[0]);

		UnifiedVectorFormat input_vdata;
		child_vec.ToUnifiedFormat(child_count, input_vdata);

		UnaryExecutor::Execute<list_entry_t, string_t>(
		    args.data[0], result, args.size(), [&](const list_entry_t &entry) {
			    const auto offset = entry.offset;
			    const auto length = entry.length;

			    // First figure out if we have Z or M
			    bool has_z = false;
			    bool has_m = false;
			    bool all_points = true;
			    bool all_lines = true;
			    bool all_polygons = true;

			    // First pass, check if we have Z or M
			    for (idx_t out_idx = offset; out_idx < offset + length; out_idx++) {
				    const auto row_idx = input_vdata.sel->get_index(out_idx);
				    if (!input_vdata.validity.RowIsValid(row_idx)) {
					    continue;
				    }
				    auto &blob = UnifiedVectorFormat::GetData<string_t>(input_vdata)[row_idx];

				    // TODO: Peek dont deserialize
				    const auto geom = lstate.Deserialize(blob);
				    has_z = has_z || geom.has_z();
				    has_m = has_m || geom.has_m();
			    }

			    sgl::geometry collection(sgl::geometry_type::INVALID, has_z, has_m);

			    for (idx_t out_idx = offset; out_idx < offset + length; out_idx++) {
				    const auto row_idx = input_vdata.sel->get_index(out_idx);
				    if (input_vdata.validity.RowIsValid(row_idx)) {
					    continue;
				    }

				    auto &blob = UnifiedVectorFormat::GetData<string_t>(input_vdata)[row_idx];
				    // TODO: Deserialize to heap immediately
				    auto geom = lstate.Deserialize(blob);

				    // TODO: Peek dont deserialize
				    if (geom.is_empty()) {
					    continue;
				    }

				    all_points = all_points && geom.get_type() == sgl::geometry_type::POINT;
				    all_lines = all_lines && geom.get_type() == sgl::geometry_type::LINESTRING;
				    all_polygons = all_polygons && geom.get_type() == sgl::geometry_type::POLYGON;

				    // Force Z and M so that the dimensions match
				    sgl::ops::force_zm(lstate.GetAllocator(), &geom, has_z, has_m, 0, 0);

				    const auto mem = lstate.GetArena().Allocate(sizeof(sgl::geometry));
				    const auto part = new (mem) sgl::geometry(geom);

				    // Append to collection
				    collection.append_part(part);
			    }

			    // Figure out the type of the collection
			    if (all_points) {
				    collection.set_type(sgl::geometry_type::MULTI_POINT);
			    } else if (all_lines) {
				    collection.set_type(sgl::geometry_type::MULTI_LINESTRING);
			    } else if (all_polygons) {
				    collection.set_type(sgl::geometry_type::MULTI_POLYGON);
			    } else {
				    collection.set_type(sgl::geometry_type::MULTI_GEOMETRY);
			    }

			    // Serialize the collection
			    return lstate.Serialize(result, collection);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
	Collects a list of geometries into a collection geometry.
	- If all geometries are `POINT`'s, a `MULTIPOINT` is returned.
	- If all geometries are `LINESTRING`'s, a `MULTILINESTRING` is returned.
	- If all geometries are `POLYGON`'s, a `MULTIPOLYGON` is returned.
	- Otherwise if the input collection contains a mix of geometry types, a `GEOMETRYCOLLECTION` is returned.

	Empty and `NULL` geometries are ignored. If all geometries are empty or `NULL`, a `GEOMETRYCOLLECTION EMPTY` is returned.
	)";

	static constexpr auto EXAMPLE = R"(
	-- With all POINT's, a MULTIPOINT is returned
	SELECT ST_Collect([ST_Point(1, 2), ST_Point(3, 4)]);
	----
	MULTIPOINT (1 2, 3 4)

	-- With mixed geometry types, a GEOMETRYCOLLECTION is returned
	SELECT ST_Collect([ST_Point(1, 2), ST_GeomFromText('LINESTRING(3 4, 5 6)')]);
	----
	GEOMETRYCOLLECTION (POINT (1 2), LINESTRING (3 4, 5 6))

	-- Note that the empty geometry is ignored, so the result is a MULTIPOINT
	SELECT ST_Collect([ST_Point(1, 2), NULL, ST_GeomFromText('GEOMETRYCOLLECTION EMPTY')]);
	----
	MULTIPOINT (1 2)

	-- If all geometries are empty or NULL, a GEOMETRYCOLLECTION EMPTY is returned
	SELECT ST_Collect([NULL, ST_GeomFromText('GEOMETRYCOLLECTION EMPTY')]);
	----
	GEOMETRYCOLLECTION EMPTY

	-- Tip: You can use the `ST_Collect` function together with the `list()` aggregate function to collect multiple rows of geometries into a single geometry collection:

	CREATE TABLE points (geom GEOMETRY);

	INSERT INTO points VALUES (ST_Point(1, 2)), (ST_Point(3, 4));

	SELECT ST_Collect(list(geom)) FROM points;
	----
	MULTIPOINT (1 2, 3 4)
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Collect", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geoms", LogicalType::LIST(GeoTypes::GEOMETRY()));
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_CollectionExtract
//----------------------------------------------------------------------------------------------------------------------
// TODO: Implement
struct ST_CollectionExtract {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, int32_t, string_t>(
		    args.data[0], args.data[2], result, args.size(), [&](const string_t &blob, int32_t requested_type) {
			    const auto geom = lstate.Deserialize(blob);
			    const auto type = geom.get_type();
			    const auto is_empty = geom.is_empty();

			    switch (requested_type) {
			    case 1:
				    switch (type) {
				    case sgl::geometry_type::MULTI_POINT:
				    case sgl::geometry_type::POINT:
					    return blob;
				    case sgl::geometry_type::MULTI_GEOMETRY: {
					    if (is_empty) {
						    return lstate.Serialize(result, sgl::multi_point::make_empty());
					    }
				    }
				    default:
					    return lstate.Serialize(result, sgl::point::make_empty());
				    }
				    break;
			    case 2:
				    switch (type) {
				    case sgl::geometry_type::MULTI_LINESTRING:
				    case sgl::geometry_type::LINESTRING:
					    return blob;
				    case sgl::geometry_type::MULTI_GEOMETRY: {
					    if (is_empty) {
						    return lstate.Serialize(result, sgl::multi_linestring::make_empty());
					    }
				    }
				    default:
					    return lstate.Serialize(result, sgl::linestring::make_empty());
				    }
				    break;
			    case 3:
				    switch (type) {
				    case sgl::geometry_type::MULTI_POLYGON:
				    case sgl::geometry_type::POLYGON:
					    return blob;
				    case sgl::geometry_type::MULTI_GEOMETRY: {
					    if (is_empty) {
						    return lstate.Serialize(result, sgl::multi_polygon::make_empty());
					    }
				    }
				    default:
					    return lstate.Serialize(result, sgl::polygon::make_empty());
				    }
				    break;
			    default:
				    throw InvalidInputException("Invalid requested type parameter for collection extract, must be 1 "
				                                "(POINT), 2 (LINESTRING) or 3 (POLYGON)");
			    }
		    });
	}

	static void Register(DatabaseInstance &db) {
	}
};

//======================================================================================================================
// ST_Contains
//======================================================================================================================

struct ST_Contains {

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D -> POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	// TODO: This should probably be revised. Im not sure if the current implementation is entirely accurate

	static void Operation(Vector &in_point, Vector &in_polygon, Vector &result, idx_t count) {
		enum class Side { LEFT, RIGHT, ON };

		in_polygon.Flatten(count);
		in_point.Flatten(count);

		// Setup point vectors
		auto &p_children = StructVector::GetEntries(in_point);
		auto p_x_data = FlatVector::GetData<double>(*p_children[0]);
		auto p_y_data = FlatVector::GetData<double>(*p_children[1]);

		// Setup polygon vectors
		auto polygon_entries = ListVector::GetData(in_polygon);
		auto &ring_vec = ListVector::GetEntry(in_polygon);
		auto ring_entries = ListVector::GetData(ring_vec);
		auto &coord_vec = ListVector::GetEntry(ring_vec);
		auto &coord_children = StructVector::GetEntries(coord_vec);
		auto x_data = FlatVector::GetData<double>(*coord_children[0]);
		auto y_data = FlatVector::GetData<double>(*coord_children[1]);

		auto result_data = FlatVector::GetData<bool>(result);

		for (idx_t polygon_idx = 0; polygon_idx < count; polygon_idx++) {
			auto polygon = polygon_entries[polygon_idx];
			auto polygon_offset = polygon.offset;
			auto polygon_length = polygon.length;
			bool first = true;

			// does the point lie inside the polygon?
			bool contains = false;

			auto x = p_x_data[polygon_idx];
			auto y = p_y_data[polygon_idx];

			for (idx_t ring_idx = polygon_offset; ring_idx < polygon_offset + polygon_length; ring_idx++) {
				auto ring = ring_entries[ring_idx];
				auto ring_offset = ring.offset;
				auto ring_length = ring.length;

				auto x1 = x_data[ring_offset];
				auto y1 = y_data[ring_offset];
				int winding_number = 0;

				for (idx_t coord_idx = ring_offset + 1; coord_idx < ring_offset + ring_length; coord_idx++) {
					// foo foo foo
					auto x2 = x_data[coord_idx];
					auto y2 = y_data[coord_idx];

					if (x1 == x2 && y1 == y2) {
						x1 = x2;
						y1 = y2;
						continue;
					}

					auto y_min = std::min(y1, y2);
					auto y_max = std::max(y1, y2);

					if (y > y_max || y < y_min) {
						x1 = x2;
						y1 = y2;
						continue;
					}

					auto side = Side::ON;
					double side_v = ((x - x1) * (y2 - y1) - (x2 - x1) * (y - y1));
					if (side_v == 0) {
						side = Side::ON;
					} else if (side_v < 0) {
						side = Side::LEFT;
					} else {
						side = Side::RIGHT;
					}

					if (side == Side::ON &&
					    (((x1 <= x && x < x2) || (x1 >= x && x > x2)) || ((y1 <= y && y < y2) || (y1 >= y && y > y2)))) {

						// return Contains::ON_EDGE;
						contains = false;
						break;
					} else if (side == Side::LEFT && (y1 < y && y <= y2)) {
						winding_number++;
					} else if (side == Side::RIGHT && (y2 <= y && y < y1)) {
						winding_number--;
					}

					x1 = x2;
					y1 = y2;
				}
				bool in_ring = winding_number != 0;
				if (first) {
					if (!in_ring) {
						// if the first ring is not inside, then the point is not inside the polygon
						contains = false;
						break;
					} else {
						// if the first ring is inside, then the point is inside the polygon
						// but might be inside a hole, so we continue
						contains = true;
					}
				} else {
					if (in_ring) {
						// if the hole is inside, then the point is not inside the polygon
						contains = false;
						break;
					} // else continue
				}
				first = false;
			}
			result_data[polygon_idx] = contains;
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		Operation(args.data[0], args.data[1], result, args.size());
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------

	// TODO: Add example
	static constexpr auto DESCRIPTION = "";
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Contains", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::POLYGON_2D());
				variant.AddParameter("geom2", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

//======================================================================================================================
// ST_Dimension
//======================================================================================================================

struct ST_Dimension {

	//------------------------------------------------------------------------------------------------------------------
	// Execute
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			return sgl::ops::max_surface_dimension(&geom);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = "Returns the dimension of a geometry.";

	static constexpr auto EXAMPLE = R"(
	select st_dimension('POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))'::geometry);
	----
	2
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Dimension", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::INTEGER);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Distance
//======================================================================================================================

struct ST_Distance {

	//------------------------------------------------------------------------------------------------------------------
	// Helpers
	//------------------------------------------------------------------------------------------------------------------
	// TODO: Move this to SGL, into VectorOperations or deprecate.
	template <class T>
	static PointXY<T> ClosestPointOnSegment(const PointXY<T> &p, const PointXY<T> &p1, const PointXY<T> &p2) {
		// If the segment is a Vertex, then return that Vertex
		if (p1.ApproxEqualTo(p2)) {
			return p1;
		}
		auto n1 = ((p.x - p1.x) * (p2.x - p1.x) + (p.y - p1.y) * (p2.y - p1.y));
		auto n2 = ((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
		auto r = n1 / n2;
		// If r is less than 0, then the Point is outside the segment in the p1 direction
		if (r <= 0) {
			return p1;
		}
		// If r is greater than 1, then the Point is outside the segment in the p2 direction
		if (r >= 1) {
			return p2;
		}
		// Interpolate between p1 and p2
		return PointXY<T>(p1.x + r * (p2.x - p1.x), p1.y + r * (p2.y - p1.y));
	}

	template <class T>
	static double DistanceToSegmentSquared(const PointXY<T> &px, const PointXY<T> &ax, const PointXY<T> &bx) {
		auto point = ClosestPointOnSegment(px, ax, bx);
		auto dx = px.x - point.x;
		auto dy = px.y - point.y;
		return dx * dx + dy * dy;
	}

	//------------------------------------------------------------------------------
	// POINT_2D/POINT_2D
	//------------------------------------------------------------------------------
	static void ExecutePointPoint(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 2);
		auto &left = args.data[0];
		auto &right = args.data[1];
		auto count = args.size();

		left.Flatten(count);
		right.Flatten(count);

		auto &left_entries = StructVector::GetEntries(left);
		auto &right_entries = StructVector::GetEntries(right);

		auto left_x = FlatVector::GetData<double>(*left_entries[0]);
		auto left_y = FlatVector::GetData<double>(*left_entries[1]);
		auto right_x = FlatVector::GetData<double>(*right_entries[0]);
		auto right_y = FlatVector::GetData<double>(*right_entries[1]);

		auto out_data = FlatVector::GetData<double>(result);
		for (idx_t i = 0; i < count; i++) {
			out_data[i] = std::sqrt(std::pow(left_x[i] - right_x[i], 2) + std::pow(left_y[i] - right_y[i], 2));
		}

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------
	// POINT_2D/LINESTRING_2D
	//------------------------------------------------------------------------------
	static void PointLineStringOperation(Vector &in_point, Vector &in_line, Vector &result, idx_t count) {

		// Set up the point vectors
		in_point.Flatten(count);
		auto &p_children = StructVector::GetEntries(in_point);
		auto &p_x = p_children[0];
		auto &p_y = p_children[1];
		auto p_x_data = FlatVector::GetData<double>(*p_x);
		auto p_y_data = FlatVector::GetData<double>(*p_y);

		// Set up the line vectors
		in_line.Flatten(count);

		auto &inner = ListVector::GetEntry(in_line);
		auto &children = StructVector::GetEntries(inner);
		auto &x = children[0];
		auto &y = children[1];
		auto x_data = FlatVector::GetData<double>(*x);
		auto y_data = FlatVector::GetData<double>(*y);
		auto lines = ListVector::GetData(in_line);

		auto result_data = FlatVector::GetData<double>(result);
		for (idx_t i = 0; i < count; i++) {
			auto offset = lines[i].offset;
			auto length = lines[i].length;

			double min_distance = std::numeric_limits<double>::max();
			auto p = PointXY<double>(p_x_data[i], p_y_data[i]);

			// Loop over the segments and find the closes one to the point
			for (idx_t j = 0; j < length - 1; j++) {
				auto a = PointXY<double>(x_data[offset + j], y_data[offset + j]);
				auto b = PointXY<double>(x_data[offset + j + 1], y_data[offset + j + 1]);

				auto distance = DistanceToSegmentSquared(p, a, b);
				if (distance < min_distance) {
					min_distance = distance;

					if (min_distance == 0) {
						break;
					}
				}
			}
			result_data[i] = std::sqrt(min_distance);
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	static void ExecutePointLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 2);
		auto &in_point = args.data[0];
		auto &in_line = args.data[1];
		auto count = args.size();
		PointLineStringOperation(in_point, in_line, result, count);
	}

	static void ExecuteLineStringPoint(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 2);
		auto &in_line = args.data[0];
		auto &in_point = args.data[1];
		auto count = args.size();
		PointLineStringOperation(in_point, in_line, result, count);
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	// TODO: add example/description
	static constexpr auto DESCRIPTION = "";
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Distance", [](ScalarFunctionBuilder &func) {

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point1", GeoTypes::POINT_2D());
				variant.AddParameter("point2", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::DOUBLE);
				variant.SetFunction(ExecutePointPoint);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(LogicalType::DOUBLE);
				variant.SetFunction(ExecutePointLineString);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::DOUBLE);
				variant.SetFunction(ExecuteLineStringPoint);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Dump
//======================================================================================================================

struct ST_Dump {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		auto &arena = lstate.GetArena();
		auto count = args.size();

		auto &geom_vec = args.data[0];
		UnifiedVectorFormat geom_format;
		geom_vec.ToUnifiedFormat(count, geom_format);

		idx_t total_geom_count = 0;
		idx_t total_path_count = 0;

		for (idx_t out_row_idx = 0; out_row_idx < count; out_row_idx++) {
			auto in_row_idx = geom_format.sel->get_index(out_row_idx);

			if (!geom_format.validity.RowIsValid(in_row_idx)) {
				FlatVector::SetNull(result, out_row_idx, true);
				continue;
			}

			auto &blob = UnifiedVectorFormat::GetData<string_t>(geom_format)[in_row_idx];
			auto geom = lstate.Deserialize(blob);

			vector<std::tuple<sgl::geometry*, vector<int32_t>>> stack;
			vector<std::tuple<sgl::geometry*, vector<int32_t>>> items;

			stack.emplace_back(geom, vector<int32_t>());

			while (!stack.empty()) {
				auto current = stack.back();
				auto current_geom = std::get<0>(current);
				auto current_path = std::get<1>(current);

				stack.pop_back();
				if (current_geom->is_collection()) {
					// Push all children

					auto head = current_geom->get_first_part();
					if(!head) {
						continue;
					}

					stack.emplace_back(head, current_path);

					while(head != current_geom->get_last_part()) {
						head = head->get_next();
						stack.emplace_back(head, current_path);
					}

				} else {
					items.push_back(current);
				}
			}

			// Finally reverse the results
			std::reverse(items.begin(), items.end());

			// Push to the result vector
			auto result_entries = ListVector::GetData(result);

			auto geom_offset = total_geom_count;
			auto geom_length = items.size();

			result_entries[out_row_idx].length = geom_length;
			result_entries[out_row_idx].offset = geom_offset;

			total_geom_count += geom_length;

			ListVector::Reserve(result, total_geom_count);
			ListVector::SetListSize(result, total_geom_count);

			auto &result_list = ListVector::GetEntry(result);
			auto &result_list_children = StructVector::GetEntries(result_list);
			auto &result_geom_vec = result_list_children[0];
			auto &result_path_vec = result_list_children[1];

			// The child geometries must share the same properties as the parent geometry
			auto geom_data = FlatVector::GetData<string_t>(*result_geom_vec);
			for (idx_t i = 0; i < geom_length; i++) {
				// Write the geometry
				auto item_blob = std::get<0>(items[i]);
				geom_data[geom_offset + i] = lstate.Serialize(*result_geom_vec, *item_blob);

				// Now write the paths
				auto &path = std::get<1>(items[i]);
				auto path_offset = total_path_count;
				auto path_length = path.size();

				total_path_count += path_length;

				ListVector::Reserve(*result_path_vec, total_path_count);
				ListVector::SetListSize(*result_path_vec, total_path_count);

				auto path_entries = ListVector::GetData(*result_path_vec);

				path_entries[geom_offset + i].offset = path_offset;
				path_entries[geom_offset + i].length = path_length;

				auto &path_data_vec = ListVector::GetEntry(*result_path_vec);
				auto path_data = FlatVector::GetData<int32_t>(path_data_vec);

				for (idx_t j = 0; j < path_length; j++) {
					path_data[path_offset + j] = path[j];
				}
			}
		}

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
	Dumps a geometry into a list of sub-geometries and their "path" in the original geometry.
	)";

	static constexpr auto EXAMPLE = R"(
	select st_dump('MULTIPOINT(1 2,3 4)'::geometry);
	----
	[{'geom': 'POINT(1 2)', 'path': [0]}, {'geom': 'POINT(3 4)', 'path': [1]}]
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {

		FunctionBuilder::RegisterScalar(db, "ST_Dump", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());

				variant.SetReturnType(LogicalType::LIST(
				LogicalType::STRUCT({
					{"geom", GeoTypes::GEOMETRY()},
					{"path", LogicalType::LIST(
						LogicalType::INTEGER)}})));

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//======================================================================================================================
// ST_Extent
//======================================================================================================================
// TODO: WKB Implementation

struct ST_Extent {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		const auto &bbox_vec = StructVector::GetEntries(result);
		const auto min_x_data = FlatVector::GetData<double>(*bbox_vec[0]);
		const auto min_y_data = FlatVector::GetData<double>(*bbox_vec[1]);
		const auto max_x_data = FlatVector::GetData<double>(*bbox_vec[2]);
		const auto max_y_data = FlatVector::GetData<double>(*bbox_vec[3]);

		UnifiedVectorFormat input_vdata;
		args.data[0].ToUnifiedFormat(args.size(), input_vdata);
		const auto input_data = UnifiedVectorFormat::GetData<geometry_t>(input_vdata);

		const auto count = args.size();

		for (idx_t out_idx = 0; out_idx < count; out_idx++) {
			const auto row_idx = input_vdata.sel->get_index(out_idx);
			if (!input_vdata.validity.RowIsValid(row_idx)) {
				// null in -> null out
				FlatVector::SetNull(result, out_idx, true);
				continue;
			}

			const auto &blob = input_data[row_idx];
			const auto geom = lstate.Deserialize(blob);

			auto bbox = sgl::box_xy::smallest();

			if (!sgl::ops::try_get_extent_xy(&geom, &bbox)) {
				// no vertices -> no extent -> return null
				FlatVector::SetNull(result, out_idx, true);
				continue;
			}

			min_x_data[out_idx] = bbox.min.x;
			min_y_data[out_idx] = bbox.min.y;
			max_x_data[out_idx] = bbox.max.x;
			max_y_data[out_idx] = bbox.max.y;
		}

		if (args.AllConstant()) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Returns the minimal bounding box enclosing the input geometry
	)";

	// TODO: Example
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Extent", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::BOX_2D());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_ExteriorRing
//======================================================================================================================

struct ST_ExteriorRing {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
		    args.data[0], result, args.size(), [&](const string_t &blob, ValidityMask &mask, const idx_t idx) {
			    // TODO: Peek dont deserialize
			    const auto geom = lstate.Deserialize(blob);

			    if (geom.get_type() != sgl::geometry_type::POLYGON) {
				    mask.SetInvalid(idx);
				    return string_t {};
			    }

			    if (geom.is_empty()) {
				    const auto empty = sgl::polygon::make_empty(geom.has_z(), geom.has_m());
				    return lstate.Serialize(result, empty);
			    }

			    const auto shell = geom.get_first_part();
			    return lstate.Serialize(result, *shell);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePolygon(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);
		auto &poly_vec = args.data[0];
		auto poly_entries = ListVector::GetData(poly_vec);
		auto &ring_vec = ListVector::GetEntry(poly_vec);
		auto ring_entries = ListVector::GetData(ring_vec);
		auto &vertex_vec = ListVector::GetEntry(ring_vec);
		auto &vertex_vec_children = StructVector::GetEntries(vertex_vec);
		auto poly_x_data = FlatVector::GetData<double>(*vertex_vec_children[0]);
		auto poly_y_data = FlatVector::GetData<double>(*vertex_vec_children[1]);

		auto count = args.size();
		UnifiedVectorFormat poly_format;
		poly_vec.ToUnifiedFormat(count, poly_format);

		// First figure out how many vertices we need
		idx_t total_vertex_count = 0;
		for (idx_t i = 0; i < count; i++) {
			auto row_idx = poly_format.sel->get_index(i);
			if (poly_format.validity.RowIsValid(row_idx)) {
				auto poly = poly_entries[row_idx];
				if (poly.length != 0) {
					// We only care about the exterior ring (first entry)
					auto &ring = ring_entries[poly.offset];
					total_vertex_count += ring.length;
				}
			}
		}

		// Now we can allocate the result vector
		auto &line_vec = result;
		ListVector::Reserve(line_vec, total_vertex_count);
		ListVector::SetListSize(line_vec, total_vertex_count);

		auto line_entries = ListVector::GetData(line_vec);
		auto &line_coord_vec = StructVector::GetEntries(ListVector::GetEntry(line_vec));
		auto line_data_x = FlatVector::GetData<double>(*line_coord_vec[0]);
		auto line_data_y = FlatVector::GetData<double>(*line_coord_vec[1]);

		// Now we can fill the result vector
		idx_t line_data_offset = 0;
		for (idx_t i = 0; i < count; i++) {
			auto row_idx = poly_format.sel->get_index(i);
			if (poly_format.validity.RowIsValid(row_idx)) {
				auto poly = poly_entries[row_idx];

				if (poly.length == 0) {
					line_entries[i].offset = 0;
					line_entries[i].length = 0;
					continue;
				}

				// We only care about the exterior ring (first entry)
				auto &ring = ring_entries[poly.offset];

				auto &line_entry = line_entries[i];
				line_entry.offset = line_data_offset;
				line_entry.length = ring.length;

				for (idx_t coord_idx = 0; coord_idx < ring.length; coord_idx++) {
					line_data_x[line_entry.offset + coord_idx] = poly_x_data[ring.offset + coord_idx];
					line_data_y[line_entry.offset + coord_idx] = poly_y_data[ring.offset + coord_idx];
				}

				line_data_offset += ring.length;
			} else {
				FlatVector::SetNull(line_vec, i, true);
			}
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = "Returns the exterior ring (shell) of a polygon geometry.";

	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_ExteriorRing", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("polygon", GeoTypes::POLYGON_2D());
				variant.SetReturnType(GeoTypes::LINESTRING_2D());

				variant.SetFunction(ExecutePolygon);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_FlipCoordinates
//======================================================================================================================

struct ST_FlipCoordinates {

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D
	//------------------------------------------------------------------------------------------------------------------
	// TODO: We should be able to optimize these and avoid the flatten
	static void ExecutePoint(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto count = args.size();

		// TODO: Avoid flatten
		input.Flatten(count);

		auto &coords_in = StructVector::GetEntries(input);
		auto x_data_in = FlatVector::GetData<double>(*coords_in[0]);
		auto y_data_in = FlatVector::GetData<double>(*coords_in[1]);

		auto &coords_out = StructVector::GetEntries(result);
		auto x_data_out = FlatVector::GetData<double>(*coords_out[0]);
		auto y_data_out = FlatVector::GetData<double>(*coords_out[1]);

		memcpy(x_data_out, y_data_in, count * sizeof(double));
		memcpy(y_data_out, x_data_in, count * sizeof(double));

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto count = args.size();

		// TODO: Avoid flatten
		input.Flatten(count);

		auto coord_vec_in = ListVector::GetEntry(input);
		auto &coords_in = StructVector::GetEntries(coord_vec_in);
		auto x_data_in = FlatVector::GetData<double>(*coords_in[0]);
		auto y_data_in = FlatVector::GetData<double>(*coords_in[1]);

		auto coord_count = ListVector::GetListSize(input);
		ListVector::Reserve(result, coord_count);
		ListVector::SetListSize(result, coord_count);

		auto line_entries_in = ListVector::GetData(input);
		auto line_entries_out = ListVector::GetData(result);
		memcpy(line_entries_out, line_entries_in, count * sizeof(list_entry_t));

		auto coord_vec_out = ListVector::GetEntry(result);
		auto &coords_out = StructVector::GetEntries(coord_vec_out);
		auto x_data_out = FlatVector::GetData<double>(*coords_out[0]);
		auto y_data_out = FlatVector::GetData<double>(*coords_out[1]);

		memcpy(x_data_out, y_data_in, coord_count * sizeof(double));
		memcpy(y_data_out, x_data_in, coord_count * sizeof(double));

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePolygon(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto count = args.size();

		// TODO: Avoid flatten
		input.Flatten(count);

		auto ring_vec_in = ListVector::GetEntry(input);
		auto ring_count = ListVector::GetListSize(input);

		auto coord_vec_in = ListVector::GetEntry(ring_vec_in);
		auto &coords_in = StructVector::GetEntries(coord_vec_in);
		auto x_data_in = FlatVector::GetData<double>(*coords_in[0]);
		auto y_data_in = FlatVector::GetData<double>(*coords_in[1]);

		auto coord_count = ListVector::GetListSize(ring_vec_in);

		ListVector::Reserve(result, ring_count);
		ListVector::SetListSize(result, ring_count);
		auto ring_vec_out = ListVector::GetEntry(result);
		ListVector::Reserve(ring_vec_out, coord_count);
		ListVector::SetListSize(ring_vec_out, coord_count);

		auto ring_entries_in = ListVector::GetData(input);
		auto ring_entries_out = ListVector::GetData(result);
		memcpy(ring_entries_out, ring_entries_in, count * sizeof(list_entry_t));

		auto coord_entries_in = ListVector::GetData(ring_vec_in);
		auto coord_entries_out = ListVector::GetData(ring_vec_out);
		memcpy(coord_entries_out, coord_entries_in, ring_count * sizeof(list_entry_t));

		auto coord_vec_out = ListVector::GetEntry(ring_vec_out);
		auto &coords_out = StructVector::GetEntries(coord_vec_out);
		auto x_data_out = FlatVector::GetData<double>(*coords_out[0]);
		auto y_data_out = FlatVector::GetData<double>(*coords_out[1]);

		memcpy(x_data_out, y_data_in, coord_count * sizeof(double));
		memcpy(y_data_out, x_data_in, coord_count * sizeof(double));

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// BOX_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteBox(DataChunk &args, ExpressionState &state, Vector &result) {

		auto input = args.data[0];
		auto count = args.size();

		// TODO: Avoid flatten
		input.Flatten(count);

		auto &children_in = StructVector::GetEntries(input);
		auto min_x_in = FlatVector::GetData<double>(*children_in[0]);
		auto min_y_in = FlatVector::GetData<double>(*children_in[1]);
		auto max_x_in = FlatVector::GetData<double>(*children_in[2]);
		auto max_y_in = FlatVector::GetData<double>(*children_in[3]);

		auto &children_out = StructVector::GetEntries(result);
		auto min_x_out = FlatVector::GetData<double>(*children_out[0]);
		auto min_y_out = FlatVector::GetData<double>(*children_out[1]);
		auto max_x_out = FlatVector::GetData<double>(*children_out[2]);
		auto max_y_out = FlatVector::GetData<double>(*children_out[3]);

		memcpy(min_x_out, min_y_in, count * sizeof(double));
		memcpy(min_y_out, min_x_in, count * sizeof(double));
		memcpy(max_x_out, max_y_in, count * sizeof(double));
		memcpy(max_y_out, max_x_in, count * sizeof(double));
	}

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	// TODO: Move this to SGL, make non-recursive
	static void FlipPoint(ArenaAllocator &alloc, sgl::geometry *geom) {
		if(!geom->is_empty()) {
			const auto vertex_count = geom->get_count();
			const auto vertex_size = geom->get_vertex_size();
			const auto vertex_data = geom->get_vertex_data();

			// Copy the vertex data
			const auto new_vertex_data = alloc.AllocateAligned(vertex_count * vertex_size);
			memcpy(new_vertex_data, vertex_data, vertex_count * vertex_size);

			// Flip the x and y coordinates
			const auto vertex_ptr = reinterpret_cast<double*>(new_vertex_data);
			std::swap(vertex_ptr[0], vertex_ptr[1]);

			// Update the vertex data
			geom->set_vertex_data(new_vertex_data, 1);
		}
	}

	static void FlipLineString(ArenaAllocator &alloc, sgl::geometry *geom) {
		if(!geom->is_empty()) {
			const auto vertex_count = geom->get_count();
			const auto vertex_size = geom->get_vertex_size();
			const auto vertex_data = geom->get_vertex_data();

			// Copy the vertex data
			const auto new_vertex_data = alloc.AllocateAligned(vertex_count * vertex_size);
			memcpy(new_vertex_data, vertex_data, vertex_count * vertex_size);

			// Flip the x and y coordinates
			for(idx_t i = 0; i < vertex_count; i++) {
				const auto x_ptr = reinterpret_cast<double*>(new_vertex_data + i * vertex_size);
				const auto y_ptr = reinterpret_cast<double*>(new_vertex_data + i * vertex_size + sizeof(double));

				std::swap(*x_ptr, *y_ptr);
			}

			// Update the vertex data
			geom->set_vertex_data(new_vertex_data, vertex_count);
		}
	}

	static void FlipPolygon(ArenaAllocator &alloc, sgl::geometry *geom) {
		const auto tail = geom->get_last_part();
		auto head = tail;
		if(head) {
			do {
				head = head->get_next();
				FlipLineString(alloc, head);
			} while(head != tail);
		}
	}

	static void FlipRecursive(ArenaAllocator &alloc, sgl::geometry *geom) {
		switch(geom->get_type()) {
			case sgl::geometry_type::POINT: FlipPoint(alloc, geom); break;
			case sgl::geometry_type::LINESTRING: FlipLineString(alloc, geom); break;
			case sgl::geometry_type::POLYGON: FlipPolygon(alloc, geom); break;
			case sgl::geometry_type::MULTI_POINT: {
				const auto tail = geom->get_last_part();
				auto head = tail;
				if(head) {
					do {
						FlipPoint(alloc, head);
						head = head->get_next();
					} while(head != tail);
				}
			} break;
			case sgl::geometry_type::MULTI_LINESTRING: {
				const auto tail = geom->get_last_part();
				auto head = tail;
				if(head) {
					do {
						FlipLineString(alloc, head);
						head = head->get_next();
					} while(head != tail);
				}
			} break;
			case sgl::geometry_type::MULTI_POLYGON: {
				const auto tail = geom->get_last_part();
				auto head = tail;
				if(head) {
					do {
						FlipPolygon(alloc, head);
						head = head->get_next();
					} while(head != tail);
				}
			} break;
			case sgl::geometry_type::MULTI_GEOMETRY: {
				const auto tail = geom->get_last_part();
				auto head = tail;
				if(head) {
					do {
						FlipRecursive(alloc, head);
						head = head->get_next();
					} while(head != tail);
				}
			} break;
			default:
				D_ASSERT(false);
				break;
		}
	}


	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {

		auto input = args.data[0];
		auto count = args.size();

		UnaryExecutor::Execute<string_t, string_t>(input, result, count, [&](const string_t &blob) {
			// This is pretty memory intensive, so reset arena after each call
			auto &lstate = LocalState::ResetAndGet(state);
			auto &arena = lstate.GetArena();

			// Deserialize the geometry
			auto geom = lstate.Deserialize(blob);
			// Flip the coordinates
			FlipRecursive(arena, &geom);

			// Serialize the result
			return lstate.Serialize(result, geom);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Description
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Returns a new geometry with the coordinates of the input geometry "flipped" so that x = y and y = x
	)";

	// TODO: Add example
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_FlipCoordinates", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.SetReturnType(GeoTypes::POINT_2D());

				variant.SetFunction(ExecutePoint);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(GeoTypes::LINESTRING_2D());

				variant.SetFunction(ExecuteLineString);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("polygon", GeoTypes::POLYGON_2D());
				variant.SetReturnType(GeoTypes::POLYGON_2D());

				variant.SetFunction(ExecutePolygon);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2D());
				variant.SetReturnType(GeoTypes::BOX_2D());

				variant.SetFunction(ExecuteBox);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//======================================================================================================================
// ST_Force 2D/3DZ/3DM/4D
//======================================================================================================================

template <class IMPL>
struct ST_ForceBase {

	//------------------------------------------------------------------------------------------------------------------
	// Execute
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		auto &alloc = lstate.GetAllocator();

		auto has_z = IMPL::HAS_Z;
		auto has_m = IMPL::HAS_M;

		auto &input = args.data[0];
		const auto count = args.size();

		// TODO: This can be optimized to avoid de/serialization if the vertex type already matches

		if (has_z && has_m) {
			auto &z_values = args.data[1];
			auto &m_values = args.data[2];

			TernaryExecutor::Execute<string_t, double, double, string_t>(
			    input, z_values, m_values, result, count, [&](const string_t &blob, double z, double m) {
				    auto geom = lstate.Deserialize(blob);
				    sgl::ops::force_zm(alloc, &geom, true, true, z, m);
				    return lstate.Serialize(result, geom);
			    });

			return;
		}

		if (has_z || has_m) {
			auto &zm_values = args.data[1];

			BinaryExecutor::Execute<string_t, double, string_t>(
			    input, zm_values, result, count, [&](const string_t &blob, double zm) {
				    const auto def_z = has_z ? zm : 0;
				    const auto def_m = has_m ? zm : 0;

				    auto geom = lstate.Deserialize(blob);
				    sgl::ops::force_zm(alloc, &geom, has_z, has_m, def_z, def_m);
				    return lstate.Serialize(result, geom);
			    });

			return;
		}

		UnaryExecutor::Execute<string_t, string_t>(input, result, count, [&](const string_t &blob) {
			auto geom = lstate.Deserialize(blob);
			sgl::ops::force_zm(alloc, &geom, false, false, 0, 0);
			return lstate.Serialize(result, geom);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar([](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(IMPL::DESCRIPTION);
			func.SetExample(IMPL::EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Force2D : ST_ForceBase<ST_Force2D> {
	static auto constexpr NAME = "ST_Force2D";
	static auto constexpr HAS_Z = false;
	static auto constexpr HAS_M = false;
	static auto constexpr EXAMPLE = "";
	static auto constexpr DESCRIPTION = R"(
		Forces the vertices of a geometry to have X and Y components

		This function will drop any Z and M values from the input geometry, if present. If the input geometry is already 2D, it will be returned as is.
		)";
};

struct ST_Force3DZ : ST_ForceBase<ST_Force3DZ> {
	static auto constexpr NAME = "ST_Force3DZ";
	static auto constexpr HAS_Z = true;
	static auto constexpr HAS_M = false;
	static auto constexpr EXAMPLE = "";
	static auto constexpr DESCRIPTION = R"(
		Forces the vertices of a geometry to have X, Y and Z components

		The following cases apply:
		- If the input geometry has a M component but no Z component, the M component will be replaced with the new Z value.
		- If the input geometry has a Z component but no M component, it will be returned as is.
		- If the input geometry has both a Z component and a M component, the M component will be removed.
		- Otherwise, if the input geometry has neither a Z or M component, the new Z value will be added to the vertices of the input geometry.
		)";
};

struct ST_Force3DM : ST_ForceBase<ST_Force3DM> {
	static auto constexpr NAME = "ST_Force3DM";
	static auto constexpr HAS_Z = false;
	static auto constexpr HAS_M = true;
	static auto constexpr EXAMPLE = "";
	static auto constexpr DESCRIPTION = R"(
		Forces the vertices of a geometry to have X, Y and M components

		The following cases apply:
		- If the input geometry has a Z component but no M component, the Z component will be replaced with the new M value.
		- If the input geometry has a M component but no Z component, it will be returned as is.
		- If the input geometry has both a Z component and a M component, the Z component will be removed.
		- Otherwise, if the input geometry has neither a Z or M component, the new M value will be added to the vertices of the input geometry.
		)";
};

struct ST_Force4D : ST_ForceBase<ST_Force4D> {
	static auto constexpr NAME = "ST_Force4D";
	static auto constexpr HAS_Z = true;
	static auto constexpr HAS_M = true;
	static auto constexpr EXAMPLE = "";
	static auto constexpr DESCRIPTION = R"(
		Forces the vertices of a geometry to have X, Y, Z and M components

		The following cases apply:
		- If the input geometry has a Z component but no M component, the new M value will be added to the vertices of the input geometry.
		- If the input geometry has a M component but no Z component, the new Z value will be added to the vertices of the input geometry.
		- If the input geometry has both a Z component and a M component, the geometry will be returned as is.
		- Otherwise, if the input geometry has neither a Z or M component, the new Z and M values will be added to the vertices of the input geometry.
		)";
};

//======================================================================================================================
// ST_GeometryType
//======================================================================================================================

struct ST_GeometryType {

	//------------------------------------------------------------------------------------------------------------------
	// Binding
	//------------------------------------------------------------------------------------------------------------------
	// This function is a bit botched, but we cant change it without breaking backwards compatability
	// therefore, we use these constants for the geometry type values, instead of the normal type enum

	static constexpr uint8_t LEGACY_POINT_TYPE = 0;
	static constexpr uint8_t LEGACY_LINESTRING_TYPE = 1;
	static constexpr uint8_t LEGACY_POLYGON_TYPE = 2;
	static constexpr uint8_t LEGACY_MULTIPOINT_TYPE = 3;
	static constexpr uint8_t LEGACY_MULTILINESTRING_TYPE = 4;
	static constexpr uint8_t LEGACY_MULTIPOLYGON_TYPE = 5;
	static constexpr uint8_t LEGACY_GEOMETRYCOLLECTION_TYPE = 6;
	static constexpr uint8_t LEGACY_UNKNOWN_TYPE = 7;

	static unique_ptr<FunctionData> Bind(ClientContext &context, ScalarFunction &bound_function,
	                                     vector<unique_ptr<Expression>> &arguments) {
		// Create an enum type for all geometry types
		// Ensure that these are in the same order as the GeometryType enum
		const vector<string> enum_values = {"POINT", "LINESTRING", "POLYGON", "MULTIPOINT", "MULTILINESTRING",
		                                    "MULTIPOLYGON", "GEOMETRYCOLLECTION",
		                                    // or...
		                                    "UNKNOWN"};

		bound_function.return_type = GeoTypes::CreateEnumType("GEOMETRY_TYPE", enum_values);
		return nullptr;
	}

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, uint8_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			// TODO: Peek dont deserialize

			const auto geom = lstate.Deserialize(blob);
			switch (geom.get_type()) {
			case sgl::geometry_type::POINT:
				return LEGACY_POINT_TYPE;
			case sgl::geometry_type::LINESTRING:
				return LEGACY_LINESTRING_TYPE;
			case sgl::geometry_type::POLYGON:
				return LEGACY_POLYGON_TYPE;
			case sgl::geometry_type::MULTI_POINT:
				return LEGACY_MULTIPOINT_TYPE;
			case sgl::geometry_type::MULTI_LINESTRING:
				return LEGACY_MULTILINESTRING_TYPE;
			case sgl::geometry_type::MULTI_POLYGON:
				return LEGACY_MULTIPOLYGON_TYPE;
			case sgl::geometry_type::MULTI_GEOMETRY:
				return LEGACY_GEOMETRYCOLLECTION_TYPE;
			default:
				return LEGACY_UNKNOWN_TYPE;
			}
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePoint(DataChunk &args, ExpressionState &state, Vector &result) {
		result.SetVectorType(VectorType::CONSTANT_VECTOR);
		*ConstantVector::GetData<uint8_t>(result) = LEGACY_POINT_TYPE;
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		result.SetVectorType(VectorType::CONSTANT_VECTOR);
		*ConstantVector::GetData<uint8_t>(result) = LEGACY_LINESTRING_TYPE;
	}

	//------------------------------------------------------------------------------------------------------------------
	// POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePolygon(DataChunk &args, ExpressionState &state, Vector &result) {
		result.SetVectorType(VectorType::CONSTANT_VECTOR);
		*ConstantVector::GetData<uint8_t>(result) = LEGACY_POLYGON_TYPE;
	}

	//------------------------------------------------------------------------------------------------------------------
	// WKB
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteWKB(DataChunk &args, ExpressionState &state, Vector &result) {

		UnaryExecutor::Execute<string_t, uint8_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			BinaryReader cursor(blob.GetData(), blob.GetSize());

			const auto le = cursor.Read<uint8_t>();
			const auto type = le ? cursor.Read<uint32_t>() : cursor.ReadBE<uint32_t>();
			const auto normalized_type = (type & 0xffff) % 1000;

			if (normalized_type == 0 || normalized_type > 7) {
				return LEGACY_UNKNOWN_TYPE;
			}

			// Return the geometry type
			// Subtract 1 since the WKB type is 1-indexed
			return static_cast<uint8_t>(normalized_type - 1);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
	Returns a 'GEOMETRY_TYPE' enum identifying the input geometry type. Possible enum return types are: `POINT`, `LINESTRING`, `POLYGON`, `MULTIPOINT`, `MULTILINESTRING`, `MULTIPOLYGON`, and `GEOMETRYCOLLECTION`.
	)";

	static constexpr auto EXAMPLE = R"(
	SELECT DISTINCT ST_GeometryType(ST_GeomFromText('POINT(1 1)'));
	----
	POINT
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_GeometryType", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalTypeId::ANY);

				variant.SetBind(Bind);
				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalTypeId::ANY);

				variant.SetBind(Bind);
				variant.SetFunction(ExecutePoint);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(LogicalTypeId::ANY);

				variant.SetBind(Bind);
				variant.SetFunction(ExecuteLineString);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("polygon", GeoTypes::POLYGON_2D());
				variant.SetReturnType(LogicalTypeId::ANY);

				variant.SetBind(Bind);
				variant.SetFunction(ExecutePolygon);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("wkb", GeoTypes::WKB_BLOB());
				variant.SetReturnType(LogicalTypeId::ANY);

				variant.SetBind(Bind);
				variant.SetFunction(ExecuteWKB);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_GeomFromHEXWKB
//======================================================================================================================

struct ST_GeomFromHEXWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_GeomFromText
//----------------------------------------------------------------------------------------------------------------------
struct ST_GeomFromText {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_GeomFromWKB
//----------------------------------------------------------------------------------------------------------------------
struct ST_GeomFromWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//======================================================================================================================
// ST_HasZ
//======================================================================================================================

struct ST_HasZ {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [&](const string_t &blob) {
			// TODO: Peek without deserializing!
			const auto geom = lstate.Deserialize(blob);
			return geom.has_z();
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// WKB
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteWKB(DataChunk &args, ExpressionState &state, Vector &result) {
		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [](const string_t &wkb) {
			BinaryReader cursor(wkb.GetData(), wkb.GetSize());

			const auto le = cursor.Read<uint8_t>();
			const auto type = le ? cursor.Read<uint32_t>() : cursor.ReadBE<uint32_t>();

			// Check for ISO WKB and EWKB Z flag;
			const auto flags = (type & 0xffff) / 1000;
			return flags == 1 || flags == 3 || ((type & 0x80000000) != 0);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = "Check if the input geometry has Z values.";

	static constexpr auto EXAMPLE = R"(
	-- HasZ for a 2D geometry
	SELECT ST_HasZ(ST_GeomFromText('POINT(1 1)'));
	----
	false

	-- HasZ for a 3DZ geometry
	SELECT ST_HasZ(ST_GeomFromText('POINT Z(1 1 1)'));
	----
	true

	-- HasZ for a 3DM geometry
	SELECT ST_HasZ(ST_GeomFromText('POINT M(1 1 1)'));
	----
	false

	-- HasZ for a 4D geometry
	SELECT ST_HasZ(ST_GeomFromText('POINT ZM(1 1 1 1)'));
	----
	true
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_HasZ", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("wkb", GeoTypes::WKB_BLOB());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetFunction(ExecuteWKB);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_HasM
//======================================================================================================================

struct ST_HasM {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [&](const string_t &blob) {
			// TODO: Peek without deserializing!
			const auto geom = lstate.Deserialize(blob);
			return geom.has_m();
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// WKB_BLOB
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteWKB(DataChunk &args, ExpressionState &state, Vector &result) {
		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [](const string_t &wkb) {
			BinaryReader cursor(wkb.GetData(), wkb.GetSize());

			const auto le = cursor.Read<uint8_t>();
			const auto type = le ? cursor.Read<uint32_t>() : cursor.ReadBE<uint32_t>();

			// Check for ISO WKB and EWKB M flag;
			const auto flags = (type & 0xffff) / 1000;
			return flags == 2 || flags == 3 || ((type & 0x40000000) != 0);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = "Check if the input geometry has M values.";

	static constexpr auto EXAMPLE = R"(
	-- HasM for a 2D geometry
	SELECT ST_HasM(ST_GeomFromText('POINT(1 1)'));
	----
	false

	-- HasM for a 3DZ geometry
	SELECT ST_HasM(ST_GeomFromText('POINT Z(1 1 1)'));
	----
	false

	-- HasM for a 3DM geometry
	SELECT ST_HasM(ST_GeomFromText('POINT M(1 1 1)'));
	----
	true

	-- HasM for a 4D geometry
	SELECT ST_HasM(ST_GeomFromText('POINT ZM(1 1 1 1)'));
	----
	true
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_HasM", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("wkb", GeoTypes::WKB_BLOB());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetFunction(ExecuteWKB);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_ZMFlag
//======================================================================================================================

struct ST_ZMFlag {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			const auto has_z = geom.has_z();
			const auto has_m = geom.has_m();

			if (has_z && has_m) {
				return 3;
			}
			if (has_z) {
				return 2;
			}
			if (has_m) {
				return 1;
			}
			return 0;
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// WKB
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteWKB(DataChunk &args, ExpressionState &state, Vector &result) {
		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(), [](const string_t &wkb) {
			BinaryReader cursor(wkb.GetData(), wkb.GetSize());

			const auto le = cursor.Read<uint8_t>();
			const auto type = le ? cursor.Read<uint32_t>() : cursor.ReadBE<uint32_t>();

			// Check for ISO WKB and EWKB Z and M flags
			const uint32_t iso_wkb_props = (type & 0xffff) / 1000;
			const auto has_z = (iso_wkb_props == 1) || (iso_wkb_props == 3) || ((type & 0x80000000) != 0);
			const auto has_m = (iso_wkb_props == 2) || (iso_wkb_props == 3) || ((type & 0x40000000) != 0);

			if (has_z && has_m) {
				return 3;
			}
			if (has_z) {
				return 2;
			}
			if (has_m) {
				return 1;
			}
			return 0;
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
	Returns a flag indicating the presence of Z and M values in the input geometry.
	0 = No Z or M values
	1 = M values only
	2 = Z values only
	3 = Z and M values
	)";

	static constexpr auto EXAMPLE = R"(
	-- ZMFlag for a 2D geometry
	SELECT ST_ZMFlag(ST_GeomFromText('POINT(1 1)'));
	----
	0

	-- ZMFlag for a 3DZ geometry
	SELECT ST_ZMFlag(ST_GeomFromText('POINT Z(1 1 1)'));
	----
	2

	-- ZMFlag for a 3DM geometry
	SELECT ST_ZMFlag(ST_GeomFromText('POINT M(1 1 1)'));
	----
	1

	-- ZMFlag for a 4D geometry
	SELECT ST_ZMFlag(ST_GeomFromText('POINT ZM(1 1 1 1)'));
	----
	3
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_ZMFlag", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::INTEGER);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("wkb", GeoTypes::WKB_BLOB());
				variant.SetReturnType(LogicalType::INTEGER);

				variant.SetFunction(ExecuteWKB);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Distance_Sphere
//======================================================================================================================

struct ST_Distance_Sphere {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, string_t, double>(
		    args.data[0], args.data[1], result, args.size(), [&](const string_t &l_blob, const string_t &r_blob) {
			    const auto lhs = lstate.Deserialize(l_blob);
			    const auto rhs = lstate.Deserialize(r_blob);

			    if (lhs.get_type() != sgl::geometry_type::POINT || rhs.get_type() != sgl::geometry_type::POINT) {
				    throw InvalidInputException("ST_Distance_Sphere only accepts POINT geometries");
			    }

			    if (lhs.is_empty() || rhs.is_empty()) {
				    throw InvalidInputException("ST_Distance_Sphere does not accept empty geometries");
			    }

			    const auto lv = lhs.get_vertex_xy(0);
			    const auto rv = rhs.get_vertex_xy(0);

			    return sgl::util::haversine_distance(lv.x, lv.y, rv.x, rv.y);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePoint(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 2);
		auto &left = args.data[0];
		auto &right = args.data[1];
		auto count = args.size();

		using POINT_TYPE = StructTypeBinary<double, double>;
		using DISTANCE_TYPE = PrimitiveType<double>;

		GenericExecutor::ExecuteBinary<POINT_TYPE, POINT_TYPE, DISTANCE_TYPE>(
		    left, right, result, count, [&](POINT_TYPE left, POINT_TYPE right) {
			    return sgl::util::haversine_distance(left.a_val, left.b_val, right.a_val, right.b_val);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Returns the haversine (great circle) distance between two geometries.

	    - Only supports POINT geometries.
	    - Returns the distance in meters.
	    - The input is expected to be in WGS84 (EPSG:4326) coordinates, using a [latitude, longitude] axis order.
	)";

	// TODO: Example
	static constexpr auto EXAMPLE = R"";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Distance_Sphere", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point1", GeoTypes::POINT_2D());
				variant.AddParameter("point2", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(ExecutePoint);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Hilbert
//======================================================================================================================
struct ST_Hilbert {

	//------------------------------------------------------------------------------------------------------------------
	// Hilbert Curve Encoding
	// From (Public Domain): https://github.com/rawrunprotected/hilbert_curves
	//------------------------------------------------------------------------------------------------------------------
	static uint32_t Interleave(uint32_t x) {
		x = (x | (x << 8)) & 0x00FF00FF;
		x = (x | (x << 4)) & 0x0F0F0F0F;
		x = (x | (x << 2)) & 0x33333333;
		x = (x | (x << 1)) & 0x55555555;
		return x;
	}

	static uint32_t HilbertEncode(uint32_t n, uint32_t x, uint32_t y) {
		x = x << (16 - n);
		y = y << (16 - n);

		// Initial prefix scan round, prime with x and y
		uint32_t a = x ^ y;
		uint32_t b = 0xFFFF ^ a;
		uint32_t c = 0xFFFF ^ (x | y);
		uint32_t d = x & (y ^ 0xFFFF);
		uint32_t A = a | (b >> 1);
		uint32_t B = (a >> 1) ^ a;
		uint32_t C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
		uint32_t D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;

		a = A;
		b = B;
		c = C;
		d = D;
		A = ((a & (a >> 2)) ^ (b & (b >> 2)));
		B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));
		C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
		D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));

		a = A;
		b = B;
		c = C;
		d = D;
		A = ((a & (a >> 4)) ^ (b & (b >> 4)));
		B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));
		C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
		D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));

		// Final round and projection
		a = A;
		b = B;
		c = C;
		d = D;
		C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
		D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));

		// Undo transformation prefix scan
		a = C ^ (C >> 1);
		b = D ^ (D >> 1);

		// Recover index bits
		uint32_t i0 = x ^ y;
		uint32_t i1 = b | (0xFFFF ^ (i0 | a));

		return ((Interleave(i1) << 1) | Interleave(i0)) >> (32 - 2 * n);
	}

	static uint32_t FloatToUint32(float f) {
		if (std::isnan(f)) {
			return 0xFFFFFFFF;
		}
		uint32_t res;
		memcpy(&res, &f, sizeof(res));
		if ((res & 0x80000000) != 0) {
			res ^= 0xFFFFFFFF;
		} else {
			res |= 0x80000000;
		}
		return res;
	}

	//------------------------------------------------------------------------------------------------------------------
	// BOX_2D / BOX_2F
	//------------------------------------------------------------------------------------------------------------------
	template <class T>
	static void ExecuteBox(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &input_vec = args.data[0];
		auto &bounds_vec = args.data[1];
		auto count = args.size();

		constexpr auto max_hilbert = std::numeric_limits<uint16_t>::max();

		using BOX_TYPE = StructTypeQuaternary<T, T, T, T>;
		using UINT32_TYPE = PrimitiveType<uint32_t>;

		GenericExecutor::ExecuteBinary<BOX_TYPE, BOX_TYPE, UINT32_TYPE>(
		    input_vec, bounds_vec, result, count, [&](BOX_TYPE &box, BOX_TYPE &bounds) {
			    const auto x = box.a_val + (box.c_val - box.a_val) / static_cast<T>(2);
			    const auto y = box.b_val + (box.d_val - box.b_val) / static_cast<T>(2);

			    const auto hilbert_width = max_hilbert / (bounds.c_val - bounds.a_val);
			    const auto hilbert_height = max_hilbert / (bounds.d_val - bounds.b_val);

			    // TODO: Check for overflow
			    const auto hilbert_x = static_cast<uint32_t>((x - bounds.a_val) * hilbert_width);
			    const auto hilbert_y = static_cast<uint32_t>((y - bounds.b_val) * hilbert_height);
			    const auto h = HilbertEncode(16, hilbert_x, hilbert_y);
			    return UINT32_TYPE {h};
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// LON/LAT
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLonlat(DataChunk &args, ExpressionState &state, Vector &result) {
		using DOUBLE_TYPE = PrimitiveType<double>;
		using UINT32_TYPE = PrimitiveType<uint32_t>;
		using BOX_TYPE = StructTypeQuaternary<double, double, double, double>;

		auto constexpr max_hilbert = std::numeric_limits<uint16_t>::max();

		GenericExecutor::ExecuteTernary<DOUBLE_TYPE, DOUBLE_TYPE, BOX_TYPE, UINT32_TYPE>(
		    args.data[0], args.data[1], args.data[3], result, args.size(),
		    [&](DOUBLE_TYPE x, DOUBLE_TYPE y, BOX_TYPE &box) {
			    const auto hilbert_width = max_hilbert / (box.c_val - box.a_val);
			    const auto hilbert_height = max_hilbert / (box.d_val - box.b_val);

			    // TODO: Check for overflow
			    const auto hilbert_x = static_cast<uint32_t>((x.val - box.a_val) * hilbert_width);
			    const auto hilbert_y = static_cast<uint32_t>((y.val - box.b_val) * hilbert_height);
			    const auto h = HilbertEncode(16, hilbert_x, hilbert_y);
			    return UINT32_TYPE {h};
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		UnaryExecutor::ExecuteWithNulls<geometry_t, uint32_t>(
		    args.data[0], result, args.size(),
		    [&](const geometry_t &geom, ValidityMask &mask, idx_t out_idx) -> uint32_t {
			    // TODO: This is shit, dont rely on cached bounds
			    Box2D<double> bounds;
			    if (!geom.TryGetCachedBounds(bounds)) {
				    mask.SetInvalid(out_idx);
				    return 0;
			    }

			    Box2D<float> bounds_f;
			    bounds_f.min.x = MathUtil::DoubleToFloatDown(bounds.min.x);
			    bounds_f.min.y = MathUtil::DoubleToFloatDown(bounds.min.y);
			    bounds_f.max.x = MathUtil::DoubleToFloatUp(bounds.max.x);
			    bounds_f.max.y = MathUtil::DoubleToFloatUp(bounds.max.y);

			    const auto dx = bounds_f.min.x + (bounds_f.max.x - bounds_f.min.x) / 2;
			    const auto dy = bounds_f.min.y + (bounds_f.max.y - bounds_f.min.y) / 2;

			    const auto hx = FloatToUint32(dx);
			    const auto hy = FloatToUint32(dy);

			    return HilbertEncode(16, hx, hy);
		    });
	}

	static void ExecuteGeometryWithBounds(DataChunk &args, ExpressionState &state, Vector &result) {

		auto constexpr max_hilbert = std::numeric_limits<uint16_t>::max();

		using BOX_TYPE = StructTypeQuaternary<double, double, double, double>;
		using GEOM_TYPE = PrimitiveType<geometry_t>;
		using UINT32_TYPE = PrimitiveType<uint32_t>;

		GenericExecutor::ExecuteBinary<GEOM_TYPE, BOX_TYPE, UINT32_TYPE>(
		    args.data[0], args.data[1], result, args.size(), [&](const GEOM_TYPE &geom_type, const BOX_TYPE &bounds) {
			    const auto geom = geom_type.val;

			    // TODO: This is shit, dont rely on cached bounds
			    Box2D<double> geom_bounds;
			    if (!geom.TryGetCachedBounds(geom_bounds)) {
				    throw InvalidInputException(
				        "ST_Hilbert(geom, bounds) requires that all geometries have a bounding box");
			    }

			    const auto dx = geom_bounds.min.x + (geom_bounds.max.x - geom_bounds.min.x) / 2;
			    const auto dy = geom_bounds.min.y + (geom_bounds.max.y - geom_bounds.min.y) / 2;

			    const auto hilbert_width = max_hilbert / (bounds.c_val - bounds.a_val);
			    const auto hilbert_height = max_hilbert / (bounds.d_val - bounds.b_val);
			    // TODO: Check for overflow
			    const auto hilbert_x = static_cast<uint32_t>((dx - bounds.a_val) * hilbert_width);
			    const auto hilbert_y = static_cast<uint32_t>((dy - bounds.b_val) * hilbert_height);

			    const auto h = HilbertEncode(16, hilbert_x, hilbert_y);
			    return UINT32_TYPE {h};
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Encodes the X and Y values as the hilbert curve index for a curve covering the given bounding box.
		If a geometry is provided, the center of the approximate bounding box is used as the point to encode.
		If no bounding box is provided, the hilbert curve index is mapped to the full range of a single-presicion float.
		For the BOX_2D and BOX_2DF variants, the center of the box is used as the point to encode.
	)";

	// TODO: example
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		// TODO: All of these needs examples and docs

		FunctionBuilder::RegisterScalar(db, "ST_Hilbert", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("x", LogicalType::DOUBLE);
				variant.AddParameter("y", LogicalType::DOUBLE);
				variant.AddParameter("bounds", GeoTypes::BOX_2D());
				variant.SetReturnType(LogicalType::UINTEGER);

				variant.SetFunction(ExecuteLonlat);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("bounds", GeoTypes::BOX_2D());
				variant.SetReturnType(LogicalType::UINTEGER);

				variant.SetFunction(ExecuteGeometryWithBounds);
				variant.SetInit(LocalState::Init);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::UINTEGER);

				variant.SetFunction(ExecuteGeometry);
				variant.SetInit(LocalState::Init);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2D());
				variant.AddParameter("bounds", GeoTypes::BOX_2D());
				variant.SetReturnType(LogicalType::UINTEGER);

				variant.SetFunction(ExecuteBox<double>);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2DF());
				variant.AddParameter("bounds", GeoTypes::BOX_2DF());
				variant.SetReturnType(LogicalType::UINTEGER);

				variant.SetFunction(ExecuteBox<float>);
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);
		});
	}
};

//======================================================================================================================
// ST_Intersects
//======================================================================================================================
// TODO: Implement
struct ST_Intersects {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//======================================================================================================================
// ST_IntersectsExtent
//======================================================================================================================
// TODO: Implement
struct ST_IntersectsExtent {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//======================================================================================================================
// ST_IsClosed
//======================================================================================================================
struct ST_IsClosed {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			switch (geom.get_type()) {
			case sgl::geometry_type::LINESTRING:
				return sgl::linestring::is_closed(&geom);
			case sgl::geometry_type::MULTI_LINESTRING:
				return sgl::multi_linestring::is_closed(&geom);
			default:
				// TODO: We should support more than just LINESTRING and MULTILINESTRING (like PostGIS does)
				throw InvalidInputException("ST_IsClosed only accepts LINESTRING and MULTILINESTRING geometries");
			}
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_IsClosed", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Check if a geometry is closed");
				// TODO: Set example
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_IsEmpty
//======================================================================================================================

struct ST_IsEmpty {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			const auto vertex_count = sgl::ops::vertex_count(&geom);
			return vertex_count == 0;
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLinestring(DataChunk &args, ExpressionState &state, Vector &result) {
		UnaryExecutor::Execute<list_entry_t, bool>(args.data[0], result, args.size(),
		                                           [&](const list_entry_t &line) { return line.length == 0; });
	}

	//------------------------------------------------------------------------------------------------------------------
	// POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePolygon(DataChunk &args, ExpressionState &state, Vector &result) {
		UnaryExecutor::Execute<list_entry_t, bool>(args.data[0], result, args.size(),
		                                           [&](const list_entry_t &poly) { return poly.length == 0; });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Returns true if the geometry is "empty".
	)";
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_IsEmpty", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetFunction(ExecuteLinestring);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("polygon", GeoTypes::POLYGON_2D());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetFunction(ExecutePolygon);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Length
//======================================================================================================================

struct ST_Length {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			return sgl::ops::length(&geom);
		});
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLinestring(DataChunk &args, ExpressionState &state, Vector &result) {
		D_ASSERT(args.data.size() == 1);

		auto &line_vec = args.data[0];
		auto count = args.size();

		auto &coord_vec = ListVector::GetEntry(line_vec);
		auto &coord_vec_children = StructVector::GetEntries(coord_vec);
		auto x_data = FlatVector::GetData<double>(*coord_vec_children[0]);
		auto y_data = FlatVector::GetData<double>(*coord_vec_children[1]);

		UnaryExecutor::Execute<list_entry_t, double>(line_vec, result, count, [&](const list_entry_t &line) {
			auto offset = line.offset;
			auto length = line.length;
			double sum = 0;
			// Loop over the segments
			for (idx_t j = offset; j < offset + length - 1; j++) {
				auto x1 = x_data[j];
				auto y1 = y_data[j];
				auto x2 = x_data[j + 1];
				auto y2 = y_data[j + 1];
				sum += std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
			}
			return sum;
		});

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Returns the length of the input line geometry
	)";

	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Length", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("linestring", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(ExecuteLinestring);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_MakeEnvelope
//----------------------------------------------------------------------------------------------------------------------
struct ST_MakeEnvelope {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		auto &min_x_vec = args.data[0];
		auto &min_y_vec = args.data[1];
		auto &max_x_vec = args.data[2];
		auto &max_y_vec = args.data[3];

		using DOUBLE_TYPE = PrimitiveType<double>;
		using STRING_TYPE = PrimitiveType<string_t>;

		GenericExecutor::ExecuteQuaternary<DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, STRING_TYPE>(
		    min_x_vec, min_y_vec, max_x_vec, max_y_vec, result, args.size(),
		    [&](const DOUBLE_TYPE vmin_x, const DOUBLE_TYPE vmin_y, const DOUBLE_TYPE vmax_x,
		        const DOUBLE_TYPE vmax_y) {
			    const auto min_x = vmin_x.val;
			    const auto min_y = vmin_y.val;
			    const auto max_x = vmax_x.val;
			    const auto max_y = vmax_y.val;

			    // This is pretty cool, we dont even need to allocate anything
			    const double buffer[10] = {min_x, min_y, min_x, max_y, max_x, max_y, max_x, min_y, min_x, min_y};

			    auto ring = sgl::linestring::make_empty(false, false);
			    ring.set_vertex_data(reinterpret_cast<const char *>(buffer), 5);

			    auto poly = sgl::polygon::make_empty();
			    poly.append_part(&ring);

			    return lstate.Serialize(result, poly);
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_MakeEnvelope", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("min_x", LogicalType::DOUBLE);
				variant.AddParameter("min_y", LogicalType::DOUBLE);
				variant.AddParameter("max_x", LogicalType::DOUBLE);
				variant.AddParameter("max_y", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Create a rectangular polygon from min/max coordinates");
				// todo: example
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_MakeLine
//----------------------------------------------------------------------------------------------------------------------
struct ST_MakeLine {
	static void ExecuteList(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		auto &child_vec = ListVector::GetEntry(args.data[0]);
		auto child_len = ListVector::GetListSize(args.data[0]);

		UnifiedVectorFormat format;
		child_vec.ToUnifiedFormat(child_len, format);

		UnaryExecutor::Execute<list_entry_t, string_t>(
		    args.data[0], result, args.size(), [&](const list_entry_t &entry) {
			    const auto offset = entry.offset;
			    const auto length = entry.length;

			    if (length < 2) {
				    // Early out if we don't have enough geometries
				    throw InvalidInputException("ST_MakeLine requires at least 2 geometries");
			    }

			    uint32_t line_length = 0;
			    // First pass, filter types, count non-null entries

			    for (idx_t i = offset; i < offset + length; i++) {
				    const auto mapped_idx = format.sel->get_index(i);
				    if (format.validity.RowIsValid(mapped_idx)) {
					    continue;
				    }
				    auto &blob = UnifiedVectorFormat::GetData<string_t>(format)[mapped_idx];

				    // TODO: Peek without deserializing
				    const auto geom = lstate.Deserialize(blob);
				    if (geom.get_type() != sgl::geometry_type::POINT) {
					    throw InvalidInputException("ST_MakeLine only accepts POINT geometries");
				    }

				    // TODO: Support Z and M
				    if (geom.has_z() || geom.has_m()) {
					    throw InvalidInputException(
					        "ST_MakeLine from list does not accept POINT geometries with Z or M values");
				    }

				    line_length++;
			    }

			    if (line_length < 2) {
				    throw InvalidInputException("ST_MakeLine requires at least 2 non-null geometries");
			    }

			    const auto line_data = lstate.GetArena().AllocateAligned(line_length * 2 * sizeof(double));

			    // Second pass, copy over the vertex data
			    uint32_t vertex_idx = 0;
			    for (idx_t i = offset; i < offset + length; i++) {
				    D_ASSERT(vertex_idx < line_length);

				    const auto mapped_idx = format.sel->get_index(i);
				    if (format.validity.RowIsValid(mapped_idx)) {
					    continue;
				    }
				    auto &blob = UnifiedVectorFormat::GetData<string_t>(format)[mapped_idx];

				    const auto point = lstate.Deserialize(blob);
				    const auto point_data = point.get_vertex_data();

				    memcpy(line_data + vertex_idx * 2 * sizeof(double), point_data, 2 * sizeof(double));
				    vertex_idx++;
			    }

			    D_ASSERT(vertex_idx == line_length);

			    auto line = sgl::linestring::make_empty(false, false);
			    line.set_vertex_data(line_data, line_length);

			    return lstate.Serialize(result, line);
		    });
	}

	static void ExecuteBinary(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, string_t, string_t>(
		    args.data[0], args.data[1], result, args.size(), [&](const string_t &l_blob, const string_t &r_blob) {
			    const auto l_geom = lstate.Deserialize(l_blob);
			    const auto r_geom = lstate.Deserialize(r_blob);

			    if (l_geom.get_type() != sgl::geometry_type::POINT || r_geom.get_type() != sgl::geometry_type::POINT) {
				    throw InvalidInputException("ST_MakeLine only accepts POINT geometries");
			    }

			    if (l_geom.is_empty() || r_geom.is_empty()) {
				    throw InvalidInputException("ST_MakeLine does not accept empty POINT geometries");
			    }

			    const auto has_z = l_geom.has_z() || r_geom.has_z();
			    const auto has_m = l_geom.has_m() || r_geom.has_m();

			    auto linestring = sgl::linestring::make_empty(has_z, has_m);

			    // Create a buffer large enough to store two vertices
			    double buffer[8] = {0};

			    const auto v1 = l_geom.get_vertex_xyzm(0);
			    const auto v2 = r_geom.get_vertex_xyzm(0);

			    // TODO: this is a bit ugly, add proper append method to sgl instead
			    idx_t idx = 0;
			    buffer[idx++] = v1.x;
			    buffer[idx++] = v1.y;
			    if (has_z) {
				    buffer[idx++] = l_geom.has_z() ? v1.zm : 0;
			    }
			    if (has_m) {
				    buffer[idx++] = l_geom.has_m() ? l_geom.has_z() ? v1.m : v1.zm : 0;
			    }
			    buffer[idx++] = v2.x;
			    buffer[idx++] = v2.y;
			    if (has_z) {
				    buffer[idx++] = r_geom.has_z() ? v2.zm : 0;
			    }
			    if (has_m) {
				    buffer[idx++] = r_geom.has_m() ? r_geom.has_z() ? v2.m : v2.zm : 0;
			    }

			    linestring.set_vertex_data(reinterpret_cast<const char *>(buffer), 2);

			    return lstate.Serialize(result, linestring);
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_MakeLine", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geoms", LogicalType::LIST(GeoTypes::GEOMETRY()));
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteList);

				variant.SetDescription("Create a LINESTRING from a list of POINT geometries");
				variant.SetExample("SELECT ST_MakeLine([ST_Point(0, 0), ST_Point(1, 1)]);");
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("start", GeoTypes::GEOMETRY());
				variant.AddParameter("end", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteBinary);

				variant.SetDescription("Create a LINESTRING from two POINT geometries");
				variant.SetExample("SELECT ST_MakeLine(ST_Point(0, 0), ST_Point(1, 1));");
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_MakePolygon
//----------------------------------------------------------------------------------------------------------------------
struct ST_MakePolygon {
	static void ExecuteFromShell(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			auto line = lstate.Deserialize(blob);

			if (line.get_type() != sgl::geometry_type::LINESTRING) {
				throw InvalidInputException("ST_MakePolygon only accepts LINESTRING geometries");
			}

			if (line.get_count() < 4) {
				throw InvalidInputException("ST_MakePolygon shell requires at least 4 vertices");
			}

			if (!sgl::linestring::is_closed(&line)) {
				throw std::runtime_error("ST_MakePolygon shell must be closed (first and last vertex must be equal)");
			}

			auto polygon = sgl::polygon::make_empty(line.has_z(), line.has_m());
			polygon.append_part(&line);

			return lstate.Serialize(result, polygon);
		});
	}

	static void ExecuteFromRings(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		auto &child_vec = ListVector::GetEntry(args.data[0]);
		auto child_len = ListVector::GetListSize(args.data[0]);

		UnifiedVectorFormat child_format;
		child_vec.ToUnifiedFormat(child_len, child_format);

		BinaryExecutor::Execute<string_t, list_entry_t, string_t>(
		    args.data[0], args.data[1], result, args.size(), [&](const string_t &blob, const list_entry_t &hole_list) {
			    // First, setup shell
			    auto shell = lstate.Deserialize(blob);
			    if (shell.get_type() != sgl::geometry_type::LINESTRING) {
				    throw InvalidInputException("ST_MakePolygon only accepts LINESTRING geometries");
			    }
			    // TODO: Support Z and M
			    if (shell.has_z() || shell.has_m()) {
				    throw InvalidInputException("ST_MakePolygon from list does not support Z or M values");
			    }
			    if (shell.get_count() < 4) {
				    throw InvalidInputException("ST_MakePolygon shell requires at least 4 vertices");
			    }
			    if (!sgl::linestring::is_closed(&shell)) {
				    throw InvalidInputException(
				        "ST_MakePolygon shell must be closed (first and last vertex must be equal)");
			    }

			    // Make a polygon!
			    auto polygon = sgl::polygon::make_empty(false, false);

			    // Append the shell
			    polygon.append_part(&shell);

			    // Now setup the rings
			    const auto holes_offset = hole_list.offset;
			    const auto holes_length = hole_list.length;

			    for (idx_t hole_idx = 0; hole_idx < holes_length; hole_idx++) {
				    const auto mapped_idx = child_format.sel->get_index(holes_offset + hole_idx);
				    if (child_format.validity.RowIsValid(mapped_idx)) {
					    continue;
				    }

				    const auto &hole_blob = UnifiedVectorFormat::GetData<string_t>(child_format)[mapped_idx];

				    // Allocate a new hole and deserialize into the memory
				    auto hole_mem = lstate.GetArena().AllocateAligned(sizeof(sgl::geometry));
				    const auto hole = new (hole_mem) sgl::geometry();

				    // TODO: Make this nicer... Add a deserialize in place method to the context
				    *hole = lstate.Deserialize(hole_blob);

				    if (hole->get_type() != sgl::geometry_type::LINESTRING) {
					    throw InvalidInputException("ST_MakePolygon hole #%lu is not a LINESTRING", hole_idx + 1);
				    }
				    if (hole->has_z() || hole->has_m()) {
					    throw InvalidInputException("ST_MakePolygon hole #%lu has Z or M values", hole_idx + 1);
				    }
				    if (hole->get_count() < 4) {
					    throw InvalidInputException("ST_MakePolygon hole requires at least 4 vertices");
				    }
				    if (!sgl::linestring::is_closed(hole)) {
					    throw InvalidInputException(
					        "ST_MakePolygon hole #%lu must be closed (first and last vertex must be equal)",
					        hole_idx + 1);
				    }

				    // Add the hole to the polygon
				    polygon.append_part(hole);
			    }

			    // Now serialize the polygon
			    return lstate.Serialize(result, polygon);
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_MakePolygon", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("shell", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteFromShell);

				// TODO: Set example & docs
				variant.SetDescription("Create a POLYGON from a LINESTRING shell");
				variant.SetExample("SELECT ST_MakePolygon(ST_LineString([ST_Point(0, 0), ST_Point(1, 0), ST_Point(1, "
				                   "1), ST_Point(0, 0)]));");
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("shell", GeoTypes::GEOMETRY());
				variant.AddParameter("holes", LogicalType::LIST(GeoTypes::GEOMETRY()));
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteFromRings);

				// TODO: Set example & docs
				variant.SetDescription("Create a POLYGON from a LINESTRING shell and a list of LINESTRING holes");
				variant.SetExample("SELECT ST_MakePolygon(ST_LineString([ST_Point(0, 0), ST_Point(1, 0), ST_Point(1, "
				                   "1), ST_Point(0, 0)]), [ST_LineString([ST_Point(0.25, 0.25), ST_Point(0.75, 0.25), "
				                   "ST_Point(0.75, 0.75), ST_Point(0.25, 0.25)])]);");
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_Multi
//----------------------------------------------------------------------------------------------------------------------
// TODO: Implement
struct ST_Multi {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			auto geom = lstate.Deserialize(blob);
			const auto has_z = geom.has_z();
			const auto has_m = geom.has_m();

			switch (geom.get_type()) {
			case sgl::geometry_type::POINT: {
				auto mpoint = sgl::multi_point::make_empty(has_z, has_m);
				mpoint.append_part(&geom);
				return lstate.Serialize(result, mpoint);
			}
			case sgl::geometry_type::LINESTRING: {
				auto mline = sgl::multi_linestring::make_empty(has_z, has_m);
				mline.append_part(&geom);
				return lstate.Serialize(result, mline);
			}
			case sgl::geometry_type::POLYGON: {
				auto mpoly = sgl::multi_polygon::make_empty(has_z, has_m);
				mpoly.append_part(&geom);
				return lstate.Serialize(result, mpoly);
			}
			default:
				// Just return the original geometry
				return blob;
			}
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Multi", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Convert a geometry to a MULTI* geometry");
				// TODO: Set example
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_NGeometries / ST_NumGeometries
//----------------------------------------------------------------------------------------------------------------------
// TODO: Implement
struct ST_NGeometries {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			if (geom.is_single_part()) {
				return geom.is_empty() ? 0 : 1;
			}
			if (geom.is_multi_part()) {
				return static_cast<int32_t>(geom.get_count());
			}
			return 0;
		});
	}

	static void Register(DatabaseInstance &db) {
		// TODO: Maybe make a macro for the aliases
		for (auto &alias : {"ST_NumGeometries", "ST_NGeometries"}) {
			FunctionBuilder::RegisterScalar(db, alias, [](ScalarFunctionBuilder &func) {
				func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
					variant.AddParameter("geom", GeoTypes::GEOMETRY());
					variant.SetReturnType(LogicalType::INTEGER);

					variant.SetInit(LocalState::Init);
					variant.SetFunction(Execute);

					variant.SetDescription("Returns the number of geometries in a geometry collection");
					// todo: Set example
				});
				func.SetTag("ext", "spatial");
				func.SetTag("category", "property");
			});
		}
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_NInteriorRings / ST_NumInteriorRings
//----------------------------------------------------------------------------------------------------------------------
// TODO: Implement
struct ST_NInteriorRings {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::ExecuteWithNulls<string_t, int32_t>(
		    args.data[0], result, args.size(), [&](const string_t &blob, ValidityMask &validity, idx_t idx) {
			    const auto geom = lstate.Deserialize(blob);

			    if (geom.get_type() != sgl::geometry_type::POLYGON) {
				    validity.SetInvalid(idx);
				    return 0;
			    }

			    const auto n_rings = static_cast<int32_t>(geom.get_count());
			    return n_rings == 0 ? 0 : n_rings - 1;
		    });
	}

	static void Register(DatabaseInstance &db) {
		// todo: maybe make a macro for the aliases
		for (auto &alias : {"ST_NumInteriorRings", "ST_NInteriorRings"}) {
			FunctionBuilder::RegisterScalar(db, alias, [](ScalarFunctionBuilder &func) {
				func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
					variant.AddParameter("geom", GeoTypes::GEOMETRY());
					variant.SetReturnType(LogicalType::INTEGER);

					variant.SetInit(LocalState::Init);
					variant.SetFunction(Execute);

					variant.SetDescription("Returns the number of interior rings in a polygon");
					// TODO: Set example
				});
				func.SetTag("ext", "spatial");
				func.SetTag("category", "property");
			});
		}
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_NPoints
//----------------------------------------------------------------------------------------------------------------------
// TODO: Implement
struct ST_NPoints {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			return sgl::ops::vertex_count(&geom);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_NPoints", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::INTEGER);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Returns the number of vertices in a geometry");
				// todo: Set example
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_Perimeter
//----------------------------------------------------------------------------------------------------------------------
// TODO: Implement
struct ST_Perimeter {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [&](const string_t &blob) {
			const auto geom = lstate.Deserialize(blob);
			return sgl::ops::perimeter(&geom);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Perimeter", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Compute the perimeter of a geometry");
				// TODO: Set example
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Point
//======================================================================================================================

struct ST_Point {

	//------------------------------------------------------------------------------------------------------------------
	// Execute
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<double, double, string_t>(
		    args.data[0], args.data[1], result, args.size(), [&](const double x, const double y) {
			    const double buffer[2] = {x, y};

			    sgl::geometry geometry;
			    geometry.set_type(sgl::geometry_type::POINT);
			    geometry.set_vertex_data(reinterpret_cast<const uint8_t *>(buffer), 1);

			    return lstate.Serialize(result, geometry);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Creates a point geometry from X and Y coordinates.
	)";

	// TODO: example
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Point", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("x", LogicalType::DOUBLE);
				variant.AddParameter("y", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetFunction(Execute);
				variant.SetInit(LocalState::Init);

				variant.SetDescription("Creates a GEOMETRY point");
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_PointN
//----------------------------------------------------------------------------------------------------------------------
struct ST_PointN {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		/*
		BinaryExecutor::ExecuteWithNulls<geometry_t, int32_t, geometry_t>(
		    args.data[0], args.data[1], result, args.size(),
		    [&](const string_t &blob, int32_t index, ValidityMask &mask, idx_t row_idx) {

		        // TODO: peek without deserialing
		        if (input.GetType() != GeometryType::LINESTRING) {
		            mask.SetInvalid(row_idx);
		            return string_t {};
		        }

		        const auto line = lstate.Deserialize(blob);
		        const auto point_count = line.get_count();;

		        if (point_count == 0 || index == 0 || index < -static_cast<int64_t>(point_count) ||
		            index > static_cast<int64_t>(point_count)) {
		            mask.SetInvalid(row_idx);
		            return geometry_t {};
		        }

		        auto actual_index = index < 0 ? point_count + index : index - 1;
		        auto point = LineString::GetPointAsReference(line, actual_index);
		        return lstate.Serialize(result, point);
		    });
		    */
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_Points
//----------------------------------------------------------------------------------------------------------------------
struct ST_Points {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &blob) {
			// Deserialize the geometry
			const auto geom = lstate.Deserialize(blob);
			const auto has_z = geom.has_z();
			const auto has_m = geom.has_m();

			// Create a new result multipoint
			auto mpoint = sgl::multi_point::make_empty(has_z, has_m);

			sgl::ops::visit_vertices(&geom, [&](const uint8_t *vertex_data) {
				// Allocate a new point
				auto point_mem = lstate.GetArena().AllocateAligned(sizeof(sgl::geometry));

				// Create a new point
				const auto point = new (point_mem) sgl::geometry(sgl::geometry_type::POINT, has_z, has_m);
				point->set_vertex_data(vertex_data, 1);

				// Append the point to the multipoint
				mpoint.append_part(point);
			});

			// Serialize the multipoint
			return lstate.Serialize(result, mpoint);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Points", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Returns a MULTIPOINT containing all the vertices of a geometry");
				// TODO: Set example
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_QuadKey
//----------------------------------------------------------------------------------------------------------------------
struct ST_QuadKey {

	static void GetQuadKey(double lon, double lat, int32_t level, char *buffer) {

		lat = std::max(-85.05112878, std::min(85.05112878, lat));
		lon = std::max(-180.0, std::min(180.0, lon));

		const auto lat_rad = lat * PI / 180.0;
		const auto x = static_cast<int32_t>((lon + 180.0) / 360.0 * (1 << level));
		const auto y = static_cast<int32_t>((1.0 - std::log(std::tan(lat_rad) + 1.0 / std::cos(lat_rad)) / PI) / 2.0 *
		                                    (1 << level));

		for (int i = level; i > 0; --i) {
			char digit = '0';
			const int32_t mask = 1 << (i - 1);
			if ((x & mask) != 0) {
				digit += 1;
			}
			if ((y & mask) != 0) {
				digit += 2;
			}
			buffer[level - i] = digit;
		}
	}
	//------------------------------------------------------------------------------------------------------------------
	// Geometry Function
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		auto &point_in = args.data[0];
		auto &level_in = args.data[1];

		BinaryExecutor::Execute<string_t, int32_t, string_t>(
		    point_in, level_in, result, args.size(), [&](const string_t &blob, const int32_t level) {
			    if (level < 1 || level > 23) {
				    throw InvalidInputException("ST_QuadKey: Level must be between 1 and 23");
			    }

			    const auto point = lstate.Deserialize(blob);
			    if (point.get_type() != sgl::geometry_type::POINT) {
				    throw InvalidInputException("ST_QuadKey: Only POINT geometries are supported");
			    }

			    if (point.is_empty()) {
				    throw InvalidInputException("ST_QuadKey: Empty geometries are not supported");
			    }

			    const auto vertex = point.get_vertex_xy(0);

			    char buffer[64];
			    GetQuadKey(vertex.x, vertex.y, level, buffer);
			    return StringVector::AddString(result, buffer, level);
		    });
	}

	// Point Function  -------------------------------------------------------------------------------------------------
	static void ExecuteLonLat(DataChunk &args, ExpressionState &state, Vector &result) {

		auto &lon_in = args.data[0];
		auto &lat_in = args.data[1];
		auto &lev_in = args.data[2];

		TernaryExecutor::Execute<double, double, int32_t, string_t>(
		    lon_in, lat_in, lev_in, result, args.size(), [&](const double lon, const double lat, const int32_t level) {
			    if (level < 1 || level > 23) {
				    throw InvalidInputException("ST_QuadKey: Level must be between 1 and 23");
			    }
			    char buffer[64];
			    GetQuadKey(lon, lat, level, buffer);
			    return StringVector::AddString(result, buffer, level);
		    });
	}

	// Register --------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_QuadKey", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("longitude", LogicalType::DOUBLE);
				variant.AddParameter("latitude", LogicalType::DOUBLE);
				variant.AddParameter("level", LogicalType::INTEGER);
				variant.SetReturnType(LogicalType::VARCHAR);
				variant.SetFunction(ExecuteLonLat);

				// TODO: Set example
				// variant.SetExample(DOC_EXAMPLE);
				// variant.SetDescription(DOC_DESCRIPTION);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::GEOMETRY());
				variant.AddParameter("level", LogicalType::INTEGER);
				variant.SetReturnType(LogicalType::VARCHAR);
				variant.SetFunction(ExecuteGeometry);
				variant.SetInit(LocalState::Init);

				// TODO: Set example
				// variant.SetExample(DOC_EXAMPLE);
				// variant.SetDescription(DOC_DESCRIPTION);
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_RemoveRepeatedPoints
//======================================================================================================================

struct ST_RemoveRepeatedPoints {

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto count = args.size();
		UnifiedVectorFormat format;
		input.ToUnifiedFormat(count, format);

		auto in_line_entries = ListVector::GetData(input);
		auto &in_line_vertex_vec = StructVector::GetEntries(ListVector::GetEntry(input));
		auto in_x_data = FlatVector::GetData<double>(*in_line_vertex_vec[0]);
		auto in_y_data = FlatVector::GetData<double>(*in_line_vertex_vec[1]);

		auto out_line_entries = ListVector::GetData(result);
		auto &out_line_vertex_vec = StructVector::GetEntries(ListVector::GetEntry(result));

		idx_t out_offset = 0;
		for (idx_t out_row_idx = 0; out_row_idx < count; out_row_idx++) {

			auto in_row_idx = format.sel->get_index(out_row_idx);
			if (!format.validity.RowIsValid(in_row_idx)) {
				FlatVector::SetNull(result, out_row_idx, true);
				continue;
			}
			auto in = in_line_entries[in_row_idx];
			auto in_offset = in.offset;
			auto in_length = in.length;

			// Special case: if the line has less than 3 points, we can't remove any points
			if (in_length < 3) {

				ListVector::Reserve(result, out_offset + in_length);
				auto out_x_data = FlatVector::GetData<double>(*out_line_vertex_vec[0]);
				auto out_y_data = FlatVector::GetData<double>(*out_line_vertex_vec[1]);

				// If the line has less than 3 points, we can't remove any points
				// so we just copy the line
				out_line_entries[out_row_idx] = list_entry_t {out_offset, in_length};
				for (idx_t coord_idx = 0; coord_idx < in_length; coord_idx++) {
					out_x_data[out_offset + coord_idx] = in_x_data[in_offset + coord_idx];
					out_y_data[out_offset + coord_idx] = in_y_data[in_offset + coord_idx];
				}
				out_offset += in_length;
				continue;
			}

			// First pass, calculate how many points we need to keep
			// We always keep the first and last point, so we start at 2
			uint32_t points_to_keep = 0;

			auto last_x = in_x_data[in_offset];
			auto last_y = in_y_data[in_offset];
			points_to_keep++;

			for (idx_t i = 1; i < in_length; i++) {
				auto curr_x = in_x_data[in_offset + i];
				auto curr_y = in_y_data[in_offset + i];

				if (curr_x != last_x || curr_y != last_y) {
					points_to_keep++;
					last_x = curr_x;
					last_y = curr_y;
				}
			}

			// Special case: there is only 1 unique point in the line, so just keep
			// the start and end points
			if (points_to_keep == 1) {
				out_line_entries[out_row_idx] = list_entry_t {out_offset, 2};
				ListVector::Reserve(result, out_offset + 2);
				auto out_x_data = FlatVector::GetData<double>(*out_line_vertex_vec[0]);
				auto out_y_data = FlatVector::GetData<double>(*out_line_vertex_vec[1]);
				out_x_data[out_offset] = in_x_data[in_offset];
				out_y_data[out_offset] = in_y_data[in_offset];
				out_x_data[out_offset + 1] = in_x_data[in_offset + in_length - 1];
				out_y_data[out_offset + 1] = in_y_data[in_offset + in_length - 1];
				out_offset += 2;
				continue;
			}

			// Set the list entry
			out_line_entries[out_row_idx] = list_entry_t {out_offset, points_to_keep};

			// Second pass, copy the points we need to keep
			ListVector::Reserve(result, out_offset + points_to_keep);
			auto out_x_data = FlatVector::GetData<double>(*out_line_vertex_vec[0]);
			auto out_y_data = FlatVector::GetData<double>(*out_line_vertex_vec[1]);

			// Copy the first point
			out_x_data[out_offset] = in_x_data[in_offset];
			out_y_data[out_offset] = in_y_data[in_offset];
			out_offset++;

			// Copy the middle points (skip the last one, we'll copy it at the end)
			last_x = in_x_data[in_offset];
			last_y = in_y_data[in_offset];

			for (idx_t i = 1; i < in_length; i++) {
				auto curr_x = in_x_data[in_offset + i];
				auto curr_y = in_y_data[in_offset + i];

				if (curr_x != last_x || curr_y != last_y) {
					out_x_data[out_offset] = curr_x;
					out_y_data[out_offset] = curr_y;
					last_x = curr_x;
					last_y = curr_y;
					out_offset++;
				}
			}
		}
		ListVector::SetListSize(result, out_offset);

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D (With Tolerance)
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLineStringWithTolerance(DataChunk &args, ExpressionState &state, Vector &result) {
		auto input = args.data[0];
		auto tolerance = args.data[1];
		auto count = args.size();
		UnifiedVectorFormat format;
		input.ToUnifiedFormat(count, format);

		UnifiedVectorFormat tolerance_format;
		tolerance.ToUnifiedFormat(count, tolerance_format);

		auto in_line_entries = ListVector::GetData(input);
		auto &in_line_vertex_vec = StructVector::GetEntries(ListVector::GetEntry(input));
		auto in_x_data = FlatVector::GetData<double>(*in_line_vertex_vec[0]);
		auto in_y_data = FlatVector::GetData<double>(*in_line_vertex_vec[1]);

		auto out_line_entries = ListVector::GetData(result);
		auto &out_line_vertex_vec = StructVector::GetEntries(ListVector::GetEntry(result));

		idx_t out_offset = 0;

		for (idx_t out_row_idx = 0; out_row_idx < count; out_row_idx++) {
			auto in_row_idx = format.sel->get_index(out_row_idx);
			auto in_tol_idx = tolerance_format.sel->get_index(out_row_idx);
			if (!format.validity.RowIsValid(in_row_idx) || !tolerance_format.validity.RowIsValid(in_tol_idx)) {
				FlatVector::SetNull(result, out_row_idx, true);
				continue;
			}

			auto in = in_line_entries[in_row_idx];
			auto in_offset = in.offset;
			auto in_length = in.length;

			auto tolerance = Load<double>(tolerance_format.data + in_tol_idx);
			auto tolerance_squared = tolerance * tolerance;

			if (in_length < 3) {

				ListVector::Reserve(result, out_offset + in_length);
				auto out_x_data = FlatVector::GetData<double>(*out_line_vertex_vec[0]);
				auto out_y_data = FlatVector::GetData<double>(*out_line_vertex_vec[1]);

				// If the line has less than 3 points, we can't remove any points
				// so we just copy the line
				out_line_entries[out_row_idx] = list_entry_t {out_offset, in_length};
				for (idx_t coord_idx = 0; coord_idx < in_length; coord_idx++) {
					out_x_data[out_offset + coord_idx] = in_x_data[in_offset + coord_idx];
					out_y_data[out_offset + coord_idx] = in_y_data[in_offset + coord_idx];
				}
				out_offset += in_length;
				continue;
			}

			// First pass, calculate how many points we need to keep
			uint32_t points_to_keep = 0;

			auto last_x = in_x_data[in_offset];
			auto last_y = in_y_data[in_offset];
			points_to_keep++;

			for (idx_t i = 1; i < in_length; i++) {
				auto curr_x = in_x_data[in_offset + i];
				auto curr_y = in_y_data[in_offset + i];

				auto dist_squared = (curr_x - last_x) * (curr_x - last_x) + (curr_y - last_y) * (curr_y - last_y);

				if (dist_squared > tolerance_squared) {
					last_x = curr_x;
					last_y = curr_y;
					points_to_keep++;
				}
			}

			// Special case: there is only 1 unique point in the line, so just keep
			// the start and end points
			if (points_to_keep == 1) {
				out_line_entries[out_row_idx] = list_entry_t {out_offset, 2};
				ListVector::Reserve(result, out_offset + 2);
				auto out_x_data = FlatVector::GetData<double>(*out_line_vertex_vec[0]);
				auto out_y_data = FlatVector::GetData<double>(*out_line_vertex_vec[1]);
				out_x_data[out_offset] = in_x_data[in_offset];
				out_y_data[out_offset] = in_y_data[in_offset];
				out_x_data[out_offset + 1] = in_x_data[in_offset + in_length - 1];
				out_y_data[out_offset + 1] = in_y_data[in_offset + in_length - 1];
				out_offset += 2;
				continue;
			}

			// Set the list entry
			out_line_entries[out_row_idx] = list_entry_t {out_offset, points_to_keep};

			// Second pass, copy the points we need to keep
			ListVector::Reserve(result, out_offset + points_to_keep);
			auto out_x_data = FlatVector::GetData<double>(*out_line_vertex_vec[0]);
			auto out_y_data = FlatVector::GetData<double>(*out_line_vertex_vec[1]);

			// Copy the first point
			out_x_data[out_offset] = in_x_data[in_offset];
			out_y_data[out_offset] = in_y_data[in_offset];
			out_offset++;

			// With tolerance its different, we always keep the first and last point
			// regardless of distance to the previous point
			// Copy the middle points
			last_x = in_x_data[in_offset];
			last_y = in_y_data[in_offset];

			for (idx_t i = 1; i < in_length - 1; i++) {

				auto curr_x = in_x_data[in_offset + i];
				auto curr_y = in_y_data[in_offset + i];

				auto dist_squared = (curr_x - last_x) * (curr_x - last_x) + (curr_y - last_y) * (curr_y - last_y);
				if (dist_squared > tolerance_squared) {
					out_x_data[out_offset] = curr_x;
					out_y_data[out_offset] = curr_y;
					last_x = curr_x;
					last_y = curr_y;
					out_offset++;
				}
			}

			// Copy the last point
			out_x_data[points_to_keep - 1] = in_x_data[in_offset + in_length - 1];
			out_y_data[points_to_keep - 1] = in_y_data[in_offset + in_length - 1];
			out_offset++;
		}
		ListVector::SetListSize(result, out_offset);

		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Remove repeated points from a LINESTRING.
	)";

	// TODO: example
	static constexpr auto EXAMPLE = R"";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_RemoveRepeatedPoints", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("line", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(GeoTypes::LINESTRING_2D());

				variant.SetFunction(ExecuteLineString);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("line", GeoTypes::LINESTRING_2D());
				variant.AddParameter("tolerance", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::LINESTRING_2D());

				variant.SetFunction(ExecuteLineStringWithTolerance);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

//======================================================================================================================
// ST_StartPoint
//======================================================================================================================

struct ST_StartPoint {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
		    args.data[0], result, args.size(), [&](const string_t &blob, ValidityMask &mask, const idx_t idx) {
			    // TODO: Peek without deserializing!
			    const auto geom = lstate.Deserialize(blob);

			    if (geom.get_type() != sgl::geometry_type::LINESTRING) {
				    mask.SetInvalid(idx);
				    return string_t {};
			    }

			    if (geom.is_empty()) {
				    mask.SetInvalid(idx);
				    return string_t {};
			    }

			    const auto vertex_data = geom.get_vertex_data();
			    auto point = sgl::geometry(sgl::geometry_type::POINT, geom.has_z(), geom.has_m());
			    point.set_vertex_data(vertex_data, 1);

			    return lstate.Serialize(result, point);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		auto geom_vec = args.data[0];
		auto count = args.size();

		UnifiedVectorFormat geom_format;
		geom_vec.ToUnifiedFormat(count, geom_format);

		auto line_vertex_entries = ListVector::GetData(geom_vec);
		auto &line_vertex_vec = ListVector::GetEntry(geom_vec);
		auto &line_vertex_vec_children = StructVector::GetEntries(line_vertex_vec);
		auto line_x_data = FlatVector::GetData<double>(*line_vertex_vec_children[0]);
		auto line_y_data = FlatVector::GetData<double>(*line_vertex_vec_children[1]);

		auto &point_vertex_children = StructVector::GetEntries(result);
		auto point_x_data = FlatVector::GetData<double>(*point_vertex_children[0]);
		auto point_y_data = FlatVector::GetData<double>(*point_vertex_children[1]);

		for (idx_t out_row_idx = 0; out_row_idx < count; out_row_idx++) {
			auto in_row_idx = geom_format.sel->get_index(out_row_idx);

			if (!geom_format.validity.RowIsValid(in_row_idx)) {
				FlatVector::SetNull(result, out_row_idx, true);
				continue;
			}

			auto line = line_vertex_entries[in_row_idx];
			auto line_offset = line.offset;
			auto line_length = line.length;

			if (line_length == 0) {
				FlatVector::SetNull(result, out_row_idx, true);
				continue;
			}

			point_x_data[out_row_idx] = line_x_data[line_offset];
			point_y_data[out_row_idx] = line_y_data[line_offset];
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Returns the start point of a LINESTRING.
	)";

	// todo: add example
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_StartPoint", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("line", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(GeoTypes::POINT_2D());

				variant.SetFunction(ExecuteLineString);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_EndPoint
//======================================================================================================================

struct ST_EndPoint {

	//------------------------------------------------------------------------------------------------------------------
	// GEOMETRY
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
		    args.data[0], result, args.size(), [&](const string_t &blob, ValidityMask &mask, const idx_t idx) {
			    // TODO: Peek without deserializing!
			    const auto geom = lstate.Deserialize(blob);

			    if (geom.get_type() != sgl::geometry_type::LINESTRING) {
				    mask.SetInvalid(idx);
				    return string_t {};
			    }

			    if (geom.is_empty()) {
				    mask.SetInvalid(idx);
				    return string_t {};
			    }

			    const auto vertex_count = geom.get_count();
			    const auto vertex_size = geom.get_vertex_size();
			    const auto vertex_data = geom.get_vertex_data();

			    const auto point_data = vertex_data + ((vertex_count - 1) * vertex_size);

			    auto point = sgl::geometry(sgl::geometry_type::POINT, geom.has_z(), geom.has_m());
			    point.set_vertex_data(point_data, 1);

			    return lstate.Serialize(result, point);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// LINESTRING_2D
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteLineString(DataChunk &args, ExpressionState &state, Vector &result) {
		auto geom_vec = args.data[0];
		auto count = args.size();

		UnifiedVectorFormat geom_format;
		geom_vec.ToUnifiedFormat(count, geom_format);

		auto line_vertex_entries = ListVector::GetData(geom_vec);
		auto &line_vertex_vec = ListVector::GetEntry(geom_vec);
		auto &line_vertex_vec_children = StructVector::GetEntries(line_vertex_vec);
		auto line_x_data = FlatVector::GetData<double>(*line_vertex_vec_children[0]);
		auto line_y_data = FlatVector::GetData<double>(*line_vertex_vec_children[1]);

		auto &point_vertex_children = StructVector::GetEntries(result);
		auto point_x_data = FlatVector::GetData<double>(*point_vertex_children[0]);
		auto point_y_data = FlatVector::GetData<double>(*point_vertex_children[1]);

		for (idx_t out_row_idx = 0; out_row_idx < count; out_row_idx++) {
			auto in_row_idx = geom_format.sel->get_index(out_row_idx);

			if (!geom_format.validity.RowIsValid(in_row_idx)) {
				FlatVector::SetNull(result, out_row_idx, true);
				continue;
			}

			auto line = line_vertex_entries[in_row_idx];
			auto line_offset = line.offset;
			auto line_length = line.length;

			if (line_length == 0) {
				FlatVector::SetNull(result, out_row_idx, true);
				continue;
			}

			point_x_data[out_row_idx] = line_x_data[line_offset + line_length - 1];
			point_y_data[out_row_idx] = line_y_data[line_offset + line_length - 1];
		}
		if (count == 1) {
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
		}
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
		Returns the end point of a LINESTRING.
	)";

	// TODO: add example
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_EndPoint", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("line", GeoTypes::LINESTRING_2D());
				variant.SetReturnType(GeoTypes::POINT_2D());

				variant.SetFunction(ExecuteLineString);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

//======================================================================================================================
// ST_Within
//======================================================================================================================

struct ST_Within {

	//------------------------------------------------------------------------------------------------------------------
	// POINT_2D -> POLYGON_2D
	//------------------------------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &point_in = args.data[0];
		auto &polygon_in = args.data[1];

		// Just execute ST_Contains, but reversed
		ST_Contains::Operation(point_in, polygon_in, result, args.size());
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	// TODO: add example
	static constexpr auto DESCRIPTION = "";
	static constexpr auto EXAMPLE = "";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		// ST_Within is the inverse of ST_Contains
		FunctionBuilder::RegisterScalar(db, "ST_Within", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::POINT_2D());
				variant.AddParameter("geom2", GeoTypes::POLYGON_2D());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetFunction(Execute);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};



enum class VertexOrdinate { X, Y, Z, M };

template <class OP>
struct PointAccessFunctionBase {
	static size_t GetOrdinateOffset(const sgl::geometry &geom) {
		switch (OP::ORDINATE) {
		case VertexOrdinate::X:
			return 0;
		case VertexOrdinate::Y:
			return 1;
		case VertexOrdinate::Z:
			return 2;
		case VertexOrdinate::M:
			return geom.has_z() ? 3 : 2;
		default:
			return 0;
		}
	}

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::ExecuteWithNulls<string_t, double>(
		    args.data[0], result, args.size(), [&](const string_t &blob, ValidityMask &mask, const idx_t idx) {
			    const auto geom = lstate.Deserialize(blob);

			    if (geom.get_type() != sgl::geometry_type::POINT) {
				    throw InvalidInputException("%s only supports POINT geometries", OP::NAME);
			    }

			    if (geom.is_empty()) {
				    mask.SetInvalid(idx);
				    return 0.0;
			    }

			    if (OP::ORDINATE == VertexOrdinate::Z && !geom.has_z()) {
				    mask.SetInvalid(idx);
				    return 0.0;
			    }

			    if (OP::ORDINATE == VertexOrdinate::M && !geom.has_m()) {
				    mask.SetInvalid(idx);
				    return 0.0;
			    }

			    const auto vertex_data = geom.get_vertex_data();
			    const auto offset = GetOrdinateOffset(geom);

			    double res = 0.0;
			    memcpy(&res, vertex_data + offset * sizeof(double), sizeof(double));
			    return res;
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, OP::NAME, [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription(OP::DESCRIPTION);
				variant.SetExample(OP::EXAMPLE);
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct VertexMinAggOp {
	static double Init() {
		return std::numeric_limits<double>::max();
	}
	static double Merge(const double a, const double b) {
		return std::min(a, b);
	}
};

struct VertexMaxAggOp {
	static double Init() {
		return std::numeric_limits<double>::lowest();
	}
	static double Merge(const double a, const double b) {
		return std::max(a, b);
	}
};

template <class OP, class AGG>
struct VertexAggFunctionBase {
	static size_t GetOrdinateOffset(const sgl::geometry &geom) {
		switch (OP::ORDINATE) {
		case VertexOrdinate::X:
			return 0;
		case VertexOrdinate::Y:
			return 1;
		case VertexOrdinate::Z:
			return 2;
		case VertexOrdinate::M:
			return geom.has_z() ? 3 : 2;
		default:
			return 0;
		}
	}

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::ExecuteWithNulls<string_t, double>(
		    args.data[0], result, args.size(), [&](const string_t &blob, ValidityMask &mask, const idx_t idx) {
			    const auto geom = lstate.Deserialize(blob);

			    if (geom.is_empty()) {
				    mask.SetInvalid(idx);
				    return 0.0;
			    }
			    if (OP::ORDINATE == VertexOrdinate::Z && !geom.has_z()) {
				    mask.SetInvalid(idx);
				    return 0.0;
			    }
			    if (OP::ORDINATE == VertexOrdinate::M && !geom.has_m()) {
				    mask.SetInvalid(idx);
				    return 0.0;
			    }

			    const auto offset = GetOrdinateOffset(geom);

			    double res = AGG::Init();

			    sgl::ops::visit_vertices(&geom, [&](const uint8_t *vertex) {
				    double val = 0.0;
				    memcpy(&val, vertex + offset * sizeof(double), sizeof(double));

				    res = AGG::Merge(res, val);
			    });

			    return res;
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, OP::NAME, [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription(OP::DESCRIPTION);
				variant.SetExample(OP::EXAMPLE);
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_X : PointAccessFunctionBase<ST_X> {
	static constexpr auto NAME = "ST_X";
	static constexpr auto DESCRIPTION = "Returns the X coordinate of a point geometry";
	static constexpr auto EXAMPLE = "SELECT ST_X(ST_Point(1, 2))";
	static constexpr auto ORDINATE = VertexOrdinate::X;
};

struct ST_XMax : VertexAggFunctionBase<ST_XMax, VertexMaxAggOp> {
	static auto constexpr NAME = "ST_XMax";
	static auto constexpr DESCRIPTION = "Returns the maximum X coordinate of a geometry";
	static auto constexpr EXAMPLE = "SELECT ST_XMax(ST_Point(1, 2))";
	static auto constexpr ORDINATE = VertexOrdinate::X;
};

struct ST_XMin : VertexAggFunctionBase<ST_XMin, VertexMinAggOp> {
	static constexpr auto NAME = "ST_XMin";
	static constexpr auto DESCRIPTION = "Returns the minimum X coordinate of a geometry";
	static constexpr auto EXAMPLE = "SELECT ST_XMin(ST_Point(1, 2))";
	static constexpr auto ORDINATE = VertexOrdinate::X;
};

struct ST_Y : PointAccessFunctionBase<ST_Y> {
	static constexpr auto NAME = "ST_Y";
	static constexpr auto DESCRIPTION = "Returns the Y coordinate of a point geometry";
	static constexpr auto EXAMPLE = "SELECT ST_Y(ST_Point(1, 2))";
	static constexpr auto ORDINATE = VertexOrdinate::Y;
};

struct ST_YMax : VertexAggFunctionBase<ST_YMax, VertexMaxAggOp> {
	static constexpr auto NAME = "ST_YMax";
	static constexpr auto DESCRIPTION = "Returns the maximum Y coordinate of a geometry";
	static constexpr auto EXAMPLE = "SELECT ST_YMax(ST_Point(1, 2))";
	static constexpr auto ORDINATE = VertexOrdinate::Y;
};

struct ST_YMin : VertexAggFunctionBase<ST_YMin, VertexMinAggOp> {
	static constexpr auto NAME = "ST_YMin";
	static constexpr auto DESCRIPTION = "Returns the minimum Y coordinate of a geometry";
	static constexpr auto EXAMPLE = "SELECT ST_YMin(ST_Point(1, 2))";
	static constexpr auto ORDINATE = VertexOrdinate::Y;
};

struct ST_Z : PointAccessFunctionBase<ST_Z> {
	static constexpr auto NAME = "ST_Z";
	static constexpr auto DESCRIPTION = "Returns the Z coordinate of a point geometry";
	static constexpr auto EXAMPLE = "SELECT ST_Z(ST_Point(1, 2, 3))";
	static constexpr auto ORDINATE = VertexOrdinate::Z;
};

struct ST_ZMax : VertexAggFunctionBase<ST_ZMax, VertexMaxAggOp> {
	static auto constexpr NAME = "ST_ZMax";
	static auto constexpr DESCRIPTION = "Returns the maximum Z coordinate of a geometry";
	static auto constexpr EXAMPLE = "SELECT ST_ZMax(ST_Point(1, 2, 3))";
	static auto constexpr ORDINATE = VertexOrdinate::Z;
};

struct ST_ZMin : VertexAggFunctionBase<ST_ZMin, VertexMinAggOp> {
	static constexpr auto NAME = "ST_ZMin";
	static constexpr auto DESCRIPTION = "Returns the minimum Z coordinate of a geometry";
	static constexpr auto EXAMPLE = "SELECT ST_ZMin(ST_Point(1, 2, 3))";
	static constexpr auto ORDINATE = VertexOrdinate::Z;
};

struct ST_M : PointAccessFunctionBase<ST_M> {
	static constexpr auto NAME = "ST_M";
	static constexpr auto DESCRIPTION = "Returns the M coordinate of a point geometry";
	static constexpr auto EXAMPLE = "SELECT ST_M(ST_Point(1, 2, 3, 4))";
	static constexpr auto ORDINATE = VertexOrdinate::M;
};

struct ST_MMax : VertexAggFunctionBase<ST_MMax, VertexMaxAggOp> {
	static constexpr auto NAME = "ST_MMax";
	static constexpr auto DESCRIPTION = "Returns the maximum M coordinate of a geometry";
	static constexpr auto EXAMPLE = "SELECT ST_MMax(ST_Point(1, 2, 3, 4))";
	static constexpr auto ORDINATE = VertexOrdinate::M;
};

struct ST_MMin : VertexAggFunctionBase<ST_MMin, VertexMinAggOp> {
	static constexpr auto NAME = "ST_MMin";
	static constexpr auto DESCRIPTION = "Returns the minimum M coordinate of a geometry";
	static constexpr auto EXAMPLE = "SELECT ST_MMin(ST_Point(1, 2, 3, 4))";
	static constexpr auto ORDINATE = VertexOrdinate::M;
};

} // namespace

//######################################################################################################################
// Register
//######################################################################################################################

void CoreModule::RegisterSpatialFunctions(DatabaseInstance &db) {
	ST_Area::Register(db);

	// 7 functions to go!
	ST_AsGeoJSON::Register(db);
	ST_AsText::Register(db);
	ST_AsWKB::Register(db);
	ST_AsHEXWKB::Register(db);
	ST_AsSVG::Register(db);
	ST_Centroid::Register(db);
	ST_Collect::Register(db);
	ST_CollectionExtract::Register(db);
	ST_Contains::Register(db);
	ST_Dimension::Register(db);
	ST_Distance::Register(db);
	ST_Dump::Register(db);
	ST_EndPoint::Register(db);
	ST_Extent::Register(db);
	ST_ExteriorRing::Register(db);
	ST_FlipCoordinates::Register(db);
	ST_Force2D::Register(db);
	ST_Force3DZ::Register(db);
	ST_Force3DM::Register(db);
	ST_Force4D::Register(db);
	ST_GeometryType::Register(db);
	/*
	ST_GeomFromHEXWKB::Register(db);
	ST_GeomFromText::Register(db);
	ST_GeomFromWKB::Register(db);
	*/
	ST_HasZ::Register(db);
	ST_HasM::Register(db);
	ST_ZMFlag::Register(db);
	ST_Distance_Sphere::Register(db);
	ST_Hilbert::Register(db);
	// ST_Intersects::Register(db); - not applicable now
	// ST_IntersectsExtent::Register(db); - not applicable now
	ST_IsClosed::Register(db);
	ST_IsEmpty::Register(db);
	ST_Length::Register(db);
	ST_MakeEnvelope::Register(db);
	ST_MakeLine::Register(db);
	ST_MakePolygon::Register(db);
	ST_Multi::Register(db);
	ST_NGeometries::Register(db);
	ST_NInteriorRings::Register(db);
	ST_NPoints::Register(db);
	ST_Perimeter::Register(db);
	ST_Point::Register(db);
	ST_PointN::Register(db);
	ST_Points::Register(db);
	ST_QuadKey::Register(db);
	ST_RemoveRepeatedPoints::Register(db);
	ST_StartPoint::Register(db);
	ST_Within::Register(db);
	ST_X::Register(db);
	ST_XMax::Register(db);
	ST_XMin::Register(db);
	ST_Y::Register(db);
	ST_YMax::Register(db);
	ST_YMin::Register(db);
	ST_Z::Register(db);
	ST_ZMax::Register(db);
	ST_ZMin::Register(db);
	ST_M::Register(db);
	ST_MMax::Register(db);
	ST_MMin::Register(db);
}

} // namespace core

} // namespace spatial