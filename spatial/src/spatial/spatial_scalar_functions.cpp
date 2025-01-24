#include "spatial/core/function_builder.hpp"
#include "spatial/core/module.hpp"
#include "spatial/core/types.hpp"
#include "spatial/core/util/math.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"

#define SGL_ASSERT(x) D_ASSERT(x)
#include "sgl/sgl.hpp"

namespace spatial {
namespace core {

//------------------------------------------------------------------------------
// Util
//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------
// LocalState
//------------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------------------------------------
// ST_AsGeoJSON
//----------------------------------------------------------------------------------------------------------------------
struct ST_AsGeoJSON {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_AsText
//----------------------------------------------------------------------------------------------------------------------
struct ST_AsText {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_AsWKB
//----------------------------------------------------------------------------------------------------------------------
struct ST_AsWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_AsHEXWKB
//----------------------------------------------------------------------------------------------------------------------
struct ST_AsHEXWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_AsSVG
//----------------------------------------------------------------------------------------------------------------------
struct ST_AsSVG {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_Centroid
//----------------------------------------------------------------------------------------------------------------------
struct ST_Centroid {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
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

//----------------------------------------------------------------------------------------------------------------------
// ST_Contains
//----------------------------------------------------------------------------------------------------------------------
struct ST_Contains {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
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
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//======================================================================================================================
// ST_Dump
//======================================================================================================================

struct ST_Dump {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
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

// TODO: Implement
struct ST_FlipCoordinates {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
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

//----------------------------------------------------------------------------------------------------------------------
// ST_RemoveRepeatedPoints
//----------------------------------------------------------------------------------------------------------------------
struct ST_RemoveRepeatedPoints {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	}

	static void Register(DatabaseInstance &db) {
	}
};

//----------------------------------------------------------------------------------------------------------------------
// ST_StartPoint
//----------------------------------------------------------------------------------------------------------------------
struct ST_StartPoint {

	// Geometry Function -----------------------------------------------------------------------------------------------
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
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

	// Register --------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_StartPoint", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Returns the start point of a LINESTRING");
				// todo: Set example
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_EndPoint {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
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

	static void Register(DatabaseInstance &db) {
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

void CoreModule::RegisterSpatialFunctions(DatabaseInstance &db) {
	ST_Area::Register(db);

	// 10 functions to go!
	/*
	ST_AsGeoJSON::Register(db);
	ST_AsText::Register(db);
	ST_AsWKB::Register(db);
	ST_AsHEXWKB::Register(db);
	ST_AsSVG::Register(db);
	// ST_Centroid::Register(db); - not applicable now
	*/
	ST_Collect::Register(db);
	ST_CollectionExtract::Register(db);
	// ST_Contains::Register(db); - not applicable now
	ST_Dimension::Register(db);
	/*
	// ST_Distance::Register(db); -- not applicable now
	ST_Dump::Register(db);
	*/
	ST_EndPoint::Register(db);
	ST_Extent::Register(db);
	ST_ExteriorRing::Register(db);
	/*
	ST_FlipCoordinates::Register(db);
	*/
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
	// ST_RemoveRepeatedPoints::Register(db); - not applicable right now
	ST_StartPoint::Register(db);
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