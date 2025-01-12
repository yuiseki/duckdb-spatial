
#include "spatial/core/function_builder.hpp"
#include "spatial/core/module.hpp"

#include "spatial/core/types.hpp"
#include "spatial/core/util/math.hpp"

#define SGL_ASSERT(x) D_ASSERT(x)
#include "sgl/sgl.hpp"

#include "duckdb/common/vector_operations/generic_executor.hpp"

namespace spatial {
namespace core {

//------------------------------------------------------------------------------
// Util
//------------------------------------------------------------------------------

class BinaryReader {
public:
	BinaryReader (const char* ptr, const char* end) : beg(ptr), end(end), ptr(ptr) { }
	BinaryReader (const char* buffer, const size_t size) : BinaryReader(buffer, buffer + size) { }

	template<class T>
	T Read() {
		static_assert(std::is_trivially_copyable<T>::value, "Type must be trivially copyable");
		CheckSize(sizeof(T));
		T value;
		memcpy(&value, ptr, sizeof(T));
		ptr += sizeof(T);
		return value;
	}

	const char* Reserve(const size_t size) {
		CheckSize(size);
		const char* result = ptr;
		ptr += size;
		return result;
	}

	void Skip(const size_t size) {
		CheckSize(size);
		ptr += size;
	}

private:
	void CheckSize(const size_t size) const {
		if(ptr + size > end) {
			throw InternalException("Buffer overflow");
		}
	}

	const char* beg;
	const char* end;
	const char* ptr;
};

class BinaryWriter {
public:
	BinaryWriter(char* ptr, char* end) : beg(ptr), end(end), ptr(ptr) { }
	BinaryWriter(char* buffer, const size_t size) : BinaryWriter(buffer, buffer + size) { }

	template<class T>
	void Write(const T &value) {
		static_assert(std::is_trivially_copyable<T>::value, "Type must be trivially copyable");
		CheckSize(sizeof(T));
		memcpy(ptr, &value, sizeof(T));
		ptr += sizeof(T);
	}

	char* Reserve(const size_t size) {
		CheckSize(size);
		char* result = ptr;
		ptr += size;
		return result;
	}

	void Skip(const size_t size, const bool zero = false) {
		CheckSize(size);
		if(zero) {
			memset(ptr, 0, size);
		}
		ptr += size;
	}

	void Copy(const char* buffer, const size_t size) {
		CheckSize(size);
		memcpy(ptr, buffer, size);
		ptr += size;
	}
private:
	void CheckSize(const size_t size) const {
		if(ptr + size > end) {
			throw InternalException("Buffer overflow");
		}
	}

	char* beg;
	char* end;
	char* ptr;
};

// todo:
struct Serde {
	static size_t GetRequiredSize(const sgl::geometry &geom);
	static void Serialize(const sgl::geometry &geom, char* buffer, size_t buffer_size);
	static void Deserialize(sgl::geometry &result, ArenaAllocator &arena, const char* buffer, size_t buffer_size);
};

static size_t GetRequiredSizeInternal(const sgl::geometry *geom) {
	const auto vertex_size = geom->get_vertex_size();
	const auto part_count = geom->get_count();

	switch(geom->get_type()) {
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
			if(!tail) {
				return size;
			}
			auto part = tail;
			do {
				part = part->get_next();
				size += 4 + part->get_count() * vertex_size;
			} while(part != tail);

			if(part_count % 2 == 1) {
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
			if(!tail) {
				return size;
			}
			auto part = tail;
			do {
				part = part->get_next();
				size += GetRequiredSizeInternal(part);
			} while(part != tail);
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


static void SerializeVertices(BinaryWriter &cursor, const sgl::geometry *geom, const uint32_t count,
							  const bool has_z, const bool has_m, const bool has_bbox, const uint32_t vsize,
							  sgl::box_xyzm &bbox) {

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

	if(type < sgl::geometry_type::POINT || type > sgl::geometry_type::MULTI_GEOMETRY) {
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

void Serde::Deserialize(sgl::geometry &result, ArenaAllocator &arena, const char* buffer, size_t buffer_size) {

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

	if(has_bbox) {
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
	class LocalState final : public FunctionLocalState {
	public:
		explicit LocalState(ClientContext &context) : arena(BufferAllocator::Get(context)) { }
		static unique_ptr<FunctionLocalState> Init(ExpressionState &state, const BoundFunctionExpression &expr, FunctionData *bind_data);
		static LocalState &ResetAndGet(ExpressionState &state);

		// De/Serialize geometries
		sgl::geometry Deserialize(const string_t &blob);
		string_t Serialize(Vector &vector, const sgl::geometry &geom);

		ArenaAllocator &GetArena() { return arena; }
	private:
		ArenaAllocator arena;
	};

	unique_ptr<FunctionLocalState> LocalState::Init(ExpressionState &state, const BoundFunctionExpression &expr, FunctionData *bind_data) {
		return make_uniq_base<FunctionLocalState, LocalState>(state.GetContext());
	}

	LocalState& LocalState::ResetAndGet(ExpressionState &state) {
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
}

//------------------------------------------------------------------------------
// Functions
//------------------------------------------------------------------------------

namespace {

struct ST_Area {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				const auto geom = lstate.Deserialize(blob);
				return sgl::ops::area(&geom);
		});
	}

	static void Register(DatabaseInstance &db) {

		FunctionBuilder::RegisterScalar(db, "ST_Area", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Compute the area of a geometry");
				// TODO: Set example
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_AsGeoJSON {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_AsText {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_AsWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_AsHEXWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_AsSVG {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Centroid {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Collect {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_CollectionExtract {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Contains {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Dimension {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				const auto geom = lstate.Deserialize(blob);
				return sgl::ops::max_surface_dimension(&geom);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Dimension", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::INTEGER);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Returns the max surface dimension of a geometry");
				// TODO: Set example
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_Distance {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Dump {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_EndPoint {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Extent {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_ExteriorRing {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				const auto geom = lstate.Deserialize(blob);

				if(geom.get_type() != sgl::geometry_type::POLYGON) {
					return lstate.Serialize(result, sgl::polygon::make_empty());
				}

				const auto shell = geom.get_first_part();

				if(!shell) {
					return lstate.Serialize(result, sgl::polygon::make_empty());
				}

				return lstate.Serialize(result, *shell);
		});
	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_FlipCoordinates {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Force {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_GeometryType {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_GeomFromHEXWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_GeomFromText {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_GeomFromWKB {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Has {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Haversine {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, string_t, double>(args.data[0], args.data[1], result, args.size(),
			[&](const string_t &l_blob, const string_t &r_blob) {
				const auto lhs = lstate.Deserialize(l_blob);
				const auto rhs = lstate.Deserialize(r_blob);

				if(lhs.get_type() != sgl::geometry_type::POINT || rhs.get_type() != sgl::geometry_type::POINT) {
					throw InvalidInputException("ST_Distance_Sphere only accepts POINT geometries");
				}

				if(lhs.is_empty() || rhs.is_empty()) {
					throw InvalidInputException("ST_Distance_Sphere does not accept empty geometries");
				}

				const auto lv = lhs.get_vertex_xy(0);
				const auto rv = rhs.get_vertex_xy(0);

				return sgl::util::haversine_distance(lv.x, lv.y, rv.x, rv.y);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Distance_Sphere", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription(R"(
				    Returns the haversine distance between two geometries.

				    - Only supports POINT geometries.
				    - Returns the distance in meters.
				    - The input is expected to be in WGS84 (EPSG:4326) coordinates, using a [latitude, longitude] axis order.
				)");
			});
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_Hilbert {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Intersects {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_IntersectsExtent {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_IsClosed {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				const auto geom = lstate.Deserialize(blob);
				switch(geom.get_type()) {
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

struct ST_IsEmpty {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				const auto geom = lstate.Deserialize(blob);
				const auto vertex_count = sgl::ops::vertex_count(&geom);
				return vertex_count == 0;
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_IsEmpty", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Check if a geometry is empty");
				// todo: Set example
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_Length {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				const auto geom = lstate.Deserialize(blob);
				return sgl::ops::length(&geom);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Length", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);

				variant.SetDescription("Compute the length of a geometry");
				// TODO: Set example
			});

			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

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
			[&](const DOUBLE_TYPE vmin_x, const DOUBLE_TYPE vmin_y, const DOUBLE_TYPE vmax_x, const DOUBLE_TYPE vmax_y) {

				const auto min_x = vmin_x.val;
				const auto min_y = vmin_y.val;
				const auto max_x = vmax_x.val;
				const auto max_y = vmax_y.val;

				// This is pretty cool, we dont even need to allocate anything
				const double buffer[10] = {min_x, min_y, min_x, max_y, max_x, max_y, max_x, min_y, min_x, min_y};

				auto ring = sgl::linestring::make_empty(false, false);
				ring.set_vertex_data(reinterpret_cast<const char*>(buffer), 5);

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

struct ST_MakeLine {
	static void ExecuteList(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		auto &child_vec = ListVector::GetEntry(args.data[0]);
		auto child_len = ListVector::GetListSize(args.data[0]);

		UnifiedVectorFormat format;
		child_vec.ToUnifiedFormat(child_len, format);

		UnaryExecutor::Execute<list_entry_t, string_t>(args.data[0], result, args.size(),
			[&](const list_entry_t &entry) {
				const auto offset = entry.offset;
				const auto length = entry.length;

				if(length < 2) {
					// Early out if we don't have enough geometries
					throw InvalidInputException("ST_MakeLine requires at least 2 geometries");
				}

				uint32_t line_length = 0;
				// First pass, filter types, count non-null entries

				for(idx_t i = offset; i < offset + length; i++) {
					const auto mapped_idx = format.sel->get_index(i);
					if(format.validity.RowIsValid(mapped_idx)) {
						continue;
					}
					auto &blob = UnifiedVectorFormat::GetData<string_t>(format)[mapped_idx];

					// TODO: Peek without deserializing
					const auto geom = lstate.Deserialize(blob);
					if(geom.get_type() != sgl::geometry_type::POINT) {
						throw InvalidInputException("ST_MakeLine only accepts POINT geometries");
					}

					// TODO: Support Z and M
					if(geom.has_z() || geom.has_m()) {
						throw InvalidInputException("ST_MakeLine from list does not accept POINT geometries with Z or M values");
					}

					line_length++;
				}

				if(line_length < 2) {
					throw InvalidInputException("ST_MakeLine requires at least 2 non-null geometries");
				}

				const auto line_data = lstate.GetArena().AllocateAligned(line_length * 2 * sizeof(double));

				// Second pass, copy over the vertex data
				uint32_t vertex_idx = 0;
				for(idx_t i = offset; i < offset + length; i++) {
					D_ASSERT(vertex_idx < line_length);

					const auto mapped_idx = format.sel->get_index(i);
					if(format.validity.RowIsValid(mapped_idx)) {
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
			args.data[0], args.data[1], result, args.size(),
			[&](const string_t &l_blob, const string_t &r_blob) {

				const auto l_geom = lstate.Deserialize(l_blob);
				const auto r_geom = lstate.Deserialize(r_blob);

				if(l_geom.get_type() != sgl::geometry_type::POINT || r_geom.get_type() != sgl::geometry_type::POINT) {
					throw InvalidInputException("ST_MakeLine only accepts POINT geometries");
				}

				if(l_geom.is_empty() || r_geom.is_empty()) {
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
				if(has_z) {
					buffer[idx++] = l_geom.has_z() ? v1.zm : 0;
				}
				if(has_m) {
					buffer[idx++] = l_geom.has_m() ? l_geom.has_z() ? v1.m : v1.zm : 0;
				}
				buffer[idx++] = v2.x;
				buffer[idx++] = v2.y;
				if(has_z) {
					buffer[idx++] = r_geom.has_z() ? v2.zm : 0;
				}
				if(has_m) {
					buffer[idx++] = r_geom.has_m() ? r_geom.has_z() ? v2.m : v2.zm : 0;
				}

				linestring.set_vertex_data(reinterpret_cast<const char*>(buffer), 2);

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

struct ST_MakePolygon {
	static void ExecuteFromShell(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				auto line = lstate.Deserialize(blob);

				if(line.get_type() != sgl::geometry_type::LINESTRING) {
					throw InvalidInputException("ST_MakePolygon only accepts LINESTRING geometries");
				}

				if(line.get_count() < 4) {
					throw InvalidInputException("ST_MakePolygon shell requires at least 4 vertices");
				}

				if(!sgl::linestring::is_closed(&line)) {
					throw std::runtime_error("ST_MakePolygon shell must be closed (first and last vertex must be equal)");
				}

				auto polygon = sgl::polygon::make_empty(line.has_z(), line.has_m());
				polygon.append_part(&line);

				return lstate.Serialize(result, polygon);
		});
	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Multi {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(),
			[&](const string_t &blob) {

				auto geom = lstate.Deserialize(blob);
				const auto has_z = geom.has_z();
				const auto has_m = geom.has_m();

				switch(geom.get_type()) {
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

struct ST_NGeometries {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
				const auto geom = lstate.Deserialize(blob);
				if(geom.is_single_part()) {
					return geom.is_empty() ? 0 : 1;
				}
				if(geom.is_multi_part()) {
					return static_cast<int32_t>(geom.get_count());
				}
				return 0;
			});
	}

	static void Register(DatabaseInstance &db) {
		// TODO: Maybe make a macro for the aliases
		for(auto &alias : {"ST_NumGeometries", "ST_NGeometries"}) {
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

struct ST_NInteriorRings {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::ExecuteWithNulls<string_t, int32_t>(args.data[0], result, args.size(),
			[&](const string_t &blob, ValidityMask &validity, idx_t idx) {
				const auto geom = lstate.Deserialize(blob);

				if(geom.get_type() != sgl::geometry_type::POLYGON) {
					validity.SetInvalid(idx);
					return 0;
				}

				const auto n_rings = static_cast<int32_t>(geom.get_count());
				return n_rings == 0 ? 0 : n_rings - 1;
		});
	}

	static void Register(DatabaseInstance &db) {
		// todo: maybe make a macro for the aliases
		for(auto &alias : {"ST_NumInteriorRings", "ST_NInteriorRings"}) {
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

struct ST_NPoints {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, int32_t>(args.data[0], result, args.size(),
			[&](const string_t &blob) {
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

struct ST_Perimeter {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Point {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<double, double, string_t>(args.data[0], args.data[1], result, args.size(),
			[&](const double x, const double y) {

				const double buffer[2] = {x, y};

				sgl::geometry geometry;
				geometry.set_type(sgl::geometry_type::POINT);
				geometry.set_vertex_data(reinterpret_cast<const uint8_t*>(buffer), 1);

				return lstate.Serialize(result, geometry);
		});
	}

	static void Register(DatabaseInstance &db) {

	}
};

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

struct ST_Points {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_QuadKey {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_RemoveRepeatedPoints {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_StartPoint {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};


enum class VertexOrdinate {
	X,
	Y,
	Z,
	M
};

template<class OP>
struct PointAccessFunctionBase {
	static size_t GetOrdinateOffset(const sgl::geometry& geom) {
		switch(OP::ORDINATE) {
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

		UnaryExecutor::ExecuteWithNulls<string_t, double>(args.data[0], result, args.size(),
			[&](const string_t &blob, ValidityMask &mask, const idx_t idx) {
				const auto geom = lstate.Deserialize(blob);

				if(geom.get_type() != sgl::geometry_type::POINT) {
					throw InvalidInputException("%s only supports POINT geometries", OP::NAME);
				}

				if(geom.is_empty()) {
					mask.SetInvalid(idx);
					return 0.0;
				}

				if(OP::ORDINATE == VertexOrdinate::Z && !geom.has_z()) {
					mask.SetInvalid(idx);
					return 0.0;
				}

				if(OP::ORDINATE == VertexOrdinate::M && !geom.has_m()) {
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


template<class OP, class AGG>
struct VertexAggFunctionBase {
	static size_t GetOrdinateOffset(const sgl::geometry& geom) {
		switch(OP::ORDINATE) {
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
		UnaryExecutor::ExecuteWithNulls<string_t, double>(args.data[0], result, args.size(),
			[&](const string_t &blob, ValidityMask &mask, const idx_t idx) {
				const auto geom = lstate.Deserialize(blob);

				if(geom.is_empty()) {
					mask.SetInvalid(idx);
					return 0.0;
				}
				if(OP::ORDINATE == VertexOrdinate::Z && !geom.has_z()) {
					mask.SetInvalid(idx);
					return 0.0;
				}
				if(OP::ORDINATE == VertexOrdinate::M && !geom.has_m()) {
					mask.SetInvalid(idx);
					return 0.0;
				}

				const auto offset = GetOrdinateOffset(geom);

				double res = AGG::Init();

				sgl::ops::visit_vertices(&geom, [&](const uint8_t* vertex) {

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

	// 25 functions to go!

	/*
	ST_AsGeoJSON::Register(db);
	ST_AsText::Register(db);
	ST_AsWKB::Register(db);
	ST_AsHEXWKB::Register(db);
	ST_AsSVG::Register(db);
	// ST_Centroid::Register(db); - not applicable now
	ST_Collect::Register(db);
	ST_CollectionExtract::Register(db);
	// ST_Contains::Register(db); - not applicable now
	*/
	ST_Dimension::Register(db);
	/*
	// ST_Distance::Register(db); -- not applicable now
	ST_Dump::Register(db);
	ST_EndPoint::Register(db);
	ST_Extent::Register(db);
	ST_ExteriorRing::Register(db);
	ST_FlipCoordinates::Register(db);
	ST_Force::Register(db);
	ST_GeometryType::Register(db);
	ST_GeomFromHEXWKB::Register(db);
	ST_GeomFromText::Register(db);
	ST_GeomFromWKB::Register(db);
	ST_Has::Register(db);
	*/
	ST_Haversine::Register(db);
	/*
	ST_Hilbert::Register(db);
	// ST_Intersects::Register(db); - not applicable now
	// ST_IntersectsExtent::Register(db); - not applicable now
	*/
	ST_IsClosed::Register(db);
	ST_IsEmpty::Register(db);
	ST_Length::Register(db);

	ST_MakeEnvelope::Register(db);
	ST_MakeLine::Register(db);
	//ST_MakePolygon::Register(db);

	ST_Multi::Register(db);
	ST_NGeometries::Register(db);

	ST_NInteriorRings::Register(db);
	ST_NPoints::Register(db);

	//ST_Perimeter::Register(db);
	ST_Point::Register(db);
	ST_PointN::Register(db);
	//ST_Points::Register(db);
	//ST_QuadKey::Register(db);
	//ST_RemoveRepeatedPoints::Register(db);
	//ST_StartPoint::Register(db);

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

}

}