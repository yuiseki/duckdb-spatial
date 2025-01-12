
#include "spatial/core/function_builder.hpp"
#include "spatial/core/module.hpp"

#include "spatial/core/types.hpp"
#include "spatial/core/util/math.hpp"

#define SGL_ASSERT(x) D_ASSERT(x)
#include "sgl/sgl.hpp"


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

	}

	static void Register(DatabaseInstance &db) {

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

	}

	static void Register(DatabaseInstance &db) {

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

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_IsEmpty {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Length {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_MakeEnvelope {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_MakeLine {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

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

				auto polygon = sgl::polygon::make_empty();
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

				switch(geom.get_type()) {
					case sgl::geometry_type::POINT: {
						auto mpoint = sgl::multi_point::make_empty();
						mpoint.append_part(&geom);
						return lstate.Serialize(result, mpoint);
					}
					case sgl::geometry_type::LINESTRING: {
						auto mline = sgl::multi_line_string::make_empty();
						mline.append_part(&geom);
						return lstate.Serialize(result, mline);
					}
					case sgl::geometry_type::POLYGON: {
						auto mpoly = sgl::multi_polygon::make_empty();
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

	}
};

struct ST_NGeometries {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_NInteriorRings {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_NPoints {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

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

struct ST_X {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_XMax {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_XMin {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Y {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_YMax {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_YMin {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_Z {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_ZMax {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_ZMin {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_M {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_MMax {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

struct ST_MMin {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

	}

	static void Register(DatabaseInstance &db) {

	}
};

} // namespace

void CoreModule::RegisterSpatialFunctions(DatabaseInstance &db) {
	ST_Area::Register(db);
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
	ST_Force::Register(db);
	ST_GeometryType::Register(db);
	ST_GeomFromHEXWKB::Register(db);
	ST_GeomFromText::Register(db);
	ST_GeomFromWKB::Register(db);
	ST_Has::Register(db);
	ST_Haversine::Register(db);
	ST_Hilbert::Register(db);
	ST_Intersects::Register(db);
	ST_IntersectsExtent::Register(db);
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