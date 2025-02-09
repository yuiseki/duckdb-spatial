#include "spatial/modules/proj/proj_module.hpp"
#include "spatial/spatial_types.hpp"
#include "spatial/util/function_builder.hpp"
#include "spatial/geometry/sgl.hpp"
#include "spatial/geometry/geometry_serialization.hpp"

#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/planner/expression/bound_function_expression.hpp"

#include "proj.h"
#include "geodesic.h"
#include "sqlite3.h"

// We embed the whole proj.db in the proj_db.c file, which we then link into the extension binary
// We can then use the sqlite3 "memvfs" (which we also statically link to) to point to the proj.db database in memory
// To genereate the proj_db.c file, we use the following command:
// `xxd -i proj.db > proj_db.c`
// Then rename the array to proj_db and the length to proj_db_len if necessary
// We link these from the proj_db.c file externally instead of #include:ing so our IDE doesnt go haywire
extern "C" unsigned char proj_db[];
extern "C" unsigned int proj_db_len;
extern "C" int sqlite3_memvfs_init(sqlite3 *, char **, const sqlite3_api_routines *);

// Specialize hash for std::pair<std::string, std::string> so we can use it as a key in an unordered_map
template <>
struct std::hash<std::pair<std::string, std::string>> {
	size_t operator()(pair<string, string> const &v) const noexcept {
		const auto lhs = std::hash<string> {}(v.first);
		const auto rhs = std::hash<string> {}(v.second);
		// Shift by one so we dont match the hash of the reversed pair
		return lhs ^ (rhs << 1);
	}
};

namespace duckdb {

namespace {

//######################################################################################################################
// PROJ Module & SQLITE VFS Registration
//######################################################################################################################

struct ProjModule {
	static void RegisterVFS(DatabaseInstance &db);
	static PJ_CONTEXT *GetThreadProjContext();
};

PJ_CONTEXT *ProjModule::GetThreadProjContext() {

	const auto ctx = proj_context_create();

	// We set the default context proj.db path to the one in the binary here
	// Otherwise GDAL will try to load the proj.db from the system
	// Any PJ_CONTEXT we create after this will inherit these settings
	const auto path = StringUtil::Format("file:/proj.db?ptr=%llu&sz=%lu&max=%lu", static_cast<void *>(proj_db),
	                                     proj_db_len, proj_db_len);

	proj_context_set_sqlite3_vfs_name(ctx, "memvfs");
	const auto ok = proj_context_set_database_path(ctx, path.c_str(), nullptr, nullptr);
	if (!ok) {
		throw InternalException("Could not set proj.db path");
	}

	// Dont log errors to stderr
	proj_log_level(ctx, PJ_LOG_NONE);

	// Dont allow network
	proj_context_set_enable_network(ctx, false);

	return ctx;
}

// IMPORTANT: Make sure this module is loaded before any other modules that use proj (like GDAL)
void ProjModule::RegisterVFS(DatabaseInstance &db) {
	// we use the sqlite "memvfs" to store the proj.db database in the extension binary itself
	// this way we don't have to worry about the user having the proj.db database installed
	// on their system. We therefore have to tell proj to use memvfs as the sqlite3 vfs and
	// point it to the segment of the binary that contains the proj.db database

	sqlite3_initialize();
	sqlite3_memvfs_init(nullptr, nullptr, nullptr);
	const auto vfs = sqlite3_vfs_find("memvfs");
	if (!vfs) {
		throw InternalException("Could not find sqlite memvfs extension");
	}
	sqlite3_vfs_register(vfs, 0);

	// We set the default context proj.db path to the one in the binary here
	// Otherwise GDAL will try to load the proj.db from the system
	// Any PJ_CONTEXT we create after this will inherit these settings (on this thread?)
	const auto path = StringUtil::Format("file:/proj.db?ptr=%llu&sz=%lu&max=%lu", static_cast<void *>(proj_db),
	                                     proj_db_len, proj_db_len);

	proj_context_set_sqlite3_vfs_name(nullptr, "memvfs");

	const auto ok = proj_context_set_database_path(nullptr, path.c_str(), nullptr, nullptr);
	if (!ok) {
		throw InternalException("Could not set proj.db path");
	}
}

//######################################################################################################################
// Coordinate Transformation Functions
//######################################################################################################################

//======================================================================================================================
// Local State
//======================================================================================================================

struct ProjCRSDelete {
	void operator()(PJ *crs) const {
		proj_destroy(crs);
	}
};

using ProjCRS = unique_ptr<PJ, ProjCRSDelete>;

struct ProjFunctionLocalState final : FunctionLocalState {

	PJ_CONTEXT *proj_ctx;
	ArenaAllocator arena;

	// Cache for PJ* objects
	unordered_map<std::pair<string, string>, ProjCRS> crs_cache;

	explicit ProjFunctionLocalState(ClientContext &context)
	    : proj_ctx(ProjModule::GetThreadProjContext()), arena(BufferAllocator::Get(context)) {
	}

	~ProjFunctionLocalState() override {
		proj_context_destroy(proj_ctx);
	}

	sgl::geometry Deserialize(const string_t &blob);
	string_t Serialize(Vector &vector, const sgl::geometry &geom);

	static unique_ptr<FunctionLocalState> Init(ExpressionState &state, const BoundFunctionExpression &expr,
	                                           FunctionData *bind_data) {
		auto result = make_uniq<ProjFunctionLocalState>(state.GetContext());
		return std::move(result);
	}

	static ProjFunctionLocalState &ResetAndGet(ExpressionState &state) {
		auto &local_state = ExecuteFunctionState::GetFunctionState(state)->Cast<ProjFunctionLocalState>();
		local_state.arena.Reset();
		return local_state;
	}

	PJ *GetOrCreateProjection(const string &source, const string &target, bool normalize) {
		const auto crs_entry = crs_cache.find({source, target});
		if (crs_entry != crs_cache.end()) {
			return crs_entry->second.get();
		}

		auto crs = proj_create_crs_to_crs(proj_ctx, source.c_str(), target.c_str(), nullptr);
		if (!crs) {
			throw InvalidInputException("Could not create projection: " + source + " -> " + target);
		}

		if (normalize) {
			const auto normalized_crs = proj_normalize_for_visualization(proj_ctx, crs);
			proj_destroy(crs);
			if (!normalized_crs) {
				throw InvalidInputException("Could not normalize projection: " + source + " -> " + target);
			}
			crs = normalized_crs;
		}

		crs_cache[{source, target}] = ProjCRS(crs);
		return crs;
	}
};

sgl::geometry ProjFunctionLocalState::Deserialize(const string_t &blob) {
	sgl::geometry geom;
	Serde::Deserialize(geom, arena, blob.GetDataUnsafe(), blob.GetSize());
	return geom;
}

string_t ProjFunctionLocalState::Serialize(Vector &vector, const sgl::geometry &geom) {
	const auto size = Serde::GetRequiredSize(geom);
	auto blob = StringVector::EmptyString(vector, size);
	Serde::Serialize(geom, blob.GetDataWriteable(), size);
	blob.Finalize();
	return blob;
}

//======================================================================================================================
// ST_Transform
//======================================================================================================================

struct ST_Transform {

	//------------------------------------------------------------------------------------------------------------------
	// Bind
	//------------------------------------------------------------------------------------------------------------------
	struct BindData final : FunctionData {
		bool normalize = false;

		unique_ptr<FunctionData> Copy() const override {
			auto result = make_uniq<BindData>();
			result->normalize = normalize;
			return std::move(result);
		}

		bool Equals(const FunctionData &other) const override {
			auto &data = other.Cast<BindData>();
			return normalize == data.normalize;
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &ctx, ScalarFunction &, vector<unique_ptr<Expression>> &args) {
		auto result = make_uniq<BindData>();
		if (args.size() == 4) {
			// Ensure the "always_xy" parameter is a constant
			const auto &arg = args[3];
			if (arg->HasParameter()) {
				throw InvalidInputException("The 'always_xy' parameter must be a constant");
			}
			if (!arg->IsFoldable()) {
				throw InvalidInputException("The 'always_xy' parameter must be a constant");
			}
			result->normalize = BooleanValue::Get(ExpressionExecutor::EvaluateScalar(ctx, *arg));
		}
		return std::move(result);
	}

	//------------------------------------------------------------------------------------------------------------------
	// Execute (POINT_2D)
	//------------------------------------------------------------------------------------------------------------------
	static void ExecutePoint(DataChunk &args, ExpressionState &state, Vector &result) {
		using POINT_TYPE = StructTypeBinary<double, double>;
		using PROJ_TYPE = PrimitiveType<string_t>;

		auto &lstate = ProjFunctionLocalState::ResetAndGet(state);
		auto &func_expr = state.expr.Cast<BoundFunctionExpression>();
		const auto &info = func_expr.bind_info->Cast<BindData>();

		GenericExecutor::ExecuteTernary<POINT_TYPE, PROJ_TYPE, PROJ_TYPE, POINT_TYPE>(
		    args.data[0], args.data[2], args.data[3], result, args.size(),
		    [&](const POINT_TYPE &point_in, const PROJ_TYPE &source, const PROJ_TYPE target) {
			    const auto source_str = source.val.GetString();
			    const auto target_str = target.val.GetString();

			    const auto crs = lstate.GetOrCreateProjection(source_str, target_str, info.normalize);

			    POINT_TYPE point_out;
			    const auto transformed = proj_trans(crs, PJ_FWD, proj_coord(point_in.a_val, point_in.b_val, 0, 0)).xy;
			    point_out.a_val = transformed.x;
			    point_out.b_val = transformed.y;

			    return point_out;
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Execute (BOX_2D)
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteBox(DataChunk &args, ExpressionState &state, Vector &result) {
		using BOX_TYPE = StructTypeQuaternary<double, double, double, double>;
		using PROJ_TYPE = PrimitiveType<string_t>;

		auto &lstate = ProjFunctionLocalState::ResetAndGet(state);
		auto &func_expr = state.expr.Cast<BoundFunctionExpression>();
		const auto &info = func_expr.bind_info->Cast<BindData>();

		GenericExecutor::ExecuteTernary<BOX_TYPE, PROJ_TYPE, PROJ_TYPE, BOX_TYPE>(
		    args.data[0], args.data[1], args.data[2], result, args.size(),
		    [&](const BOX_TYPE &box_in, const PROJ_TYPE source, const PROJ_TYPE &target) {
			    const auto source_str = source.val.GetString();
			    const auto target_str = target.val.GetString();

			    const auto crs = lstate.GetOrCreateProjection(source_str, target_str, info.normalize);

			    // TODO: this may be interesting to use, but at that point we can only return a BOX_TYPE
			    constexpr int densify_pts = 0;
			    BOX_TYPE box_out;
			    proj_trans_bounds(lstate.proj_ctx, crs, PJ_FWD, box_in.a_val, box_in.b_val, box_in.c_val, box_in.d_val,
			                      &box_out.a_val, &box_out.b_val, &box_out.c_val, &box_out.d_val, densify_pts);
			    return box_out;
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Execute (GEOMETRY)
	//------------------------------------------------------------------------------------------------------------------
	static void ExecuteGeometry(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = ProjFunctionLocalState::ResetAndGet(state);

		auto &func_expr = state.expr.Cast<BoundFunctionExpression>();
		const auto &info = func_expr.bind_info->Cast<BindData>();

		TernaryExecutor::Execute<string_t, string_t, string_t, string_t>(
		    args.data[0], args.data[1], args.data[2], result, args.size(),
		    [&](const string_t &input_geom, const string_t &source, const string_t &target) {
			    const auto source_str = source.GetString();
			    const auto target_str = target.GetString();

			    const auto crs = lstate.GetOrCreateProjection(source_str, target_str, info.normalize);

			    auto geom = lstate.Deserialize(input_geom);

			    sgl::ops::replace_vertices_xy(&geom, crs, [](void *arg, sgl::vertex_xy *vertex) {
				    const auto crs_ptr = static_cast<PJ *>(arg);
				    const auto transformed = proj_trans(crs_ptr, PJ_FWD, proj_coord(vertex->x, vertex->y, 0, 0)).xy;
				    vertex->x = transformed.x;
				    vertex->y = transformed.y;
			    });

			    return lstate.Serialize(result, geom);
		    });
	}

	//------------------------------------------------------------------------------------------------------------------
	// Documentation
	//------------------------------------------------------------------------------------------------------------------
	static constexpr auto DESCRIPTION = R"(
	Transforms a geometry between two coordinate systems

	The source and target coordinate systems can be specified using any format that the [PROJ library](https://proj.org) supports.

	The third optional `always_xy` parameter can be used to force the input and output geometries to be interpreted as having a [easting, northing] coordinate axis order regardless of what the source and target coordinate system definition says. This is particularly useful when transforming to/from the [WGS84/EPSG:4326](https://en.wikipedia.org/wiki/World_Geodetic_System) coordinate system (what most people think of when they hear "longitude"/"latitude" or "GPS coordinates"), which is defined as having a [latitude, longitude] axis order even though [longitude, latitude] is commonly used in practice (e.g. in [GeoJSON](https://tools.ietf.org/html/rfc7946)). More details available in the [PROJ documentation](https://proj.org/en/9.3/faq.html#why-is-the-axis-ordering-in-proj-not-consistent).

	DuckDB spatial vendors its own static copy of the PROJ database of coordinate systems, so if you have your own installation of PROJ on your system the available coordinate systems may differ to what's available in other GIS software.
	)";

	static constexpr auto EXAMPLE = R"(
	-- Transform a geometry from EPSG:4326 to EPSG:3857 (WGS84 to WebMercator)
	-- Note that since WGS84 is defined as having a [latitude, longitude] axis order
	-- we follow the standard and provide the input geometry using that axis order,
	-- but the output will be [easting, northing] because that is what's defined by
	-- WebMercator.

	SELECT ST_AsText(
	    ST_Transform(
	        st_point(52.373123, 4.892360),
	        'EPSG:4326',
	        'EPSG:3857'
	    )
	);
	----
	POINT (544615.0239773799 6867874.103539125)

	-- Alternatively, let's say we got our input point from e.g. a GeoJSON file,
	-- which uses WGS84 but with [longitude, latitude] axis order. We can use the
	-- `always_xy` parameter to force the input geometry to be interpreted as having
	-- a [northing, easting] axis order instead, even though the source coordinate
	-- reference system definition (WGS84) says otherwise.

	SELECT ST_AsText(
	    ST_Transform(
	        -- note the axis order is reversed here
	        st_point(4.892360, 52.373123),
	        'EPSG:4326',
	        'EPSG:3857',
	        always_xy := true
	    )
	);
	----
	POINT (544615.0239773799 6867874.103539125)
	)";

	//------------------------------------------------------------------------------------------------------------------
	// Register
	//------------------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Transform", [](ScalarFunctionBuilder &func) {
			func.AddVariant([&](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2D());
				variant.AddParameter("source_crs", LogicalType::VARCHAR);
				variant.AddParameter("target_crs", LogicalType::VARCHAR);
				variant.SetReturnType(GeoTypes::BOX_2D());

				variant.SetInit(ProjFunctionLocalState::Init);
				variant.SetBind(Bind);
				variant.SetFunction(ExecuteBox);
			});

			func.AddVariant([&](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("box", GeoTypes::BOX_2D());
				variant.AddParameter("source_crs", LogicalType::VARCHAR);
				variant.AddParameter("target_crs", LogicalType::VARCHAR);
				variant.AddParameter("always_xy", LogicalType::BOOLEAN);
				variant.SetReturnType(GeoTypes::BOX_2D());

				variant.SetInit(ProjFunctionLocalState::Init);
				variant.SetBind(Bind);
				variant.SetFunction(ExecuteBox);
			});

			func.AddVariant([&](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.AddParameter("source_crs", LogicalType::VARCHAR);
				variant.AddParameter("target_crs", LogicalType::VARCHAR);
				variant.SetReturnType(GeoTypes::POINT_2D());

				variant.SetInit(ProjFunctionLocalState::Init);
				variant.SetBind(Bind);
				variant.SetFunction(ExecutePoint);
			});

			func.AddVariant([&](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("point", GeoTypes::POINT_2D());
				variant.AddParameter("source_crs", LogicalType::VARCHAR);
				variant.AddParameter("target_crs", LogicalType::VARCHAR);
				variant.AddParameter("always_xy", LogicalType::BOOLEAN);
				variant.SetReturnType(GeoTypes::POINT_2D());

				variant.SetInit(ProjFunctionLocalState::Init);
				variant.SetBind(Bind);
				variant.SetFunction(ExecutePoint);
			});

			func.AddVariant([&](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("source_crs", LogicalType::VARCHAR);
				variant.AddParameter("target_crs", LogicalType::VARCHAR);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(ProjFunctionLocalState::Init);
				variant.SetBind(Bind);
				variant.SetFunction(ExecuteGeometry);
			});

			func.AddVariant([&](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("source_crs", LogicalType::VARCHAR);
				variant.AddParameter("target_crs", LogicalType::VARCHAR);
				variant.AddParameter("always_xy", LogicalType::BOOLEAN);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(ProjFunctionLocalState::Init);
				variant.SetBind(Bind);
				variant.SetFunction(ExecuteGeometry);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "conversion");
		});
	}
};

//######################################################################################################################
// Geodesic Functions
//######################################################################################################################

constexpr auto EARTH_A = 6378137;
constexpr auto EARTH_F = 1 / 298.257223563;

//======================================================================================================================
// Local State
//======================================================================================================================

struct GeodesicLocalState final : FunctionLocalState {

	ArenaAllocator arena;
	GeometryAllocator alloc;
	geod_geodesic geod = {};
	geod_polygon poly = {};
	double accum = 0;

	explicit GeodesicLocalState(ClientContext &context, bool is_line)
	    : arena(BufferAllocator::Get(context)), alloc(arena) {

		// Initialize the geodesic object for earth
		geod_init(&geod, EARTH_A, EARTH_F);
		geod_polygon_init(&poly, is_line ? 1 : 0);
	}

	static unique_ptr<FunctionLocalState> InitPolygon(ExpressionState &state, const BoundFunctionExpression &expr,
	                                                  FunctionData *bind_data) {
		return make_uniq<GeodesicLocalState>(state.GetContext(), false);
	}

	static unique_ptr<FunctionLocalState> InitLine(ExpressionState &state, const BoundFunctionExpression &expr,
	                                               FunctionData *bind_data) {
		return make_uniq<GeodesicLocalState>(state.GetContext(), true);
	}

	static GeodesicLocalState &ResetAndGet(ExpressionState &state) {
		auto &local_state = ExecuteFunctionState::GetFunctionState(state)->Cast<GeodesicLocalState>();
		local_state.arena.Reset();
		return local_state;
	}

	sgl::geometry Deserialize(const string_t &blob) {
		sgl::geometry geom;
		Serde::Deserialize(geom, arena, blob.GetDataUnsafe(), blob.GetSize());
		return geom;
	}
};

//======================================================================================================================
// ST_Area_Spheroid
//======================================================================================================================

struct ST_Area_Spheroid {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

		auto &lstate = GeodesicLocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [&](const string_t &input) {
			const auto geom = lstate.Deserialize(input);

			// Reset the state
			lstate.accum = 0;

			// Visit all polygons
			sgl::ops::visit_by_dimension(&geom, 2, &lstate, [](void *arg, const sgl::geometry *part) {
				if (part->get_type() != sgl::geometry_type::POLYGON) {
					return;
				}

				auto &sstate = *static_cast<GeodesicLocalState *>(arg);

				// Calculate the area of the polygon
				const auto tail = part->get_last_part();
				auto ring = tail;
				if (!ring) {
					return;
				}

				const auto head = ring->get_next();

				do {
					ring = ring->get_next();

					const auto vertex_count = ring->get_count();
					if (vertex_count < 4) {
						continue;
					}

					geod_polygon_clear(&sstate.poly);

					// Dont add the last vertex
					for (auto i = 0; i < vertex_count - 1; i++) {
						const auto vertex = ring->get_vertex_xy(i);
						geod_polygon_addpoint(&sstate.geod, &sstate.poly, vertex.x, vertex.y);
					}

					double area = 0;
					geod_polygon_compute(&sstate.geod, &sstate.poly, 0, 1, &area, nullptr);

					if (ring == head) {
						sstate.accum += std::abs(area);
					} else {
						sstate.accum -= std::abs(area);
					}
				} while (ring != tail);
			});

			return lstate.accum;
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Area_Spheroid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(GeodesicLocalState::InitPolygon);
				variant.SetFunction(Execute);
			});
		});
	}
};

//======================================================================================================================
// ST_Perimeter_Spheroid
//======================================================================================================================
struct ST_Perimeter_Spheroid {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

		auto &lstate = GeodesicLocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [&](const string_t &input) {
			const auto geom = lstate.Deserialize(input);

			// Reset the state
			lstate.accum = 0;

			// Visit all polygons
			sgl::ops::visit_by_dimension(&geom, 2, &lstate, [](void *arg, const sgl::geometry *part) {
				if (part->get_type() != sgl::geometry_type::POLYGON) {
					return;
				}

				auto &sstate = *static_cast<GeodesicLocalState *>(arg);

				// Calculate the perimeter of the polygon
				const auto tail = part->get_last_part();
				auto ring = tail;
				if (!ring) {
					return;
				}
				do {
					ring = ring->get_next();

					const auto vertex_count = ring->get_count();
					if (vertex_count < 4) {
						continue;
					}

					geod_polygon_clear(&sstate.poly);

					// Dont add the last vertex
					for (auto i = 0; i < vertex_count - 1; i++) {
						const auto vertex = ring->get_vertex_xy(i);
						geod_polygon_addpoint(&sstate.geod, &sstate.poly, vertex.x, vertex.y);
					}

					double perimeter = 0;
					geod_polygon_compute(&sstate.geod, &sstate.poly, 0, 1, nullptr, &perimeter);
					// Add the perimeter of the ring
					sstate.accum += perimeter;

				} while (ring != tail);
			});

			return lstate.accum;
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Perimeter_Spheroid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(GeodesicLocalState::InitPolygon);
				variant.SetFunction(Execute);
			});
		});
	}
};

//======================================================================================================================
// ST_Length_Spheroid
//======================================================================================================================

struct ST_Length_Spheroid {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {

		auto &lstate = GeodesicLocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [&](const string_t &input) {
			const auto geom = lstate.Deserialize(input);

			// Reset the state
			lstate.accum = 0;

			// Visit all polygons
			sgl::ops::visit_by_dimension(&geom, 1, &lstate, [](void *arg, const sgl::geometry *part) {
				if (part->get_type() != sgl::geometry_type::LINESTRING) {
					return;
				}

				auto &sstate = *static_cast<GeodesicLocalState *>(arg);

				const auto vertex_count = part->get_count();
				if (vertex_count < 2) {
					return;
				}

				geod_polygon_clear(&sstate.poly);

				for (auto i = 0; i < vertex_count; i++) {
					const auto vertex = part->get_vertex_xy(i);
					geod_polygon_addpoint(&sstate.geod, &sstate.poly, vertex.x, vertex.y);
				}

				// Calculate the length of the linestring
				double length = 0;
				geod_polygon_compute(&sstate.geod, &sstate.poly, 0, 1, nullptr, &length);

				sstate.accum += length;
			});

			return lstate.accum;
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Length_Spheroid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(GeodesicLocalState::InitLine);
				variant.SetFunction(Execute);
			});
		});
	}
};

//======================================================================================================================
// ST_Distance_Spheroid
//======================================================================================================================

struct ST_Distance_Spheroid {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		using POINT_TYPE = StructTypeBinary<double, double>;
		using DISTANCE_TYPE = PrimitiveType<double>;

		geod_geodesic geod = {};
		geod_init(&geod, EARTH_A, EARTH_F);

		GenericExecutor::ExecuteBinary<POINT_TYPE, POINT_TYPE, DISTANCE_TYPE>(
		    args.data[0], args.data[1], result, args.size(), [&](const POINT_TYPE &p1, const POINT_TYPE &p2) {
			    double distance;
			    geod_inverse(&geod, p1.a_val, p1.b_val, p2.a_val, p2.b_val, &distance, nullptr, nullptr);
			    return distance;
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Distance_Spheroid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("p1", GeoTypes::POINT_2D());
				variant.AddParameter("p2", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(Execute);
			});
		});
	}
};

//======================================================================================================================
// ST_DWithin_Spheroid
//======================================================================================================================

struct ST_DWithin_Spheroid {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		using POINT_TYPE = StructTypeBinary<double, double>;
		using DISTANCE_TYPE = PrimitiveType<double>;
		using BOOL_TYPE = PrimitiveType<bool>;

		geod_geodesic geod = {};
		geod_init(&geod, EARTH_A, EARTH_F);

		GenericExecutor::ExecuteTernary<POINT_TYPE, POINT_TYPE, DISTANCE_TYPE, BOOL_TYPE>(
		    args.data[0], args.data[1], args.data[2], result, args.size(),
		    [&](const POINT_TYPE &p1, const POINT_TYPE &p2, const DISTANCE_TYPE &limit) {
			    double distance;
			    geod_inverse(&geod, p1.a_val, p1.b_val, p2.a_val, p2.b_val, &distance, nullptr, nullptr);
			    return distance <= limit.val;
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_DWithin_Spheroid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("p1", GeoTypes::POINT_2D());
				variant.AddParameter("p2", GeoTypes::POINT_2D());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetFunction(Execute);
			});
		});
	}
};

} // namespace

//######################################################################################################################
// Module Registration
//######################################################################################################################
void RegisterProjModule(DatabaseInstance &db) {

	// Register the VFS for the proj.db database
	ProjModule::RegisterVFS(db);

	// Coordinate Transform Function
	ST_Transform::Register(db);

	// Geodesic Functions
	ST_Area_Spheroid::Register(db);
	ST_Perimeter_Spheroid::Register(db);
	ST_Length_Spheroid::Register(db);
	ST_Distance_Spheroid::Register(db);
	ST_DWithin_Spheroid::Register(db);
}

} // namespace duckdb