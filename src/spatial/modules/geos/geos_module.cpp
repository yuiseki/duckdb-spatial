#include "spatial/modules/geos/geos_module.hpp"
#include "spatial/modules/geos/geos_geometry.hpp"
#include "spatial/modules/geos/geos_serde.hpp"
#include "spatial/spatial_types.hpp"
#include "spatial/util/function_builder.hpp"

#include "duckdb/common/vector_operations/senary_executor.hpp"
#include "duckdb/main/extension_util.hpp"

namespace duckdb {

//------------------------------------------------------------------------------
// Local State
//------------------------------------------------------------------------------

namespace {

class LocalState final : public FunctionLocalState {
public:
	static unique_ptr<FunctionLocalState> Init(ExpressionState &state, const BoundFunctionExpression &expr,
	                                           FunctionData *bind_data) {
		return make_uniq<LocalState>(state.GetContext());
	}

	static LocalState &ResetAndGet(ExpressionState &state) {
		auto &local_state = ExecuteFunctionState::GetFunctionState(state)->Cast<LocalState>();
		return local_state;
	}

	GEOSContextHandle_t GetContext() const {
		return ctx;
	}

	GeosGeometry Deserialize(const string_t &blob) const;
	string_t Serialize(Vector &result, const GeosGeometry &geom) const;

	// Most GEOS functions do not use an arena, so just use the default allocator
	explicit LocalState(ClientContext &context) {
		ctx = GEOS_init_r();

		GEOSContext_setErrorMessageHandler_r(
		    ctx, [](const char *message, void *) { throw InvalidInputException(message); }, nullptr);
	}

	~LocalState() override {
		GEOS_finish_r(ctx);
	}

private:
	GEOSContextHandle_t ctx;
};

string_t LocalState::Serialize(Vector &result, const GeosGeometry &geom) const {
	// Get the size of the serialized geometry
	const auto raw = geom.get_raw();
	const auto size = GeosSerde::GetRequiredSize(ctx, raw);

	// Allocate a blob of the correct size
	auto blob = StringVector::EmptyString(result, size);
	const auto ptr = blob.GetDataWriteable();

	// Serialize the geometry into the blob
	GeosSerde::Serialize(ctx, raw, ptr, size);

	// Finalize and return the blob
	blob.Finalize();
	return blob;
}

GeosGeometry LocalState::Deserialize(const string_t &blob) const {
	const auto blob_ptr = blob.GetData();
	const auto blob_len = blob.GetSize();

	const auto geom = GeosSerde::Deserialize(ctx, blob_ptr, blob_len);

	if (geom == nullptr) {
		throw InvalidInputException("Could not deserialize geometry");
	}

	return GeosGeometry(ctx, geom);
}

} // namespace

//------------------------------------------------------------------------------
// Base Functions
//------------------------------------------------------------------------------

namespace {

template <class IMPL, class RETURN_TYPE = bool>
class SymmetricPreparedBinaryFunction {
public:
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		auto &lhs_vec = args.data[0];
		auto &rhs_vec = args.data[1];

		const auto lhs_is_const =
		    lhs_vec.GetVectorType() == VectorType::CONSTANT_VECTOR && !ConstantVector::IsNull(lhs_vec);
		const auto rhs_is_const =
		    rhs_vec.GetVectorType() == VectorType::CONSTANT_VECTOR && !ConstantVector::IsNull(rhs_vec);

		if (lhs_is_const && rhs_is_const) {
			// Both are const, just execute once
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
			const auto &lhs_blob = ConstantVector::GetData<string_t>(lhs_vec)[0];
			const auto &rhs_blob = ConstantVector::GetData<string_t>(rhs_vec)[0];
			const auto lhs_geom = lstate.Deserialize(lhs_blob);
			const auto rhs_geom = lstate.Deserialize(rhs_blob);
			ConstantVector::GetData<RETURN_TYPE>(result)[0] = IMPL::ExecutePredicateNormal(lhs_geom, rhs_geom);

		} else if (lhs_is_const != rhs_is_const) {
			// One of the two is const, prepare the const one and execute on the non-const one
			auto &const_vec = lhs_is_const ? lhs_vec : rhs_vec;
			auto &probe_vec = lhs_is_const ? rhs_vec : lhs_vec;

			const auto &const_blob = ConstantVector::GetData<string_t>(const_vec)[0];
			const auto const_geom = lstate.Deserialize(const_blob);
			const auto const_prep = const_geom.get_prepared();

			UnaryExecutor::Execute<string_t, RETURN_TYPE>(
			    probe_vec, result, args.size(), [&](const string_t &probe_blob) {
				    const auto probe_geom = lstate.Deserialize(probe_blob);
				    return IMPL::ExecutePredicatePrepared(const_prep, probe_geom);
			    });
		} else {
			// Both are non-const, just execute normally
			BinaryExecutor::Execute<string_t, string_t, RETURN_TYPE>(
			    lhs_vec, rhs_vec, result, args.size(), [&](const string_t &lhs_blob, const string_t &rhs_blob) {
				    const auto lhs = lstate.Deserialize(lhs_blob);
				    const auto rhs = lstate.Deserialize(rhs_blob);
				    return IMPL::ExecutePredicateNormal(lhs, rhs);
			    });
		}
	}
};

template <class IMPL, class RETURN_TYPE = bool>
class AsymmetricPreparedBinaryFunction {
public:
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		auto &lhs_vec = args.data[0];
		auto &rhs_vec = args.data[1];

		const auto lhs_is_const =
		    lhs_vec.GetVectorType() == VectorType::CONSTANT_VECTOR && !ConstantVector::IsNull(lhs_vec);
		const auto rhs_is_const =
		    rhs_vec.GetVectorType() == VectorType::CONSTANT_VECTOR && !ConstantVector::IsNull(rhs_vec);

		if (lhs_is_const && rhs_is_const) {
			// Both are const, just execute once
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
			const auto &lhs_blob = ConstantVector::GetData<string_t>(lhs_vec)[0];
			const auto &rhs_blob = ConstantVector::GetData<string_t>(rhs_vec)[0];
			const auto lhs_geom = lstate.Deserialize(lhs_blob);
			const auto rhs_geom = lstate.Deserialize(rhs_blob);
			ConstantVector::GetData<RETURN_TYPE>(result)[0] = IMPL::ExecutePredicateNormal(lhs_geom, rhs_geom);

		} else if (lhs_is_const) {
			// Prepare the left const and run on the non-const right
			// Because this predicate is not symmetric, we can't just swap the two, so we only prepare the left
			const auto lhs_blob = ConstantVector::GetData<string_t>(lhs_vec)[0];
			const auto lhs_geom = lstate.Deserialize(lhs_blob);
			const auto lhs_prep = lhs_geom.get_prepared();

			UnaryExecutor::Execute<string_t, RETURN_TYPE>(rhs_vec, result, args.size(), [&](const string_t &rhs_blob) {
				const auto rhs_geom = lstate.Deserialize(rhs_blob);
				return IMPL::ExecutePredicatePrepared(lhs_prep, rhs_geom);
			});
		} else {
			// Both are non-const, just execute normally
			BinaryExecutor::Execute<string_t, string_t, RETURN_TYPE>(
			    lhs_vec, rhs_vec, result, args.size(), [&](const string_t &lhs_blob, const string_t &rhs_blob) {
				    const auto lhs = lstate.Deserialize(lhs_blob);
				    const auto rhs = lstate.Deserialize(rhs_blob);
				    return IMPL::ExecutePredicateNormal(lhs, rhs);
			    });
		}
	}
};

} // namespace

//------------------------------------------------------------------------------
// Functions
//------------------------------------------------------------------------------

namespace {

struct ST_Boundary {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
		    args.data[0], result, args.size(), [&](const string_t &geom_blob, ValidityMask &mask, idx_t row_idx) {
			    const auto geom = lstate.Deserialize(geom_blob);
			    if (geom.type() == GEOS_GEOMETRYCOLLECTION) {
				    mask.SetInvalid(row_idx);
				    return string_t {};
			    }
			    const auto boundary = geom.get_boundary();

			    return lstate.Serialize(result, boundary);
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Boundary", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the boundary of a geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Buffer {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, double, string_t>(args.data[0], args.data[1], result, args.size(),
		                                                    [&](const string_t &blob, double radius) {
			                                                    const auto geom = lstate.Deserialize(blob);
			                                                    const auto buffer = geom.get_buffer(radius, 8);
			                                                    return lstate.Serialize(result, buffer);
		                                                    });
	}

	static void ExecuteWithSegments(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		TernaryExecutor::Execute<string_t, double, int32_t, string_t>(
		    args.data[0], args.data[1], args.data[2], result, args.size(),
		    [&](const string_t &blob, double radius, int32_t segments) {
			    const auto geom = lstate.Deserialize(blob);
			    const auto buffer = geom.get_buffer(radius, segments);
			    return lstate.Serialize(result, buffer);
		    });
	}

	template <class T>
	static T TryParseStringArgument(const char *name, const vector<string> &keys, const vector<T> &values,
	                                const string_t &arg) {
		D_ASSERT(keys.size() == values.size());
		for (idx_t i = 0; i < keys.size(); i++) {
			if (StringUtil::CIEquals(keys[i], arg.GetString())) {
				return values[i];
			}
		}

		auto candidates = StringUtil::Join(keys, ", ");
		throw InvalidInputException("Unknown %s: '%s', accepted inputs: %s", name, arg.GetString().c_str(),
		                            candidates.c_str());
	}

	static void ExecuteWithStyle(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		SenaryExecutor::Execute<string_t, double, int32_t, string_t, string_t, double, string_t>(
		    args, result,
		    [&](const string_t &blob, double radius, int32_t segments, const string_t &cap_style_str,
		        const string_t &join_style_str, double mitre_limit) {
			    const auto geom = lstate.Deserialize(blob);
			    const auto cap_style = TryParseStringArgument<GEOSBufCapStyles>(
			        "cap style", {"CAP_ROUND", "CAP_FLAT", "CAP_SQUARE"},
			        {GEOSBUF_CAP_ROUND, GEOSBUF_CAP_FLAT, GEOSBUF_CAP_SQUARE}, cap_style_str);

			    const auto join_style = TryParseStringArgument<GEOSBufJoinStyles>(
			        "join style", {"JOIN_ROUND", "JOIN_MITRE", "JOIN_BEVEL"},
			        {GEOSBUF_JOIN_ROUND, GEOSBUF_JOIN_MITRE, GEOSBUF_JOIN_BEVEL}, join_style_str);

			    const auto buffer = geom.get_buffer_style(radius, segments, cap_style, join_style, mitre_limit);
			    return lstate.Serialize(result, buffer);
			    ;
		    });
	}

	static constexpr auto DESCRIPTION = R"(
	    Returns a buffer around the input geometry at the target distance

	    `geom` is the input geometry.

	    `distance` is the target distance for the buffer, using the same units as the input geometry.

	    `num_triangles` represents how many triangles that will be produced to approximate a quarter circle. The larger the number, the smoother the resulting geometry. The default value is 8.

	    `join_style` must be one of "JOIN_ROUND", "JOIN_MITRE", "JOIN_BEVEL". This parameter is case-insensitive.

	    `cap_style` must be one of "CAP_ROUND", "CAP_FLAT", "CAP_SQUARE". This parameter is case-insensitive.

	    `mitre_limit` only applies when `join_style` is "JOIN_MITRE". It is the ratio of the distance from the corner to the mitre point to the corner radius. The default value is 1.0.

	    This is a planar operation and will not take into account the curvature of the earth.
	)";
	static constexpr auto EXAMPLE = "";

	static void Register(DatabaseInstance &db) {

		FunctionBuilder::RegisterScalar(db, "ST_Buffer", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("distance", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("distance", LogicalType::DOUBLE);
				variant.AddParameter("num_triangles", LogicalType::INTEGER);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteWithSegments);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("distance", LogicalType::DOUBLE);
				variant.AddParameter("num_triangles", LogicalType::INTEGER);
				variant.AddParameter("cap_style", LogicalType::VARCHAR);
				variant.AddParameter("join_style", LogicalType::VARCHAR);
				variant.AddParameter("mitre_limit", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteWithStyle);
			});

			func.SetDescription(DESCRIPTION);
			func.SetExample(EXAMPLE);

			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Centroid {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			const auto centroid = geom.get_centroid();
			return lstate.Serialize(result, centroid);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Centroid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the centroid of a geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Contains : AsymmetricPreparedBinaryFunction<ST_Contains> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.contains(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.contains(rhs);
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Contains", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the first geometry contains the second geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_ContainsProperly : AsymmetricPreparedBinaryFunction<ST_ContainsProperly> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		// We have no choice but to prepare the left geometry
		const auto lhs_prep = lhs.get_prepared();
		return lhs_prep.contains_properly(rhs);
	}

	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.contains_properly(rhs);
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_ContainsProperly", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the first geometry contains the second geometry properly");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_ConvexHull {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			const auto hull = geom.get_convex_hull();
			return lstate.Serialize(result, hull);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_ConvexHull", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the convex hull enclosing the geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_CoveredBy : AsymmetricPreparedBinaryFunction<ST_CoveredBy> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.covered_by(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.covered_by(rhs);
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_CoveredBy", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the first geometry is covered by the second geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_Covers : AsymmetricPreparedBinaryFunction<ST_Covers> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.covers(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.covers(rhs);
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Covers", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the first geometry covers the second geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_Crosses : SymmetricPreparedBinaryFunction<ST_Crosses> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.crosses(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.crosses(rhs);
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Crosses", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometries cross each other");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_Difference {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, string_t, string_t>(args.data[0], args.data[1], result, args.size(),
		                                                      [&](const string_t &lhs_blob, const string_t &rhs_blob) {
			                                                      const auto lhs = lstate.Deserialize(lhs_blob);
			                                                      const auto rhs = lstate.Deserialize(rhs_blob);
			                                                      const auto difference = lhs.get_difference(rhs);
			                                                      return lstate.Serialize(result, difference);
		                                                      });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Difference", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the difference between two geometries");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Disjoint : SymmetricPreparedBinaryFunction<ST_Disjoint> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.disjoint(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.disjoint(rhs);
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Disjoint", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometries are disjoint");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_Distance : SymmetricPreparedBinaryFunction<ST_Distance, double> {
	static double ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.distance_to(rhs);
	}
	static double ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.distance_to(rhs);
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Distance", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::DOUBLE);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the planar distance between two geometries");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "measurement");
		});
	}
};

struct ST_DistanceWithin {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		// Because this takes an extra argument, we cant reuse the SymmetricPreparedBinary...

		const auto &lstate = LocalState::ResetAndGet(state);

		auto &lhs_vec = args.data[0];
		auto &rhs_vec = args.data[1];
		auto &arg_vec = args.data[2];

		const auto lhs_is_const =
		    lhs_vec.GetVectorType() == VectorType::CONSTANT_VECTOR && !ConstantVector::IsNull(lhs_vec);
		const auto rhs_is_const =
		    rhs_vec.GetVectorType() == VectorType::CONSTANT_VECTOR && !ConstantVector::IsNull(rhs_vec);
		const auto arg_is_const =
		    arg_vec.GetVectorType() == VectorType::CONSTANT_VECTOR && !ConstantVector::IsNull(arg_vec);

		if (lhs_is_const && rhs_is_const && arg_is_const) {
			// Both geometries (and the argument) are constant, so only execute it once
			result.SetVectorType(VectorType::CONSTANT_VECTOR);
			const auto &lhs_blob = ConstantVector::GetData<string_t>(lhs_vec)[0];
			const auto &rhs_blob = ConstantVector::GetData<string_t>(rhs_vec)[0];
			const auto &arg_dist = ConstantVector::GetData<double>(arg_vec)[0];
			const auto lhs_geom = lstate.Deserialize(lhs_blob);
			const auto rhs_geom = lstate.Deserialize(rhs_blob);

			ConstantVector::GetData<bool>(result)[0] = lhs_geom.distance_within(rhs_geom, arg_dist);
		} else if (lhs_is_const && rhs_is_const && !arg_is_const) {
			// The geometries are constant, but the distance is not, prepare the larger one and execute unary

			const auto &lhs_blob = ConstantVector::GetData<string_t>(lhs_vec)[0];
			const auto &rhs_blob = ConstantVector::GetData<string_t>(rhs_vec)[0];

			const auto rhs_bigger = rhs_blob.GetSize() > lhs_blob.GetSize();

			const auto large_geom = rhs_bigger ? lstate.Deserialize(rhs_blob) : lstate.Deserialize(lhs_blob);
			const auto probe_geom = rhs_bigger ? lstate.Deserialize(lhs_blob) : lstate.Deserialize(rhs_blob);

			// PreparedDistanceWithin only works if one is prepared. so just choose the larger one
			const auto prep_geom = large_geom.get_prepared();

			UnaryExecutor::Execute<double, bool>(arg_vec, result, args.size(), [&](const double arg_dist) {
				return prep_geom.distance_within(probe_geom, arg_dist);
			});

		} else if (lhs_is_const != rhs_is_const) {
			// One of the two is const, prepare the const one and execute on the non-const one
			auto &const_vec = lhs_is_const ? lhs_vec : rhs_vec;
			auto &probe_vec = lhs_is_const ? rhs_vec : lhs_vec;

			const auto &const_blob = ConstantVector::GetData<string_t>(const_vec)[0];
			const auto const_geom = lstate.Deserialize(const_blob);
			const auto const_prep = const_geom.get_prepared();

			BinaryExecutor::Execute<string_t, double, bool>(probe_vec, arg_vec, result, args.size(),
			                                                [&](const string_t &probe_blob, double distance) {
				                                                const auto probe_geom = lstate.Deserialize(probe_blob);
				                                                return const_prep.distance_within(probe_geom, distance);
			                                                });
		} else {
			// Both are non-const, just execute normally
			TernaryExecutor::Execute<string_t, string_t, double, bool>(
			    lhs_vec, rhs_vec, arg_vec, result, args.size(),
			    [&](const string_t &lhs_blob, const string_t &rhs_blob, double distance) {
				    const auto lhs = lstate.Deserialize(lhs_blob);
				    const auto rhs = lstate.Deserialize(rhs_blob);
				    return lhs.distance_within(rhs, distance);
			    });
		}
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_DWithin", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.AddParameter("distance", LogicalType::DOUBLE);
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription(R"(
				Returns if two geometries are within a target distance of each-other
			)");

			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_Equals {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, string_t, bool>(args.data[0], args.data[1], result, args.size(),
		                                                  [&](const string_t &lhs_blob, const string_t &rhs_blob) {
			                                                  const auto lhs = lstate.Deserialize(lhs_blob);
			                                                  const auto rhs = lstate.Deserialize(rhs_blob);
			                                                  return lhs.equals(rhs);
		                                                  });
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Equals", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometries are equal");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_Envelope {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			const auto intersection = geom.get_envelope();
			return lstate.Serialize(result, intersection);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Envelope", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the envelope of a geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Intersection {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, string_t, string_t>(args.data[0], args.data[1], result, args.size(),
		                                                      [&](const string_t &lhs_blob, const string_t &rhs_blob) {
			                                                      const auto lhs = lstate.Deserialize(lhs_blob);
			                                                      const auto rhs = lstate.Deserialize(rhs_blob);
			                                                      const auto intersection = lhs.get_intersection(rhs);
			                                                      return lstate.Serialize(result, intersection);
		                                                      });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Intersection", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the intersection of two geometries");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Intersects : SymmetricPreparedBinaryFunction<ST_Intersects> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.intersects(rhs);
	}

	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.intersects(rhs);
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Intersects", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometries intersect");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_IsRing {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			return geom.is_ring();
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_IsRing", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometry is a ring (both ST_IsClosed and ST_IsSimple).");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_IsSimple {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			return geom.is_simple();
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_IsSimple", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometry is simple");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_IsValid {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, bool>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			// GEOS can only construct geometries with a valid amount of vertices.
			// So if deserialization fails, it cant be valid
			try {
				const auto geom = lstate.Deserialize(geom_blob);
				return geom.is_valid();
			} catch (...) {
				return false;
			}
		});
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_IsValid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometry is valid");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "property");
		});
	}
};

struct ST_LineMerge {

	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(),
		                                           [&](const string_t &geometry_blob) {
			                                           const auto geometry = lstate.Deserialize(geometry_blob);
			                                           const auto merged = geometry.get_linemerged(false);
			                                           return lstate.Serialize(result, merged);
		                                           });
	}

	static void ExecuteWithDirection(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		BinaryExecutor::Execute<string_t, bool, string_t>(args.data[0], args.data[1], result, args.size(),
		                                                  [&](const string_t &geometry_blob, bool preserve_direction) {
			                                                  const auto geometry = lstate.Deserialize(geometry_blob);
			                                                  const auto merged =
			                                                      geometry.get_linemerged(preserve_direction);
			                                                  return lstate.Serialize(result, merged);
		                                                  });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_LineMerge", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());
				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("preserve_direction", LogicalType::BOOLEAN);
				variant.SetReturnType(GeoTypes::GEOMETRY());
				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteWithDirection);
			});

			func.SetDescription(R"("Merges" the input line geometry, optionally taking direction into account.)");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_MakeValid {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			const auto valid = geom.get_made_valid();
			return lstate.Serialize(result, valid);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_MakeValid", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns a valid representation of the geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Normalize {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			geom.normalize_in_place();
			return lstate.Serialize(result, geom);
		});
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Normalize", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns a normalized representation of the geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Overlaps : SymmetricPreparedBinaryFunction<ST_Overlaps> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.overlaps(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.overlaps(rhs);
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Overlaps", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometries overlap");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_PointOnSurface {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			const auto point = geom.get_point_on_surface();
			return lstate.Serialize(result, point);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_PointOnSurface", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns a point guaranteed to lie on the surface of the geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_ReducePrecision {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, double, string_t>(
		    args.data[0], args.data[1], result, args.size(), [&](const string_t &geom_blob, double precision) {
			    const auto geom = lstate.Deserialize(geom_blob);
			    const auto reduced = geom.get_reduced_precision(precision);
			    return lstate.Serialize(result, reduced);
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_ReducePrecision", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("precision", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the geometry with all vertices reduced to the given precision");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_RemoveRepeatedPoints {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			const auto reduced = geom.get_without_repeated_points(0);
			return lstate.Serialize(result, reduced);
		});
	}

	static void ExecuteWithTolerance(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, double, string_t>(
		    args.data[0], args.data[1], result, args.size(), [&](const string_t &geom_blob, double tolerance) {
			    const auto geom = lstate.Deserialize(geom_blob);
			    const auto reduced = geom.get_without_repeated_points(tolerance);
			    return lstate.Serialize(result, reduced);
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_RemoveRepeatedPoints", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("tolerance", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(ExecuteWithTolerance);
			});

			func.SetDescription("Returns the geometry with repeated points removed");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Reverse {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), [&](const string_t &geom_blob) {
			const auto geom = lstate.Deserialize(geom_blob);
			const auto reversed = geom.get_reversed();
			return lstate.Serialize(result, reversed);
		});
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Reverse", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the geometry with the order of its vertices reversed");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_ShortestLine {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		BinaryExecutor::Execute<string_t, string_t, string_t>(args.data[0], args.data[1], result, args.size(),
		                                                      [&](const string_t &lhs_blob, const string_t &rhs_blob) {
			                                                      const auto lhs = lstate.Deserialize(lhs_blob);
			                                                      const auto rhs = lstate.Deserialize(rhs_blob);
			                                                      const auto line = lhs.get_shortest_line(rhs);
			                                                      return lstate.Serialize(result, line);
		                                                      });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_ShortestLine", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the shortest line between two geometries");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "measurement");
		});
	}
};

struct ST_Simplify {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		BinaryExecutor::Execute<string_t, double, string_t>(args.data[0], args.data[1], result, args.size(),
		                                                    [&](const string_t &geom_blob, double tolerance) {
			                                                    const auto geom = lstate.Deserialize(geom_blob);
			                                                    const auto simplified = geom.get_simplified(tolerance);
			                                                    return lstate.Serialize(result, simplified);
		                                                    });
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Simplify", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("tolerance", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns a simplified version of the geometry");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_SimplifyPreserveTopology {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);
		BinaryExecutor::Execute<string_t, double, string_t>(
		    args.data[0], args.data[1], result, args.size(), [&](const string_t &geom_blob, double tolerance) {
			    const auto geom = lstate.Deserialize(geom_blob);
			    const auto simplified = geom.get_simplified_topo(tolerance);
			    return lstate.Serialize(result, simplified);
		    });
	}

	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_SimplifyPreserveTopology", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom", GeoTypes::GEOMETRY());
				variant.AddParameter("tolerance", LogicalType::DOUBLE);
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns a simplified version of the geometry that preserves topology");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Touches : SymmetricPreparedBinaryFunction<ST_Touches> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.touches(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.touches(rhs);
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Touches", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the geometries touch");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

struct ST_Union {
	static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
		const auto &lstate = LocalState::ResetAndGet(state);

		BinaryExecutor::Execute<string_t, string_t, string_t>(args.data[0], args.data[1], result, args.size(),
		                                                      [&](const string_t &lhs_blob, const string_t &rhs_blob) {
			                                                      const auto lhs = lstate.Deserialize(lhs_blob);
			                                                      const auto rhs = lstate.Deserialize(rhs_blob);
			                                                      const auto unioned = lhs.get_union(rhs);
			                                                      return lstate.Serialize(result, unioned);
		                                                      });
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Union", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(GeoTypes::GEOMETRY());

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns the union of two geometries");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "construction");
		});
	}
};

struct ST_Within : AsymmetricPreparedBinaryFunction<ST_Within> {
	static bool ExecutePredicateNormal(const GeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.within(rhs);
	}
	static bool ExecutePredicatePrepared(const PreparedGeosGeometry &lhs, const GeosGeometry &rhs) {
		return lhs.within(rhs);
	}
	static void Register(DatabaseInstance &db) {
		FunctionBuilder::RegisterScalar(db, "ST_Within", [](ScalarFunctionBuilder &func) {
			func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
				variant.AddParameter("geom1", GeoTypes::GEOMETRY());
				variant.AddParameter("geom2", GeoTypes::GEOMETRY());
				variant.SetReturnType(LogicalType::BOOLEAN);

				variant.SetInit(LocalState::Init);
				variant.SetFunction(Execute);
			});

			func.SetDescription("Returns true if the first geometry is within the second");
			func.SetTag("ext", "spatial");
			func.SetTag("category", "relation");
		});
	}
};

//######################################################################################################################
// Aggregate Functions
//######################################################################################################################

//======================================================================================================================
// Base GEOS-based unary aggregate
//======================================================================================================================
struct GeosUnaryAggState {
	GEOSGeometry *geom = nullptr;
	GEOSContextHandle_t context = nullptr;
};

struct GeosUnaryAggFunction {

	// Serialize a GEOS geometry
	static string_t Serialize(const GEOSContextHandle_t context, Vector &result, const GEOSGeometry *geom) {
		D_ASSERT(geom);
		const auto size = GeosSerde::GetRequiredSize(context, geom);
		auto blob = StringVector::EmptyString(result, size);
		const auto ptr = blob.GetDataWriteable();

		// Serialize the geometry
		GeosSerde::Serialize(context, geom, ptr, size);

		blob.Finalize();
		return blob;
	}

	// Deserialize a GEOS geometry
	static GEOSGeometry *Deserialize(const GEOSContextHandle_t context, const string_t &blob) {
		const auto ptr = blob.GetData();
		const auto size = blob.GetSize();

		return GeosSerde::Deserialize(context, ptr, size);
	}

	template <class STATE>
	static void Initialize(STATE &state) {
		state.geom = nullptr;
		state.context = GEOS_init_r();
	}

	template <class STATE, class OP>
	static void Combine(const STATE &source, STATE &target, AggregateInputData &data) {
		if (!source.geom) {
			return;
		}
		if (!target.geom) {
			target.geom = GEOSGeom_clone_r(target.context, source.geom);
			return;
		}
		auto curr = target.geom;
		target.geom = OP::Merge(target.context, curr, source.geom);
		GEOSGeom_destroy_r(target.context, curr);
	}

	template <class INPUT_TYPE, class STATE, class OP>
	static void Operation(STATE &state, const INPUT_TYPE &input, AggregateUnaryInput &) {
		if (!state.geom) {
			state.geom = Deserialize(state.context, input);
		} else {
			auto next = Deserialize(state.context, input);
			auto curr = state.geom;
			state.geom = OP::Merge(state.context, curr, next);
			GEOSGeom_destroy_r(state.context, next);
			GEOSGeom_destroy_r(state.context, curr);
		}
	}

	template <class INPUT_TYPE, class STATE, class OP>
	static void ConstantOperation(STATE &state, const INPUT_TYPE &input, AggregateUnaryInput &, idx_t) {
		// There is no point in doing anything else, intersection and union is idempotent
		if (!state.geom) {
			state.geom = Deserialize(state.context, input);
		}
	}

	template <class T, class STATE>
	static void Finalize(STATE &state, T &target, AggregateFinalizeData &finalize_data) {
		if (!state.geom) {
			finalize_data.ReturnNull();
		} else {
			target = Serialize(state.context, finalize_data.result, state.geom);
		}
	}

	template <class STATE>
	static void Destroy(STATE &state, AggregateInputData &) {
		if (state.geom) {
			GEOSGeom_destroy_r(state.context, state.geom);
			state.geom = nullptr;
		}
		if (state.context) {
			GEOS_finish_r(state.context);
			state.context = nullptr;
		}
	}

	static bool IgnoreNull() {
		return true;
	}
};

//======================================================================================================================
// ST_Union_Agg
//======================================================================================================================

struct ST_Union_Agg : GeosUnaryAggFunction {
	static GEOSGeometry *Merge(const GEOSContextHandle_t context, const GEOSGeometry *curr, const GEOSGeometry *next) {
		return GEOSUnion_r(context, curr, next);
	}

	static void Register(DatabaseInstance &db) {
		auto func = AggregateFunction::UnaryAggregateDestructor<GeosUnaryAggState, string_t, string_t, ST_Union_Agg>(
		    GeoTypes::GEOMETRY(), GeoTypes::GEOMETRY());
		func.name = "ST_Union_Agg";

		ExtensionUtil::RegisterFunction(db, func);
	}
};

//======================================================================================================================
// ST_Intersection_Agg
//======================================================================================================================

struct ST_Intersection_Agg : GeosUnaryAggFunction {
	static GEOSGeometry *Merge(const GEOSContextHandle_t context, const GEOSGeometry *curr, const GEOSGeometry *next) {
		return GEOSIntersection_r(context, curr, next);
	}

	static void Register(DatabaseInstance &db) {
		auto func =
		    AggregateFunction::UnaryAggregateDestructor<GeosUnaryAggState, string_t, string_t, ST_Intersection_Agg>(
		        GeoTypes::GEOMETRY(), GeoTypes::GEOMETRY());
		func.name = "ST_Intersection_Agg";

		ExtensionUtil::RegisterFunction(db, func);
	}
};

} // namespace

//######################################################################################################################
// Register Module
//######################################################################################################################

void RegisterGEOSModule(DatabaseInstance &db) {

	// Scalar Functions
	ST_Boundary::Register(db);
	ST_Buffer::Register(db);
	ST_Centroid::Register(db);
	ST_Contains::Register(db);
	ST_ContainsProperly::Register(db);
	ST_ConvexHull::Register(db);
	ST_CoveredBy::Register(db);
	ST_Covers::Register(db);
	ST_Crosses::Register(db);
	ST_Difference::Register(db);
	ST_Disjoint::Register(db);
	ST_Distance::Register(db);
	ST_DistanceWithin::Register(db);
	ST_Equals::Register(db);
	ST_Envelope::Register(db);
	ST_Intersection::Register(db);
	ST_Intersects::Register(db);
	ST_IsRing::Register(db);
	ST_IsSimple::Register(db);
	ST_IsValid::Register(db);
	ST_LineMerge::Register(db);
	ST_MakeValid::Register(db);
	ST_Normalize::Register(db);
	ST_Overlaps::Register(db);
	ST_PointOnSurface::Register(db);
	ST_ReducePrecision::Register(db);
	ST_RemoveRepeatedPoints::Register(db);
	ST_Reverse::Register(db);
	ST_ShortestLine::Register(db);
	ST_Simplify::Register(db);
	ST_SimplifyPreserveTopology::Register(db);
	ST_Touches::Register(db);
	ST_Union::Register(db);
	ST_Within::Register(db);

	// Aggregate Functions
	ST_Union_Agg::Register(db);
	ST_Intersection_Agg::Register(db);
}

} // namespace duckdb
