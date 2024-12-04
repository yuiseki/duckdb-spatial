#include "duckdb/common/vector_operations/binary_executor.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "spatial/common.hpp"
#include "spatial/core/types.hpp"
#include "spatial/geos/functions/common.hpp"
#include "spatial/geos/functions/scalar.hpp"
#include "spatial/geos/geos_executor.hpp"
#include "spatial/geos/geos_wrappers.hpp"

#include "spatial/core/function_builder.hpp"

namespace spatial {

namespace geos {

using namespace spatial::core;

static void ExecuteContainsProperlyPrepared(GEOSFunctionLocalState &lstate, Vector &left, Vector &right, idx_t count,
                                            Vector &result) {
	auto &ctx = lstate.ctx.GetCtx();

	if (left.GetVectorType() == VectorType::CONSTANT_VECTOR && right.GetVectorType() != VectorType::CONSTANT_VECTOR) {
		auto &left_blob = FlatVector::GetData<geometry_t>(left)[0];
		auto left_geom = lstate.ctx.Deserialize(left_blob);
		auto left_prepared = make_uniq_geos(ctx, GEOSPrepare_r(ctx, left_geom.get()));

		UnaryExecutor::Execute<geometry_t, bool>(right, result, count, [&](geometry_t &right_blob) {
			auto right_geometry = lstate.ctx.Deserialize(right_blob);
			auto ok = GEOSPreparedContainsProperly_r(ctx, left_prepared.get(), right_geometry.get());
			return ok == 1;
		});
	} else {
		// ContainsProperly only has a prepared version, so we just prepare the left one always
		BinaryExecutor::Execute<geometry_t, geometry_t, bool>(
		    left, right, result, count, [&](geometry_t &left_blob, geometry_t &right_blob) {
			    auto left_geometry = lstate.ctx.Deserialize(left_blob);
			    auto right_geometry = lstate.ctx.Deserialize(right_blob);

			    auto left_prepared = make_uniq_geos(ctx, GEOSPrepare_r(ctx, left_geometry.get()));

			    auto ok = GEOSPreparedContainsProperly_r(ctx, left_prepared.get(), right_geometry.get());
			    return ok == 1;
		    });
	}
}

static void ContainsProperlyFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = GEOSFunctionLocalState::ResetAndGet(state);
	auto &left = args.data[0];
	auto &right = args.data[1];
	auto count = args.size();
	ExecuteContainsProperlyPrepared(lstate, left, right, count, result);
}

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
static constexpr const char *DOC_DESCRIPTION = R"(
    Returns true if geom1 "properly contains" geom2
)";

static constexpr const char *DOC_EXAMPLE = R"(

)";
//------------------------------------------------------------------------------
// Register functions
//------------------------------------------------------------------------------
void GEOSScalarFunctions::RegisterStContainsProperly(DatabaseInstance &db) {

	FunctionBuilder::RegisterScalar(db, "ST_ContainsProperly", [](ScalarFunctionBuilder &func) {
		func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
			variant.AddParameter("geom1", GeoTypes::GEOMETRY());
			variant.AddParameter("geom2", GeoTypes::GEOMETRY());
			variant.SetReturnType(LogicalType::BOOLEAN);
			variant.SetFunction(ContainsProperlyFunction);
			variant.SetInit(GEOSFunctionLocalState::Init);

			variant.SetExample(DOC_EXAMPLE);
			variant.SetDescription(DOC_DESCRIPTION);
		});

		func.SetTag("ext", "spatial");
		func.SetTag("category", "relation");
	});
}

} // namespace geos

} // namespace spatial
