#include "duckdb/common/vector_operations/binary_executor.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "spatial/common.hpp"
#include "spatial/core/types.hpp"
#include "spatial/geos/functions/common.hpp"
#include "spatial/geos/functions/scalar.hpp"
#include "spatial/geos/geos_wrappers.hpp"
#include "spatial/core/function_builder.hpp"

namespace spatial {

namespace geos {

using namespace spatial::core;

static void IntersectionFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = GEOSFunctionLocalState::ResetAndGet(state);
	auto &ctx = lstate.ctx.GetCtx();
	BinaryExecutor::Execute<geometry_t, geometry_t, geometry_t>(
	    args.data[0], args.data[1], result, args.size(), [&](geometry_t left, geometry_t right) {
		    auto left_geom = lstate.ctx.Deserialize(left);
		    auto right_geom = lstate.ctx.Deserialize(right);

		    auto result_geom = make_uniq_geos(ctx, GEOSIntersection_r(ctx, left_geom.get(), right_geom.get()));
		    return lstate.ctx.Serialize(result, result_geom);
	    });
}

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
static constexpr const char *DOC_DESCRIPTION = R"(
    Returns the "intersection" of geom1 and geom2
)";

static constexpr const char *DOC_EXAMPLE = R"()";
//------------------------------------------------------------------------------
// Register Functions
//------------------------------------------------------------------------------
void GEOSScalarFunctions::RegisterStIntersection(DatabaseInstance &db) {

	FunctionBuilder::RegisterScalar(db, "ST_Intersection", [](ScalarFunctionBuilder &func) {
		func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
			variant.AddParameter("geom1", GeoTypes::GEOMETRY());
			variant.AddParameter("geom2", GeoTypes::GEOMETRY());
			variant.SetReturnType(GeoTypes::GEOMETRY());
			variant.SetFunction(IntersectionFunction);
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
