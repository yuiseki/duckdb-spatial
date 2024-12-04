#include "spatial/common.hpp"
#include "spatial/core/types.hpp"
#include "spatial/geos/functions/scalar.hpp"
#include "spatial/geos/functions/common.hpp"
#include "spatial/geos/geos_wrappers.hpp"
#include "spatial/geos/geos_executor.hpp"
#include "spatial/core/function_builder.hpp"

namespace spatial {

namespace geos {

using namespace spatial::core;

static void OverlapsFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = GEOSFunctionLocalState::ResetAndGet(state);
	auto &left = args.data[0];
	auto &right = args.data[1];
	auto count = args.size();
	GEOSExecutor::ExecuteSymmetricPreparedBinary(lstate, left, right, count, result, GEOSOverlaps_r,
	                                             GEOSPreparedOverlaps_r);
}

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
static constexpr const char *DOC_DESCRIPTION = R"(
    Returns true if geom1 "overlaps" geom2
)";

static constexpr const char *DOC_EXAMPLE = R"()";
//------------------------------------------------------------------------------
// Register Functions
//------------------------------------------------------------------------------

void GEOSScalarFunctions::RegisterStOverlaps(DatabaseInstance &db) {
	FunctionBuilder::RegisterScalar(db, "ST_Overlaps", [](ScalarFunctionBuilder &func) {
		func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
			variant.AddParameter("geom1", GeoTypes::GEOMETRY());
			variant.AddParameter("geom2", GeoTypes::GEOMETRY());
			variant.SetReturnType(LogicalType::BOOLEAN);
			variant.SetFunction(OverlapsFunction);
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
