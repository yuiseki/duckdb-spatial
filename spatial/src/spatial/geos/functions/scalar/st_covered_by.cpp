#include "spatial/common.hpp"
#include "spatial/core/types.hpp"
#include "spatial/geos/functions/scalar.hpp"
#include "spatial/geos/functions/common.hpp"
#include "spatial/geos/geos_wrappers.hpp"
#include "spatial/geos/geos_executor.hpp"
#include "spatial/core/function_builder.hpp"

#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"

namespace spatial {

namespace geos {

using namespace spatial::core;

static void CoveredByFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = GEOSFunctionLocalState::ResetAndGet(state);
	auto &left = args.data[0];
	auto &right = args.data[1];
	auto count = args.size();
	GEOSExecutor::ExecuteNonSymmetricPreparedBinary(lstate, left, right, count, result, GEOSCoveredBy_r,
	                                                GEOSPreparedCoveredBy_r);
}

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
static constexpr const char *DOC_DESCRIPTION = R"(
    Returns true if geom1 is "covered" by geom2
)";

static constexpr const char *DOC_EXAMPLE = R"()";
//------------------------------------------------------------------------------
// Register functions
//------------------------------------------------------------------------------
void GEOSScalarFunctions::RegisterStCoveredBy(DatabaseInstance &db) {

	FunctionBuilder::RegisterScalar(db, "ST_CoveredBy", [](ScalarFunctionBuilder &func) {
		func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
			variant.AddParameter("geom1", GeoTypes::GEOMETRY());
			variant.AddParameter("geom2", GeoTypes::GEOMETRY());
			variant.SetReturnType(LogicalType::BOOLEAN);
			variant.SetFunction(CoveredByFunction);
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
