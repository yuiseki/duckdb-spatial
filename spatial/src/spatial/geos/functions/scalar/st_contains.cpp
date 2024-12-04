#include "spatial/common.hpp"
#include "spatial/core/types.hpp"
#include "spatial/geos/functions/scalar.hpp"
#include "spatial/geos/functions/common.hpp"
#include "spatial/geos/geos_wrappers.hpp"
#include "spatial/geos/geos_executor.hpp"

#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"

#include "spatial/core/function_builder.hpp"

namespace spatial {

namespace geos {

using namespace spatial::core;

static void ContainsFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = GEOSFunctionLocalState::ResetAndGet(state);
	auto &left = args.data[0];
	auto &right = args.data[1];
	auto count = args.size();
	GEOSExecutor::ExecuteNonSymmetricPreparedBinary(lstate, left, right, count, result, GEOSContains_r,
	                                                GEOSPreparedContains_r);
}

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
static constexpr const char *DOC_DESCRIPTION = R"(
Returns true if geom1 contains geom2.
)";

static constexpr const char *DOC_EXAMPLE = R"(
select st_contains('POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))'::geometry, 'POINT(0.5 0.5)'::geometry);
----
true
)";

//------------------------------------------------------------------------------
// Register functions
//------------------------------------------------------------------------------
void GEOSScalarFunctions::RegisterStContains(DatabaseInstance &db) {

	FunctionBuilder::RegisterScalar(db, "ST_Contains", [](ScalarFunctionBuilder &func) {
		func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
			variant.AddParameter("geom1", GeoTypes::GEOMETRY());
			variant.AddParameter("geom2", GeoTypes::GEOMETRY());
			variant.SetReturnType(LogicalType::BOOLEAN);
			variant.SetFunction(ContainsFunction);
			variant.SetInit(GEOSFunctionLocalState::Init);

			variant.SetExample(DOC_EXAMPLE);
			variant.SetDescription(DOC_DESCRIPTION);
		});

		func.SetDescription(DOC_DESCRIPTION);

		func.SetTag("ext", "spatial");
		func.SetTag("category", "relation");
	});
}

} // namespace geos

} // namespace spatial
