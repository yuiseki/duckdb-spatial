#include "spatial/common.hpp"
#include "spatial/core/types.hpp"
#include "spatial/core/functions/scalar.hpp"
#include "spatial/core/functions/common.hpp"
#include "spatial/core/geometry/geometry.hpp"
#include "spatial/core/function_builder.hpp"

#include "duckdb/common/vector_operations/generic_executor.hpp"

namespace spatial {

namespace core {

static void MakeEnvelopeFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = GeometryFunctionLocalState::ResetAndGet(state);
	auto count = args.size();

	auto &min_x_vec = args.data[0];
	auto &min_y_vec = args.data[1];
	auto &max_x_vec = args.data[2];
	auto &max_y_vec = args.data[3];

	using DOUBLE_TYPE = PrimitiveType<double>;
	using GEOMETRY_TYPE = PrimitiveType<geometry_t>;

	GenericExecutor::ExecuteQuaternary<DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, DOUBLE_TYPE, GEOMETRY_TYPE>(
	    min_x_vec, min_y_vec, max_x_vec, max_y_vec, result, count,
	    [&](DOUBLE_TYPE x_min, DOUBLE_TYPE y_min, DOUBLE_TYPE x_max, DOUBLE_TYPE y_max) {
		    auto box = Polygon::CreateFromBox(lstate.arena, x_min.val, y_min.val, x_max.val, y_max.val);
		    return Geometry::Serialize(box, result);
	    });
}

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
static constexpr const char *DOC_DESCRIPTION = R"(
    Returns a minimal bounding box polygon enclosing the input geometry
)";

static constexpr const char *DOC_EXAMPLE = R"(

)";
//------------------------------------------------------------------------------
// Register Functions
//------------------------------------------------------------------------------
void CoreScalarFunctions::RegisterStMakeEnvelope(DatabaseInstance &db) {
	FunctionBuilder::RegisterScalar(db, "ST_MakeEnvelope", [](ScalarFunctionBuilder &func) {
		func.AddVariant([](ScalarFunctionVariantBuilder &variant) {
			variant.AddParameter("min_x", LogicalType::DOUBLE);
			variant.AddParameter("min_y", LogicalType::DOUBLE);
			variant.AddParameter("max_x", LogicalType::DOUBLE);
			variant.AddParameter("max_y", LogicalType::DOUBLE);
			variant.SetReturnType(GeoTypes::GEOMETRY());
			variant.SetFunction(MakeEnvelopeFunction);
			variant.SetInit(GeometryFunctionLocalState::Init);

			variant.SetExample(DOC_EXAMPLE);
			variant.SetDescription(DOC_DESCRIPTION);
		});

		func.SetTag("ext", "spatial");
		func.SetTag("category", "construction");
	});
}

} // namespace core

} // namespace spatial
