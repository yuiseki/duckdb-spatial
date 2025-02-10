#include "spatial/geometry/bbox.hpp"
#include "spatial/geometry/geometry_serialization.hpp"
#include "spatial/geometry/sgl.hpp"
#include "spatial/geometry/geometry_type.hpp"
#include "spatial/modules/main/spatial_functions.hpp"
#include "spatial/spatial_types.hpp"
#include "spatial/util/function_builder.hpp"

namespace duckdb {

namespace {

struct ExtentAggState {
	bool is_set;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
};

//------------------------------------------------------------------------
// ENVELOPE AGG
//------------------------------------------------------------------------
struct ExtentAggFunction {
	template <class STATE>
	static void Initialize(STATE &state) {
		state.is_set = false;
	}

	template <class STATE, class OP>
	static void Combine(const STATE &source, STATE &target, AggregateInputData &) {
		if (!source.is_set) {
			return;
		}
		if (!target.is_set) {
			target = source;
			return;
		}
		target.xmin = std::min(target.xmin, source.xmin);
		target.xmax = std::max(target.xmax, source.xmax);
		target.ymin = std::min(target.ymin, source.ymin);
		target.ymax = std::max(target.ymax, source.ymax);
	}

	template <class INPUT_TYPE, class STATE, class OP>
	static void Operation(STATE &state, const INPUT_TYPE &input, AggregateUnaryInput &) {
		Box2D<double> bbox;
		if (!state.is_set) {
			if (input.TryGetCachedBounds(bbox)) {
				state.is_set = true;
				state.xmin = bbox.min.x;
				state.xmax = bbox.max.x;
				state.ymin = bbox.min.y;
				state.ymax = bbox.max.y;
			}
		} else {
			if (input.TryGetCachedBounds(bbox)) {
				state.xmin = std::min(state.xmin, bbox.min.x);
				state.xmax = std::max(state.xmax, bbox.max.x);
				state.ymin = std::min(state.ymin, bbox.min.y);
				state.ymax = std::max(state.ymax, bbox.max.y);
			}
		}
	}

	template <class INPUT_TYPE, class STATE, class OP>
	static void ConstantOperation(STATE &state, const INPUT_TYPE &input, AggregateUnaryInput &agg, idx_t) {
		Operation<INPUT_TYPE, STATE, OP>(state, input, agg);
	}

	template <class T, class STATE>
	static void Finalize(STATE &state, T &target, AggregateFinalizeData &finalize_data) {
		if (!state.is_set) {
			finalize_data.ReturnNull();
		} else {
			// We can create the bounding box polygon directly on the stack
			double buf[10];
			buf[0] = state.xmin;
			buf[1] = state.ymin;

			buf[2] = state.xmin;
			buf[3] = state.ymax;

			buf[4] = state.xmax;
			buf[5] = state.ymax;

			buf[6] = state.xmax;
			buf[7] = state.ymin;

			buf[8] = state.xmin;
			buf[9] = state.ymin;

			auto ring = sgl::linestring::make_empty();
			ring.set_vertex_data(reinterpret_cast<const char*>(buf), 5);

			auto bbox = sgl::polygon::make_empty();
			bbox.append_part(&ring);

			const auto size = Serde::GetRequiredSize(bbox);
			auto blob = StringVector::EmptyString(finalize_data.result, size);
			Serde::Serialize(bbox, blob.GetDataWriteable(), size);
			blob.Finalize();

			// TODO: dont use geometry_t here
			target = geometry_t(blob);
		}
	}

	static bool IgnoreNull() {
		return true;
	}
};

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
// static constexpr DocTag DOC_TAGS[] = {{"ext", "spatial"}, {"category", "construction"}};
static constexpr const char *DOC_DESCRIPTION = R"(
    Computes the minimal-bounding-box polygon containing the set of input geometries
)";
static constexpr const char *DOC_EXAMPLE = R"(
	SELECT ST_Extent_Agg(geom) FROM UNNEST([ST_Point(1,1), ST_Point(5,5)]) AS _(geom);
	-- POLYGON ((1 1, 1 5, 5 5, 5 1, 1 1))
)";

static constexpr const char *DOC_ALIAS_DESCRIPTION = R"(
	Alias for [ST_Extent_Agg](#st_extent_agg).

	Computes the minimal-bounding-box polygon containing the set of input geometries.
)";


} // namespace

//------------------------------------------------------------------------
// Register
//------------------------------------------------------------------------
void RegisterSpatialAggregateFunctions(DatabaseInstance &db) {

	// TODO: Dont use geometry_t here
	const auto agg = AggregateFunction::UnaryAggregate<ExtentAggState, geometry_t, geometry_t, ExtentAggFunction>(
		GeoTypes::GEOMETRY(), GeoTypes::GEOMETRY());

	FunctionBuilder::RegisterAggregate(db, "ST_Extent_Agg", [&](AggregateFunctionBuilder &func) {

		func.SetFunction(agg);
		func.SetDescription(DOC_DESCRIPTION);
		func.SetExample(DOC_EXAMPLE);

		func.SetTag("ext", "spatial");
		func.SetTag("category", "construction");
	});

	FunctionBuilder::RegisterAggregate(db, "ST_Envelope_Agg", [&](AggregateFunctionBuilder &func) {

		func.SetFunction(agg);
		func.SetDescription(DOC_ALIAS_DESCRIPTION);
		func.SetExample(DOC_EXAMPLE);

		func.SetTag("ext", "spatial");
		func.SetTag("category", "construction");
	});
}

} // namespace duckdb