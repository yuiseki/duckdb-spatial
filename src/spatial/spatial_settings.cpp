#include "spatial/spatial_settings.hpp"

#include "duckdb/main/database.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/client_context.hpp"

namespace duckdb {

static constexpr auto GEOMETRY_ALWAYS_XY = "geometry_always_xy";

bool SpatialSettings::AlwaysXY(ClientContext &context, bool &is_set) {
	Value value;
	if (context.TryGetCurrentSetting(GEOMETRY_ALWAYS_XY, value)) {
		is_set = true;
		return value.GetValue<bool>();
	}

	is_set = false;
	return false;
}

void SpatialSettings::Register(ExtensionLoader &loader) {

	auto &db = loader.GetDatabaseInstance();

	db.config.AddExtensionOption(GEOMETRY_ALWAYS_XY,
	                             "ignore the defined axis order of coordinate reference systems and always treat them "
	                             "as (easting, northing) instead of e.g. (latitude, longitude)",
	                             LogicalType::BOOLEAN, Value(LogicalTypeId::BOOLEAN));
}

} // namespace duckdb
