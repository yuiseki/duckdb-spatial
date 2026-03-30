#include "spatial/spatial_types.hpp"

#include "duckdb/common/type_visitor.hpp"
#include "duckdb/common/types/geometry_crs.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

LogicalType GeoTypes::POINT_2D() {
	auto type = LogicalType::STRUCT({{"x", LogicalType::DOUBLE}, {"y", LogicalType::DOUBLE}});
	type.SetAlias("POINT_2D");
	return type;
}

LogicalType GeoTypes::POINT_3D() {
	auto type =
	    LogicalType::STRUCT({{"x", LogicalType::DOUBLE}, {"y", LogicalType::DOUBLE}, {"z", LogicalType::DOUBLE}});
	type.SetAlias("POINT_3D");
	return type;
}

LogicalType GeoTypes::POINT_4D() {
	auto type = LogicalType::STRUCT({{"x", LogicalType::DOUBLE},
	                                 {"y", LogicalType::DOUBLE},
	                                 {"z", LogicalType::DOUBLE},
	                                 {"m", LogicalType::DOUBLE}});
	type.SetAlias("POINT_4D");
	return type;
}

LogicalType GeoTypes::BOX_2D() {
	auto type = LogicalType::STRUCT({{"min_x", LogicalType::DOUBLE},
	                                 {"min_y", LogicalType::DOUBLE},
	                                 {"max_x", LogicalType::DOUBLE},
	                                 {"max_y", LogicalType::DOUBLE}});
	type.SetAlias("BOX_2D");
	return type;
}

LogicalType GeoTypes::BOX_2DF() {
	auto type = LogicalType::STRUCT({{"min_x", LogicalType::FLOAT},
	                                 {"min_y", LogicalType::FLOAT},
	                                 {"max_x", LogicalType::FLOAT},
	                                 {"max_y", LogicalType::FLOAT}});
	type.SetAlias("BOX_2DF");
	return type;
}

LogicalType GeoTypes::LINESTRING_2D() {
	auto type = LogicalType::LIST(LogicalType::STRUCT({{"x", LogicalType::DOUBLE}, {"y", LogicalType::DOUBLE}}));
	type.SetAlias("LINESTRING_2D");
	return type;
}

LogicalType GeoTypes::LINESTRING_3D() {
	auto type = LogicalType::LIST(
	    LogicalType::STRUCT({{"x", LogicalType::DOUBLE}, {"y", LogicalType::DOUBLE}, {"z", LogicalType::DOUBLE}}));
	type.SetAlias("LINESTRING_3D");
	return type;
}

LogicalType GeoTypes::POLYGON_2D() {
	auto type = LogicalType::LIST(
	    LogicalType::LIST(LogicalType::STRUCT({{"x", LogicalType::DOUBLE}, {"y", LogicalType::DOUBLE}})));
	type.SetAlias("POLYGON_2D");
	return type;
}

LogicalType GeoTypes::POLYGON_3D() {
	auto type = LogicalType::LIST(LogicalType::LIST(
	    LogicalType::STRUCT({{"x", LogicalType::DOUBLE}, {"y", LogicalType::DOUBLE}, {"z", LogicalType::DOUBLE}})));
	type.SetAlias("POLYGON_3D");
	return type;
}

LogicalType GeoTypes::CreateEnumType(const string &name, const vector<string> &members) {
	auto varchar_vector = Vector(LogicalType::VARCHAR, members.size());
	auto varchar_data = FlatVector::GetData<string_t>(varchar_vector);
	for (idx_t i = 0; i < members.size(); i++) {
		auto str = string_t(members[i]);
		varchar_data[i] = str.IsInlined() ? str : StringVector::AddString(varchar_vector, str);
	}
	auto enum_type = LogicalType::ENUM(name, varchar_vector, members.size());
	enum_type.SetAlias(name);
	return enum_type;
}

static unique_ptr<FunctionData> PropagateTypesInternal(ClientContext &context, BaseScalarFunction &bound_function,
                                                       vector<unique_ptr<Expression>> &arguments) {

	CoordinateReferenceSystem crs;
	auto found_crs = false;

	for (idx_t arg_idx = 0; arg_idx < arguments.size(); arg_idx++) {
		const auto &arg = arguments[arg_idx];

		auto has_crs = false;

		TypeVisitor::Contains(arg->return_type, [&](const LogicalType &type) {
			if (type.id() == LogicalTypeId::GEOMETRY && GeoType::HasCRS(type)) {
				has_crs = true;
				if (!found_crs) {
					found_crs = true;
					crs = GeoType::GetCRS(type);
				} else {
					auto &type_crs = GeoType::GetCRS(type);
					if (!crs.Equals(type_crs)) {
						throw BinderException(
						    arg->query_location,
						    "Cannot call function '%s' with geometries of different coordinate reference systems "
						    "(CRS).\n"
						    "First geometry type is in '%s' which is not compatible with '%s'.\n"
						    " * Use 'ST_Transform' to convert geometries to the same CRS before passing them to this "
						    "function.\n"
						    " * Use 'ST_SetCRS' to explicitly override the CRS of a geometry expression, without "
						    "performing a "
						    "transformation.\n",
						    bound_function.name, crs.GetIdentifier(), type_crs.GetIdentifier());
					}
				}
			}
			return false;
		});

		if (has_crs) {
			// Override the type so that we set the CRS
			bound_function.arguments[arg_idx] = arg->return_type;
		}
	}

	const auto return_has_geom = TypeVisitor::Contains(bound_function.return_type, LogicalTypeId::GEOMETRY);
	if (found_crs && return_has_geom) {
		// If the return type is geometry, we need to set the CRS on it as well
		bound_function.return_type =
		    TypeVisitor::VisitReplace(bound_function.return_type, [&](const LogicalType &type) {
			    if (type.id() == LogicalTypeId::GEOMETRY) {
				    return LogicalType::GEOMETRY(crs);
			    }
			    return type;
		    });
	}

	return nullptr;
}
unique_ptr<FunctionData> GeoTypes::PropagateCRS(ClientContext &context, ScalarFunction &bound_function,
                                                vector<unique_ptr<Expression>> &arguments) {
	return PropagateTypesInternal(context, bound_function, arguments);
}

unique_ptr<FunctionData> GeoTypes::PropagateCRS(ClientContext &context, AggregateFunction &bound_function,
                                                vector<unique_ptr<Expression>> &arguments) {
	return PropagateTypesInternal(context, bound_function, arguments);
}

void GeoTypes::Register(ExtensionLoader &loader) {

	// POINT_2D
	loader.RegisterType("POINT_2D", GeoTypes::POINT_2D());

	// POINT_3D
	loader.RegisterType("POINT_3D", GeoTypes::POINT_3D());

	// POINT_4D
	loader.RegisterType("POINT_4D", GeoTypes::POINT_4D());

	// LineString2D
	loader.RegisterType("LINESTRING_2D", GeoTypes::LINESTRING_2D());

	// LineString3D
	loader.RegisterType("LINESTRING_3D", GeoTypes::LINESTRING_3D());

	// Polygon2D
	loader.RegisterType("POLYGON_2D", GeoTypes::POLYGON_2D());

	// Polygon3D
	loader.RegisterType("POLYGON_3D", GeoTypes::POLYGON_3D());

	// Box2D
	loader.RegisterType("BOX_2D", GeoTypes::BOX_2D());

	// Box2DF
	loader.RegisterType("BOX_2DF", GeoTypes::BOX_2DF());
}

} // namespace duckdb
