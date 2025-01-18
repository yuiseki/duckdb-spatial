
#include "spatial/core/geoarrow.hpp"

#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/common/arrow/schema_metadata.hpp"
#include "duckdb/function/table/arrow/arrow_duck_schema.hpp"
#include "geos/vend/json.hpp"
#include "spatial/core/types.hpp"

namespace spatial {

namespace core {

struct GeoArrowWKB {
	static shared_ptr<ArrowType> GetType(const ArrowSchema &schema, const ArrowSchemaMetadata &schema_metadata) {
		string extension_metadata = schema_metadata.GetOption(ArrowSchemaMetadata::ARROW_METADATA_KEY);
		if (!extension_metadata.empty()) {
			nlohmann::json parsed(extension_metadata);
			if (parsed.is_object() && parsed.contains("edges")) {
				auto edges = parsed["edges"];
				if (!edges.is_string() || edges != "planar") {
					throw InternalException("Unsupported GeoArrow metadata: " + extension_metadata);
				}
			}
		}

		const auto format = string(schema.format);
		if (format == "z") {
			return make_shared_ptr<ArrowType>(GeoTypes::GEOMETRY(),
			                                  make_uniq<ArrowStringInfo>(ArrowVariableSizeType::NORMAL));
		} else if (format == "Z") {
			return make_shared_ptr<ArrowType>(GeoTypes::GEOMETRY(),
			                                  make_uniq<ArrowStringInfo>(ArrowVariableSizeType::NORMAL));
		} else if (format == "vz") {
			return make_shared_ptr<ArrowType>(GeoTypes::GEOMETRY(),
			                                  make_uniq<ArrowStringInfo>(ArrowVariableSizeType::NORMAL));
		}
		throw InvalidInputException("Arrow extension type \"%s\" not supported for geoarrow.wkb", format.c_str());
	}

	static void PopulateSchema(DuckDBArrowSchemaHolder &root_holder, ArrowSchema &schema, const LogicalType &type,
	                           ClientContext &context, const ArrowTypeExtension &extension) {
		ArrowSchemaMetadata schema_metadata;
		schema_metadata.AddOption(ArrowSchemaMetadata::ARROW_EXTENSION_NAME, "geoarrow.wkb");
		schema_metadata.AddOption(ArrowSchemaMetadata::ARROW_METADATA_KEY, "{}");
		root_holder.metadata_info.emplace_back(schema_metadata.SerializeMetadata());
		schema.metadata = root_holder.metadata_info.back().get();

		const auto options = context.GetClientProperties();
		if (options.arrow_offset_size == ArrowOffsetSize::LARGE) {
			schema.format = "Z";
		} else {
			schema.format = "z";
		}
	}

	static void ArrowToDuck(ClientContext &context, Vector &source, Vector &result, idx_t count) {
		// Call geomfromwkb()?
	}

	static void DuckToArrow(ClientContext &context, Vector &source, Vector &result, idx_t count) {
		// Call aswkb()?
	}
};

void GeoArrow::Register(DatabaseInstance &db) {
	// Arrow Extensions
	db.config.RegisterArrowExtension(
	    {"geoarrow.wkb", GeoArrowWKB::PopulateSchema, GeoArrowWKB::GetType,
	     make_shared_ptr<ArrowTypeExtensionData>(GeoTypes::GEOMETRY(), LogicalType::BLOB, GeoArrowWKB::ArrowToDuck,
	                                             GeoArrowWKB::DuckToArrow)});
}

} // namespace core

} // namespace spatial
