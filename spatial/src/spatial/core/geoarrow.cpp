
#include "spatial/core/geoarrow.hpp"

#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/common/arrow/schema_metadata.hpp"
#include "duckdb/function/table/arrow/arrow_duck_schema.hpp"
#include "spatial/core/geometry/wkb_reader.hpp"
#include "spatial/core/geometry/wkb_writer.hpp"
#include "spatial/core/types.hpp"
#include "yyjson.h"

namespace spatial {

namespace core {

struct GeoArrowWKB {
	static shared_ptr<ArrowType> GetType(const ArrowSchema &schema, const ArrowSchemaMetadata &schema_metadata) {
		// Validate extension metadata. This metadata also contains a CRS, which we drop
		// because the GEOMETRY type does not implement a CRS at the type level.
		string extension_metadata = schema_metadata.GetOption(ArrowSchemaMetadata::ARROW_METADATA_KEY);
		if (!extension_metadata.empty()) {
			using namespace duckdb_yyjson_spatial;

			unique_ptr<yyjson_doc, void (*)(yyjson_doc *)> doc(
			    yyjson_read(extension_metadata.data(), extension_metadata.size(), YYJSON_READ_NOFLAG), yyjson_doc_free);
			if (!doc) {
				throw SerializationException("Invalid JSON in GeoArrow metadata");
			}

			yyjson_val *val = yyjson_doc_get_root(doc.get());
			if (!yyjson_is_obj(val)) {
				throw SerializationException("Invalid GeoArrow metadata: not a JSON object");
			}

			yyjson_val *edges = yyjson_obj_get(val, "edges");
			if (edges && yyjson_is_str(edges) && std::strcmp(yyjson_get_str(edges), "planar") != 0) {
				throw NotImplementedException("Can't import non-planar edges");
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
		ArenaAllocator arena(BufferAllocator::Get(context));
		WKBReader reader(arena);

		UnaryExecutor::ExecuteWithNulls<string_t, geometry_t>(source, result, count,
		                                                      [&](string_t input, ValidityMask &mask, idx_t idx) {
			                                                      auto geom = reader.Deserialize(input);
			                                                      return Geometry::Serialize(geom, result);
		                                                      });
	}

	static void DuckToArrow(ClientContext &context, Vector &source, Vector &result, idx_t count) {
		WKBWriter writer;
		UnaryExecutor::Execute<geometry_t, string_t>(source, result, count,
		                                             [&](geometry_t input) { return writer.Write(input, result); });
	}
};

void RegisterArrowExtensions(DBConfig &config) {
	config.RegisterArrowExtension(
	    {"geoarrow.wkb", GeoArrowWKB::PopulateSchema, GeoArrowWKB::GetType,
	     make_shared_ptr<ArrowTypeExtensionData>(GeoTypes::GEOMETRY(), LogicalType::BLOB, GeoArrowWKB::ArrowToDuck,
	                                             GeoArrowWKB::DuckToArrow)});
}

class GeoArrowRegisterFunctionData : public TableFunctionData {
public:
	GeoArrowRegisterFunctionData() : finished(false) {
	}
	bool finished {false};
};

static inline duckdb::unique_ptr<FunctionData> GeoArrowRegisterBind(ClientContext &context,
                                                                    TableFunctionBindInput &input,
                                                                    vector<LogicalType> &return_types,
                                                                    vector<string> &names) {
	names.push_back("registered");
	return_types.push_back(LogicalType::BOOLEAN);
	return make_uniq<GeoArrowRegisterFunctionData>();
}

void GeoArrowRegisterScan(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->CastNoConst<GeoArrowRegisterFunctionData>();
	if (data.finished) {
		return;
	}

	DBConfig &config = DatabaseInstance::GetDatabase(context).config;
	if (config.HasArrowExtension(GeoTypes::GEOMETRY())) {
		output.SetValue(0, 0, false);
	} else {
		RegisterArrowExtensions(config);
		output.SetValue(0, 0, true);
	}

	output.SetCardinality(1);
	data.finished = true;
}

void GeoArrow::Register(DatabaseInstance &db) {
	TableFunction register_func("register_geoarrow_extensions", {}, GeoArrowRegisterScan, GeoArrowRegisterBind);
	ExtensionUtil::RegisterFunction(db, register_func);
}

} // namespace core

} // namespace spatial
