#include "duckdb/parser/parsed_data/create_table_function_info.hpp"

#include "spatial/common.hpp"
#include "spatial/gdal/functions.hpp"

#include "ogrsf_frmts.h"

namespace spatial {

namespace gdal {

//----------------------------------------------------------------------------------------------------------------------
// ST_Read
//----------------------------------------------------------------------------------------------------------------------
struct STRead {

	// Bind ------------------------------------------------------------------------------------------------------------
	struct STReadData : public FunctionData {};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<string> &names) {
		return nullptr;
	}

	// Global State ----------------------------------------------------------------------------------------------------
	struct GlobalState : public GlobalTableFunctionState {};

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
		return nullptr;
	}

	// Local State -----------------------------------------------------------------------------------------------------
	struct LocalState : public LocalTableFunctionState {};

	static unique_ptr<LocalTableFunctionState> InitLocal(ClientContext &context, LocalTableFunctionState &state) {
		return nullptr;
	}

	// Execute ---------------------------------------------------------------------------------------------------------
	static void Execute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	}

	static void Finalize(ClientContext &context, unique_ptr<GlobalTableFunctionState> state) {
	}

	// Docs ------------------------------------------------------------------------------------------------------------
	static constexpr const char *DOC_DESCRIPTION = R"(
    Read and import a variety of geospatial file formats using the GDAL library.

    The `ST_Read` table function is based on the [GDAL](https://gdal.org/index.html) translator library and enables reading spatial data from a variety of geospatial vector file formats as if they were DuckDB tables.

    > See [ST_Drivers](#st_drivers) for a list of supported file formats and drivers.

    Except for the `path` parameter, all parameters are optional.

    | Parameter | Type | Description |
    | --------- | -----| ----------- |
    | `path` | VARCHAR | The path to the file to read. Mandatory |
    | `sequential_layer_scan` | BOOLEAN | If set to true, the table function will scan through all layers sequentially and return the first layer that matches the given layer name. This is required for some drivers to work properly, e.g., the OSM driver. |
    | `spatial_filter` | WKB_BLOB | If set to a WKB blob, the table function will only return rows that intersect with the given WKB geometry. Some drivers may support efficient spatial filtering natively, in which case it will be pushed down. Otherwise the filtering is done by GDAL which may be much slower. |
    | `open_options` | VARCHAR[] | A list of key-value pairs that are passed to the GDAL driver to control the opening of the file. E.g., the GeoJSON driver supports a FLATTEN_NESTED_ATTRIBUTES=YES option to flatten nested attributes. |
    | `layer` | VARCHAR | The name of the layer to read from the file. If NULL, the first layer is returned. Can also be a layer index (starting at 0). |
    | `allowed_drivers` | VARCHAR[] | A list of GDAL driver names that are allowed to be used to open the file. If empty, all drivers are allowed. |
    | `sibling_files` | VARCHAR[] | A list of sibling files that are required to open the file. E.g., the ESRI Shapefile driver requires a .shx file to be present. Although most of the time these can be discovered automatically. |
    | `spatial_filter_box` | BOX_2D | If set to a BOX_2D, the table function will only return rows that intersect with the given bounding box. Similar to spatial_filter. |
    | `keep_wkb` | BOOLEAN | If set, the table function will return geometries in a wkb_geometry column with the type WKB_BLOB (which can be cast to BLOB) instead of GEOMETRY. This is useful if you want to use DuckDB with more exotic geometry subtypes that DuckDB spatial doesnt support representing in the GEOMETRY type yet. |

    Note that GDAL is single-threaded, so this table function will not be able to make full use of parallelism.

    By using `ST_Read`, the spatial extension also provides “replacement scans” for common geospatial file formats, allowing you to query files of these formats as if they were tables directly.

    ```sql
    SELECT * FROM './path/to/some/shapefile/dataset.shp';
    ```

    In practice this is just syntax-sugar for calling ST_Read, so there is no difference in performance. If you want to pass additional options, you should use the ST_Read table function directly.

    The following formats are currently recognized by their file extension:

    | Format | Extension |
    | ------ | --------- |
    | ESRI ShapeFile | .shp |
    | GeoPackage | .gpkg |
    | FlatGeoBuf | .fgb |
	)";

	// Register --------------------------------------------------------------------------------------------------------
	static void Register(DatabaseInstance &db) {
		TableFunction func("ST_Read", {}, Execute, Bind, InitGlobal, InitLocal, Finalize);

		ExtensionUtil::RegisterFunction(db, func);
	}
};

// Simple table function to list all the drivers available
unique_ptr<FunctionData> GdalDriversTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                        vector<LogicalType> &return_types, vector<string> &names) {
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::BOOLEAN);
	return_types.emplace_back(LogicalType::BOOLEAN);
	return_types.emplace_back(LogicalType::BOOLEAN);
	return_types.emplace_back(LogicalType::VARCHAR);
	names.emplace_back("short_name");
	names.emplace_back("long_name");
	names.emplace_back("can_create");
	names.emplace_back("can_copy");
	names.emplace_back("can_open");
	names.emplace_back("help_url");

	auto driver_count = GDALGetDriverCount();
	auto result = make_uniq<BindData>(driver_count);
	return std::move(result);
}

unique_ptr<GlobalTableFunctionState> GdalDriversTableFunction::Init(ClientContext &context,
                                                                    TableFunctionInitInput &input) {
	auto result = make_uniq<State>();
	return std::move(result);
}

void GdalDriversTableFunction::Execute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &state = (State &)*input.global_state;
	auto &bind_data = (BindData &)*input.bind_data;

	idx_t count = 0;
	auto next_idx = MinValue<idx_t>(state.current_idx + STANDARD_VECTOR_SIZE, bind_data.driver_count);

	for (; state.current_idx < next_idx; state.current_idx++) {
		auto driver = GDALGetDriver((int)state.current_idx);

		// Check if the driver is a vector driver
		if (GDALGetMetadataItem(driver, GDAL_DCAP_VECTOR, nullptr) == nullptr) {
			continue;
		}

		auto short_name = Value::CreateValue(GDALGetDriverShortName(driver));
		auto long_name = Value::CreateValue(GDALGetDriverLongName(driver));

		const char *create_flag = GDALGetMetadataItem(driver, GDAL_DCAP_CREATE, nullptr);
		auto create_value = Value::CreateValue(create_flag != nullptr);

		const char *copy_flag = GDALGetMetadataItem(driver, GDAL_DCAP_CREATECOPY, nullptr);
		auto copy_value = Value::CreateValue(copy_flag != nullptr);
		const char *open_flag = GDALGetMetadataItem(driver, GDAL_DCAP_OPEN, nullptr);
		auto open_value = Value::CreateValue(open_flag != nullptr);

		auto help_topic_flag = GDALGetDriverHelpTopic(driver);
		auto help_topic_value = help_topic_flag == nullptr
		                            ? Value(LogicalType::VARCHAR)
		                            : Value(StringUtil::Format("https://gdal.org/%s", help_topic_flag));

		output.data[0].SetValue(count, short_name);
		output.data[1].SetValue(count, long_name);
		output.data[2].SetValue(count, create_value);
		output.data[3].SetValue(count, copy_value);
		output.data[4].SetValue(count, open_value);
		output.data[5].SetValue(count, help_topic_value);
		count++;
	}
	output.SetCardinality(count);
}

//------------------------------------------------------------------------------
// Documentation
//------------------------------------------------------------------------------
static constexpr DocTag DOC_TAGS[] = {{"ext", "spatial"}};

static constexpr const char *DOC_DESCRIPTION = R"(
    Returns the list of supported GDAL drivers and file formats

    Note that far from all of these drivers have been tested properly, and some may require additional options to be passed to work as expected. If you run into any issues please first consult the [consult the GDAL docs](https://gdal.org/drivers/vector/index.html).
)";

static constexpr const char *DOC_EXAMPLE = R"(
    SELECT * FROM ST_Drivers();
)";

//------------------------------------------------------------------------------
// Register
//------------------------------------------------------------------------------
void GdalDriversTableFunction::Register(DatabaseInstance &db) {
	TableFunction func("ST_Drivers", {}, Execute, Bind, Init);

	ExtensionUtil::RegisterFunction(db, func);
	DocUtil::AddDocumentation(db, "ST_Drivers", DOC_DESCRIPTION, DOC_EXAMPLE, DOC_TAGS);
}

} // namespace gdal

} // namespace spatial