#pragma once

#include "duckdb/common/string.hpp"
#include "duckdb/common/vector.hpp"
#include "duckdb/common/helper.hpp"

namespace duckdb {

class ExtensionLoader;
struct LogicalType;
struct FunctionData;
class ScalarFunction;
class AggregateFunction;
class ClientContext;
class Expression;
class BindScalarFunctionInput;
class BindAggregateFunctionInput;

typedef unique_ptr<FunctionData> (*bind_scalar_function_t)(BindScalarFunctionInput &input);
typedef unique_ptr<FunctionData> (*bind_aggregate_function_t)(BindAggregateFunctionInput &input);

struct GeoTypes {
	static LogicalType POINT_2D();
	static LogicalType POINT_3D();
	static LogicalType POINT_4D();
	static LogicalType LINESTRING_2D();
	static LogicalType LINESTRING_3D();
	static LogicalType POLYGON_2D();
	static LogicalType POLYGON_3D();
	static LogicalType BOX_2D();
	static LogicalType BOX_2DF();

	static void Register(ExtensionLoader &loader);

	static LogicalType CreateEnumType(const string &name, const vector<string> &members);

	static unique_ptr<FunctionData> PropagateCRS(BindScalarFunctionInput &input);
	static unique_ptr<FunctionData> PropagateCRS(BindAggregateFunctionInput &input);


	template <bind_scalar_function_t BIND>
	static unique_ptr<FunctionData> PropagateCRS(BindScalarFunctionInput &input) {
		PropagateCRS(input);
		return BIND(input);
	}

	template <bind_aggregate_function_t BIND>
	static unique_ptr<FunctionData> PropagateCRS(BindAggregateFunctionInput &input) {
		PropagateCRS(input);
		return BIND(input);
	}

};

} // namespace duckdb
