#pragma once

#include "duckdb.hpp"

#include "spatial/common.hpp"
#include "duckdb/function/function_set.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/parser/parsed_data/create_function_info.hpp"

namespace spatial {

namespace core {
//------------------------------------------------------------------------------
// Scalar Function Variant Builder
//------------------------------------------------------------------------------

class ScalarFunctionVariantBuilder {
	friend class ScalarFunctionBuilder;

public:
	void AddParameter(const char *name, const LogicalType &type);
	void SetReturnType(LogicalType type);
	void SetFunction(scalar_function_t fn);
	void SetInit(init_local_state_t init);
	void SetBind(bind_scalar_function_t bind);
	void SetDescription(const string &desc);
	void SetExample(const string &ex);

private:
	explicit ScalarFunctionVariantBuilder() : function({}, LogicalTypeId::INVALID, nullptr) {
	}

	ScalarFunction function;
	FunctionDescription description = {};
};

inline void ScalarFunctionVariantBuilder::AddParameter(const char *name, const LogicalType &type) {
	function.arguments.emplace_back(type);
	description.parameter_names.emplace_back(name);
	description.parameter_types.emplace_back(type);
}

inline void ScalarFunctionVariantBuilder::SetReturnType(LogicalType type) {
	function.return_type = std::move(type);
}

inline void ScalarFunctionVariantBuilder::SetFunction(scalar_function_t fn) {
	function.function = fn;
}

inline void ScalarFunctionVariantBuilder::SetInit(init_local_state_t init) {
	function.init_local_state = init;
}

inline void ScalarFunctionVariantBuilder::SetBind(bind_scalar_function_t bind) {
	function.bind = bind;
}

inline void ScalarFunctionVariantBuilder::SetDescription(const string &desc) {
	description.description = desc;
}

inline void ScalarFunctionVariantBuilder::SetExample(const string &ex) {
	description.examples.emplace_back(ex);
}

//------------------------------------------------------------------------------
// Scalar Function Builder
//------------------------------------------------------------------------------

class ScalarFunctionBuilder {
	friend class FunctionBuilder;

public:
	template <class CALLBACK>
	void AddVariant(CALLBACK &&callback);
	void SetTag(const string &key, const string &value);
	void SetDescription(const string &desc);
	void SetExample(const string &ex);

private:
	explicit ScalarFunctionBuilder(const char *name) : set(name) {
	}

	ScalarFunctionSet set;
	vector<FunctionDescription> descriptions = {};
	unordered_map<string, string> tags = {};

	// If not set by a variant
	string default_description;
	string default_example;
};

inline void ScalarFunctionBuilder::SetDescription(const string &desc) {
	default_description = desc;
}

inline void ScalarFunctionBuilder::SetExample(const string &ex) {
	default_example = ex;
}

inline void ScalarFunctionBuilder::SetTag(const string &key, const string &value) {
	tags[key] = value;
}

template <class CALLBACK>
void ScalarFunctionBuilder::AddVariant(CALLBACK &&callback) {
	ScalarFunctionVariantBuilder builder;

	callback(builder);

	// A return type is required
	if (builder.function.return_type.id() == LogicalTypeId::INVALID) {
		throw InternalException("Return type not set in ScalarFunctionBuilder::AddVariant");
	}

	// Add the new variant to the set
	set.AddFunction(std::move(builder.function));

	// Add the default description if not set by the variant
	if (builder.description.description.empty()) {
		builder.description.description = default_description;
	}

	if (builder.description.examples.empty()) {
		builder.description.examples.emplace_back(default_example);
	}

	// Add the description
	descriptions.emplace_back(std::move(builder.description));
}

//------------------------------------------------------------------------------
// Function Builder
//------------------------------------------------------------------------------

class FunctionBuilder {
public:
	template <class CALLBACK>
	static void RegisterScalar(DatabaseInstance &db, const char *name, CALLBACK &&callback);

private:
	static void Register(DatabaseInstance &db, const char *name, ScalarFunctionBuilder &builder);
};

template <class CALLBACK>
void FunctionBuilder::RegisterScalar(DatabaseInstance &db, const char *name, CALLBACK &&callback) {
	ScalarFunctionBuilder builder(name);
	callback(builder);

	Register(db, name, builder);
}

} // namespace core

} // namespace spatial