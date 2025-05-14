#pragma once
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/parser/parsed_data/create_index_info.hpp"
#include "duckdb/planner/operator/logical_extension_operator.hpp"
#include "duckdb/common/serializer/serializer.hpp"
#include "duckdb/common/serializer/deserializer.hpp"

namespace duckdb {

class LogicalCreateRTreeIndex final : public LogicalExtensionOperator {
public:
	static constexpr auto OPERATOR_TYPE_NAME = "logical_rtree_create_index";

	// Info for index creation
	unique_ptr<CreateIndexInfo> info;

	//! The table to create the index for
	TableCatalogEntry &table;

	//! Unbound expressions to be used in the optimizer
	vector<unique_ptr<Expression>> unbound_expressions;

public:
	LogicalCreateRTreeIndex(unique_ptr<CreateIndexInfo> info_p, vector<unique_ptr<Expression>> expressions_p,
	                        TableCatalogEntry &table_p);
	void ResolveTypes() override;
	void ResolveColumnBindings(ColumnBindingResolver &res, vector<ColumnBinding> &bindings) override;

	// Actually create and plan the index creation
	PhysicalOperator &CreatePlan(ClientContext &context, PhysicalPlanGenerator &planner) override;

public:
	void Serialize(Serializer &writer) const override;
	static unique_ptr<LogicalExtensionOperator> Deserialize(Deserializer &reader);
public:
	string GetName() const override {
		return "CREATE_RTREE_INDEX";
	}

	string GetExtensionName() const override {
		return "duckdb_spatial";
	}
};

} // namespace duckdb