#include "duckdb/parser/parsed_data/create_index_info.hpp"
#include "spatial/index/rtree/rtree_index.hpp"
#include "spatial/index/rtree/rtree_index_create_logical.hpp"
#include "spatial/index/rtree/rtree_module.hpp"

#include "duckdb/main/database.hpp"

namespace duckdb {
//-------------------------------------------------------------
// Register
//-------------------------------------------------------------
void RTreeModule::RegisterIndexPlanCreate(DatabaseInstance &db) {

	db.config.operator_extensions.push_back(make_uniq<LogicalCreateRTreeIndexOperatorExtension>());
}

} // namespace duckdb