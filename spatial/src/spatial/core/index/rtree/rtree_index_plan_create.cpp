#include "duckdb/parser/parsed_data/create_index_info.hpp"
#include "spatial/core/index/rtree/rtree_module.hpp"
#include "spatial/core/index/rtree/rtree_index.hpp"
#include "spatial/core/index/rtree/rtree_index_create_logical.hpp"

namespace spatial {

namespace core {

//-------------------------------------------------------------
// Register
//-------------------------------------------------------------
void RTreeModule::RegisterIndexPlanCreate(DatabaseInstance &db) {

	db.config.operator_extensions.push_back(make_uniq<LogicalCreateRTreeIndexOperatorExtension>());
}

} // namespace core

} // namespace spatial
