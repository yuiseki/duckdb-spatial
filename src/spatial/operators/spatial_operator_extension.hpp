#pragma once

namespace duckdb {

class DatabaseInstance;
void RegisterSpatialOperatorExtension(DatabaseInstance &db);

}