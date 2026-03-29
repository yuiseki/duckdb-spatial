#pragma once
#include "duckdb/common/string.hpp"

namespace duckdb {

class ExtensionLoader;

bool IdentifyProjCRS(const char *crs, string &auth_name, string &auth_code);
void RegisterProjModule(ExtensionLoader &db);

} // namespace duckdb
