#pragma once

#include "duckdb/common/string.hpp"
#include "duckdb/common/typedefs.hpp"

namespace duckdb {

class ExtensionLoader;
class ClientContext;

struct SpatialSettings {

	static bool AlwaysXY(ClientContext &context, bool &is_set);

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
