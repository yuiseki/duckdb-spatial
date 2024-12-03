#pragma once
#include "duckdb.hpp"
#include "duckdb/main/extension_util.hpp"

namespace spatial {

struct DocTag {
	const char *key;
	const char *value;
};

struct DocUtil {
	static void AddDocumentation(duckdb::DatabaseInstance &db, const char *function_name, const char *description,
	                             const char *example,
	                             const duckdb::unordered_map<duckdb::string, duckdb::string> &tags, duckdb::vector<duckdb::string> parameter_names = {});

	// Abuse adding tags as a comment
	template <size_t N>
	static void AddDocumentation(duckdb::DatabaseInstance &db, const char *function_name, const char *description,
	                             const char *example, const DocTag (&tags)[N], duckdb::vector<duckdb::string> parameter_names = {}) {
		duckdb::unordered_map<duckdb::string, duckdb::string> tag_map;
		for (size_t i = 0; i < N; i++) {
			tag_map[tags[i].key] = tags[i].value;
		}
		AddDocumentation(db, function_name, description, example, tag_map, parameter_names);
	}
};

} // namespace spatial
