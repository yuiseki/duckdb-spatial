# This file is included by DuckDB's build system. It specifies which extension to load

# Extension from this repo
duckdb_extension_load(spatial
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}
    INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/spatial/include
    LOAD_TESTS
    DONT_LINK
    LINKED_LIBS "../../deps/local/lib/*.a"
)
