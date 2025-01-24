PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Configuration of extension
EXT_NAME=excel
EXT_CONFIG=${PROJ_DIR}extension_config.cmake

CORE_EXTENSIONS='httpfs'

# Include the Makefile from extension-ci-tools
include extension-ci-tools/makefiles/duckdb_extension.Makefile

#### Override the included format target because we have different source tree layout
format:
	find spatial/src -iname *.hpp -o -iname *.cpp | xargs clang-format --sort-includes=0 -style=file -i
	find spatial/include -iname *.hpp -o -iname *.cpp | xargs clang-format --sort-includes=0 -style=file -i
	cmake-format -i CMakeLists.txt