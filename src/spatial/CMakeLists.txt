add_subdirectory(util)
add_subdirectory(modules)
add_subdirectory(geometry)
add_subdirectory(index)

set(EXTENSION_SOURCES
    ${EXTENSION_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_extension.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_types.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_optimizers.cpp
PARENT_SCOPE)