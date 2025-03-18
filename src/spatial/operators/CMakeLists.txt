set(EXTENSION_SOURCES
    ${EXTENSION_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_join_logical.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_join_physical.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/spatial_join_optimizer.cpp
    PARENT_SCOPE)