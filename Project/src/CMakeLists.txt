# src/CMakeLists.txt

# Aggiungi i file sorgente
list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp")

list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/Import_Export.cpp")

list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/Finding_Traces.cpp")


#headers
list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp")

list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/namespace.hpp")

list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/inline.hpp")


list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/FracturesTraces.hpp")

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})


set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
