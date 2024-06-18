list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/TestFractures.hpp")
#list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/TestFractures_partTwo.hpp")


set(src_test_sources ${src_sources} PARENT_SCOPE)
set(src_test_headers ${src_headers} PARENT_SCOPE)
set(src_test_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
