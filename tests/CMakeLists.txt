cmake_minimum_required(VERSION 2.8)

add_executable(test_api api_tests.cpp)

target_link_libraries(test_api fastahack_core pthread ${ZLIB_LIBRARIES})

target_link_libraries(test_api ${LIBS})
