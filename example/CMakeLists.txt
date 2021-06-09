add_executable(composite_example composite_example.cpp)

target_include_directories(composite_example PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/AudioFile ${CMAKE_SOURCE_DIR}/api)

target_link_libraries(composite_example composite-measure ${OpenMP_CXX_LIBRARIES})
