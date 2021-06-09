add_cython_target(composite_measure CXX)

add_library(composite_measure MODULE ${composite_measure})

python_extension_module(composite_measure)

target_link_libraries(composite_measure composite-measure ${PYTHON_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

target_include_directories(composite_measure PUBLIC ${CMAKE_SOURCE_DIR}/api ${PYTHON_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS} )

install(TARGETS composite_measure LIBRARY DESTINATION composite)
