add_library(composite-measure abs.cpp composite.cpp composite_cpp.cpp composite_cpp_data.cpp composite_cpp_initialize.cpp composite_cpp_rtwutil.cpp composite_cpp_terminate.cpp 
  fft.cpp FFTImplementationCallback.cpp filter.cpp find.cpp gencoswin.cpp hann.cpp ifft.cpp interp1.cpp  
  log2.cpp minOrMax.cpp pesq_original_cpp.cpp relop.cpp rtGetInf.cpp rtGetNaN.cpp rt_nonfinite.cpp sort.cpp sortIdx.cpp sprintf.cpp)

target_include_directories(composite-measure PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/;${CMAKE_SOURCE_DIR}/api/>" "$<INSTALL_INTERFACE:>")

install(FILES ${CMAKE_SOURCE_DIR}/api/composite.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS composite-measure EXPORT composite_export DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT composite_export FILE composite-config.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR})
