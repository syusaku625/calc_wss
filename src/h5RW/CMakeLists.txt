add_library(H5 STATIC)
target_sources(H5
  PRIVATE
    hdf5_rw.cpp
)

target_link_libraries(H5 PRIVATE
  ${HDF5_LIBRARIES}
  ${HDF5_CXX_LIBRARIES}
)

target_include_directories(H5
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)