add_library(Base STATIC)
target_sources(Base
  PRIVATE
    fem_base.cpp
    math_tools.cpp
    fileIO.cpp
)

target_include_directories(Base
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)