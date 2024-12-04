cmake_minimum_required(VERSION 3.3)

project(GenDisj)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scipgendisj
   src/cmain.c)

# link to math library if it is available
find_library(LIBM m)
if(NOT LIBM)
  set(LIBM "")
endif()

target_link_libraries(scipgendisj ${SCIP_LIBRARIES} ${LIBM})

if( TARGET examples )
    add_dependencies( examples scipgendisj )
endif()
