cmake_minimum_required(VERSION 3.16)
project(balrog)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

file(GLOB globbed
        "*.h"
        "*.cpp"
        )

include_directories(include)
include_directories(include/gzip)

set(Torch_DIR libtorch/share/cmake/Torch)
find_package(Torch REQUIRED)

set(CMAKE_INSTALL_RPATH "${CMAKE_CURRENT_LIST_DIR}/libtorch/lib")

add_executable(balrog main.cpp ${globbed} include/tqdm.h include/FastaReader.cpp include/FastaReader.h include/GeneFinder.cpp include/GeneFinder.h)

target_link_libraries(balrog ${TORCH_LIBRARIES})
set_property(TARGET balrog PROPERTY CXX_STANDARD 14)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( balrog ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

include(CMakeRC)
cmrc_add_resource_library(
        cmakeresources
        NAMESPACE cmakeresources
        WHENCE data
        data/gene_model_v1.0.pt
        data/TIS_model_v1.0.pt
        data/reference_genes.tar.gz
        )

target_link_libraries(${PROJECT_NAME} cmakeresources)

install(
  TARGETS balrog
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)
