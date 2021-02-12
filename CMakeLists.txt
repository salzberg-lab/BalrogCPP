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


function(download_file url filename hash_type hash)

if(NOT EXISTS ${filename})
  file(DOWNLOAD ${url} ${filename}
       TIMEOUT 60  # seconds
       EXPECTED_HASH ${hash_type}=${hash}
       TLS_VERIFY ON)
endif()

endfunction(download_file)

if (APPLE)
    download_file(
        https://download.pytorch.org/libtorch/cpu/libtorch-macos-1.7.1.zip
        ./libtorch-macos-1.7.1.zip
        SHA256 75e032a3a65b3dc97feffdc8dc7c2d7c6504e80e63cdefb7ae2ec7458b3ad7e4
        )
    execute_process(COMMAND unzip libtorch-macos-1.7.1.zip)
    execute_process(COMMAND mv libtorch ..)
else()
    download_file(
        #https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.5.0%2Bcpu.zip
        #./libtorch-cxx11-abi-shared-with-deps-1.5.0+cpu.zip
        #SHA256 3e438237a08099a4bf014335cd0da88708da3a1678aec12a46c67305792b5fa4
        
        #https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.5.0%2Bcpu.zip
        #./libtorch-shared-with-deps-1.5.0+cpu.zip
        #SHA256 db3545b0d2b144db4292c2f0bec236febec44aa658dd54f6b3532f2848c50c8a
        
        #https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.6.0%2Bcpu.zip
        #./libtorch-shared-with-deps-1.6.0+cpu.zip
        #SHA256 https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.5.0%2Bcpu.zip
        
        https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.7.1%2Bcpu.zip
        ./libtorch-cxx11-abi-shared-with-deps-1.7.1+cpu.zip
        SHA256 cc2386c30603c0f145904450a2a28f04b8047f3e624e40f859eebc858ff8f2e9
        
        
    )
    #execute_process(COMMAND unzip libtorch-cxx11-abi-shared-with-deps-1.5.0+cpu.zip)
    #execute_process(COMMAND unzip libtorch-shared-with-deps-1.5.0+cpu.zip)
    execute_process(COMMAND unzip libtorch-cxx11-abi-shared-with-deps-1.7.1+cpu.zip)
    
    
    execute_process(COMMAND mv libtorch ..)
endif()


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
        data/reference_genes.fasta
        )

target_link_libraries(${PROJECT_NAME} cmakeresources)

install(
  TARGETS balrog
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)
