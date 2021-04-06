cmake_minimum_required (VERSION 3.10)
project (GarbledPoints)

set(CMAKE_BUILD_TYPE Release)

add_definitions(-Wfatal-errors)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-pthread -Wall -march=native -O3 -maes -mrdseed -DNDEBUG")

if (APPLE)
  set(OPENSSL_ROOT_DIR "/usr/local/opt/openssl")
endif()
find_package(OpenSSL REQUIRED)

file(GLOB sources "src/*.cc")

add_executable(GarbledPoints ${sources})
target_include_directories(GarbledPoints PUBLIC src ${OPENSSL_INCLUDE_DIR} emp-ot/src)
target_link_libraries(GarbledPoints ${OPENSSL_LIBRARIES})

