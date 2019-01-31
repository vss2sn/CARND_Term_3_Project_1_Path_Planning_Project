project(Path_Planning)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

option(CONSERVATIVE "CONSERVATIVE" ON)
if(CONSERVATIVE)
    add_definitions(-DCONSERVATIVE)
    message(STATUS "Building with conservative settings")
endif(CONSERVATIVE)

set(CXX_FLAGS "-w")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp)

set(uWS_CUSTOM_INSTALL_DIR "YES")  # NO if system-wide install
# NOTE:
# Minor addition for convenience, I have not installed
# the uWebSockets library system wide.
if(${uWS_CUSTOM_INSTALL_DIR} MATCHES "YES")

include_directories(../uWebSockets/bin/include)
link_directories(../uWebSockets/bin/lib64)

include_directories(../Ipopt-3.12.7/install/include)
link_directories(../Ipopt-3.12.7/install/lib)

endif(${uWS_CUSTOM_INSTALL_DIR} MATCHES "YES")


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


add_executable(path_planning ${sources})

target_link_libraries(path_planning z ssl uv uWS)
