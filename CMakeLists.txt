cmake_minimum_required(VERSION 3.10)
project(PBRPBD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_CLANG_TIDY clang-tidy -checks=-*,readability-*)

find_package(glfw3 CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(TBB CONFIG REQUIRED)
find_package(assimp CONFIG REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/3rdparty/flatbuffers)
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE SRC_LIST ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(PBRPBD ${SRC_LIST})

target_link_libraries(PBRPBD PRIVATE glad::glad) 
target_link_libraries(PBRPBD PRIVATE glfw)
target_link_libraries(PBRPBD PRIVATE glm::glm)
target_link_libraries(PBRPBD PRIVATE assimp::assimp)
target_link_libraries(PBRPBD PRIVATE TBB::tbb TBB::tbbmalloc TBB::tbbmalloc_proxy)