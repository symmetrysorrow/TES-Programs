cmake_minimum_required(VERSION 3.20)
project(tes_cpp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Both native tools are distributed by the one ``tes_cpp`` Python package.
# Each component keeps its own CMakeLists.txt so it can also be built alone.
set(TES_CPP_DUMP2JSON_INSTALL_DIR "tes_cpp/bin" CACHE STRING
    "Install directory for dump2json executable" FORCE)
add_subdirectory(dump2json)
add_subdirectory(posi2pulse)
