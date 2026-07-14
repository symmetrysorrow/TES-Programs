cmake_minimum_required(VERSION 3.20)
project(tes_cpp_dump2json VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(tes_cpp_dump2json src/converter.cpp)
target_include_directories(tes_cpp_dump2json
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_compile_features(tes_cpp_dump2json PUBLIC cxx_std_17)

add_executable(dump2json src/main.cpp)
target_link_libraries(dump2json PRIVATE tes_cpp_dump2json)

include(GNUInstallDirs)
install(TARGETS tes_cpp_dump2json
    EXPORT tes_cpp_dump2jsonTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# scikit-build-core uses this install location when producing a Python wheel.
set(TES_CPP_DUMP2JSON_INSTALL_DIR "tes_cpp_dump2json/bin" CACHE STRING "Install directory for dump2json executable")
install(TARGETS dump2json RUNTIME DESTINATION ${TES_CPP_DUMP2JSON_INSTALL_DIR})
