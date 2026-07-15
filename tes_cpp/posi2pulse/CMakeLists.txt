cmake_minimum_required(VERSION 3.20)
project(tes_cpp_posi2pulse VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
find_package(Eigen3 3.3 CONFIG QUIET)
if (NOT Eigen3_FOUND)
    FetchContent_Declare(eigen URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz)
    FetchContent_MakeAvailable(eigen)
endif()
find_package(nlohmann_json 3.11 CONFIG QUIET)
if (NOT nlohmann_json_FOUND)
    FetchContent_Declare(nlohmann_json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# Use a static library so the executable can link on MSVC without requiring
# DLL export annotations for the pulse implementation.
add_library(tes_cpp_pulse STATIC src/pulse.cpp)
find_package(Threads REQUIRED)
target_include_directories(tes_cpp_pulse
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(tes_cpp_pulse PUBLIC Eigen3::Eigen nlohmann_json::nlohmann_json Threads::Threads)

add_executable(posi2pulse src/pulse_main.cpp)
target_link_libraries(posi2pulse PRIVATE tes_cpp_pulse)

install(TARGETS posi2pulse RUNTIME DESTINATION tes_cpp/bin)
install(DIRECTORY include/ DESTINATION include)
