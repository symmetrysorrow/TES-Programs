cmake_minimum_required(VERSION 3.25)
project(positron_10000_raw_audit LANGUAGES CXX)
add_executable(extract_selected_dumpall
    extract_selected_dumpall.cpp
    ../../../src/dumpall_record_parser.cpp)
target_include_directories(extract_selected_dumpall PRIVATE ../../../include)
target_compile_features(extract_selected_dumpall PRIVATE cxx_std_17)

add_executable(extract_event_tracks
    extract_event_tracks.cpp
    ../../../src/dumpall_record_parser.cpp)
target_include_directories(extract_event_tracks PRIVATE ../../../include)
target_compile_features(extract_event_tracks PRIVATE cxx_std_17)
