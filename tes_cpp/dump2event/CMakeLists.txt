cmake_minimum_required(VERSION 3.20)
project(tes_cpp_dump2event VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# PHITS Dumpall structural reader + physically reconstructed Pb EGS5 ledger.
# The obsolete reaction-difference converter is no longer built.
add_library(tes_cpp_dump2event STATIC
    src/dumpall_record_parser.cpp
    src/physical_deposit.cpp)
target_include_directories(tes_cpp_dump2event
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
target_compile_features(tes_cpp_dump2event PUBLIC cxx_std_17)

add_executable(dump2event src/main.cpp)
target_link_libraries(dump2event PRIVATE tes_cpp_dump2event)

add_executable(dumpall_record_inspect src/dumpall_inspect.cpp)
target_link_libraries(dumpall_record_inspect PRIVATE tes_cpp_dump2event)

add_executable(dump2event_physical_reconstruct src/physical_dump_reconstruct.cpp)
target_link_libraries(dump2event_physical_reconstruct PRIVATE tes_cpp_dump2event)

add_executable(dumpall_transport_export src/dumpall_transport_export.cpp)
target_link_libraries(dumpall_transport_export PRIVATE tes_cpp_dump2event)

include(CTest)
option(TES_BUILD_DUMP2EVENT_TESTS "Build dump2event unit tests" ON)
if (TES_BUILD_DUMP2EVENT_TESTS)
    enable_testing()
    add_executable(test_physical_deposit_engine tests/test_physical_deposit_engine.cpp)
    target_link_libraries(test_physical_deposit_engine PRIVATE tes_cpp_dump2event)
    add_test(NAME physical_deposit_engine COMMAND test_physical_deposit_engine)
endif()

include(GNUInstallDirs)
install(TARGETS tes_cpp_dump2event
    EXPORT tes_cpp_dump2eventTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# scikit-build-core uses this install location when producing a Python wheel.
set(TES_CPP_DUMP2EVENT_INSTALL_DIR "tes_cpp_dump2event/bin" CACHE STRING "Install directory for dump2event executable")
install(TARGETS dump2event RUNTIME DESTINATION ${TES_CPP_DUMP2EVENT_INSTALL_DIR})
install(TARGETS dump2event_physical_reconstruct dumpall_transport_export
    RUNTIME DESTINATION ${TES_CPP_DUMP2EVENT_INSTALL_DIR})
install(FILES validation/build_physical_hdf5_v2.py
    DESTINATION ${CMAKE_INSTALL_DATADIR}/tes_cpp_dump2event)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../python/tes_cpp"
    DESTINATION python
    FILES_MATCHING PATTERN "*.py"
    PATTERN "__pycache__" EXCLUDE)

# Package a strictly byte-verified Pb HATCH resource when its applicable
# repository storage/use permission and access scope are confirmed.
# Do not fabricate or download unrelated stopping-power coefficients.
set(TES_VERIFIED_PB_TABLE "${CMAKE_CURRENT_SOURCE_DIR}/data/pb_egs5_phits3341.csv")
if (EXISTS "${TES_VERIFIED_PB_TABLE}")
    install(FILES "${TES_VERIFIED_PB_TABLE}" DESTINATION "${TES_CPP_DUMP2EVENT_INSTALL_DIR}")
endif()
