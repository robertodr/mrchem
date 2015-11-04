include(DefaultPathsMRCPP)
set_property(GLOBAL PROPERTY TestSources)
add_subdirectory(mwrepr)

get_property(test_sources GLOBAL PROPERTY TestSources)

include_directories(${PROJECT_SOURCE_DIR}/external/catch)
include_directories(${CMAKE_SOURCE_DIR}/src/mwreprlib)

add_executable(unit_tests.x unit_tests.cpp ${test_sources})

target_link_libraries(unit_tests.x mwrepr)
target_link_libraries(unit_tests.x ${BLAS_LIBRARIES})
target_link_libraries(unit_tests.x ${Boost_LIBRARIES})