add_executable(efp_control_all_test test.cpp)

target_link_libraries(efp_control_all_test
    PRIVATE
    Catch2::Catch2WithMain
    efp-control)

catch_discover_tests(efp_control_all_test)

add_custom_command(
    TARGET efp_control_all_test
    COMMENT "Run tests"
    POST_BUILD
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> "all_test" --output-on-failure
)
