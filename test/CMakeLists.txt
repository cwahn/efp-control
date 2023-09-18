add_executable(efp_math_test efp_math_test.cpp)

target_link_libraries(efp_math_test
    PRIVATE
    Catch2::Catch2WithMain
    efp_math)

catch_discover_tests(efp_math_test)

# add_custom_command(
#     TARGET efp_math_test
#     COMMENT "Run tests"
#     POST_BUILD
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#     COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> "efp_math_test" --output-on-failure
# )
