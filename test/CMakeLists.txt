add_executable(LCContentSmokeTest LCContentSmokeTest.cc)
target_link_libraries(LCContentSmokeTest PRIVATE ${PROJECT_NAME}::${PROJECT_NAME})

add_test(NAME LCContentSmokeTest COMMAND LCContentSmokeTest)
