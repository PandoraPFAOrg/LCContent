add_executable(LCContentSmokeTest LCContentSmokeTest.cc)
target_link_libraries(LCContentSmokeTest PRIVATE ${PROJECT_NAME}::${PROJECT_NAME})

add_test(NAME LCContentSmokeTest COMMAND LCContentSmokeTest)

add_executable(LCThetaEnergyCorrectionTest LCThetaEnergyCorrectionTest.cc)
target_link_libraries(LCThetaEnergyCorrectionTest PRIVATE ${PROJECT_NAME}::${PROJECT_NAME})

# The theta-energy lookup registry is process-global and cannot be cleared, so each
# case must run in its own process.
foreach(case no-table lookup out-of-range direction validation type-isolation
             name-ambiguity multi-instance)
    add_test(NAME LCThetaEnergyCorrection.${case}
             COMMAND LCThetaEnergyCorrectionTest ${case})
endforeach()
