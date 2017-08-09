include(CTest)

#
# instances of varying size in the Binpacking data set.
#
set(instances
    "u20_00\;9"
    "u40_00\;17"
    "u60_00\;27"
    "u80_00\;33"
#    "u100_00\;40"
#    "u120_00\;48"
# the remaining 2 from this set take too long.
)

#
# add a test for every instance
#
foreach(instance ${instances})
    list(GET instance 0 basename)
    list(GET instance 1 optval)
    #
    # call the Binpacking binary and validate the solve with the given objective value
    #
    add_test(NAME "examples-binpacking-${basename}" COMMAND $<TARGET_FILE:binpacking> -f ${CMAKE_CURRENT_SOURCE_DIR}/../data/${basename}.bpa -o ${optval} ${optval})
    set_tests_properties("examples-binpacking-${basename}" PROPERTIES PASS_REGULAR_EXPRESSION "Validation         : Success")
endforeach()