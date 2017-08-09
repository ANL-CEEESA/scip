project(GMI)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

add_executable(gmi
   src/cmain.c
   src/sepa_gmi.c)

target_link_libraries(gmi libscip m)
add_dependencies(examples gmi)

include(CTest)

#
# define the instance sets
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)
#

set(instances
    "instances/MIP/dcmulti.mps\;188182"
    "instances/MIP/misc03.mps\;3360"
    "instances/MIP/stein27.fzn\;18"
    "instances/MINLP/parincQuadratic.osil\;49920.5564"
    "instances/MINLP/tltr.mps\;48.0666666667"
    "instances/PseudoBoolean/factor-mod-size=9-P0=67-P1=349-P2=67-P3=499-P4=79-P5=347-P6=307-B.opb\;3"
    "instances/PseudoBoolean/normalized-bsg_10_4_5.opb\;-4"
    "instances/PseudoBoolean/normalized-mds_10_4_3.opb\;2"
    )
set(settings
    "gmionly"
    "gmirep"
    "gmi"
    "scipdefault")

#
# loop over the instances
#
foreach(instance ${instances})
    list(GET instance 0 path)
    list(GET instance 1 optval)
    get_filename_component(basename ${path} NAME)
    foreach(setting ${settings})
        #
        # treat the instance as a tuple (list) of two values
        #
        add_test(NAME examples-gmi-${setting}-${basename} COMMAND $<TARGET_FILE:gmi> ${CMAKE_CURRENT_SOURCE_DIR}/../../check/${path} ${CMAKE_CURRENT_SOURCE_DIR}/settings/${setting}.set -o ${optval} ${optval})
    endforeach(setting)
endforeach(instance)