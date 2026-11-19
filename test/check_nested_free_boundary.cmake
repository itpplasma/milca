# Free-boundary nested_equilibrium on a synthetic mgrid (test_mgrid_field
# tokamak): the solve with coupling_margin 0.6 must converge, and the
# solve with coupling_margin 0.3 must give the same equilibrium within the
# tolerances of test_nested_equilibrium_app (axis 1.5 mm, wb 2e-5, volume
# 1e-4, edge iota 1e-3). The exterior Steklov-Poincare coupling is exact,
# so the equilibrium does not depend on where Gamma_c lies.
file(MAKE_DIRECTORY "${WORK}")
execute_process(COMMAND "${MGRID}" tokamak "${WORK}/mgrid_tokamak.nc"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "synthetic mgrid failed: ${output}${errors}")
endif()
file(READ "${CONFIG}" deck)
set(values "")
foreach(margin 0.6 0.3)
    file(MAKE_DIRECTORY "${WORK}/margin-${margin}/output")
    string(REPLACE "extcur(1)" "coupling_margin = ${margin}, extcur(1)" text "${deck}")
    file(WRITE "${WORK}/margin-${margin}/input.nml" "${text}")
    file(CREATE_LINK "${WORK}/mgrid_tokamak.nc" "${WORK}/margin-${margin}/mgrid_tokamak.nc"
        COPY_ON_ERROR)
    execute_process(COMMAND "${APP}" --config input.nml
        WORKING_DIRECTORY "${WORK}/margin-${margin}"
        RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "nested_equilibrium failed: ${output}${errors}")
    endif()
    if(values STREQUAL "")
        # first solve: its raxis, wb, volume and edge iota become the reference
        execute_process(COMMAND "${CHECK}"
            "${WORK}/margin-${margin}/output/nested_equilibrium.nc" 0 1 1 0
            OUTPUT_VARIABLE check_output ERROR_VARIABLE errors)
        string(REGEX MATCH "iota edge: +([^ ]+) +([^ ]+) +([^ ]+) +([^ \n]+)" found
            "${check_output}")
        set(values "${CMAKE_MATCH_2} ${CMAKE_MATCH_3} ${CMAKE_MATCH_1} ${CMAKE_MATCH_4}")
        if(NOT found OR errors MATCHES "not converged")
            message(FATAL_ERROR "no converged result: ${check_output}${errors}")
        endif()
        message(STATUS "margin ${margin}: raxis, wb, volume, iota edge ${values}")
    else()
        separate_arguments(expected UNIX_COMMAND "${values}")
        execute_process(COMMAND "${CHECK}"
            "${WORK}/margin-${margin}/output/nested_equilibrium.nc" ${expected}
            RESULT_VARIABLE status OUTPUT_VARIABLE check_output ERROR_VARIABLE errors)
        if(NOT status EQUAL 0)
            message(FATAL_ERROR "result check failed: ${check_output}${errors}")
        endif()
        message(STATUS "${check_output}")
    endif()
endforeach()
