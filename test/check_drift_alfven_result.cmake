# Drift-Alfven application on a small case near the Hasegawa-Wakatani limit
# (beta = mu = 0.01, c = 1, kpar = 1, so alpha = kpar^2/c = 1), run by CTest
# with -P. Inputs: APP (drift_alfven), WORK. The saturated state must carry a
# positive particle flux and close its energy and enstrophy balances.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
file(WRITE "${WORK}/case.nml" "&drift_alfven
 beta = 0.01, mu = 0.01, c = 1.0, kpar = 1.0, kappa = 1.0, d = 2.0e-3, order = 2,
 nx = 32, ny = 32, k0 = 0.3, dt = 0.025, time_end = 500.0, average_from = 100.0,
 sample_steps = 16, segment = 256, modes = 10, amplitude = 1.0e-2, seed = 3,
 rho_s_cm = 0.2, ln_over_cs_s = 1.0e-6, te_over_e_statvolt = 2.5 /
")
execute_process(COMMAND "${APP}" --config case.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "${out}")
if(NOT code EQUAL 0)
    message(FATAL_ERROR "drift_alfven failed:\n${out}${err}")
endif()
string(REGEX MATCH "Gamma_n = +([-0-9.E+]+)" m "${out}")
set(flux "${CMAKE_MATCH_1}")
string(REGEX MATCH "E += +([-0-9.E+]+)" m "${out}")
set(energy "${CMAKE_MATCH_1}")
string(REGEX MATCH "balance residuals \\(E, W\\): +([-0-9.E+]+) +([-0-9.E+]+)" m "${out}")
set(res_e "${CMAKE_MATCH_1}")
set(res_w "${CMAKE_MATCH_2}")
if(flux STREQUAL "" OR res_e STREQUAL "" OR res_w STREQUAL "")
    message(FATAL_ERROR "output lacks Gamma_n or the balance residuals:\n${out}")
endif()
if(NOT flux GREATER 0)
    message(FATAL_ERROR "no positive particle flux: ${flux}")
endif()
foreach(r ${res_e} ${res_w})
    string(REGEX REPLACE "^-" "" a "${r}")
    if(a GREATER 0.05)
        message(FATAL_ERROR "balance residual ${r} exceeds 0.05")
    endif()
endforeach()
