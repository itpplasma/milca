# Diverted free-boundary solve with mesh refit through gs_free.
# Fixture (test/fixtures/gs_free_xpoint): a synthetic tokamak taken from no
# device (12 rectangular coil groups of 40 turns, Gamma_c on its limiter
# polygon), the profile law of p' and F F' as natural cubic splines through 10
# knots in normalized flux, and a 65 x 65 seed flux, all from a FreeGSNKE 3.1.1
# forward solve of the same case (I_P = 800 kA, lower X-point at
# (1.4447, -1.0059) m, axis at (1.6999, 0.0485) m on its 257 grid).
# Before the refit handled X-point corners (chain spline, saddle on a vertex,
# rays grazing the X-point, readback on the separatrix) this run stopped with
# refit status 9, then 3, 4 and readback status 18.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
find_program(NCDUMP ncdump REQUIRED)
foreach(name boundary.dat coils.txt profile.dat seed.dat)
    file(COPY "${FIXTURES}/${name}" DESTINATION "${WORK}")
endforeach()
file(WRITE "${WORK}/input.nml" "&gs_free
 boundary_file = 'boundary.dat', boundary_nodes = 96, degree = 2, coil_file = 'coils.txt',
 profile_file = 'profile.dat', plasma_ellipse = 1.70, 0.5, 0.9, chain_nodes = 48,
 lcfs_mode = 'xpoint', axis_guess = 1.70, 0.05, xpoint_guess = 1.47, -0.99, search_radius = 0.3,
 currents = -3359.2305625, -1926.3881174, -799.80921591, 5875.5717228, 13410.545564,
 9811.3651916, 15321.140628, -2627.821568, 17994.523203, 3139.9711821, -7504.9622108,
 -732.62318045, error_estimate = .false., seed_file = 'seed.dat'
/
")
execute_process(COMMAND "${APP}" --config input.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 600)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "gs_free failed: ${out}${err}")
endif()

function(nc_value name result)
    execute_process(COMMAND "${NCDUMP}" -v ${name} "${WORK}/output/gs_free.nc"
        OUTPUT_VARIABLE data RESULT_VARIABLE rc)
    string(REGEX MATCH "[\n ]${name} = ([^;]+);" match "${data}")
    if(NOT rc EQUAL 0 OR NOT match)
        message(FATAL_ERROR "missing result variable ${name}")
    endif()
    string(REGEX REPLACE "[ \n,]+" ";" values "${CMAKE_MATCH_1}")
    string(REGEX REPLACE "^;|;$" "" values "${values}")
    set(${result} "${values}" PARENT_SCOPE)
endfunction()

find_package(Python3 REQUIRED COMPONENTS Interpreter)
function(require_near what a b tolerance)
    execute_process(COMMAND "${Python3_EXECUTABLE}" -I -c
        "import sys; sys.exit(0 if abs(float(sys.argv[1])-float(sys.argv[2])) <= float(sys.argv[3]) else 1)"
        "${a}" "${b}" "${tolerance}" RESULT_VARIABLE near_status)
    if(NOT near_status EQUAL 0)
        message(FATAL_ERROR "${what}: ${a} differs from ${b} by more than ${tolerance}")
    endif()
endfunction()

nc_value(refits refits)
if(refits LESS 1)
    message(FATAL_ERROR "no mesh refit onto the diverted boundary")
endif()
nc_value(q_status q_status)
if(NOT q_status EQUAL 0)
    message(FATAL_ERROR "flux-surface readback failed, status ${q_status}")
endif()
nc_value(xpoint xpoint)
list(GET xpoint 0 xr)
list(GET xpoint 1 xz)
require_near("X-point R" "${xr}" 1.4447 0.005)
require_near("X-point Z" "${xz}" -1.0059 0.005)
nc_value(axis axis)
list(GET axis 0 ar)
list(GET axis 1 az)
require_near("axis R" "${ar}" 1.6999 0.02)
require_near("axis Z" "${az}" 0.0485 0.02)
# FreeGSNKE gives 800 kA with these profiles; the 1.6 % difference is an open
# cross-code question, so the gate is 3 %.
nc_value(plasma_current ip)
require_near("plasma current" "${ip}" 800000 24000)
