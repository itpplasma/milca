# gs_phase1 readback near a steep edge: on a D-shaped boundary (kappa 1.8, delta 0.4)
# with p' and FF' proportional to 1 - s^16, the contour toroidal flux must agree with
# the region integral (flux_consistency) to 8e-6. Levels packed towards the edge give
# 2.9e-6, the former uniform rho levels 2.0e-5.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}")
find_program(NCDUMP ncdump REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
execute_process(COMMAND "${Python3_EXECUTABLE}" -I -c "
import math
n = 256
with open('boundary.dat', 'w') as f:
    f.write(f'{n}\\n')
    for k in range(n):
        t = 2*math.pi*k/n
        f.write(f'{1.7 + 0.6*math.cos(t + 0.4*math.sin(t)):.15g} {1.08*math.sin(t):.15g}\\n')
m, pieces = 16, 16
knots = [i/pieces for i in range(pieces + 1)]
def hermite(scale, x0, x1):
    # C1 cubic through the values and slopes of -scale (1 - s^m) at the knots
    f = lambda x: -scale*(1 - x**m)
    d = lambda x: scale*m*x**(m - 1)
    h = x1 - x0
    slope = (f(x1) - f(x0))/h
    return [f(x0), d(x0), (3*slope - 2*d(x0) - d(x1))/h, (d(x0) + d(x1) - 2*slope)/h**2]
with open('edge.law', 'w') as f:
    f.write(f'{len(knots)} 1\\n0 3 -0.4\\n')
    f.write(' '.join(repr(x) for x in knots) + '\\n')
    for x0, x1 in zip(knots[:-1], knots[1:]):
        f.write(' '.join(repr(v) for v in hermite(4.0e4, x0, x1) + hermite(1.0, x0, x1)) + '\\n')
" WORKING_DIRECTORY "${WORK}" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "fixture generation failed")
endif()
file(WRITE "${WORK}/input.nml" "&gs_phase1 boundary_file='boundary.dat', profile_file='edge.law', degree=2, boundary_nodes=64, curved=.true. /\n")
execute_process(COMMAND "${APP}" --config input.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 120)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "gs_phase1 failed: ${out}${err}")
endif()
execute_process(COMMAND "${NCDUMP}" -v flux_consistency "${WORK}/output/gs_phase1.nc"
    OUTPUT_VARIABLE data RESULT_VARIABLE status)
string(REGEX MATCH "[\n ]flux_consistency = ([^;]+);" match "${data}")
if(NOT match)
    message(FATAL_ERROR "missing flux_consistency")
endif()
string(STRIP "${CMAKE_MATCH_1}" value)
execute_process(COMMAND "${Python3_EXECUTABLE}" -I -c
    "import sys; v = abs(float(sys.argv[1])); print('flux_consistency', v); sys.exit(0 if v <= 8e-6 else 1)"
    "${value}" RESULT_VARIABLE near)
if(NOT near EQUAL 0)
    message(FATAL_ERROR "contour and region toroidal flux differ by ${value} (limit 8e-6)")
endif()
