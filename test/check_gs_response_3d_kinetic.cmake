# Kinetic closure of the rational surface q = 2/1 of the A = 10 circular
# torus from synthetic profiles: the layer response is computed and the
# solution meets the matching condition jump of d Phi/d psin = Delta A.
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}/prof")
file(WRITE "${WORK}/coils.txt"
    "1\nloop 4\n8 -1 -0.5\n8 1 -0.5\n8 1 0.5\n8 -1 0.5\n")
file(WRITE "${WORK}/currents.txt" "1000.0\n")
file(WRITE "${WORK}/prof/n.dat" "2 3.877842551e+13\n3 3.75071175e+13\n4 3.627748793e+13\n5 3.508817041e+13\n6 3.393784336e+13\n7 3.282522851e+13\n8 3.174908952e+13\n9 3.070823056e+13\n10 2.970149501e+13\n11 2.872776418e+13\n12 2.778595605e+13\n13 2.687502406e+13\n14 2.599395597e+13\n15 2.514177273e+13\n16 2.431752738e+13\n")
file(WRITE "${WORK}/prof/Te.dat" "2 440.8842964\n3 419.3821156\n4 398.9286084\n5 379.4726306\n6 360.965532\n7 343.3610353\n8 326.61512\n9 310.6859126\n10 295.5335819\n11 281.120239\n12 267.4098432\n13 254.3681112\n14 241.9624321\n15 230.161785\n16 218.9366623\n")
file(WRITE "${WORK}/prof/Ti.dat" "2 440.8842964\n3 419.3821156\n4 398.9286084\n5 379.4726306\n6 360.965532\n7 343.3610353\n8 326.61512\n9 310.6859126\n10 295.5335819\n11 281.120239\n12 267.4098432\n13 254.3681112\n14 241.9624321\n15 230.161785\n16 218.9366623\n")
file(WRITE "${WORK}/prof/q.dat" "2 1.1915\n3 1.2965\n4 1.4015\n5 1.5065\n6 1.6115\n7 1.7165\n8 1.8215\n9 1.9265\n10 2.0315\n11 2.1365\n12 2.2415\n13 2.3465\n14 2.4515\n15 2.5565\n16 2.6615\n")
file(WRITE "${WORK}/prof/Er.dat" "2 -0.1\n3 -0.1\n4 -0.1\n5 -0.1\n6 -0.1\n7 -0.1\n8 -0.1\n9 -0.1\n10 -0.1\n11 -0.1\n12 -0.1\n13 -0.1\n14 -0.1\n15 -0.1\n16 -0.1\n")
file(WRITE "${WORK}/input.nml"
    "&l4 eqdsk='${EQDSK}', coils='coils.txt', currents(1)='currents.txt', prefix='t', n=1, rat_num(1)=2, rat_den(1)=1, nrad=12, nth_max=64, order=2, nphi=32, nfine=256, np_sample=64, model='ideal', gmres_tol=1e-10, gmres_maxit=40, rat_auto=.false., closure='kinetic', kinetic_profiles='prof', kinetic_btor=2e4, kinetic_r0=164.64, kinetic_core=1.0, kinetic_band=2.0, kinetic_nwin=128, kinetic_mmax=48 /\n")
execute_process(COMMAND "${APP}" input.nml WORKING_DIRECTORY "${WORK}"
    RESULT_VARIABLE status OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 200)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "kinetic closure failed: ${out}${err}")
endif()
if(NOT out MATCHES "kinetic layer q  -2.0000: r_s\\[cm\\]    9.7000")
    message(FATAL_ERROR "rational not placed at |q| = 2 of the profiles: ${out}")
endif()
if(NOT out MATCHES "matching \\|jump - Delta A\\|/\\|jump\\| +[0-9.]+E-(0[89]|1[0-9])")
    message(FATAL_ERROR "matching condition not met: ${out}")
endif()
