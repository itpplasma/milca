cmake_minimum_required(VERSION 3.22)

# Layering test (cmake -P). Checks the Fortran use graph of every layered
# source against the partition in PARTITION ("layer=source" lines written at
# configure time): a source uses only modules of its own layer or of a lower
# layer; among the L2 model classes only mhd uses equilibrium; no layer uses an
# application module.
set(rank_field 0)
set(rank_core 1)
set(rank_fem2d 2)
set(rank_pairs 2)
foreach(class equilibrium mhd classical quantum continuum)
    set(rank_${class} 3)
endforeach()
set(rank_apps 4)
set(sideways_mhd equilibrium)

file(STRINGS "${PARTITION}" entries)
foreach(entry IN LISTS entries)
    string(REGEX REPLACE "^([a-z0-9]+)=.*$" "\\1" layer "${entry}")
    string(REGEX REPLACE "^[a-z0-9]+=" "" source "${entry}")
    list(APPEND sources_${layer} "${source}")
    list(APPEND layers ${layer})
    file(STRINGS "${source}" lines REGEX
        "^[ \t]*[Mm][Oo][Dd][Uu][Ll][Ee][ \t]+[A-Za-z0-9_]+[ \t]*(!.*)?$")
    foreach(line IN LISTS lines)
        string(REGEX REPLACE "^[ \t]*[A-Za-z]+[ \t]+([A-Za-z0-9_]+).*$" "\\1"
            module "${line}")
        string(TOLOWER "${module}" module)
        set(module_layer_${module} ${layer})
        set(module_source_${module} "${source}")
    endforeach()
endforeach()
list(REMOVE_DUPLICATES layers)

set(violations 0)
foreach(layer IN LISTS layers)
    if(layer STREQUAL "apps")
        continue()
    endif()
    foreach(source IN LISTS sources_${layer})
        # Number lines exactly: neutralise list separators and brackets first.
        file(READ "${source}" text)
        string(REGEX REPLACE "[];[]" " " text "${text}")
        string(REPLACE "\n" ";" lines "${text}")
        set(number 0)
        foreach(line IN LISTS lines)
            math(EXPR number "${number} + 1")
            if(NOT line MATCHES
                    "^[ \t]*([Uu][Ss][Ee]|[Ss][Uu][Bb][Mm][Oo][Dd][Uu][Ll][Ee])[ \t,:(]")
                continue()
            endif()
            string(REGEX REPLACE
                "^[ \t]*[A-Za-z]+[ \t]*(,[ \t]*[A-Za-z_]+[ \t]*)?(::)?[ \t(]*([A-Za-z0-9_]+).*$"
                "\\3" module "${line}")
            string(TOLOWER "${module}" module)
            if(NOT DEFINED module_layer_${module})
                continue()
            endif()
            set(used ${module_layer_${module}})
            if(used STREQUAL layer OR rank_${used} LESS rank_${layer})
                continue()
            endif()
            if(used IN_LIST sideways_${layer})
                continue()
            endif()
            message(SEND_ERROR "${source}:${number}: ${layer} uses ${module} "
                "of ${used} (${module_source_${module}})")
            math(EXPR violations "${violations} + 1")
        endforeach()
    endforeach()
endforeach()

list(LENGTH entries nsources)
if(violations GREATER 0)
    message(FATAL_ERROR "${violations} layering violation(s)")
endif()
message(STATUS "milca layering: ${nsources} sources")
