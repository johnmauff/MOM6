set(TURBO_INFRA "FMS2" CACHE STRING "Infrastructure backend: FMS2 or TIM")
set_property(CACHE TURBO_INFRA PROPERTY STRINGS FMS2 TIM)

set(TURBO_MEMORY_MODE "dynamic_symmetric" CACHE STRING
    "MOM6 memory layout: dynamic_symmetric or dynamic_nonsymmetric")
set_property(CACHE TURBO_MEMORY_MODE PROPERTY STRINGS
    dynamic_symmetric dynamic_nonsymmetric)

if(NOT TURBO_INFRA MATCHES "^(FMS2|TIM)$")
    message(FATAL_ERROR
        "Unknown TURBO_INFRA='${TURBO_INFRA}'. Valid values: FMS2, TIM")
endif()
if(NOT TURBO_MEMORY_MODE MATCHES "^(dynamic_symmetric|dynamic_nonsymmetric)$")
    message(FATAL_ERROR
        "Unknown TURBO_MEMORY_MODE='${TURBO_MEMORY_MODE}'. "
        "Valid values: dynamic_symmetric, dynamic_nonsymmetric")
endif()

