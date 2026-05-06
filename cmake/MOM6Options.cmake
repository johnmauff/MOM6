set(TURBO_INFRA "FMS2" CACHE STRING "Infrastructure backend: FMS2 or TIM")
set_property(CACHE TURBO_INFRA PROPERTY STRINGS FMS2 TIM)

set(TURBO_MEMORY_MODE "dynamic_symmetric" CACHE STRING
    "MOM6 memory layout: dynamic_symmetric or dynamic_nonsymmetric")
set_property(CACHE TURBO_MEMORY_MODE PROPERTY STRINGS
    dynamic_symmetric dynamic_nonsymmetric)

