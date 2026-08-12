set(MOM6_INFRA "FMS2" CACHE STRING "Infrastructure backend: FMS2 or TIM")
set_property(CACHE MOM6_INFRA PROPERTY STRINGS FMS2 TIM)

set(MOM6_MEMORY_MODE "dynamic_symmetric" CACHE STRING
    "MOM6 memory layout: dynamic_symmetric or dynamic_nonsymmetric")
set_property(CACHE MOM6_MEMORY_MODE PROPERTY STRINGS
    dynamic_symmetric dynamic_nonsymmetric)

if(NOT MOM6_INFRA MATCHES "^(FMS2|TIM)$")
    message(FATAL_ERROR
        "Unknown MOM6_INFRA='${MOM6_INFRA}'. Valid values: FMS2, TIM")
endif()
if(NOT MOM6_MEMORY_MODE MATCHES "^(dynamic_symmetric|dynamic_nonsymmetric)$")
    message(FATAL_ERROR
        "Unknown MOM6_MEMORY_MODE='${MOM6_MEMORY_MODE}'. "
        "Valid values: dynamic_symmetric, dynamic_nonsymmetric")
endif()

# Resolved path to the selected memory-mode sources. Interpreting the memory
# mode once here lets every target just reference ${MOM6_MEMORY_DIR}.
set(MOM6_MEMORY_DIR "${MOM6_SOURCE_DIR}/config_src/memory/${MOM6_MEMORY_MODE}")

