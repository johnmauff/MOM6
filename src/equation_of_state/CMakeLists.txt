target_sources(mom6_ocean PRIVATE
    MOM_EOS.F90
    MOM_EOS_Jackett06.F90
    MOM_EOS_Roquet_SpV.F90
    MOM_EOS_Roquet_rho.F90
    MOM_EOS_TEOS10.F90
    MOM_EOS_UNESCO.F90
    MOM_EOS_Wright.F90
    MOM_EOS_Wright_full.F90
    MOM_EOS_Wright_red.F90
    MOM_EOS_base_type.F90
    MOM_EOS_linear.F90
    MOM_TFreeze.F90
    MOM_temperature_convert.F90
)
