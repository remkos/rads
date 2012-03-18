find_package (NETCDF REQUIRED)

set (RADS_SRCS rads.f90)

add_library (rads ${RADS_SRCS})

