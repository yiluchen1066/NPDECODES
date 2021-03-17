#if SOLUTION
# Dependencies of mastersolution:
#else
# Add your custom dependencies here:
#endif

# DIR will be provided by the calling file.

set(SOURCES
  ${DIR}/nlmatode_main.cc
  ${DIR}/nlmatode.h
  ${DIR}/nlmatode.cc
  ${DIR}/ode45.h
  ${DIR}/polyfit.h
)

set(LIBRARIES
  Eigen3::Eigen
)
