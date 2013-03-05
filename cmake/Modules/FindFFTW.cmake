if (FFTW_INCLUDES AND FFTW_LIBRARIES)
  set(FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES AND FFTW_LIBRARIES)

find_path(FFTW_INCLUDES
  NAMES
  fftw3.h
  PATHS
  $ENV{FFTWDIR}
  ${INCLUDE_INSTALL_DIR}
)

find_library(FFTW_THREAD_LIB NAMES fftw3_threads PATHS $ENV{FFTWDIR} ${LIB_INSTALL_DIR})
find_library(FFTW_LIB NAMES fftw3 PATHS $ENV{FFTWDIR} ${LIB_INSTALL_DIR})
find_library(FFTW_M_LIB NAMES m PATHS $ENV{FFTWDIR} ${LIB_INSTALL_DIR})
set(FFTW_LIBRARIES ${FFTW_LIB} ${FFTW_M_LIB})
set(FFTW_THREAD_LIBRARIES ${FFTW_THREAD_LIB} ${FFTW_LIBRARIES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
   FFTW_LIB FFTW_INCLUDES FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDES FFTW_LIBRARIES)
