if (FFTW2_INCLUDES AND FFTW2_LIBRARIES)
  set(FFTW2_FIND_QUIETLY TRUE)
endif (FFTW2_INCLUDES AND FFTW2_LIBRARIES)

find_path(FFTW2_INCLUDES
  NAMES
  fftw.h dfftw.h sfftw.h
  PATHS
  $ENV{FFTW2DIR}
  ${INCLUDE_INSTALL_DIR}
)

find_library(FFTW2_LIB NAMES fftw dfftw sfftw PATHS $ENV{FFTW2DIR} ${LIB_INSTALL_DIR})
find_library(FFTW2_M_LIB NAMES m PATHS $ENV{FFTW2DIR} ${LIB_INSTALL_DIR})
set(FFTW2_LIBRARIES ${FFTW2_LIB} ${FFTW2_M_LIB})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW2 DEFAULT_MSG
                                  FFTW2_INCLUDES FFTW2_LIBRARIES)

mark_as_advanced(FFTW2_INCLUDES FFTW2_LIBRARIES)
