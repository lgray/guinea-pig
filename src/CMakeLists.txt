
file(GLOB MODULES besselCPP.cc beamParametersCPP.cc beamCPP.cc backgroundCPP.cc 
pairsCPP.cc minijetsCPP.cc mathematicalTools.cc lumiCPP.cc 
jetParameterCPP.cc guineapigCPP.cc grvCPP.cc gridCPP.cc
fourierCPP.cc fieldCPP.cc bhabhaSamplesCPP.cc timerCPP.cc 
switchesCPP.cc splineCPP.cc rndmCPP.cc resultsCPP.cc 
principal.cc physicalTools.cc particlesCPP.cc particleBeamCPP.cc 
readData.cc memory.c tridentCPP.cc option_args.cc)

set(WARN_FLAGS "-Wall -Wundef -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math -finline -fPIC ${WARN_FLAGS}")
# Set up header macro generations
include(guinea_config)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)


add_executable(guinea ${MODULES})

if(FFTW2)
   if(FFTW2_FOUND)
      target_link_libraries(guinea ${FFTW2_LIBRARIES})
      include_directories(${FFTW2_INCLUDES})
   else()
      target_link_libraries(guinea fftw m)
   endif()
elseif(FFTW3)
   if(FFTW_FOUND)
      target_link_libraries(guinea ${FFTW_LIBRARIES})
      include_directories(${FFTW_INCLUDES})
   else()
      target_link_libraries(guinea fftw3 m)
   endif()
endif()

# Installing targets and headers:
install(TARGETS guinea
  ARCHIVE
  DESTINATION lib
  COMPONENT Libraries
  RUNTIME
  DESTINATION bin
  COMPONENT Runtime)

file(GLOB headers *.h)
install(FILES ${headers}
  DESTINATION include/guinea-pig
  COMPONENT Headers)
