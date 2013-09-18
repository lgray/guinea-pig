

set(WARN_FLAGS "-Wall -Wundef -Wextra")

foreach(LANG CXX C)

  # general flags for both compilers/languages:
  set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} -finline -fPIC ${WARN_FLAGS}")

  # only for intel compilers:
  if(CMAKE_${LANG}_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} -static-intel -mavx")
  endif()

  # only for GNU compilers:
  if(CMAKE_${LANG}_COMPILER_ID STREQUAL "GNU")
   set(CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS} -ffast-math")
  endif()

endforeach()
 
