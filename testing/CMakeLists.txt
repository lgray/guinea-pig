
file(COPY acc.dat electron.ini.gz result.ref result.cfg DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(GUNZIP NAMES gunzip)
execute_process(COMMAND ${GUNZIP} -f electron.ini.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Suppress policy CMP0026 warning (as of cmake 3.0.0)
if (${CMAKE_MAJOR_VERSION} GREATER 2)
    cmake_policy(SET CMP0026 OLD)
endif()

get_target_property(binaryname guinea LOCATION)
find_program(NDIFF NAMES numdiff ndiff)

add_test(CLIC-default ${binaryname} --pos_file=electron.ini default default result.out)

if(NDIFF)
   add_test(CLIC-default_ndiff ${NDIFF} -b -l -t CLIC-default result)
   set_tests_properties(CLIC-default_ndiff PROPERTIES
      FAIL_REGULAR_EXPRESSION "FAIL"
      DEPENDS CLIC-default)
endif()
