
add_executable(simple_sim simple_sim.f90)
target_link_libraries(simple_sim pcmc ${FFT_LIBS})
set_property(TARGET simple_sim APPEND_STRING PROPERTY Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
if(OpenMP_FOUND)
target_link_libraries(simple_sim "${OpenMP_Fortran_LIBRARIES}")
endif()

add_executable(tvf30 tvf30.f90)
target_link_libraries(tvf30 pcmc ${FFT_LIBS})
set_property(TARGET tvf30 APPEND_STRING PROPERTY Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
if(OpenMP_FOUND)
target_link_libraries(tvf30 "${OpenMP_Fortran_LIBRARIES}")
endif()


add_executable(gibbst27 gibbst27.f90)
target_link_libraries(gibbst27 pcmc ${FFT_LIBS})
set_property(TARGET gibbst27 APPEND_STRING PROPERTY Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
if(OpenMP_FOUND)
target_link_libraries(gibbst27 "${OpenMP_Fortran_LIBRARIES}")
endif()

add_executable(tparalist tparalist.f90)
target_link_libraries(tparalist pcmc ${FFT_LIBS})
set_property(TARGET tparalist APPEND_STRING PROPERTY Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
if(OpenMP_FOUND)
target_link_libraries(tparalist "${OpenMP_Fortran_LIBRARIES}")
endif()

add_executable(usernrg usernrg.f90)
target_link_libraries(usernrg pcmc ${FFT_LIBS})
set_property(TARGET usernrg APPEND_STRING PROPERTY Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules")
if(OpenMP_FOUND)
target_link_libraries(usernrg "${OpenMP_Fortran_LIBRARIES}")
endif()


add_test(NAME simple WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" COMMAND simple_sim)

add_test(NAME vf30 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" COMMAND tvf30 "${PROJECT_BINARY_DIR}/${POLCOLMC_EXE}")

add_test(NAME gibbst27 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" COMMAND gibbst27)

add_test(NAME para_list WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" COMMAND tparalist)

add_test(NAME user_defined_potential WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tests" COMMAND usernrg)
