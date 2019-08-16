#px4_add_library(RatioSolver
#        ../EmergencyAttitudeControl/solver.cpp
#)
#
#target_include_directories(RatioSolver
#        PUBLIC
#        ${CMAKE_CURRENT_SOURCE_DIR}
#)
#
#
##cmake_minimum_required(VERSION 3.12)
##project(solver)
##
##set(CMAKE_CXX_STANDARD 17)
##
##find_package(PythonLibs 2.7 REQUIRED)
##set(NUMPY_DIRS /usr/local/lib/python2.7/dist-packages/numpy/core/include)
##
##add_executable(polynomial_test plynomial_tests.cpp polynomial.h)
##add_executable(main main.cpp polynomial.h solver.cpp solver.h matplotlibcpp.h)
##
##target_include_directories(main PRIVATE ${PYTHON_INCLUDE_DIRS} ${NUMPY_DIRS})
##target_link_libraries(main ${PYTHON_LIBRARIES})
