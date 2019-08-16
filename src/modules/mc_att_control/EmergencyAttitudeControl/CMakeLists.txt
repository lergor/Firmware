find_package(Eigen3 REQUIRED)

px4_add_library(
        EmergencyAttitudeControl
        EmergencyAttitudeControl.cpp
        solver.cpp
        LQRController.cpp
)

include_directories(
        ${EIGEN3_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIRS}/eigen3
)

target_include_directories(
        EmergencyAttitudeControl
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(EmergencyAttitudeControl ${EIGEN3_LIBRARIES})