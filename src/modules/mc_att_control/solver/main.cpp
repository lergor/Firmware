#include <iostream>
#include <algorithm>
#include "modules/mc_att_control/EmergencyControl/solver.h"

int main() {
    ratio_solver::plot_dependencies(ratio_solver::drone_params(), 2.);
    double optimal_rho = ratio_solver::calculate_optimal_ratio();
    std::cout << "optimal ratio: " << optimal_rho << std::endl;
    return 0;
}