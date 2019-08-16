#pragma once

#include "lqr.h"
#include "solver.h"
#include <matrix/matrix/math.hpp>
#include <mathlib/math/Limits.hpp>
#include <mathlib/math/Functions.hpp>

#define MYPI 3.14159265358979323846f

// gazebo::solo

#define I_XX_B 0.011
#define I_ZZ_B 0.021

#define I_XX_P 9.75E-07
//#define I_ZZ_P  9.75E-07
#define I_ZZ_P  0.000274004

#define I_XX_T (I_XX_B + 4 * I_XX_P)
#define I_ZZ_T (I_ZZ_B + 4 * I_ZZ_P)

#define I_R  0.000274004
#define L_l 0.206976001

#define MOTOR_CONST 8.54858e-06
#define MOMENT_CONST 0.06
#define ROLL_MOMENT_CONST 1e-06
#define DRAG_CONST 0.000806428

#define DRONE_MASS (0.5 + 0.015 + 4*0.005)
#define W_MAX 1500
#define F_MAX ((W_MAX * W_MAX) * MOTOR_CONST)

#define ROT_CONST 0.015


class EmergencyAttitudeControl {
    enum {Xp=0, Xq, Xnx, Xny, Xr, STATE_DIM};
    enum {Utx = 0, Utz, CONTROL_DIM};
    enum {Y1 = 0, Y2, Y_DIM};
    enum {Xu1 = STATE_DIM, EXTENDED_STATE_DIM};

    static const int FULL_STATE = STATE_DIM;

    using LQR_t = LQR<double, FULL_STATE, CONTROL_DIM>;
public:
    EmergencyAttitudeControl();

    ~EmergencyAttitudeControl() = default;

    void update_state(matrix::Vector3f const &rates, matrix::Vector3f const& n_axis, double outputs_sum=0.);

    matrix::Vector3f control_attitude(matrix::Vector3f const &rates_sp, matrix::Vector3f const& n_axis_sp, float u1);

    matrix::Vector3f get_rates_sp() const;

    void set_rate_limit(const matrix::Vector3f &rate_limit);

private:
    ratio_solver::drone_params _drone;
    ratio_solver::params _stable_params;

    matrix::Matrix<float, CONTROL_DIM, FULL_STATE> _K;
    matrix::Matrix<float, CONTROL_DIM, 1> _U;
    matrix::Matrix<float, FULL_STATE, 1> _X;
    matrix::Matrix<float, Y_DIM, 1> _Y;
    const double a0 = 2 * I_ZZ_T * std::sqrt(DRONE_MASS * 9.81 / (2 * MOTOR_CONST)) / I_XX_B;

    typename LQR_t::control_feedback_t _Keigen;
    typename LQR_t::state_matrix_t _A, _Q, _P;
    typename LQR_t::control_matrix_t _R;
    typename LQR_t::control_gain_matrix_t _B;

    Eigen::Matrix<double, STATE_DIM, STATE_DIM> _A_x;
    Eigen::Matrix<double, Y_DIM, STATE_DIM> _C_x;
    Eigen::Matrix<double, STATE_DIM, CONTROL_DIM> _B_x;

    matrix::Vector3f _rates;
    matrix::Vector3f _n_axis;
    matrix::Vector3f _rate_limit;
    matrix::Vector3f _rate_sp;
    double _r_des;

    void updateA(matrix::Vector3f const &rates, matrix::Vector3f const &n_axis, double outputs_sum);

    void calcGainK();

    void setQR();

    void setB();

    void initialize_stable_params();

};
