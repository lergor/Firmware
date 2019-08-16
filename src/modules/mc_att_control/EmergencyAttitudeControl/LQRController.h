#pragma once

#include "lqr.h"
#include <matrix/matrix/math.hpp>

class LQRController {
    enum {Xp=0, Xq, Xr, STATE_DIM};
    enum {Utx = 0, Uty, Utz, CONTROL_DIM};
    enum {Y1 = 0, Y2, Y3, Y_DIM};
    enum {Xy1 = STATE_DIM, Xy2, Xy3, EXTENDED_STATE_DIM};

    static const int FULL_STATE = STATE_DIM;

    using LQR_t = LQR<double, FULL_STATE, CONTROL_DIM>;
public:
    LQRController();

    ~LQRController() = default;

    void update_state(matrix::Vector3f const &rates, matrix::Vector3f const& rates_int, double outputs_sum);

    matrix::Vector3f control_attitude(matrix::Vector3f const &rates_sp);

private:
    matrix::Matrix<float, CONTROL_DIM, FULL_STATE> _K;
    matrix::Matrix<float, CONTROL_DIM, 1> _U;
    matrix::Matrix<float, FULL_STATE, 1> _X;
    matrix::Matrix<float, Y_DIM, 1> _Y;

    typename LQR_t::control_feedback_t _Keigen;
    typename LQR_t::state_matrix_t _A, _Q, _P;
    typename LQR_t::control_matrix_t _R;
    typename LQR_t::control_gain_matrix_t _B;

    Eigen::Matrix<double, STATE_DIM, STATE_DIM> _A_x;
    Eigen::Matrix<double, Y_DIM, STATE_DIM> _C_x;
    Eigen::Matrix<double, STATE_DIM, CONTROL_DIM> _B_x;

    matrix::Vector3f _rates;
    matrix::Vector3f _rates_int;

    void updateA(matrix::Vector3f const &rates, double outputs_sum);

    void calcGainK();

    void setQR();

    void setB();

};
