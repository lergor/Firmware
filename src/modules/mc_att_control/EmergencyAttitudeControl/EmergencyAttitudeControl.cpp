#include "EmergencyAttitudeControl.h"
#include <iostream>

EmergencyAttitudeControl::EmergencyAttitudeControl() {

    initialize_stable_params();
    _r_des = _stable_params.r;
    _rate_sp = matrix::Vector3f(0.f, 0.f, -_r_des);

    _K.setZero();
    setB();
    setQR();
    updateA(matrix::Vector3f(0.f, 0.f, 0.f), matrix::Vector3f(0.f, 0.f, 1.f), 0.);
    std::cout << "B:\n" << _B << std::endl;
    std::cout << "A:\n" << _A << std::endl;
    std::cout << "stable_params:\n" << _stable_params << std::endl;
    std::cout << "rates_sp: {" << _rate_sp(0) << ", " << _rate_sp(1) << ", " << _rate_sp(2) << "}" << std::endl;

}

void EmergencyAttitudeControl::update_state(matrix::Vector3f const& rates, matrix::Vector3f const& n_axis, double outputs_sum) {
    _rates = rates;
    _n_axis = n_axis;
    updateA(_rates, _n_axis, outputs_sum);
}

void EmergencyAttitudeControl::calcGainK() {
    if(EmergencyAttitudeControl::LQR_t::compute(_A, _B, _Q, _R, _Keigen, _P)) {
        for(int i(0); i < EmergencyAttitudeControl::CONTROL_DIM; i++){
            for(int j(0); j < EmergencyAttitudeControl::FULL_STATE; j++){
                _K(i,j) = static_cast<float>(_Keigen(i,j));
            }
        }
    } else {
        std::cout << "not controllable!" << std::endl;
    }
}

matrix::Vector3f EmergencyAttitudeControl::control_attitude(matrix::Vector3f const& rates_sp, matrix::Vector3f const& n_axis_sp, float u1) {
    matrix::Vector3f rates_err = _rate_sp - _rates;
    matrix::Vector3f axis_err = n_axis_sp - _n_axis;

    _X(Xp, 0) = rates_err(0);
    _X(Xq, 0) = rates_err(1);
    _X(Xnx, 0) = axis_err(0);
    _X(Xny, 0) = axis_err(1);
    _X(Xr, 0) = rates_err(2);

    calcGainK();
    _U = _K * _X;

    matrix::Vector3f att_ctrl_lqr;
    att_ctrl_lqr(1) = _U(Utx, 0);
    att_ctrl_lqr(2) = _U(Utz, 0);
    att_ctrl_lqr(0) = att_ctrl_lqr(0) * 10 / (2.f * float(F_MAX));
    att_ctrl_lqr(1) = att_ctrl_lqr(1) * 10 / (2.f * float(F_MAX));
    att_ctrl_lqr(2) = att_ctrl_lqr(2) / (2.f * float(F_MAX));

    return att_ctrl_lqr;
}

void EmergencyAttitudeControl::setQR() {
    _Q.setZero();
    _Q.setIdentity();
    _Q(Xp, Xp) = 20.;
    _Q(Xq, Xq) = 20.;
    _Q(Xnx, Xnx) = 1000.;
    _Q(Xny, Xny) = 1000.;

    _R.setIdentity();
    _R *= 135;
}

void EmergencyAttitudeControl::updateA(matrix::Vector3f const& rates, matrix::Vector3f const& n_axis, double outputs_sum) {
    _A_x.setZero();

//    double a0 = I_ZZ_T / I_XX_B * outputs_sum;
    double a = -a0 + (I_XX_T - I_ZZ_T) * (_r_des - double(rates(2))) / I_XX_B;

    _A_x(Xp, Xq) = a;
    _A_x(Xq, Xp) = -a;
    _A_x(Xnx, Xq) = -_stable_params.nz;
    _A_x(Xnx, Xny) = _r_des;
    _A_x(Xny, Xp) = _stable_params.nz;
    _A_x(Xny, Xnx) = -_r_des;

    _A_x(Xr, Xr) = DRAG_CONST/I_XX_B;

//    Eigen::Matrix<double, CONTROL_DIM, CONTROL_DIM> Ir;
//    Ir.setIdentity();
//    Ir *= -1./ROT_CONST;

    _A.block<STATE_DIM, STATE_DIM>(0, 0) = _A_x;
//    _A.block<STATE_DIM, CONTROL_DIM>(0, STATE_DIM) = _B_x;
//    _A.block<CONTROL_DIM, CONTROL_DIM>(STATE_DIM, STATE_DIM) = Ir;
}

void EmergencyAttitudeControl::setB() {
    _B_x.setZero();
    _B_x(Xq, Utx) = L_l/I_XX_B;
    _B_x(Xr, Utz) = MOMENT_CONST/I_ZZ_B;

    _B.setZero();
    _B.block<STATE_DIM, CONTROL_DIM>(0, 0) = _B_x;

//    Eigen::Matrix<double, CONTROL_DIM, CONTROL_DIM> Ir;
//    Ir.setIdentity();
//    Ir *= -1./ROT_CONST;
//    _B.block<CONTROL_DIM, CONTROL_DIM>(STATE_DIM,0) = Ir;
}

void EmergencyAttitudeControl::initialize_stable_params() {
    _drone.MASS = DRONE_MASS;
    _drone.L = L_l;
    _drone.Ixx_P = I_XX_P;
    _drone.Izz_P = I_ZZ_P;
    _drone.Ixx_T = I_XX_B + 4 * I_XX_P;
    _drone.Izz_T = I_ZZ_B + 4 * I_ZZ_P;
    _drone.Kf = MOTOR_CONST;
    _drone.Kt = MOMENT_CONST;
    _drone.GAMMA = DRAG_CONST;
    _drone.Fmax = F_MAX;

    // calculate for case when two rotors are broken
    _stable_params = ratio_solver::calculate_params_for_ratio(0., _drone);
    _stable_params.r *= -1.;
}

matrix::Vector3f EmergencyAttitudeControl::get_rates_sp() const {
    return _rate_sp;
}

void EmergencyAttitudeControl::set_rate_limit(const matrix::Vector3f &rate_limit) {
    _rate_limit = rate_limit;
    for (int i = 0; i < 3; ++i) {
        _rate_sp(i) = math::constrain(_rate_sp(i), -_rate_limit(i), _rate_limit(i));// * 0.75f;
    }
    std::cout << "rates_sp: {" << _rate_sp(0) << ", " << _rate_sp(1) << ", " << _rate_sp(2) << "}" << std::endl;
    std::cout << "rate_limit: {" << _rate_limit(0) << ", " << _rate_limit(1) << ", " << _rate_limit(2) << "}" << std::endl;
}
