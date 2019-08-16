#include "LQRController.h"

#define MYPI 3.14159265358979323846f

// gazebo::solo

#define I_XX_B 0.011
#define I_YY_B 0.011
#define I_ZZ_B 0.021

#define I_XX_P 0
#define I_ZZ_P  9.75E-07
#define I_R  0.000274004
#define L_l 0.206976001

#define MOTOR_CONST 8.54858e-06
#define MOMENT_CONST 0.06
#define ROLL_MOMENT_CONST 1e-06
#define DRAG_CONST 0.000806428

#define DRONE_MASS (0.5 + 0.015 + 4*0.005)
#define W_MAX 1500
#define F_MAX ((W_MAX * W_MAX) * MOTOR_CONST)

LQRController::LQRController() {
    _K.setZero();
    setB();
    setQR();
    updateA(matrix::Vector3f(0.f, 0.f, 0.f), 0.);
}


void LQRController::updateA(matrix::Vector3f const& rates, double outputs_sum) {
    _A_x.setZero();
    _A_x(Xp, Xq) = (MOTOR_CONST/I_XX_B) * ((I_YY_B - I_ZZ_B) * double(rates(2)) + I_ZZ_P * outputs_sum);
    _A_x(Xq, Xp) = (MOTOR_CONST/I_YY_B) * ((I_ZZ_B - I_XX_B) * double(rates(2)) + I_ZZ_P * outputs_sum);
    _A_x(Xr, Xq) = (MOTOR_CONST/I_ZZ_B) * (I_XX_B - I_YY_B) * double(rates(0));
    _A_x(Xr, Xr) = -DRAG_CONST/I_ZZ_B;

//    _C_x(Y1, Xp) = 1;
//    _C_x(Y2, Xq) = 1;
//    _C_x(Y3, Xr) = 1;

    _A.block<STATE_DIM, STATE_DIM>(0,0) = _A_x;
//    _A.block<Y_DIM, STATE_DIM>(STATE_DIM,0) = _C_x;
}

void LQRController::update_state(matrix::Vector3f const& rates, matrix::Vector3f const& rates_int, double outputs_sum) {
    _rates = rates;
    _rates_int = rates_int;
    updateA(_rates, outputs_sum);
}

void LQRController::calcGainK() {
    if(LQRController::LQR_t::compute(_A, _B, _Q, _R, _Keigen, _P)) {
        for(int i(0); i < LQRController::CONTROL_DIM; i++){
            for(int j(0); j < LQRController::FULL_STATE; j++){
                _K(i,j) = static_cast<float>(_Keigen(i,j));
            }
        }
    }
}

matrix::Vector3f LQRController::control_attitude(matrix::Vector3f const& rates_sp) {
    matrix::Vector3f rates_err = rates_sp - _rates;

    _Y(Y1, 0) = _rates(0);
    _Y(Y2, 0) = _rates(1);
    _Y(Y3, 0) = _rates(2);

    _X(Xp, 0) = rates_err(0);
    _X(Xq, 0) = rates_err(1);
    _X(Xr, 0) = rates_err(2);

//    _X(Xy1, 0) = _rates_int(0);
//    _X(Xy2, 0) = _rates_int(1);
//    _X(Xy3, 0) = _rates_int(2);

    calcGainK();

    _U = _K * _X;
    matrix::Vector3f att_ctrl_lqr(_U);
    att_ctrl_lqr(0) = att_ctrl_lqr(0) * 10 / (2.f * float(F_MAX));
    att_ctrl_lqr(1) = att_ctrl_lqr(1) * 10 / (2.f * float(F_MAX)); //* 0.95f;
    att_ctrl_lqr(2) = att_ctrl_lqr(2) / (2.f * float(F_MAX));
    return att_ctrl_lqr;
}

void LQRController::setQR() {
    _Q.setIdentity();
    _Q.block<STATE_DIM, STATE_DIM>(0, 0) *= 10;
    _Q.block<Y_DIM, Y_DIM>(STATE_DIM, STATE_DIM) *= 5;

    _R.setIdentity();
    _R *= 1035;
}

void LQRController::setB() {
    _B_x.setZero();
    _B_x(Xp, Utx) = L_l/I_XX_B;
    _B_x(Xq, Uty) = L_l/I_YY_B;
    _B_x(Xr, Utz) = MOMENT_CONST/I_ZZ_B;

    _B.block<STATE_DIM, CONTROL_DIM>(0,0) = _B_x;
}
