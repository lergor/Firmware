#pragma once

#include "modules/mc_att_control/lqr/lqr.h"

#define m_assert(expr, msg)  if(!(expr)) { fprintf(stderr, "test failed on %s:%d\n%s\n", __FILE__, __LINE__, msg); exit(1);}

namespace tests {
    static const double EPS = 1E-3;

    // test arguments and correct control feedback matrix
    template<typename Type, int state_dim, int control_dim>
    struct TestSystem {
        typedef LQR<Type, state_dim, control_dim> LQR_N_K;

        typename LQR_N_K::state_matrix_t A;
        typename LQR_N_K::control_gain_matrix_t B;
        typename LQR_N_K::state_matrix_t Q;
        typename LQR_N_K::control_matrix_t R;
        typename LQR_N_K::control_feedback_t correct_K;
    };

    template<typename Type, int state_dim, int control_dim>
    void run_test(TestSystem<Type, state_dim, control_dim> const& test) {
        typedef LQR<Type, state_dim, control_dim> LQR_N_K;
        typename LQR_N_K::control_feedback_t K;

        m_assert(LQR_N_K::compute(test.A, test.B, test.Q, test.R, K), "cannot compute K matrix");

        for (int i = 0; i < state_dim; ++i) {
            m_assert(std::abs(test.correct_K(0, i) - K(0, i)) <= EPS, "incorrect K matrix");
        }
    }

    void run_tests();
}
