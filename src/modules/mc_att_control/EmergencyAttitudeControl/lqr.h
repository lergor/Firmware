#pragma once

#include <iostream>
#include <cmath>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#define _GLIBCXX_USE_C99_FP_MACROS_DYNAMIC 1
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Cholesky>
#pragma GCC diagnostic pop


template<typename Type, int N, int K>
class RiccatiSolver {
public:
    typedef Eigen::Matrix<Type, N, N> N_N_matrix_t; // A, Q
    typedef Eigen::Matrix<Type, K, K> K_K_matrix_t; // R
    typedef Eigen::Matrix<Type, N, K> N_K_matrix_t; // B
    typedef Eigen::Matrix<Type, 2 * N, 2 * N> hamilton_matrix_t; // H
    typedef Eigen::Matrix<std::complex<Type>, 2 * N, N> complex_eigenvectors_matrix_t;
    typedef Eigen::Matrix<std::complex<Type>, N, N> complex_N_N_matrix_t;

    static bool solve(
            const N_N_matrix_t &A,
            const N_K_matrix_t &B,
            const N_N_matrix_t &Q,
            const K_K_matrix_t &R,
            N_N_matrix_t &P
    ) {
        // Arimoto-Potter method
        // Hamilton matrix
        hamilton_matrix_t H;
        H << A, -B * R.inverse() * B.transpose(), -Q, -A.transpose();

        // eigenvalues and eigenvectors
        Eigen::EigenSolver<hamilton_matrix_t> eigen_solver(H);

        // extract stable eigenvectors
        complex_eigenvectors_matrix_t eigenvectors;
        int j = 0;
        for (int i = 0; i < 2 * N; ++i) {
            if (eigen_solver.eigenvalues()[i].real() < Type(0)) {
                eigenvectors.col(j) = eigen_solver.eigenvectors().block(0, i, 2 * N, 1);
                ++j;
            }
        }

        // calculate P with stable eigen vector matrix
        complex_N_N_matrix_t U_1, U_2;
        U_1 = eigenvectors.block(0, 0, N, N);
        U_2 = eigenvectors.block(N, 0, N, N);
        P = (U_2 * U_1.inverse()).real();
        return true;
    }

};

template<typename Type = double, int state_dim = 2, int control_dim = 1>
class LQR {
public:
    typedef Eigen::Matrix<Type, state_dim, state_dim> state_matrix_t; // A, Q
    typedef Eigen::Matrix<Type, control_dim, control_dim> control_matrix_t; // R
    typedef Eigen::Matrix<Type, state_dim, control_dim> control_gain_matrix_t; // B
    typedef Eigen::Matrix<Type, control_dim, state_dim> control_feedback_t; // K
    typedef Eigen::Matrix<Type, state_dim, state_dim * control_dim> controllability_matrix_t; // C
    typedef Eigen::Matrix<Type, state_dim * state_dim, state_dim> observability_matrix_t; // O

    LQR() = default;

    ~LQR() = default;

    static bool compute(
        const state_matrix_t &A,
        const control_gain_matrix_t &B,
        const state_matrix_t &Q,
        const control_matrix_t &R,
        control_feedback_t &K,
        state_matrix_t &P
    ) {
        if (!check_system(A, B, Q, R)) {
            return false;
        }

        if (RiccatiSolver<Type, state_dim, control_dim>::solve(A, B, Q, R, P)) {
            return calculate_control_feedback_matrix(R, B, P, K);
        }
        return false;
    }

private:

    static bool check_system(
        const state_matrix_t &A,
        const control_gain_matrix_t &B,
        const state_matrix_t &Q,
        const control_matrix_t &R
    ) {
        return check_definiteness(Q, R) && pair_is_controllable(A, B) && pair_is_observable(A, Q);
    }

    static bool calculate_control_feedback_matrix(
        const control_matrix_t &R,
        const control_gain_matrix_t &B,
        const state_matrix_t &P,
        control_feedback_t &K
    ) {
        if (Eigen::FullPivLU<control_matrix_t>(R).isInvertible()) {
            // K = R^-1 * B^T * P;
            // Eigen: if R is not invertible, values are +-inf
            K = R.inverse() * B.transpose() * P;
            return true;
        }
        return false;
    }

    // check (A, B) is controllable
    static bool pair_is_controllable(const state_matrix_t &A, const control_gain_matrix_t &B) {
        controllability_matrix_t C;
        for (int i = 0; i < state_dim; ++i) {
            C.block(0, i * control_dim, state_dim, control_dim) = A.pow(i) * B;
        }
        return (Eigen::FullPivLU<controllability_matrix_t>(C).rank() == state_dim);
    }

    // (A, E) is observable, E: Q = transpose(E) * E
    static bool pair_is_observable(const state_matrix_t &A, const state_matrix_t &Q) {
//        state_matrix_t E = Q.sqrt();
        state_matrix_t E = Eigen::LLT<state_matrix_t>(Q).matrixL().transpose();
        observability_matrix_t O;
        for (int i = 0; i < state_dim; ++i) {
            O.block(i * state_dim, 0, state_dim, state_dim) = E * A.pow(i);
        }
        return (Eigen::FullPivLU<observability_matrix_t>(O).rank() == state_dim);
    }

    // check Q is positive semi-definite and R is positive definite
    static bool check_definiteness(
        const state_matrix_t &Q,
        const control_matrix_t &R
    ) {
        const Eigen::Matrix<std::complex<Type>, state_dim, 1>& Q_eigenvalues = Q.eigenvalues();
        for (int i = 0; i < state_dim; ++i) {
            if(Q_eigenvalues(i, 0).real() < 0) {
                return false;
            }
        }
        const Eigen::Matrix<std::complex<Type>, control_dim, 1>& R_eigenvalues = R.eigenvalues();
        for (int i = 0; i < control_dim; ++i) {
            if(R_eigenvalues(i, 0).real() <= 0) {
                return false;
            }
        }
        return true;
    }

};
