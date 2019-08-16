#include "tests.h"

namespace tests {

    void run_tests() {
        static const int state_dim = 4;
        static const int control_dim = 1;

        // test 1.1
        TestSystem<double, state_dim, control_dim> test1_1;
        test1_1.A << 0, 1, 0, 0,
            0, 0, -1, 0,
            0, 0, 0, 1,
            0, 0, 9, 0;

        test1_1.Q << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 10, 0,
            0, 0, 0, 10;

        test1_1.B << 0, 0.1, 0, -0.1;
        test1_1.R << 0.1;
        test1_1.correct_K << -3.1623, -11.1724, -235.2402, -80.1039;

        run_test(test1_1);
        std::cout << "test1.1 passed" << std::endl;

        // test 1.2
        TestSystem<double, state_dim, control_dim> test1_2 = test1_1;
        test1_2.R << 0.01;
        test1_2.correct_K << -10.0000, -25.4097, -308.2620, -109.4647;

        run_test(test1_2);
        std::cout << "test1.2 passed" << std::endl;

        // test 2
        TestSystem<double, state_dim, control_dim> test2;
        test2.A <<  -0.301590, -1.768978, -0.419289, -0.468806,
                    1.571917, 0.290738, 1.296124, -0.786383,
                    -0.147690, 0.270982, 0.246432, 1.012728,
                    0.619349, 1.329404, -0.096675, 0.825942;

        test2.Q <<  1.85985, -1.04525, -0.78216, -0.63415,
                    -1.04525, 3.86754, -0.58204, 2.27457,
                    -0.78216, -0.58204, 5.76405, 1.24019,
                    -0.63415, 2.27457, 1.24019, 2.20766;

        test2.B << -1.778368, -0.996262, 0.092342, -0.624538;
        test2.R << 0.085335;
        test2.correct_K << -0.32886, -8.01441, -14.68269, -12.86050;

        run_test(test2);
        std::cout << "test2 passed" << std::endl;

        // test 3
        TestSystem<double, state_dim, control_dim> test3;
        test3.A <<  12.17021, 4.80395, 16.55317, -1.33837,
                    -1.66242, 0.23943, -2.96318, -15.65328,
                    -12.78237, -5.86784, -9.68983, 6.29172,
                    13.20191, 8.43120, -1.46325, -3.08573;

        test3.Q <<  11.1182, 5.7079, 1.4477, 7.3258,
                    5.7079, 5.9309, 1.7055, 3.2598,
                    1.4477, 1.7055, 3.9226, 2.6810,
                    7.3258, 3.2598, 2.6810, 8.2324;

        test3.B << 5.9947, 8.6622, -9.3791, 9.9874;
        test3.R << 0.019543;
        test3.correct_K << 39.8808, 6.0974, 29.7095, 39.3417;

        run_test(test3);
        std::cout << "test3 passed" << std::endl;

        std::cout << "==== success =====" << std::endl;
    }

}
