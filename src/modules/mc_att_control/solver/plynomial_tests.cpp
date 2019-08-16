#include "modules/mc_att_control/EmergencyControl/polynomial.h"
#include <iostream>
#include <cstdio>

#define CHECK(x) (assert(x))
#define CHECK_CLOSE(x, y, eps) (assert( (std::abs((x) - (y))) <= (eps)))

const double EPS = 1E-9;

using Type = double;
using polynomial = math::polynomial<Type>;

void testConstructors() {
    polynomial p0;
    CHECK(0 == p0.get_degree());
    CHECK_CLOSE(0, p0[0], EPS);

    std::vector<Type> v{1, 2, 3};
    polynomial p1(v);
    CHECK(2 == p1.get_degree());
    CHECK(v == p1.get_coeffs());
}

void testAddition() {
    polynomial p1({1, 4, 10}, 3);
    polynomial p2({10, 12, 0, 15}, 4);
    polynomial p3({1, 5, 8, 0, 5}, 5);

    polynomial correct1({12, 21, 18, 15, 5}, 4);
    polynomial correct2({11, 16, 10, 15}, 3);
    polynomial correct3({2, 9, 18, 0, 5}, 4);

    CHECK(p1 + p2 + p3 == correct1);
    CHECK(p1 + p2 == correct2);
    CHECK(p1 + p3 == correct3);
}

void testMultiplication() {
    polynomial p1({1, 5}, 1);
    polynomial p2({3, 10}, 1);
    polynomial p3({2, 1, 5}, 2);

    polynomial correct1({3, 25, 50}, 2);
    polynomial correct2({6, 23, 25, 50}, 3);
    polynomial correct3({6, 53, 140, 175, 250}, 4);

    CHECK(p1 * p2 == correct1);
    CHECK(p2 * p3 == correct2);
    CHECK(p1 * p2 * p3 == correct3);
}

void testDivision() {
    polynomial p1({-1, 1}, 1);
    polynomial p2({1, 5, 4, 10}, 3);
    polynomial p3({7, -3, 15}, 2);

    polynomial correct1({-0.5, 0.5}, 1);
    polynomial correct2({0.25, 1.25, 1, 2.5}, 3);
    polynomial correct3({-1.4, 0.6, -3}, 2);

    CHECK(p1 / 2 == correct1);
    CHECK(p2 / 4 == correct2);
    CHECK(p3 / -5 == correct3);
}

void testSyntheticDivision() {
    polynomial p1({-1, 1}, 1);
    polynomial p2({1, 5, 4, 10}, 3);
    polynomial p3({1, 1}, 1);

    CHECK((p1 * p2).synthetic_div(1) == p2);
    CHECK((p1 * p3).synthetic_div(1) == p3);
    CHECK((p2 * p3).synthetic_div(-1) == p2);
}

void testIsZero() {
    std::vector<Type> v1(10, 0), v2(5, 0);
    polynomial p1(v1, 9);
    polynomial p2(v2, 4);
    polynomial p3({1}, 0);
    polynomial p4({0, 0, 0, 25}, 3);

    CHECK(p1.is_zero());
    CHECK(p2.is_zero());
    CHECK(!p3.is_zero());
    CHECK(!p4.is_zero());
}


void testNormalized() {
    polynomial p1({1, -1}, 1);
    polynomial p2({1, 5, 4, 10}, 3);
    polynomial p3({7, -3, 15}, 2);

    polynomial correct1({-1, 1}, 1);
    polynomial correct2({0.1, 0.5, 0.4, 1}, 3);
    polynomial correct3({7. / 15, -0.2, 1}, 2);

    CHECK(p1.normalized() == correct1);
    CHECK(p2.normalized() == correct2);
    CHECK(p3.normalized() == correct3);
}

void testDerivative() {
    polynomial p1({1}, 0);
    polynomial p2({1, 5}, 1);
    polynomial p3({0, 0, 2, 4}, 3);
    polynomial p4({0, 1, 0, 2}, 3);

    polynomial q1({0}, 0);
    polynomial q2({5}, 0);
    polynomial q3({0, 4, 12}, 2);
    polynomial q4({1, 0, 6}, 2);

    CHECK(p1.derivative() == q1);
    CHECK(p2.derivative() == q2);
    CHECK(p3.derivative() == q3);
    CHECK(p4.derivative() == q4);
}

void testEvaluate() {
    polynomial p1({10}, 0);
    polynomial p2({1, 3}, 1);
    polynomial p3({0, 0, 5}, 2);
    polynomial p4({1, 0, 3, 2}, 3);

    CHECK_CLOSE(10, p1.evaluate(5), EPS);
    CHECK_CLOSE(10, p1.evaluate(1000), EPS);

    CHECK_CLOSE(19, p2.evaluate(6), EPS);
    CHECK_CLOSE(301, p2.evaluate(100), EPS);
    CHECK_CLOSE(11.2, p2.evaluate(3.4), EPS);

    CHECK_CLOSE(125, p3.evaluate(5), EPS);
    CHECK_CLOSE(245, p3.evaluate(7), EPS);
    CHECK_CLOSE(500, p3.evaluate(10), EPS);

    CHECK_CLOSE(29, p4.evaluate(2), EPS);
    CHECK_CLOSE(82, p4.evaluate(3), EPS);
}

void testSolveDegree1() {
    polynomial p1({1, 10}, 1);
    polynomial p2({-50, 10}, 1);
    polynomial p3({-13, 100}, 1);
    polynomial p4({1.0986, 139.1434}, 1);

    CHECK_CLOSE(-(1.0 / 10.0), p1.solve()[0], EPS);
    CHECK_CLOSE(5, p2.solve()[0], EPS);
    CHECK_CLOSE(13.0 / 100.0, p3.solve()[0], EPS);
    CHECK_CLOSE(-(1.0986 / 139.1434), p4.solve()[0], EPS);
}

void testSolveDegree2() {
    polynomial p1({2, 1}, 1);
    polynomial p2({4, 2}, 1);
    polynomial p3({5, 1}, 1);
    polynomial p4({10, 5}, 1);
    polynomial p5({4.8, 1}, 1);
    polynomial p6({8.12, 12}, 1);
    polynomial p7({12, 1}, 1);
    polynomial p8({16, 4}, 1);

    CHECK_CLOSE(p1.solve()[0], (p1 * p2).solve()[0], EPS);
    CHECK_CLOSE(p2.solve()[0], (p1 * p2).solve()[0], EPS);
}

void testSolveDegree3() {
    polynomial p1({-1, 1}, 1);
    polynomial p2({-5, 1}, 1);
    polynomial p3({7, 1}, 1);

    std::vector<Type> res = (p1 * p2 * p3).solve();
    CHECK((res == std::vector<Type>{5, -7, 1}));
}

void testSolveDegree4() {
    polynomial p1({-1, 1}, 1);
    polynomial p2({-5, 1}, 1);
    polynomial p3({7, 1}, 1);
    polynomial p4({5, -2}, 1);

    std::vector<Type> res0 = (p1 * p2 * p3 * p4).solve();
    CHECK_CLOSE(-7, res0[0], EPS);
    CHECK_CLOSE(5, res0[1], EPS);
    CHECK_CLOSE(2.5, res0[2], EPS);
    CHECK_CLOSE(1, res0[3], EPS);

    polynomial p5({-12, 7, 0, 0, 1});
    std::vector<Type> res1 = p5.solve();
    CHECK((res1.size() == 2));
    CHECK_CLOSE(res1[0], -std::sqrt(13) / 2 - 0.5, EPS);
    CHECK_CLOSE(res1[1], std::sqrt(13) / 2 - 0.5, EPS);
}

void testNewtonRaphson() {
    polynomial p1({-10, 1}, 1);
    polynomial p2({-4, 0, 1}, 2);
    polynomial p3({-8, 0, 0, 1}, 3);
    polynomial p4({-8, -2, 1}, 2);

    CHECK_CLOSE(10, p1.newton_raphson(), EPS);
    CHECK_CLOSE(2, p2.newton_raphson(), EPS);
    CHECK_CLOSE(2, p3.newton_raphson(), EPS);
    CHECK_CLOSE(4, p4.newton_raphson(), EPS);
}

void testSolveDegree5_6_7() {
    polynomial p1({-12, 7, 0, 0, 0, 1}, 5);
    std::vector<Type> res1 = p1.solve();
    for (auto &root : res1) {
        CHECK_CLOSE(p1.evaluate(root), 0.0, EPS);
    }

    polynomial p2({-12, -7, 0, 0, 0, 1, 11}, 6);
    std::vector<Type> res2 = p2.solve();
    for (auto &root : res2) {
        CHECK_CLOSE(p2.evaluate(root), 0.0, EPS);
    }

    polynomial p3({-17, -7, 0, 0, 0, 1, 0, 17}, 7);
    std::vector<Type> res3 = p3.solve();
    for (auto &root : res3) {
        CHECK_CLOSE(p3.evaluate(root), 0.0, EPS);
    }

    polynomial p4({0.21381376, 0, 0.11165, 0, -0.15248, 0, -0.04362, 0, 0.03263}, 8);
    std::vector<Type> res4 = p4.solve(EPS);
    for (auto &root : res4) {
        CHECK_CLOSE(p4.evaluate(root), 0.0, EPS);
    }
}

int main() {
    testConstructors();
    testAddition();
    testMultiplication();
    testDivision();
    testSyntheticDivision();
    testIsZero();
    testNormalized();
    testDerivative();
    testEvaluate();
    testSolveDegree1();
    testSolveDegree2();
    testSolveDegree3();
    testSolveDegree4();
    testNewtonRaphson();
    testSolveDegree5_6_7();
    return 0;
}