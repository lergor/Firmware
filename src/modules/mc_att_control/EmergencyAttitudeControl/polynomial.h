#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <random>

namespace math {

    const double EPS = 1E-9;

    template<typename T>
    bool is_equal(const T &l, const T &r, const T &eps = std::numeric_limits<T>::epsilon()) {
        return std::abs(l - r) <= eps;
    }

    template<typename T>
    bool is_less(const T &l, const T &r, const T &eps = std::numeric_limits<T>::epsilon()) {
        return r - l > eps;
    }


    template<typename Type>
    class polynomial {
    public:
        polynomial();

        explicit polynomial(std::vector<Type> coeffs, int degree = -1);

        polynomial(const polynomial &) = default;

        polynomial(polynomial &&) noexcept = default;

        polynomial &operator=(const polynomial &rhs);

        polynomial &operator=(polynomial &&) noexcept;

        int get_degree() const;

        const std::vector<Type> &get_coeffs() const;

        Type operator[](size_t i) const;

        Type &operator[](size_t i);

        bool is_zero() const;

        polynomial operator+(const polynomial &rhs) const;

        polynomial operator*(const polynomial &rhs) const;

        polynomial operator/(const Type &rhs) const;

        bool operator==(const polynomial &rhs) const;

        polynomial synthetic_div(Type root) const;

        polynomial derivative() const;

        Type evaluate(Type val) const;

        polynomial &normalized();

        std::vector<Type> solve(Type eps = EPS);

        Type newton_raphson(Type guess = 1.0, Type eps = EPS, int max_iters = 300, Type step = 1.0,
                            Type minimal = 1.0) const;

        friend std::ostream &operator<<(std::ostream &os, const polynomial &p) {
            for (int i = p.degree_; i >= 0; --i) {
                if (is_equal(p.coeffs_[i], Type(0))) {
                    if (i == 0)
                        os << " = 0";
                    continue;
                }
                Type coeff_mod = std::abs(p.coeffs_[i]);
                char sign = is_less(p.coeffs_[i], Type(0)) ? '-' : '+';
                if (i != p.get_degree()) {
                    os << " " << sign << " ";
                }
                if (i == 0) {
                    os << coeff_mod << " = 0";
                    continue;
                }
                if (!is_equal(coeff_mod, Type(1))) {
                    os << coeff_mod << " * ";
                }
                os << "x";
                if (i != 1) {
                    os << "^" << i;
                }
            }
            return os;
        }

    private:
        int degree_;
        std::vector<double> coeffs_;

        std::vector<Type> solve_linear() const;

        std::vector<Type> solve_quadratic() const;

        std::vector<Type> solve_cubic() const;

        std::vector<Type> solve_quartic() const;

    };

    template<typename Type>
    polynomial<Type>::polynomial() : degree_(0), coeffs_(1, 0) {}

    template<typename Type>
    polynomial<Type>::polynomial(std::vector<Type> coeffs, int degree)
        : degree_(degree == -1 ? static_cast<int>(coeffs.size()) - 1 : degree),
          coeffs_(std::move(coeffs)) {}

    template<typename Type>
    int polynomial<Type>::get_degree() const {
        return degree_;
    }

    template<typename Type>
    polynomial<Type> &polynomial<Type>::operator=(const polynomial &rhs) {
        degree_ = rhs.degree_;
        coeffs_ = rhs.coeffs_;
        return *this;
    }

    template<typename Type>
    polynomial<Type> &polynomial<Type>::operator=(polynomial &&rhs) noexcept {
        degree_ = rhs.degree_;
        coeffs_ = std::move(rhs.coeffs_);
        return *this;
    }

    template<typename Type>
    const std::vector<Type> &polynomial<Type>::get_coeffs() const {
        return coeffs_;
    }

    template<typename Type>
    Type polynomial<Type>::operator[](size_t i) const {
        return coeffs_[i];
    }

    template<typename Type>
    Type &polynomial<Type>::operator[](size_t i) {
        return coeffs_[i];
    }

    template<typename Type>
    bool polynomial<Type>::is_zero() const {
        for (auto &coeff : get_coeffs()) {
            if (!is_equal(coeff, Type(0))) {
                return false;
            }
        }
        return true;
    }

    template<typename Type>
    polynomial<Type> polynomial<Type>::operator+(const polynomial &rhs) const {
        const std::vector<Type> *coeffs_less = &rhs.get_coeffs();
        const std::vector<Type> *coeffs_more = &get_coeffs();
        if (degree_ < rhs.get_degree()) {
            std::swap(coeffs_less, coeffs_more);
        }

        std::vector<Type> temp(coeffs_more->begin(), coeffs_more->end());
        for (size_t i = 0; i < coeffs_less->size(); ++i) {
            temp[i] += coeffs_less->at(i);
        }
        return polynomial(temp);
    }

    template<typename Type>
    polynomial<Type> polynomial<Type>::operator*(const polynomial &rhs) const {
        const int prod_degree = get_degree() + rhs.get_degree();
        std::vector<double> temp(prod_degree + 1, 0);
        for (int i = 0; i <= degree_; ++i) {
            for (int j = 0; j <= rhs.get_degree(); ++j) {
                temp[i + j] += coeffs_[i] * rhs[j];
            }
        }
        return polynomial(temp);
    }

    template<typename Type>
    polynomial<Type> polynomial<Type>::operator/(const Type &rhs) const {
        std::vector<Type> coeffs(coeffs_);
        for (int i = 0; i < coeffs.size(); ++i) {
            coeffs[i] /= rhs;
        }
        return polynomial(coeffs, get_degree());
    }

    template<typename Type>
    bool polynomial<Type>::operator==(const polynomial &rhs) const {
        if (degree_ != rhs.get_degree())
            return false;

        for (int i = 0; i <= degree_; ++i) {
            if (!is_equal(coeffs_[i], rhs[i], EPS)) {
                return false;
            }
        }
        return true;
    }

    template<typename Type>
    polynomial<Type> polynomial<Type>::synthetic_div(Type root) const {
        int div_degree = degree_ - 1;
        std::vector<Type> temp(div_degree + 1, 0);
        temp[div_degree] = coeffs_[degree_];
        for (int i = div_degree - 1; i >= 0; --i) {
            temp[i] = (root * temp[i + 1]) + coeffs_[i + 1];
        }
        return polynomial(temp, div_degree);
    }


    template<typename Type>
    polynomial<Type> polynomial<Type>::derivative() const {
        if (!degree_)
            return polynomial();

        int deriv_degree = degree_ - 1;
        std::vector<Type> temp(deriv_degree + 1, 0);
        for (int i = 0; i <= deriv_degree; ++i) {
            temp[i] = (i + 1) * coeffs_[i + 1];
        }
        return polynomial(temp, deriv_degree);
    }

    template<typename Type>
    Type polynomial<Type>::evaluate(Type val) const {
        Type temp = coeffs_[0];
        for (int i = 1; i <= degree_; ++i) {
            temp += coeffs_[i] * std::pow(val, i);
        }
        return temp;
    }

    template<typename Type>
    std::vector<Type> polynomial<Type>::solve(Type eps) {
        switch (degree_) {
            case 0:
                return coeffs_;
            case 1:
                return solve_linear();
            case 2:
                return solve_quadratic();
            case 3:
                return solve_cubic();
            case 4:
                return solve_quartic();
            default: {
                polynomial poly = (*this);
                std::vector<Type> res;
                while (poly.get_degree() > 3) {
                    Type root = poly.newton_raphson();
                    if (!res.empty() && is_equal(root, res.back()))
                        break;

                    if (!is_equal(poly.evaluate(root), Type(0), eps))
                        break;

                    res.push_back(root);
                    poly = poly.synthetic_div(root);
                }
                if (poly.get_degree() <= 3) {
                    std::vector<Type> roots = poly.solve(eps);
                    for (auto &r : roots) {
                        if (is_equal(poly.evaluate(r), Type(0), eps)) {
                            res.push_back(r);
                        }
                    }
                }
                return res;

            };
        }
    }

    template<typename Type>
    std::vector<Type> polynomial<Type>::solve_linear() const {
        if (is_equal(coeffs_[1], Type(0)))
            return std::vector<Type>();
        return std::vector<Type>(1, -coeffs_[0] / coeffs_[1]);
    }

    template<typename Type>
    std::vector<Type> polynomial<Type>::solve_quadratic() const {
        if (is_equal(coeffs_[2], Type(0)))
            return solve_linear();

        Type a = coeffs_[2], b = coeffs_[1], c = coeffs_[0];
        Type D = b * b - 4 * a * c;
        if (is_less(D, Type(0)))
            return std::vector<Type>();

        if (is_equal(D, Type(0)))
            return std::vector<Type>(1, -b / (2 * a));

        Type sqrt_D = std::sqrt(D);
        return {(-b + sqrt_D) / (2 * a), (-b - sqrt_D) / (2 * a)};
    }

    template<typename Type>
    std::vector<Type> polynomial<Type>::solve_cubic() const {
        if (is_equal(coeffs_[3], Type(0)))
            return solve_quadratic();

        if (is_equal(coeffs_[1], Type(0)) && is_equal(coeffs_[2], Type(0))) {
            return std::vector<Type>(1, std::cbrt(-coeffs_[0] / coeffs_[3]));
        }
        if (is_equal(coeffs_[0], Type(0))) {
            std::vector<Type> temp(coeffs_.begin() + 1, coeffs_.begin() + 4);
            polynomial p(temp, 2);
            std::vector<Type> res = p.solve_quadratic();
            res.push_back(0);
            return res;
        }
        Type root1 = newton_raphson();
        polynomial quad = synthetic_div(root1);
        std::vector<Type> temp = quad.solve_quadratic();
        temp.push_back(root1);
        return temp;
    }

    template<typename Type>
    std::vector<Type> polynomial<Type>::solve_quartic() const {
        Type root1 = newton_raphson();
        polynomial cubic_poly = synthetic_div(root1);
        std::vector<Type> temp = cubic_poly.solve_cubic();
        temp.push_back(root1);
        return temp;
    }

    template<typename Type>
    Type polynomial<Type>::newton_raphson(Type guess, Type eps, int max_iters, Type step, Type minimal) const {
        polynomial deriv = derivative();
        Type next_guess = guess;
        Type value = evaluate(next_guess);
        int iter = 0;
        while (is_less(Type(0), std::abs(value), eps) && iter < max_iters) {
            Type deriv_at_point = deriv.evaluate(next_guess);
            if (is_equal(deriv_at_point, Type(0))) {
                next_guess += minimal;
                deriv_at_point = deriv.evaluate(next_guess);
            }
            next_guess -= step * value / deriv_at_point;
            value = evaluate(next_guess);
            ++iter;
        }
        return next_guess;
    }

    template<typename Type>
    polynomial<Type> &polynomial<Type>::normalized() {
        Type coeff0 = coeffs_.back();
        for (int i = 0; i < coeffs_.size(); ++i) {
            coeffs_[i] /= coeff0;
        }
        return (*this);
    }

} // namespace math
