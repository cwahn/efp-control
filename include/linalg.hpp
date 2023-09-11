#ifndef LINALG_HPP_
#define LINALG_HPP_

#include "efp.hpp"
#include "Eigen/Dense"

namespace efp
{
    using namespace Eigen;

    template <typename A, typename B>
    B solve(const A &a, const B &b)
    {
        return a.completeOrthogonalDecomposition().solve(b);
    }

    template <typename Scalar, int degree>
    Matrix<Scalar, degree + 1, 1> poly_coefficients(const Matrix<Scalar, degree, 1> &roots)
    {
        Matrix<Scalar, degree + 1, 1> coeffs;
        coeffs.setZero();

        coeffs[0] = 1.0;

        for (int i = 0; i < degree; ++i)
        {
            for (int j = degree; j >= 1; --j)
            {
                coeffs[j] = coeffs[j] - roots[i] * coeffs[j - 1];
            }
        }

        return coeffs;
    }
}

#endif