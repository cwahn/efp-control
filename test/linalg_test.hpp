#ifndef LINALG_TEST_HPP_
#define LINALG_TEST_HPP_

#include "catch2/catch_test_macros.hpp"
#include "test_common.hpp"

#include "linalg.hpp"

using namespace efp;
using namespace Eigen;

TEST_CASE("solve")
{
    const Matrix3d a{{1., 2., 3.},
                     {4., 5., 6.},
                     {7., 8., 10.}};
    const Vector3d b{3., 3., 4.};
    const Vector3d x_ref{-2., 1., 1.};

    // CHECK((solve(a, b) - x_ref).norm() < double_tol);
    CHECK(floating_matrix_eq(solve(a, b), x_ref));
}

TEST_CASE("eigenvalues")
{
    const Matrix3d a{{2, 1, 0},
                     {1, 3, 1},
                     {0, 1, 2}};
    const Vector3<std::complex<double>> ref{std::complex<double>{1., 0.},
                                            std::complex<double>{2., 0.},
                                            std::complex<double>{4., 0.}};

    // CHECK((eigenvalues(a) - ref).norm() < double_tol);
    CHECK(floating_matrix_eq(eigenvalues(a), ref));
}

TEST_CASE("poly_coefficients")
{
    SECTION("real")
    {
        const Vector3d roots{1., 1., 1.};
        const Vector4d ref{1., -3., 3., -1};

        CHECK(floating_matrix_eq(poly_coefficients(roots), ref));
        // CHECK(poly_coefficients(roots) == ref);
    }

    SECTION("complex")
    {
        const Vector2<Complex<double>> roots{
            Complex<double>{1., 1.},
            Complex<double>{1., -1.}};

        const Vector3<Complex<double>> ref{
            Complex<double>{1., 0.},
            Complex<double>{-2., 0.},
            Complex<double>{2., 0.},
        };

        CHECK(floating_matrix_eq(poly_coefficients(roots), ref));
        // CHECK(poly_coefficients(roots) == ref);
    }
}

// TEST_CASE("poly_coefficients")
// {
//     const Vector3d a{-1. / 2, 0, 1. / 2};
//     const Vector4d ref{1., 0., -0.25, 0.};

//     CHECK((poly_coefficients(a) - ref).norm() < double_tol);
// }

#endif