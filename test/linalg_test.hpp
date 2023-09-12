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

TEST_CASE("poly_from_roots")
{
    SECTION("real")
    {
        const Vector3d roots{1., 1., 1.};
        const Vector4d ref{1., -3., 3., -1};

        CHECK(floating_matrix_eq(poly_from_roots(roots), ref));
        // CHECK(poly_from_roots(roots) == ref);
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

        CHECK(floating_matrix_eq(poly_from_roots(roots), ref));
        // CHECK(poly_from_roots(roots) == ref);
    }
}

TEST_CASE("poly_from_matrix")
{
    const Matrix2d a{{0., 1. / 3.},
                     {-1. / 2., 0}};

    const Vector3<Complex<double>> ref{
        Complex<double>{1., 0.},
        Complex<double>{0., 0.},
        Complex<double>{0.16666667, 0.}};

    CHECK(floating_matrix_eq(poly_from_matrix(a), ref));
    // CHECK(poly_from_matrix(a) == ref);
}

#endif