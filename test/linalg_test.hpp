#ifndef LINALG_TEST_HPP_
#define LINALG_TEST_HPP_

#include "catch2/catch_test_macros.hpp"
#include "test_common.hpp"

#include "linalg.hpp"

using namespace efp;
using namespace Eigen;

TEST_CASE("dot_product")
{
    const Vector3d a{1., 2., 3.};

    CHECK(dot_product(a, a) == 14.);
}

TEST_CASE("mat_map")
{
    const Matrix2d am{{-2., -1.},
                      {1., 0.}};

    const Matrix2d refm{{4., 1.},
                        {1., 0.}};

    const auto square_ = [](auto x)
    { return square<double>(x); };

    CHECK(is_mat_approx(
        mat_map(square_, am),
        refm));
}

TEST_CASE("solve")
{
    const Matrix3d a{{1., 2., 3.},
                     {4., 5., 6.},
                     {7., 8., 10.}};
    const Vector3d b{3., 3., 4.};
    const Vector3d x_ref{-2., 1., 1.};

    // CHECK((solve(a, b) - x_ref).norm() < double_tol);
    CHECK(is_mat_approx(solve(a, b), x_ref));
}

TEST_CASE("eigenvalues")
{
    const Matrix3d a{{2, 1, 0},
                     {1, 3, 1},
                     {0, 1, 2}};
    const Vector3<std::complex<double>> ref{std::complex<double>{1., 0.},
                                            std::complex<double>{2., 0.},
                                            std::complex<double>{4., 0.}};

    // CHECK(eigenvalues(a) == ref);
    // CHECK((eigenvalues(a) - ref).norm() < double_tol);
    CHECK(is_mat_approx(eigenvalues(a), ref));
}

TEST_CASE("poly_from_roots")
{
    SECTION("real")
    {
        const Vector3d roots{1., 1., 1.};
        const Vector4d ref{1., -3., 3., -1};

        CHECK(is_mat_approx(poly_from_roots(roots), ref));
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

        CHECK(is_mat_approx(poly_from_roots(roots), ref));
    }
}

TEST_CASE("characteristic_poly")
{
    SECTION("real")
    {
        const Matrix2d a{{0., 1. / 3.},
                         {-1. / 2., 0}};

        const Vector3d ref{
            1.,
            0.,
            0.1666666666666666667};

        CHECK(is_mat_approx(characteristic_poly(a), ref));
        CHECK(characteristic_poly(a) == ref);
    }

    // todo Complex test
}

// TEST_CASE("tf_from_ss_nm")
// {
//     const Matrix2d am{{-2., -1.},
//                       {1., 0.}};
//     const Vector2d bm{1., 0.};
//     const RowVector2d cm{1., 2.};
//     const Matrix<double, 1, 1> dm{{1.}};

//     const auto tf_00 = tf_from_ss_nm(am, bm, cm, dm, 0, 0);
//     const auto num = std::get<0>(tf_00);
//     const auto den = std::get<1>(tf_00);

//     CHECK(num == Matrix<double, 3, 1>{1., 3., 3.});
//     CHECK(den == Matrix<double, 3, 1>{1., 2., 1.});
// }

#endif