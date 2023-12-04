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
        // CHECK(characteristic_poly(a) == ref);
    }

    // todo Complex test
}

TEST_CASE("tf_from_ss_nm")
{
    const Matrix2d am{{-2., -1.},
                      {1., 0.}};
    const Vector2d bm{1., 0.};
    const Matrix<double, 1, 2> cm{1., 2.};
    const Matrix<double, 1, 1> dm{{1.}};

    const auto tf_00 = tf_from_ss_nm(am, bm, cm, dm, 0, 0);
    const auto num = std::get<0>(tf_00);
    const auto den = std::get<1>(tf_00);

    // CHECK(num == Matrix<double, 3, 1>{1., 3., 3.});
    // CHECK(den == Matrix<double, 3, 1>{1., 2., 1.});
    CHECK(is_mat_approx(num, Matrix<double, 3, 1>{1., 3., 3.}));
    CHECK(is_mat_approx(den, Matrix<double, 3, 1>{1., 2., 1.}));
}

TEST_CASE("tf_from_ss")
{
    SECTION("static")
    {
        const Matrix2d am{{-2., -1.},
                          {1., 0.}};

        const Matrix2d bm{{1., 1.},
                          {0., 0.}};

        const Matrix2d cm{{1., 2.},
                          {1., 2.}};

        const Matrix2d dm{{1., 0.},
                          {0., 1.}};

        const auto tfm = tf_from_ss(am, bm, cm, dm);
        const auto tf_00 = tfm(0, 0);
        const auto num_00 = std::get<0>(tf_00);
        const auto den_00 = std::get<1>(tf_00);

        const auto tf_11 = tfm(1, 1);
        const auto num_11 = std::get<0>(tf_11);
        const auto den_11 = std::get<1>(tf_11);

        // CHECK(num == Matrix<double, 3, 1>{1., 3., 3.});
        // CHECK(den == Matrix<double, 3, 1>{1., 2., 1.});
        CHECK(is_mat_approx(num_00, Matrix<double, 3, 1>{1., 3., 3.}));
        CHECK(is_mat_approx(den_00, Matrix<double, 3, 1>{1., 2., 1.}));
        CHECK(is_mat_approx(num_11, Matrix<double, 3, 1>{1., 3., 3.}));
        CHECK(is_mat_approx(den_11, Matrix<double, 3, 1>{1., 2., 1.}));
    }

    SECTION("dynamic")
    {
        const MatrixXd am{{-2., -1.},
                          {1., 0.}};

        const MatrixXd bm{{1., 1.},
                          {0., 0.}};

        const MatrixXd cm{{1., 2.},
                          {1., 2.}};

        const MatrixXd dm{{1., 0.},
                          {0., 1.}};

        const auto tfm = tf_from_ss(am, bm, cm, dm);
        const auto tf_00 = tfm(0, 0);
        const auto num_00 = std::get<0>(tf_00);
        const auto den_00 = std::get<1>(tf_00);

        const auto tf_11 = tfm(1, 1);
        const auto num_11 = std::get<0>(tf_11);
        const auto den_11 = std::get<1>(tf_11);

        // CHECK(num == Matrix<double, 3, 1>{1., 3., 3.});
        // CHECK(den == Matrix<double, 3, 1>{1., 2., 1.});
        CHECK(is_mat_approx(num_00, Matrix<double, 3, 1>{1., 3., 3.}));
        CHECK(is_mat_approx(den_00, Matrix<double, 3, 1>{1., 2., 1.}));
        CHECK(is_mat_approx(num_11, Matrix<double, 3, 1>{1., 3., 3.}));
        CHECK(is_mat_approx(den_11, Matrix<double, 3, 1>{1., 2., 1.}));
    }
}
TEST_CASE("c2d")
{
    SECTION("static_or_dyn_0.5")
    {
        const MatrixXd am{{0., 1.},
                          {-10., -3.}};

        const MatrixXd bm{{0., 0.},
                          {10., 0.}};

        const MatrixXd cm{{1., 0.},
                          {0., 0.}};

        const MatrixXd dm{{0., 0.},
                          {0., 0.}};

        const auto dt = 0.1;
        const auto alpha = 0.5;
        const auto X = gbt(am, bm, cm, dm, 0.1, 0.5);

        const MatrixXd dam{{0.95744681, 0.08510638},
                           {-0.85106383, 0.70212766}};
        const MatrixXd dbm{{0.04255319, 0.},
                           {0.85106383, 0.}};
        const MatrixXd dcm{{0.9787234, 0.04255319},
                           {0., 0.}};
        const MatrixXd ddm{{0.0212766, 0.},
                           {0., 0.}};

        // std::cout << std::get<0>(X) << std::endl;
        // std::cout << std::get<1>(X) << std::endl;
        // std::cout << std::get<2>(X) << std::endl;
        // std::cout << std::get<3>(X) << std::endl;

        const double eps = 1e-5;
        CHECK(is_mat_approx(std::get<0>(X), dam, eps));
        CHECK(is_mat_approx(std::get<1>(X), dbm, eps));
        CHECK(is_mat_approx(std::get<2>(X), dcm, eps));
        CHECK(is_mat_approx(std::get<3>(X), ddm, eps));
    }
}

TEST_CASE("KalmanFilter Test", "[KalmanFilter]")
{
    // Define the system matrices for a known example
    Matrix2d a;             // State transition matrix
    Matrix<double, 2, 1> b; // Control input matrix
    Matrix<double, 1, 2> c; // Measurement matrix
    Matrix<double, 1, 1> d; // Feedforward matrix
    Matrix2d q;             // Process noise covariance
    Matrix<double, 1, 1> r; // Measurement noise covariance

    // Initialize matrices with values (example values here, use real ones for your system)
    a << 1.0, 1.0,
        0.0, 1.0;
    b << 0.0,
        1.0;
    c << 1.0, 0.0;
    d << 0.0;
    q << 0.1, 0.0,
        0.0, 0.1;
    r << 0.1;

    // Create an instance of the KalmanFilter
    KalmanFilter<Matrix2d, Matrix<double, 2, 1>, Matrix<double, 1, 2>, Matrix<double, 1, 1>, Matrix2d, Matrix<double, 1, 1>> filter(a, b, c, d, q, r);

    // Define the initial state and input
    Matrix<double, 1, 1> z; // Sensor measurement
    Matrix<double, 1, 1> u; // Control input

    // Initialize state and input (example values, adjust as necessary)
    z << 1.0;
    u << 1.0;

    // Expected state after one update (example value, calculate based on your system)
    Matrix<double, 2, 1> expected_state;
    expected_state << 1.0, 1.0;

    // Perform the update with the filter
    Matrix<double, 2, 1> state_estimate = filter(z, u);

    // Check if the state estimate matches the expected state
    // CHECK(is_mat_approx(state_estimate, expected_state));
    CHECK(state_estimate == expected_state);
}

#endif