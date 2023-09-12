#ifndef LINALG_HPP_
#define LINALG_HPP_

#include "complex"

#include "efp.hpp"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

namespace efp
{
    using namespace Eigen;

    template <typename A>
    using Complex = typename std::complex<A>;

    template <typename A, typename B>
    bool floating_matrix_eq(const A &a, const B &b)
    {
        // todo float
        return (a - b).norm() < std::sqrt(std::numeric_limits<double>::epsilon());
    }

    template <typename A, typename B>
    B solve(const A &a, const B &b)
    {
        return a.completeOrthogonalDecomposition().solve(b);
    }

    // todo What if input is complex matrix
    template <typename Scalar, int n>
    auto eigenvalues(const Matrix<Scalar, n, n> &a)
        -> Matrix<Complex<Scalar>, n, 1>
    {
        return EigenSolver<Matrix<Scalar, n, n>>{a}.eigenvalues();
    }

    template <typename Scalar, int n>
    Matrix<Scalar, n + 1, 1> poly_from_roots(const Matrix<Scalar, n, 1> &roots)
    {
        Matrix<Scalar, n + 1, 1> p;
        p.setZero();
        p[0] = (Scalar)1.;

        const auto add_root = [&](int i, const Scalar &root)
        {
            const auto temp = from_function(i + 2,
                                            [&](int i)
                                            { return i == 0 ? 0. : root * p[i - 1]; });

            for_each_with_index([&](int i, const Scalar &x)
                                { p[i] -= x; },
                                temp);
        };

        for_each_with_index(add_root, roots);

        return p;
    }

    template <typename Scalar, int n>
    auto poly_from_matrix(const Matrix<Scalar, n, n> &am)
        -> Matrix<Complex<Scalar>, n + 1, 1>
    {
        return poly_from_roots(eigenvalues(am));
    }

    // template <typename A, typename B, typename C, typename D, int input = 0>
    // auto tf_from_ss(const A &am, const B &bm, const C &cm, const D &dm)
    // {
    //     // A, B, C, D = abcd_normalize(A, B, C, D)

    //     // nout, nin = D.shape
    //     // if input >= nin:
    //     //     raise ValueError("System does not have the input specified.")

    //     // # make SIMO from possibly MIMO system.
    //     // B = B[:, input:input + 1]
    //     // D = D[:, input:input + 1]

    //     // try:
    //     //     den = poly(A)
    //     // except ValueError:
    //     //     den = 1

    //     // if (prod(B.shape, axis=0) == 0) and (prod(C.shape, axis=0) == 0):
    //     //     num = numpy.ravel(D)
    //     //     if (prod(D.shape, axis=0) == 0) and (prod(A.shape, axis=0) == 0):
    //     //         den = []
    //     //     return num, den

    //     // num_states = A.shape[0]
    //     // type_test = A[:, 0] + B[:, 0] + C[0, :] + D + 0.0
    //     // num = numpy.empty((nout, num_states + 1), type_test.dtype)
    //     // for k in range(nout):
    //     //     Ck = atleast_2d(C[k, :])
    //     //     num[k] = poly(A - dot(B, Ck)) + (D[k] - 1) * den

    //     // return num, den

    //     // todo input dimension check

    //     using NumType = typename D::Scalar;
    //     using DenType = NumType;

    //     constexpr int out_n = D::RowsAtCompileTime;
    //     constexpr int in_n = D::ColsAtCompileTime;
    //     static_assert(out_n < in_n, "System does not have the input specified.");

    //     Matrix<typename B::Scalar, B::RowsAtCompileTime, 1> bv = bm.col(input);
    //     Matrix<typename D::Scalar, D::RowsAtCompileTime, 1> dv = dm.col(input);

    //     DenType den;
    //     // try
    //     // {
    //     den = poly_from_roots(am);
    //     // }
    //     // catch (const std::exception &e)
    //     // {
    //     //     den = DenType(1);
    //     // }

    //     if (bv.prod() == 0 && cm.prod() == 0)
    //     {
    //         Matrix<NumType, nout, 1> num = dv.array();
    //         if (dv.prod() == 0 && am.prod() == 0)
    //         {
    //             den = DenType(0);
    //         }
    //         return std::make_tuple(num, den);
    //     }

    //     constexpr int num_states = A::RowsAtCompileTime;
    //     constexpr int num_outputs = nout;
    //     Matrix<typename A::Scalar, num_states + 1, 1> num;

    //     for (int k = 0; k < num_outputs; ++k)
    //     {
    //         Matrix<typename C::Scalar, 1, C::ColsAtCompileTime> Ck_row = cm.row(k);
    //         Matrix<typename A::Scalar, num_states, 1> poly_arg = am - bm * Ck_row;
    //         Matrix<NumType, num_states + 1, 1> num_k = poly_from_roots(poly_arg) + (dv(k) - NumType(1)) * den;
    //         num.block(k * (num_states + 1), 0, num_states + 1, 1) = num_k;
    //     }

    //     return std::make_tuple(num, den);
    // }
}

#endif