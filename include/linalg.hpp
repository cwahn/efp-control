#ifndef LINALG_HPP_
#define LINALG_HPP_

#include "complex"
#include "tuple"

#include "efp.hpp"
#include "./sfinae_math.hpp"

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

namespace efp
{
    using namespace Eigen;

    template <typename DerivedA, typename DerivedB>
    bool floating_matrix_eq(const MatrixBase<DerivedA> &a, const MatrixBase<DerivedB> &b)
    {
        // todo float
        return (a - b).norm() < std::sqrt(std::numeric_limits<double>::epsilon());
    }

    // template <typename DerivedA, typename DerivedB>
    // auto solve(const MatrixBase<DerivedA> &am, const MatrixBase<DerivedB> &bm)
    // // -> decltype(am.colPivHouseholderQr().solve(bm))
    // {
    //     return am.colPivHouseholderQr().solve(bm);
    // }

    template <typename DerivedA, typename DerivedB>
    Matrix<typename DerivedA::Scalar, DerivedA::RowsAtCompileTime, DerivedB::ColsAtCompileTime>
    solve(const MatrixBase<DerivedA> &am, const MatrixBase<DerivedB> &bm)
    {
        return am.colPivHouseholderQr().solve(bm);
    }

    // // ? What if input is complex matrix
    // template <typename Derived>
    // auto eigenvalues(const MatrixBase<Derived> &am)
    // // -> decltype(EigenSolver<Matrix<typename Derived::Scalar,
    // //                                Derived::RowsAtCompileTime,
    // //                                Derived::ColsAtCompileTime>>{am}
    // //                 .eigenvalues())
    // {
    //     return EigenSolver<Matrix<typename Derived::Scalar,
    //                               Derived::RowsAtCompileTime,
    //                               Derived::ColsAtCompileTime>>{am}
    //         .eigenvalues();
    // }

    // template <typename Derived>
    // using PolyFromRoots_t = EnableIf_t<Derived::ColsAtCompileTime == 1,
    //                                    Matrix<AssertComplex_t<Scalar_t<Derived>>,
    //                                           Derived::RowsAtCompileTime == Dynamic ? Dynamic : Derived::RowsAtCompileTime + 1,
    //                                           1>>;

    // template <typename Derived>
    // auto poly_from_roots(const MatrixBase<Derived> &roots)
    //     -> PolyFromRoots_t<Derived>
    // {
    //     PolyFromRoots_t<Derived> poly;
    //     poly.setZero();
    //     poly[0] = (AssertComplex_t<Scalar_t<Derived>>)1.;

    //     const auto add_root = [&](int i, auto root)
    //     { const auto diff = from_function(i + 2,
    //                                         [&](int i)
    //                                         { return i == 0 ? 0. : root * poly[i - 1]; });

    //       for_each_with_index([&](int i, auto x)
    //                             { poly[i] -= x; },
    //                             diff); };

    //     for_each_with_index(add_root, roots);

    //     return poly;
    // }

    // template <typename Derived>
    // auto poly_from_matrix(const MatrixBase<Derived> &am)
    // {
    //     return poly_from_roots(eigenvalues(am));
    // }

    // todo MIMO system -> Matrix of transfer function
    // todo Could it be actual matrix
    // template <typename MatA, typename MatB, typename MatC, typename MatD>
    // auto tf_from_ss(const MatA &am, const MatB &bm, const MatC &cm, const MatD &dm, int n, int m)
    // {
    //     // using Scalar = typename MatD::Scalar;

    //     // constexpr int out_n = MatD::RowsAtCompileTime;
    //     // constexpr int in_n = MatD::ColsAtCompileTime;

    //     const auto bm_n = bm.col(n);
    //     const auto dm_n = dm.col(n);

    //     const auto den = poly_from_matrix(am);

    //     const auto cm_m = cm.row(m);
    //     const auto dm_mn = dm_n[m];
    //     // const auto num = poly_from_matrix(am - bm_n.dot(cm_m)) + (dm_mn - 1) * den;
    //     const auto num = poly_from_matrix((am.array() - bm_n.dot(cm_m)).matrix()) + (dm_mn - 1.) * den.array();

    //     return std::make_tuple(num, den);
    // }

    // template <typename A, typename B, typename C, typename D, int input = 0>
    // auto tf_from_ss(const A &am, const B &bm, const C &cm, const D &dm)
    // {
    // A, B, C, D = abcd_normalize(A, B, C, D)

    // nout, nin = D.shape
    // if input >= nin:
    //     raise ValueError("System does not have the input specified.")

    // # make SIMO from possibly MIMO system.
    // B = B[:, input:input + 1]
    // D = D[:, input:input + 1]

    // try:
    //     den = poly(A)
    // except ValueError:
    //     den = 1

    // if (prod(B.shape, axis=0) == 0) and (prod(C.shape, axis=0) == 0):
    //     num = numpy.ravel(D)
    //     if (prod(D.shape, axis=0) == 0) and (prod(A.shape, axis=0) == 0):
    //         den = []
    //     return num, den

    // num_states = A.shape[0]
    // type_test = A[:, 0] + B[:, 0] + C[0, :] + D + 0.0
    // num = numpy.empty((nout, num_states + 1), type_test.dtype)
    // for k in range(nout):
    //     Ck = atleast_2d(C[k, :])
    //     num[k] = poly(A - dot(B, Ck)) + (D[k] - 1) * den

    // return num, den

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