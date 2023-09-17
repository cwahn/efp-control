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

    template <typename Derive>
    using Mat = MatrixBase<Derive>;

    template <typename A, typename B>
    bool is_approx(const Mat<A> &am, const Mat<B> &bm)
    {
        return am.isApprox(bm);
    }

    template <typename A, typename B>
    auto dot_product(const Mat<A> &am, const Mat<B> &bm)
    {
        return am.dot(bm);
    }

    template <typename A, typename B>
    auto solve(const Mat<A> &am, const Mat<B> &bm)
    {
        return am.colPivHouseholderQr().solve(bm).eval();
    }

    template <typename A>
    auto eigenvalues(const Mat<A> &am)
    {
        return EigenSolver<PlainObject_t<A>>{am}.eigenvalues().eval();
    }

    template <typename A>
    using PolyFromRoots_t = Matrix<AssertComplex_t<Scalar_t<A>>,
                                   A::RowsAtCompileTime == Dynamic
                                       ? Dynamic
                                       : A::RowsAtCompileTime + 1,
                                   1>;

    // ! Partial function. Wrong shape will abort the function.
    template <typename A>
    auto poly_from_roots(const Mat<A> &roots)
        -> PolyFromRoots_t<A>
    {
        if (A::ColsAtCompileTime != 1 && roots.cols() != 1)
        {
            abort();
        }
        else
        {
            PolyFromRoots_t<A> poly;
            poly.setZero();
            poly[0] = (AssertComplex_t<Scalar_t<A>>)1.;

            const auto add_root = [&](int i, auto root)
            { const auto diff = from_function(i + 2,
                                            [&](int i)
                                            { return i == 0 ? 0. : root * poly[i - 1]; });
                                            
            for_each_with_index([&](int i, auto x)
                                { poly[i] -= x; },
                                diff); };

            for_each_with_index(add_root, roots);

            return poly;
        }
    }

    // ! Partial function. Wrong shape will abort the function.
    template <typename A>
    auto poly_from_matrix(const Mat<A> &am)
    {
        return poly_from_roots(eigenvalues(am));
    }

    // template <typename DerivedNum, typename A, typename B, typename C, typename D>
    // auto tf_num_from_ss_nm(
    //     const Mat<DerivedNum> &num,
    //     const Mat<A> &am,
    //     const Mat<B> &bm,
    //     const Mat<C> &cm,
    //     const Mat<D> &dm,
    //     int n,
    //     int m)
    // {
    //     return poly_from_matrix((am.array() - dot_product(bm.col(n), cm.row(m))).matrix()) + (dm(m, n) - 1.) * den.array();
    // }

    // // todo MIMO system -> Matrix of transfer function
    // // todo Could it be actual matrix
    // template <typename A, typename B, typename C, typename D>
    // auto tf_from_ss(
    //     const Mat<A> &am,
    //     const Mat<B> &bm,
    //     const Mat<C> &cm,
    //     const Mat<D> &dm,
    //     int n,
    //     int m)
    // {
    //     // todo Shape check
    //     // if ()
    //     // {
    //     //     abort();
    //     // }

    //     const auto den = poly_from_matrix(am);

    //     const auto bm_n = bm.col(n);
    //     const auto dm_n = dm.col(n);
    //     const auto cm_m = cm.row(m);
    //     const auto dm_mn = dm_n[m];

    //     const auto num = poly_from_matrix((am.array() - bm_n.dot(cm_m)).matrix()) + (dm_mn - 1.) * den.array();

    //     return std::make_tuple(num, den);
    // }

    // template <typename A, typename B, typename C, typename D, int input = 0>
    // auto tf_from_ss(const A &am, const B &bm, const C &cm, const D &dm)
    // {
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