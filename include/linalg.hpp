#ifndef LINALG_HPP_
#define LINALG_HPP_

#include "complex"
#include "tuple"

#include "efp.hpp"
#include "./sfinae_math.hpp"

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

// Size error -> abort
// Value error -> Maybe

namespace efp
{
    using namespace Eigen;

    template <typename Derive>
    using Mat = MatrixBase<Derive>;

    template <typename A, typename B>
    bool is_mat_approx(const Mat<A> &am, const Mat<B> &bm)
    {
        return am.isApprox(bm);
    }

    template <typename A, typename B>
    auto dot_product(const Mat<A> &am, const Mat<B> &bm)
    {
        return am.dot(bm);
    }

    // todo Accept function argument
    template <typename F, typename A>
    auto mat_map(const F &f, const Mat<A> &am)
    {
        return am.unaryExpr(f).eval();
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
    using PolyFromRoots_t = Matrix<
        AssertComplex_t<Scalar_t<A>>,
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
            // Size dynamic size can't initiated like this
            PolyFromRoots_t<A> poly;
            if (A::RowsAtCompileTime == Dynamic || A::ColsAtCompileTime == Dynamic)
            {
                poly.resize(roots.rows() + 1, 1);
            }
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
    auto characteristic_poly(const Mat<A> &am)
    {
        constexpr bool to_complex = IsComplex<Scalar_t<A>>::value;
        const auto e_values = eigenvalues(am);
        const auto polys = poly_from_roots(e_values);

        return mat_map([](auto x)
                       { return complex_cast<to_complex>(x); },
                       poly_from_roots(eigenvalues(am)));
    }

    template <typename Den, typename A, typename B, typename C, typename D>
    auto tf_num_from_ss_nm(const Mat<Den> &den, const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm, int n, int m)
    {
        // todo Shape check
        const auto no_feed_through =
            characteristic_poly((am - bm.col(n) * cm.row(m)).matrix().eval());
        const auto feed_through = ((dm(m, n) - 1.) * den.array()).matrix();
        return (no_feed_through + feed_through).eval();
    }

    // todo MIMO system -> Matrix of transfer function
    // todo Could it be actual matrix
    template <typename A, typename B, typename C, typename D>
    auto tf_from_ss_nm(const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm, int n, int m)
    {
        // todo Shape check
        const auto den = characteristic_poly(am);
        const auto num = tf_num_from_ss_nm(den, am, bm, cm, dm, n, m);

        return std::make_tuple(num, den);
    }

    // todo MIMO system -> Matrix of transfer function
    // todo Could it be actual matrix
    template <typename A, typename B, typename C, typename D>
    auto tf_from_ss(const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm)
    {
        // todo Shape check
        const auto den = characteristic_poly(am);

        Matrix<
            std::tuple<
                decltype(characteristic_poly(am)),
                decltype(tf_num_from_ss_nm(den, am, bm, cm, dm, 0, 0))>,
            C::RowsAtCompileTime,
            B::ColsAtCompileTime>
            tfm;

        if (C::RowsAtCompileTime == Dynamic || B::ColsAtCompileTime == Dynamic)
        {
            tfm.resize(cm.rows(), bm.cols());
        }

        const auto fill_tfm = [&](int i, int j)
        {  const auto num = tf_num_from_ss_nm(den, am, bm, cm, dm, i, j);
        tfm(j,i) = std::make_tuple(num, den) ; };

        const int n = bm.cols();
        const int m = cm.rows();

        cartesian_for_index(fill_tfm, n, m);

        return tfm;
    }
}

#endif