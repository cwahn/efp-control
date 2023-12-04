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

using namespace Eigen;

namespace efp
{
    template <typename Derive>
    using Mat = MatrixBase<Derive>;

    template <typename A, typename B>
    bool is_mat_approx(const Mat<A> &am, const Mat<B> &bm)
    {
        return am.isApprox(bm);
    }
    
    template <typename A, typename B>
    bool is_mat_approx(const Mat<A> &am, const Mat<B> &bm, const double &precision)
    {
        return am.isApprox(bm, precision);
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
    using PolyFromRoots = Matrix<
        AssertComplex_t<Scalar_t<A>>,
        A::RowsAtCompileTime == Dynamic
            ? Dynamic
            : A::RowsAtCompileTime + 1,
        1>;

    // ! Partial function. Wrong shape will abort the function.
    template <typename A>
    auto poly_from_roots(const Mat<A> &roots)
        -> PolyFromRoots<A>
    {
        if (A::ColsAtCompileTime != 1 && roots.cols() != 1)
        {
            abort();
        }
        else
        {
            // Size dynamic size can't initiated like this
            PolyFromRoots<A> poly;
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

            const auto tmp = roots.eval();

            for_each_with_index(add_root, VectorView<const Scalar_t<A>>(tmp.data(), tmp.size()));

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
        const auto no_feed_through =
            characteristic_poly((am - bm.col(n) * cm.row(m)).matrix().eval());
        const auto feed_through = ((dm(m, n) - 1.) * den.array()).matrix();
        return (no_feed_through + feed_through).eval();
    }

    template <typename A, typename B, typename C, typename D>
    auto tf_from_ss_nm(const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm, int n, int m)
    {
        const auto den = characteristic_poly(am);
        const auto num = tf_num_from_ss_nm(den, am, bm, cm, dm, n, m);

        return std::make_tuple(num, den);
    }

    template <typename A, typename B, typename C, typename D>
    auto tf_from_ss(const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm)
    {
        Matrix<std::tuple<decltype(characteristic_poly(am)),
                          decltype(tf_num_from_ss_nm(characteristic_poly(am), am, bm, cm, dm, 0, 0))>,
               C::RowsAtCompileTime,
               B::ColsAtCompileTime>
            tfm;

        const int n = bm.cols();
        const int m = cm.rows();

        if (C::RowsAtCompileTime == Dynamic || B::ColsAtCompileTime == Dynamic)
        {
            tfm.resize(m, n);
        }

        const auto den = characteristic_poly(am);

        const auto fill_tfm = [&](int i, int j)
        {  const auto num = tf_num_from_ss_nm(den, am, bm, cm, dm, i, j);
        tfm(j,i) = std::make_tuple(num, den) ; };

        cartesian_for_index(fill_tfm, n, m);

        return tfm;
    }

    // template <typename A, typename B>
    // auto ss_from_tf(const Mat<A> &am, const Mat<B> &bm)
    // {
    //     const auto m = num.rows();
    //     const auto k = den.rows();

    //     if (m > k)
    //     {
    //         abort();
    //     }
    // }

    template <typename A, typename B, typename C, typename D>
    auto gbt(
        const Mat<A> &am,
        const Mat<B> &bm,
        const Mat<C> &cm,
        const Mat<D> &dm,
        const double &dt,
        const double &alpha)
    {
        // todo static size check
        const auto ima = MatrixXd::Identity(am.rows(), am.rows());
        const auto ima_alpha_dt_am = ima - alpha * dt * am;

        const auto dam = solve(ima_alpha_dt_am, ima + (1.0 - alpha) * dt * am);
        const auto dbm = solve(ima_alpha_dt_am, dt * bm);
        const auto dcm = solve(ima_alpha_dt_am.transpose(), cm.transpose()).transpose().eval();
        const auto ddm = (dm + alpha * cm * dbm).eval();

        return std::make_tuple(dam, dbm, dcm, ddm);
    }

    template <typename A, typename B, typename C, typename D>
    auto c2d_ss_euler(const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm, const double &dt)
    {
        return gbt(am, bm, cm, dm, dt, 0.);
    }

    template <typename A, typename B, typename C, typename D>
    auto c2d_ss_tustin(const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm, const double &dt)
    {
        return gbt(am, bm, cm, dm, dt, 0.5);
    }

    template <typename A, typename B, typename C, typename D>
    auto c2d_ss_backward_diff(const Mat<A> &am, const Mat<B> &bm, const Mat<C> &cm, const Mat<D> &dm, const double &dt)
    {
        return gbt(am, bm, cm, dm, dt, 1.);
    }

    template <typename A, typename B, typename C, typename D, typename Q, typename R>
    class KalmanFilter
    {
    public:
        using State = Matrix<Scalar_t<A>,
                             A::RowsAtCompileTime, 1>;
        explicit KalmanFilter(
            const Mat<A> &a,
            const Mat<B> &b,
            const Mat<C> &c,
            const Mat<D> &d,
            const Mat<Q> &q,
            const Mat<R> &r)
            : a_(a),
              b_(b),
              c_(c),
              d_(d),
              q_(q),
              r_(r),
              x_(State::Zero(a.rows(), 1)),
              p_(A::Zero(a.rows(), a.rows()))
        {
        }

        template <typename Z, typename U>
        auto operator()(const Mat<Z> &z, const Mat<U> &u)
            -> State
        {
            // Predict
            x_ = a_ * x_ + b_ * u;
            p_ = a_ * p_ * a_.transpose() + q_;

            // Update
            const auto y = z - (c_ * x_ + d_ * u);
            const auto s = c_ * p_ * c_.transpose() + r_;
            const auto k = p_ * c_.transpose() * s.inverse();
            x_ = x_ + k * y;
            p_ = (A::Identity(p_.rows(), p_.cols()) - k * c_) * p_;

            return x_;
        }

        const State state_estimate() const
        {
            return x_;
        }
        const A estimate_covariance() const
        {
            return p_;
        }

    private:
        A a_;     // State transition matrix
        B b_;     // Control input matrix
        C c_;     // Measurement matrix
        D d_;     // Feedforward (direct transmission) matrix
        Q q_;     // Process noise covariance
        R r_;     // Measurement noise covariance
        State x_; // State estimate
        A p_;     // Estimate covariance
    };

}

#endif