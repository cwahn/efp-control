#ifndef CONTROL_HPP_
#define CONTROL_HPP_

#include "efp.hpp"
#include "Eigen/Dense"

namespace efp {
using namespace Eigen;

// using Scalar = double;

// template <typename A, typename B, typename C, typename D, typename E = Scalar, typename F = Scalar>
// auto c2d_gbt(const A &am, const B &bm, const C &cm, const D &dm, const E &dt, const F &alpha)
// {
//     constexpr int dim = am.RowsAtCompileTime;

//     const auto i_mat = Matrix<Scalar, dim, dim>::Identity() - alpha * dt * am;
//     const A ad = i_mat.template bdcSvd<ComputeThinU | ComputeThinV>().solve(Matrix<Scalar, dim, dim>::Identity() - (1. - alpha) * dt * am);
//     const B bd = i_mat.template bdcSvd<ComputeThinU | ComputeThinV>().solve(dt * bm);
//     const C cd = i_mat.transpose().template bdcSvd<ComputeThinU | ComputeThinV>().solve(cm.transpose()).transpose();
//     const D dd = dm + alpha * (cm * bd);

//     return std::make_tuple(am, bd, cd, dd);
// }

// template <typename Scalar, int degree>
// Matrix<Scalar, degree + 1, 1> poly_from_roots(const Matrix<Scalar, degree, 1> &roots)
// {
//     Matrix<Scalar, degree + 1, 1> coeffs;
//     coeffs.setZero();

//     coeffs[0] = 1.0;

//     for (int i = 0; i < degree; ++i)
//     {
//         for (int j = degree; j >= 1; --j)
//         {
//             coeffs[j] = coeffs[j] - roots[i] * coeffs[j - 1];
//         }
//     }

//     return coeffs;
// }

// template <typename A, typename B, typename C, typename D, int input = 0>
// auto ss_to_tf(const A &am, const B &bm, const C &cm, const D &dm)
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

//     constexpr int nout = D::RowsAtCompileTime;
//     constexpr int nin = D::ColsAtCompileTime;
//     static_assert(input < nin, "System does not have the input specified.");

//     Matrix<typename B::Scalar, B::RowsAtCompileTime, 1> bv = bm.col(input);
//     Matrix<typename D::Scalar, D::RowsAtCompileTime, 1> dv = dm.col(input);

//     DenType den;
//     try
//     {
//         den = poly_from_roots(am);
//     }
//     catch (const std::exception &e)
//     {
//         den = DenType(1);
//     }

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

}  // namespace efp

#endif
