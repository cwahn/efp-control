#ifndef SFIANE_MATH_HPP_
#define SFIANE_MATH_HPP_

#include "complex"

#include "efp.hpp"
#include "Eigen/Dense"

namespace efp
{
    using namespace Eigen;

    template <typename... Ts>
    struct void_t_impl
    {
        typedef void type;
    };

    template <typename... Ts>
    using Void_t = typename void_t_impl<Ts...>::type;

    // todo existance of Scalar to indicate matrix type. May need some more
    template <typename A, typename = void>
    struct IsMatrixLike : FalseType
    {
    };

    template <typename A>
    struct IsMatrixLike<A, Void_t<typename A::Scalar>> : TrueType
    {
    };

    template <typename MatA>
    using Scalar_t = typename MatA::Scalar;

    template <typename A>
    using Complex = typename std::complex<A>;

    template <typename A>
    struct IsComplex : FalseType
    {
    };

    template <typename A>
    struct IsComplex<Complex<A>> : TrueType
    {
    };

    template <typename A>
    struct AssertComplex
    {
        using Type = Complex<A>;
    };

    template <typename A>
    struct AssertComplex<Complex<A>>
    {
        using Type = Complex<A>;
    };

    template <typename A>
    using AssertComplex_t = typename AssertComplex<A>::Type;
}

#endif