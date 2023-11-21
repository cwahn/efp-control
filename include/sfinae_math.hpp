#ifndef SFIANE_MATH_HPP_
#define SFIANE_MATH_HPP_

#include "complex"

#include "efp.hpp"
#include "Eigen/Dense"

namespace efp
{
    using namespace Eigen;

    // todo existance of Scalar to indicate matrix type. May need some more
    template <typename A, typename = void>
    struct IsMatrixLike : False
    {
    };

    template <typename A>
    struct IsMatrixLike<A, Void<typename A::Scalar>> : True
    {
    };

    template <typename MatA>
    using Scalar_t = typename MatA::Scalar;

    template <typename Derived>
    using PlainObject_t = typename Eigen::PlainObjectBase<Derived>::PlainObject;
}

#endif