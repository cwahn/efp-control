#ifndef LINALG_HPP_
#define LINALG_HPP_

#include "efp.hpp"
#include "Eigen/Dense"

namespace efp
{
    using namespace Eigen;

    template <typename A, typename B>
    B solve(const A &a, const B &b)
    {
        return a.completeOrthogonalDecomposition().solve(b);
    }
}

#endif