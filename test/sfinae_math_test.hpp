#ifndef SFINAE_MATH_TEST_HPP_
#define SFINAE_MATH_TEST_HPP_

#include "catch2/catch_test_macros.hpp"
#include "test_common.hpp"

#include "sfinae_math.hpp"

using namespace efp;
using namespace Eigen;

TEST_CASE("Element_t_Eigen") {
    CHECK(IsSame<Element_t<Matrix3d>, double>::value);
    CHECK(IsSame<Element_t<Vector3d>, double>::value);
}

#endif