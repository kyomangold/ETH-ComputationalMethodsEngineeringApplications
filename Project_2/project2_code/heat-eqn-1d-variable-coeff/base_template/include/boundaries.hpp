#include <cmath>
#include <functional>
#include <tuple>

namespace boundaries {
//! @return 0 for all x
double zero(double);

//! @return 1 for all x
double one(double);

//! @return exp(-10 x)
double exponential(double x);

//! Convience boundary condition type. first component is left boundary,
//! second component is right boundary.
using BoundaryCondition = std::pair<std::function<double(double)>,
      std::function<double(double)>>;
}
