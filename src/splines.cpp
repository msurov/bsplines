#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

#include <bspline/splines.h>
#include <bspline/bisect.h>
#include <bspline/eigen.h>


template <typename T, int Dim>
using SplineParameters = std::tuple<
  Eigen::Matrix<T, 1, Eigen::Dynamic>,
  Eigen::Matrix<T, Dim, Eigen::Dynamic>,
  int
>;

/**
 * @brief get knot index of the given spline
 */
template <typename T, int N>
inline int spline_knot_index(Eigen::Matrix<T, 1, N> const& knots, int k, T const& x)
{
  const int n = knots.size();
  int i = bisect(knots, x);
  i = std::clamp<int>(i, k, n - k - 2);
  return i;
}

/**
 * @brief compute knots array for a given array of argumnents
 */
template <typename T, int Nx>
inline Eigen::Matrix<T, 1, Eigen::Dynamic> get_avg_knots(Eigen::Matrix<T, 1, Nx> const& x, int degree)
{
  const int nx = x.size();
  const int nknots = std::max(nx + degree + 1, 2 * degree + 2);
  Eigen::Matrix<T, 1, Eigen::Dynamic> knots(1, nknots);

  for (int j = 1; j < nx - degree; ++j)
    knots(j + degree) = x.segment(j, degree).mean();

  knots.segment(0, degree + 1).fill(x(0));
  knots.segment(nknots - degree - 1, degree + 1).fill(x(nx - 1));
  return knots;
}

/**
 * @brief c-style function
 * Evaluate spline basis functions $B_{i-k..i,k}(x)$
 * or derivatives of order nder
 *
 *  @param pti  pointer to i-th knot s.t. pknot[0] <= x < pknot[1]
 *  @param k    degree of the spline
 *  @param x    argument to evaluate the basis functions at
 *  @param d    pointer to a buffer of size (k+1) to write basis functions values to
 *  @param nder order of derivative to evaluate (0 if functions needed)
 */
template <typename T>
static void spline_eval_basis(T const* pti, int k, T const& x, T* basis, int nder = 0)
{
  for (int j = 0; j < k; ++j)
    basis[j] = T(0);
  basis[k] = T(1);

  for (int r = 1; r < k + 1 - nder; ++r) {
    {
      int a = r;
      auto beta = (pti[-a + r + 1] - x) / (pti[-a + r + 1] - pti[-a + 1]);
      basis[k - a] = beta * basis[k - a + 1];
    }

    for (int a = r - 1; a > 0; --a) {
      auto alpha = (x - pti[-a]) / (pti[-a + r] - pti[-a]);
      auto beta = (pti[-a + r + 1] - x) / (pti[-a + r + 1] - pti[-a + 1]);
      basis[k - a] = alpha * basis[k - a] + beta * basis[k - a + 1];
    }

    {
      int a = 0;
      auto alpha = (x - pti[-a]) / (pti[-a + r] - pti[-a]);
      basis[k - a] = alpha * basis[k - a];
    }
  }

  for (int r = k + 1 - nder; r < k + 1; ++r) {
    {
      int a = r;
      auto beta = r / (pti[-a + r + 1] - pti[-a + 1]);
      basis[k - a] = -beta * basis[k - a + 1];
    }

    for (int a = r - 1; a > 0; --a) {
      auto alpha = r / (pti[-a + r] - pti[-a]);
      auto beta = r / (pti[-a + r + 1] - pti[-a + 1]);
      basis[k - a] = alpha * basis[k - a] - beta * basis[k - a + 1];
    }

    {
      int a = 0;
      auto alpha = r / (pti[-a + r] - pti[-a]);
      basis[k - a] = alpha * basis[k - a];
    }
  }
}

/**
 * @brief evaluate spline value
 * @param knots is the array of spline knots
 * @param ctrls is the array of spline control points
 * @param k is the degree of the spline
 * @param arg is the argument to evaluate spline at
 * @param nder derivative order to evaluate
 * @return value of the spline of vector or scalar type
 */
template <typename T, int Dim, int NKnots, int NCtrls>
inline auto spline_eval(Eigen::Matrix<T, 1, NKnots> const& knots, Eigen::Matrix<T, Dim, NCtrls> const& ctrls, int k,
                        T const& arg, int nder = 0)
{
  const int i = spline_knot_index(knots, k, arg);

  T basis[k + 1];
  spline_eval_basis(&knots(i), k, arg, basis, nder);

  Eigen::Matrix<T, Dim, 1> result = basis[0] * ctrls.col(i - k);
  for (int j = 1; j <= k; ++j)
    result += basis[j] * ctrls.col(i - k + j);

  if constexpr (Dim == 1)
    return result(0);
  else
    return result;
}

template <typename T, int Nx>
static bool is_increasing(Eigen::Matrix<T, 1, Nx> const& x)
{
  if (Nx <= 1)
    return true;

  T const* p = &x(0);
  T prev = p[0];

  for (int i = 1; i < Nx; ++ i) {
    if (p[i] <= prev)
      return false;
    prev = p[i];
  }
  return true;
}

/**
 * @brief Generate spline based on sparse matrices
 *
 * @note This function generates a spline that should be faster than standard Eigen's spline
 *
 * @param x_scaled Scaled function's arguments
 * @param y Function's value
 */
template <typename T, int Dim, int Nx, int Ny>
static SplineParameters<T, Dim> spline_interp(Eigen::Matrix<T, 1, Nx> const& x, Eigen::Matrix<T, Dim, Ny> const& y, int deg)
{
  using SparseMatrixType = Eigen::SparseMatrix<T, Eigen::RowMajor>;

  assert(deg > 0);
  assert(is_increasing(x));
  const auto npts = static_cast<int>(x.size());
  assert(npts > deg);
 
  auto knots = get_avg_knots(x, deg);
  SparseMatrixType A(npts, npts);
  A.reserve((deg + 1) * npts);

  int j = deg;
  T basis[deg + 1];

  for (int i = 0; i < npts; ++i) {
    const T v = x(i);
    for (; j < npts - 1 && knots(j + 1) <= v; ++j) {}
    spline_eval_basis(&knots(j), deg, v, basis);
    for (int l = 0; l <= deg; ++l) {
      A.insert(i, j - deg + l) = basis[l];
    }
  }

  A.makeCompressed();

  Eigen::SparseLU<SparseMatrixType, Eigen::NaturalOrdering<int>> solver(A);
  Eigen::Matrix<T, Eigen::Dynamic, Dim> ctrls = solver.solve(y.transpose());
  return {knots, ctrls.transpose(), deg};
}

template <typename T, int Dim>
Spline<T, Dim>::Spline(std::vector<T> const& knots, std::vector<ElemType> const& ctrls, int degree)
{
  assert(knots.size() == ctrls.size() + degree + 1);
  copy(knots, knots_);
  copy(ctrls, ctrls_);
  degree_ = degree;
}

template <typename T, int Dim>
typename Spline<T, Dim>::ElemType Spline<T, Dim>::evaluate(T x, int nder) const
{
  /// @todo: use optional extrapolation?
  // x = std::clamp(x, get_min_arg(), get_max_arg());
  return spline_eval(knots_, ctrls_, degree_, x, nder);
}

template <typename T, int Dim>
void Spline<T, Dim>::interpolate(std::vector<T> const& x, std::vector<ElemType> const& y, int degree)
{
  Eigen::Matrix<T, 1, Eigen::Dynamic> x_;
  Eigen::Matrix<T, Dim, Eigen::Dynamic> y_;
  copy(x, x_); /// @todo: make without copying
  copy(y, y_);
  std::tie(knots_, ctrls_, degree_) = spline_interp(x_, y_, degree);
}

template class Spline<double, 1>;
template class Spline<double, 2>;
template class Spline<double, 3>;
