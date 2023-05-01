#pragma once

#include <Eigen/Dense>


/**
 * @brief Class for B-Spline interpolation
 */
template <typename T, int Dim>
class Spline {
public:
  using ElemType = typename std::conditional_t<Dim == 1, T, Eigen::Matrix<T, Dim, 1>>;

  Spline() = default;
  Spline(std::vector<T> const& knots, std::vector<ElemType> const& ctrls, int degree);

  ElemType evaluate(T x, int nder = 0) const;
  inline ElemType operator()(const T& x, int nder = 0) const { return evaluate(x, nder); }
  inline T argmin() const {
    return knots_(0);
  }
  inline T argmax() const { 
    const int n = knots_.cols();
    return knots_(n - 1);
  }
  inline auto const& knots() const { return knots_; }
  inline auto const& ctrls() const { return ctrls_; }
  inline int degree() const { return degree_; }
  void interpolate(std::vector<T> const& x, std::vector<ElemType> const& y, int degree);

private:
  static_assert(Dim > 0, "dimension must be positive");

  Eigen::Matrix<T, 1, Eigen::Dynamic> knots_;
  Eigen::Matrix<T, Dim, Eigen::Dynamic> ctrls_;
  int degree_ = 1;
};


using Spline1d = Spline<double, 1>;
using Spline2d = Spline<double, 2>;
using Spline3d = Spline<double, 3>;

template <typename T>
inline Spline<T, 1> spline_interp(std::vector<T> const& x, std::vector<T> const& y, int degree)
{
  Spline<T, 1> sp;
  sp.interpolate(x, y, degree);
  return sp;
}

template <typename T, int Dim>
inline Spline<T, Dim> spline_interp(std::vector<T> const& x, std::vector<Eigen::Matrix<T, Dim, 1>> const& y, int degree)
{
  Spline<T, 1> sp;
  sp.interpolate(x, y, degree);
  return sp;
}
