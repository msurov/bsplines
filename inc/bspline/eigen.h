#include <Eigen/Dense>
#include <Eigen/Core>


template <typename T, int Nx>
void copy(std::vector<T> const& from, Eigen::Matrix<T, 1, Nx>& to)
{
  const int sz = from.size();
  to.resize(1, sz);
  T const* pfrom = &from[0];
  T* pto = &to(0);
  for (int i = 0; i < sz; ++ i) {
    *pto = *pfrom;
    ++ pfrom;
    ++ pto;
  }
}

template <typename T, int Dim, int Nx>
void copy(std::vector<Eigen::Matrix<T, Dim, 1>> const& from, Eigen::Matrix<T, Dim, Nx>& to)
{
  using ElemType = Eigen::Matrix<T, Dim, 1>;

  const int sz = from.size();
  to.resize(Dim, sz);
  ElemType const* pfrom = &from[0];
  for (int i = 0; i < sz; ++ i) {
    to.col(i) = *pfrom;
    ++ pfrom;
  }
}
