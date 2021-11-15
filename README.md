# BSplines library

The library provides functions to compute value and derivatives of a BSpline of degree `k` specified by knots array `t`, control points array `c`.


## API

The function 
```cxx
  template <typename KnotsArr, typename CtrlsArr>
  inline typename CtrlsArr::value_type
  spline_eval(
      KnotsArr const& t,
      CtrlsArr const& c,
      int k,
      typename KnotsArr::value_type const& x,
      int nder = 0
  )
```
evaluates derivative of order `nder` at point `x`. The function works with raw pointers, 
so the datatypes `KnotsArr` and `CtrlsArr` are assumed to be of concept `std::random_access_iterator`. 

The class 
```cxx
  template <typename ArgType, typename ValType>
  class SplineT
```
provides the same functionality, but in object-oriented manner. 
The passed knots and control points arrays are copied to std::vector when constructing a new instance of `SplineT`.


## Example
1. Let us define parameters of the spline:
    ```cxx
      vector<double> t {-1.48585135, -1.48585135, -1.48585135, -1.48585135, -1.13146015, -0.94764448,
          -0.27285909, 0.00700573, 0.10728107, 0.13703105, 0.7289507, 0.7289507, 0.7289507, 0.7289507 };
      vector<double> c {-0.99639434, -0.98639558, -0.93904056, -0.73356079, -0.40509636, -0.05280605,
          0.08375401, 0.32264514, 0.51920575, 0.66608736 };
      int k = 3;
    ```
2. Compute value of the spline at point `x`
    ```cxx
      double x = 0.27;
      double val = spline_eval(t, c, k, x);
      double der = spline_eval(t, c, k, x, 1);
      double der2 = spline_eval(t, c, k, x, 2);
    ```
More examples can be found in the `tests/test_bspline.cpp` file.
