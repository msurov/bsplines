#pragma once


namespace bspline
{
    /**
     * The template argument Arr is an arbitrary array with methods
     *  int size();
     * and
     *  operator [] (int index)
     **/

    /**
     * @brief for given val, and sorted array arr find i s.t.
     * 
     *      arr[i] <= val < arr[i+1]
     * or 
     *      0 if val < arr[0]
     * or 
     *      size-1 of val >= arr[size-1]
     */
    template <typename Arr, typename Val>
    inline int bisect(Arr const& arr, Val const& val)
    {
        int from = 0;
        int to = arr.size();

        while (true)
        {
            if (to - from <= 1)
                return from;

            int i = (from + to) / 2;
            
            if (arr[i] > val)
                to = i;
            else
                from = i;
        }
    }

    template <typename T>
    inline T clamp(T const& v, T const& low, T const& hi)
    {
        return 
            v < low ? low : 
            v > hi ? hi : 
            v;
    }

    template <typename T>
    inline T max(T const& a, T const& b)
    {
        return a > b ? a : b;
    }

    /**
     * @brief return i such that
     *  t[i] <= x < t[i+1]
     */
    template <typename KnotsArray>
    inline int spline_knot_index(KnotsArray const& t, int k, typename KnotsArray::value_type const& x)
    {
        int i = bisect(t, x);
        int n = t.size();
        i = clamp<int>(i, k, n - k - 2);
        return i;
    }


    /**
     * @brief De Boor's algorithm
     *
     *  t   knots
     *  c   controls
     *  k   degree of spline
     *  i   knots index x; can be found by get_knot_index(t, k, x)
     *  x   x
     *  nder order of derivative
     */
    template <class CtrlsArr, class KnotsArr>
    typename CtrlsArr::value_type spline_deboor(KnotsArr const& t, CtrlsArr const& c, int k,
        int i, typename KnotsArr::value_type const& x, int nder = 0)
    {
        using Coef = typename KnotsArr::value_type;
        using Result = typename CtrlsArr::value_type;

        Result d[k + 1];

        for (int j = 0; j <= k; ++ j)
        {
            d[j] = c[i - k + j];
        }

        for (int der = 1; der <= nder; ++ der)
        {
            for (int j = k; j >= der; --j)
            {
                d[j] = Coef(k - der + 1) * (d[j] - d[j-1]) / (t[i+j-der+1] - t[i-k+j]);
            }
        }

        const Coef _1(1);

        for (int r = 1 + nder; r <= k; ++ r)
        {
            for (int j = k; j >= r; --j)
            {
                Coef alpha = (x - t[i-k+j]) / (t[i+j+1-r] - t[i-k+j]);
                d[j] = (_1 - alpha) * d[j-1] + alpha * d[j];
            }
        }

        return d[k];
    }


    /**
     * @brief c-style function
     * Evaluate spline basis functions B_{i-k..i,k}(x)
     * or derivative of order nder of them
     *
     *  @param pti  pointer to i-th knot s.t. pknot[0] <= x < pknot[1]
     *  @param k    degree of the spline
     *  @param x    argument to evaluate the basis functions at
     *  @param d    pointer to a buffer of size (k+1) to write basis functions values to
     *  @param nder order of derivative to evaluate (0 if functions needed)
     */
    template <typename Arg>
    static void __spline_eval_basis(
        Arg const* pti, int k, Arg const& x,
        Arg* d, int nder=0
        )
    {
        for (int j = 0; j < k; ++ j)
            d[j] = Arg(0);
        d[k] = Arg(1);

        for (int r = 1; r < k + 1 - nder; ++ r)
        {
            {
                int a = r;
                auto beta = (pti[-a + r + 1] - x) / (pti[-a + r + 1] - pti[-a + 1]);
                d[k-a] = beta * d[k-a + 1];
            }

            for (int a = r - 1; a > 0; -- a)
            {
                auto alpha = (x - pti[-a]) / (pti[-a+r] - pti[-a]);
                auto beta = (pti[-a + r + 1] - x) / (pti[-a + r + 1] - pti[-a + 1]);
                d[k-a] = alpha * d[k-a] + beta * d[k-a + 1];
            }

            {
                int a = 0;
                auto alpha = (x - pti[-a]) / (pti[-a+r] - pti[-a]);
                d[k-a] = alpha * d[k-a];
            }
        }

        for (int r = k + 1 - nder; r < k + 1; ++ r)
        {
            {
                int a = r;
                auto beta = r / (pti[-a + r + 1] - pti[-a + 1]);
                d[k-a] = - beta * d[k-a+1];
            }

            for (int a = r - 1; a > 0; -- a)
            {
                auto alpha = r / (pti[-a + r] - pti[-a]);
                auto beta = r / (pti[-a + r + 1] - pti[-a + 1]);
                d[k-a] = alpha * d[k-a] - beta * d[k-a+1];
            }

            {
                int a = 0;
                auto alpha = r / (pti[-a + r] - pti[-a]);
                d[k-a] = alpha * d[k-a];
            }
        }
    }

    /**
     * @brief Evaluate spline basis functions B_{i-k..i,k}(x)
     * or derivative of order nder of them
     *  @param t    knots array
     *  @param k    degree of spline
     *  @param i    index of the knot s.t. t[i] <= x < t[i+1]
     *  @param x    argument to evaluate the basis functions at
     *  @param d    pointer to a buffer of size (k+1) to write basis functions values
     *  @param nder order of derivative to evaluate (0 if no derivatives needed)
     */
    template <typename KnotsArray>
    static void spline_eval_basis(KnotsArray const& t, int k, int i, typename KnotsArray::value_type const& x, typename KnotsArray::value_type* d, int nder=0)
    {
        using Arg = typename KnotsArray::value_type;
        Arg knots[max(2*k, 1)];
        // spline valuedepends on the knots i-k+1...i+k
        // copy them into the temporarily buffer, and pass 
        // the pointer to t_i location
        for (int j = -k + 1; j <= k; ++ j)
            knots[k - 1 + j] = t[i + j];
        __spline_eval_basis(&knots[k - 1], k, x, d, nder);
    }

    /**
     * @brief evaluate value of the spline function using 
     *  De Boor's algorithm
     * @param t spline knots
     * @param c spline coefficients (controls)
     * @param k spline order
     * @param x spline argument
     * @param nder order of derivative
     * @result spline value at x (or its derivative)
     */
    template <typename KnotsArr, typename CtrlsArr>
    inline typename CtrlsArr::value_type spline_eval_v1(
        KnotsArr const& t_, CtrlsArr const& c_, int k,
        typename KnotsArr::value_type const& x, int nder = 0)
    {
        int i = spline_knot_index(t_, k, x);
        return spline_deboor(t_, c_, k, i, x, nder);
    }

    /**
     * @brief evaluate value of the spline function
     *  by evaluating basis function
     * @param t spline knots
     * @param c spline coefficients (controls)
     * @param k spline order
     * @param x spline argument
     * @param nder order of derivative
     * @result spline value at x (or its derivative)
     */
    template <typename KnotsArr, typename CtrlsArr>
    inline typename CtrlsArr::value_type
    spline_eval(
        KnotsArr const& t,
        CtrlsArr const& c,
        int k,
        typename KnotsArr::value_type const& x,
        int nder = 0
        )
    {
        using Result = typename CtrlsArr::value_type;
        using Coef = typename KnotsArr::value_type;
        int i = spline_knot_index(t, k, x);
        Coef d[k + 1];
        spline_eval_basis(t, k, i, x, d, nder);
        Result result = d[0] * c[i - k];
        for (int j = 1; j <= k; ++ j)
            result += d[j] * c[i - k + j];
        return result;
    }
}
