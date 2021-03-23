#pragma once
#include "bspline.h"


namespace bspline
{
    template <typename ArgType, typename ValType>
    class SplineT
    {
    private:
        const int _k;
        std::vector<ArgType> _t;
        std::vector<ValType> _c;

    public:
        SplineT(std::vector<ArgType> const& knots, std::vector<ValType> const& ctrls, int k) :
            _k(k), _t(knots), _c(ctrls) {}

        inline ValType operator() (ArgType arg, int nder=0) const
        {
            return spline_eval(_t, _c, _k, arg, nder);
        }

        inline ArgType minarg() const
        {
            return *_t.begin();
        }

        inline ArgType maxarg() const
        {
            return *_t.rbegin();
        }
    };
}
