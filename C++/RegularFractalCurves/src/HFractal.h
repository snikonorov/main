#ifndef HFRACTAL_H_INCLUDED
#define HFRACTAL_H_INCLUDED

#include <complex>

#include "SGL.h"

/// cos(2*pi*x), approximated by an order 6 polynomial
///
inline double pCos(double x)
{
    const double M_PI2 = M_PI*M_PI;

    ///double r = x - (int)x + int(x < 0);
    double r = x - floor(x);

    r -= 0.5;
    r *= r;

    return r*(0.5*(48 - M_PI2) + r*(4*(M_PI2 - 24) - 8*(M_PI2 - 16)*r)) - 1;
}

/// harmonic fractal

namespace HFractal
{
    struct HFractal
    {
        double a0;      /// base amplitude
        double f0;      /// base frequency [Hz]

        double ma;      /// amplitude multiplication factor
        double mf;      /// frequency multiplication factor

        std::function<std::complex<double>(double)> F;        /// base function

        HFractal
        (
            double a0 = 1,
            double f0 = 10,
            double ma = 0.5,
            double mf = 7,
            std::function<std::complex<double>(double)> F = [](double x)
            {
                return std::complex<double>(pCos(x), pCos(x - 0.25));
            }
        )
        : a0(a0), f0(f0), ma(ma), mf(mf), F(F)
        {
        }

        decltype(auto) operator()(unsigned N)
        {
            return [*this, N](double x)
            {
                std::complex<double> res = 0;

                double a = 1;
                double q = 1;

                for(unsigned k = 0; k < N; k++)
                {
                    //double a = std::pow(ma, k);
                    //double phase = k/(2*M_PI);

                    //double q = std::pow(mf, k);

                    auto A = a; //*std::complex<double>(pCos(phase), pCos(phase - 0.25));

                    res += A*F(q*f0*x);

                    a *= ma;
                    q *= mf;
                }

                return a0*res;
            };
        }
    };

    struct FractalGraph: public Drawable
    {

    };
}

#endif // HFRACTAL_H_INCLUDED
