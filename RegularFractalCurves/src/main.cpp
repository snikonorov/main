#include <complex>

#include "RFractal.h"
#include "HFractal.h"

using namespace std::complex_literals;

template<class TFunction>
decltype(auto) HarmonicSignal(const TFunction& A, double f0, unsigned N)
{
    return [A, f0, N](double x)
    {
        std::complex<double> res = 0;

        for(unsigned k = 0; k <= N; k++)
        {
            double f = k*f0;

            res += A(k)*std::complex<double>(pCos(f*x), pCos(f*x - 0.25));
        }

        return res;
    };
}

inline double rnd(double x)
{
    return Rand();
}

///----------------------------------------------------------------------

int main()
{
    ///srand((unsigned)time(0));
    srand(-5);

    auto* w = SCreateWindow(1000, 800, L"Regular Fractal Curves");
    w->style.DynamicRedraw = false;

    {
        auto color_table = [](unsigned k)
        {
            static std::vector<pixel> C =
            {
                pixel(0, 0, 190),
                pixel(0, 130, 0),
                pixel(200, 0, 0),
                pixel(50, 50, 50),
                pixel(10, 50, 200),
                pixel(190, 60, 0),
                pixel(50, 0, 200)
            };

            return C[k % C.size()];
        };

        ///-------------------------------------------------------------

        auto* chart = new TChart();

        chart->Scale(100);
        //chart->Scale(100*Point(4, 1));

        ///-------------------------------------------------------------

        /// guideline points
        ///
        static std::vector<Point> data =
        {
            Point(0, 0), Point(0.25, 1), Point(1, 0.9), Point(1.4, 0.2),
            Point(2, -0.1), Point(2.5, -1.5), Point(3.5, -0.3), Point(5, 0)
        };

        //static auto data = zipmap(rnd, range(0, 5));

        ///-------------------------------------------------------------

        /// Weierstrass-like fractal
        ///
        if(1)
        {
            /// the guiding curve
            auto s = spline(data, spline::Akima);

            auto x = range(s.x0(), s.x1(), 1e6);

            /// plot a polyline approximation of `s`
            ///
            if(1)
            {
                chart->push_back
                (
                    new Graph
                    (
                        zipmap(s, range(x.front(), x.back(), 1000)),
                        0.5*color_table(0)
                    )
                );
            }

            ///----------------------------------------------------------------------------

            unsigned N_harmonics = 3;
            unsigned N_layers = 1;

            /// an array of complex amplitude functions {A[0], ...}
            /// where A[n](k) returns a complex amplitide for a harmonic `k`
            /// (independent of fundamental frequency)
            ///
            auto vA = map
            (
                [N_harmonics](unsigned k) -> std::function<std::complex<double>(double)>
                {
                    auto f = range(0, N_harmonics);

                    auto a = map
                    (
                        [N_harmonics](double f)
                        {
                            return (1 + 0.5*Rand(-1, 1))*(1 - exp(-10*f))/N_harmonics;
                        },
                        f
                    );

                    auto p = map
                    (
                        [](double f)
                        {
                            return Rand();
                        },
                        f
                    );

                    return [a, p](unsigned k)
                    {
                        auto phase = p[k];

                        return a[k]*std::complex<double>(pCos(phase), pCos(phase - 0.25));
                        //return a[k]*std::complex<double>(cos(2*M_PI*phase), cos(2*M_PI*(phase - 0.25)));
                    };
                },
                range(N_layers)
            );

            /// generating a harmonic signal using the first layer
            /// of harmonic amplitude functions
            ///
            auto A = vA[0];
            auto H = HarmonicSignal(A, 1.0, N_harmonics);

            /// draw components of `H`
            ///
            if(0)
            {
                auto x = range(0, 1, 1000);

                chart->push_back(new Graph(zipmap([&H](double x){ return H(x).real(); }, x), color_table(0)));
                chart->push_back(new Graph(zipmap([&H](double x){ return H(x).imag(); }, x), color_table(1)));
            }

            /// draw the magnitude of the harmonic amplitude function `A`
            ///
            if(0)
            {
                chart->push_back
                (
                    new Graph
                    (
                        map
                        (
                            [&A](unsigned k){ return Point(k, std::abs(A(k))); },
                            range(0, N_harmonics)
                        ),
                        color_table(2)
                    )
                );
            }

            ///----------------------------------------------------------------------------

            unsigned N = 9;

            auto F = HFractal::HFractal
            (
                max_amplitude(s)/20,
                1/(s.x1() - s.x0()),
                0.4,
                4.0,
                H
            )(N);

            auto f = [F, s](double x){ return std::real(F(x)) + s(x); };

            /// drawing the fractal curve approximation
            chart->push_back(new Graph(zipmap(f, x), color_table(3)));
        }

        ///-------------------------------------------------------------

        /// zig-zag fractal
        ///
        if(0)
        {
            /// draw the guidline
            if(1)
                chart->push_back(new Graph(data, 0.5*color_table(0)));

            RFractal::x_min = data.front().x;
            RFractal::x_max = data.back().x;

            const unsigned N = 9;      /// number of iterations
            double a = 0.5;

            /// draw the fractal
            chart->push_back(new RFractal::FractalGraph(data, a, N, color_table(3)));

            /// draw the fractal border sections
            if(0)
                chart->push_back(new RFractal::FBGraph(data, a));
        }

        *w += chart;
    }

    w->Show();

    return 0;
}
