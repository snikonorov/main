#ifndef RFRACTAL_H_INCLUDED
#define RFRACTAL_H_INCLUDED

#include <deque>
#include <cmath>

#include "SGL.h"

#ifndef M_PI
const double M_PI = 3.141592653589793;
#endif

namespace RFractal
{
    double x_min = 0;
    double x_max = 1;

    /// fractal 'amplitude'
    ///
    /// (x is supposed to be from [0, 1])
    ///
    double A(double x)
    {
        double q = (x - x_min)/(x_max - x_min) + 2;

        return 1.0 + 0.5*sin(200*M_PI*q*q);
    }

    /// fractal 'sign' function
    double S(double x)
    {
        double q = (x - x_min)/(x_max - x_min) + 2;

        return sign(sin(300*M_PI*q*q));
    }

    std::pair<Point, Point> Split(const Point& p1, const Point& p2, double a)
    {
        //auto p0 = (p1 + p2)/2;

        auto l1 = (3*p1 + p2)/4;
        auto l2 = (3*p2 + p1)/4;

        auto c1 = l1 + a*A(l1.x)*S(l1.x)*(p2 - p1).ortho()/4;
        auto c2 = l2 + a*A(l2.x)*S(l2.x)*(p1 - p2).ortho()/4;

        ///auto c1 = l1 + a*(p2 - p1).ortho()/4;
        ///auto c2 = l2 + a*(p1 - p2).ortho()/4;

        return {c1, c2};
    }

    decltype(auto) CreateFractal(const std::vector<Point>& data, double a, unsigned N)
    {
        auto p = data;
        decltype(p) _p;

        for(unsigned k = 0; k < N; k++)
        {
            _p.clear();
            _p.reserve(4*(p.size() - 1) + 1);

            unsigned Q = p.size() - 1;
            for(unsigned q = 0; q < Q; q++)
            {
                auto& p1 = p[q];
                auto& p2 = p[q+1];

                auto p0 = (p1 + p2)/2;

                auto c = Split(p1, p2, a);

                _p.push_back(p1);
                _p.push_back(c.first);
                _p.push_back(p0);
                _p.push_back(c.second);
            }

            _p.push_back(p.back());

            p = _p;
        }

        return p;
    }

    struct Node
    {
        Point p1, p2;
        unsigned level;

        Node(const Point& p1, const Point& p2, unsigned level = 0)
            : p1(p1), p2(p2), level(level)
        {
        }

        /// one 'fractal' iteration
        decltype(auto) split(double a) const
        {
            auto p0 = (p1 + p2)/2;
            auto c = Split(p1, p2, a);

            return std::array<Node, 4>
            {
                Node{p1, c.first, level+1}, Node{c.first, p0, level+1},
                Node{p0, c.second, level+1}, Node{c.second, p2, level+1}
            };
        }
    };

    ///--------------------------------------------------------------------------

    /// [rectangle]
    decltype(auto) GetFractalBorder(const std::vector<Point>& data, double a)
    {
        constexpr double SQRT_2_4 = sqrt(2)/4;

        auto q = a*SQRT_2_4;

        double ma = q < 1 ? q/(1 - q) : 0.0;
        ///double ma = q < 1 ? 0.25*a/(1 - q) : 0.0;    /// smaller

        Rect R;

        if(ma > 0 && data.size())
        {
            R = data.front();

            for(unsigned k = 0; k < data.size()-1; k++)
            {
                const auto& p1 = data[k];
                const auto& p2 = data[k+1];

                auto l = ma*(p2 - p1);
                auto n = l.ortho();

                R.Merge(p1+n);
                R.Merge(p2+n);
                R.Merge(p1-n);
                R.Merge(p2-n);

                R.Merge(p1-l);
                R.Merge(p2+l);
            }
        }

        return R;
    }

    /// [circle]
    /// (assuming a*sqrt(2)/4 < 1)
    /// returns border radius squared
    ///
    decltype(auto) GetSegmentBorder(const Point& p1, const Point& p2, double a)
    {
        constexpr double SQRT_2_4 = sqrt(2)/4;

        auto ca = 0.5 + 0.25*a/(1 - a*SQRT_2_4);      /// smaller

        return ca*ca*(p1 - p2).Sqr();
    };

    ///--------------------------------------------------------------------------

    struct FractalGraph: public Drawable
    {
        std::vector<Point> fdata;   /// all fractal points

        std::vector<Point> data;
        double a;
        unsigned N;

        pixel color;

        FractalGraph
        (
            const std::vector<Point>& data,
            double a = 1.0,
            unsigned N = 8,
            const pixel& color = pixel(0, 0, 200)
        )
        : data(data), a(a), N(N), color(color)
        {

            fdata = CreateFractal(data, a, N);
        }

        /// get a 'node' at position 'idx' at a given 'level'
        ///
        std::pair<const Point&, const Point&> GetNode(unsigned idx, unsigned level) const
        {
            /// n -> 4*(n - 1) + 1  =>
            /// idx -> 4*idx

            /// n << (2*p) == n*4^p
            ///
            unsigned q1 = idx     << (2*(N - level));
            unsigned q2 = (idx+1) << (2*(N - level));

            return {fdata[q1], fdata[q2]};
        }

        /// 'Method 1'
        ///
        int TestPixel(const TChart* chart, double x, double y) const
        {
            auto s2r = [chart](auto x, auto y)
            {
               return chart->ScreenToReal(Point(x, y));
            };

            /*auto r2s = [chart](const Point& p)
            {
                return chart->RealToScreen(p);
            };*/

            static std::vector<Node> nodes;

            const double M = std::max(chart->K.x, -chart->K.y); /// magnitude
            const double M2 = M*M;

            auto p = s2r(x, y);

            /// initialization

            nodes.clear();

            for(unsigned k = 0; k < data.size()-1; k++)
            {
                nodes.emplace_back(data[k], data[k+1]);
            }

            /// search

            while(nodes.size())
            {
                auto node = nodes.back();
                nodes.pop_back();

                auto R2 = GetSegmentBorder(node.p1, node.p2, a);
                auto p0 = (node.p1 + node.p2)/2;

                double d2 = (p - p0).Sqr2();

                if(d2 < R2)    /// point is inside the border
                {
                    /// draw current border
                    //DrawCircle(r2s(p0), sqrt(M2*R2), pixel(200, 100, 100, 100));

                    if(0.5*M2*R2 <= 1)   /// point is inside (~)
                                         /// a smaller than a pixel segment
                    {
                        /// draw this pixel and finish the search

                        //PutWuPixel(r2s(p0), color);
                        //PutWuPixel(x+0.5, y+0.5, color);

                        return 1;

                        //break;
                    }
                    else            /// subdivision needed
                    {
                        auto q = node.split(a);

                        nodes.reserve(nodes.size() + q.size());
                        nodes.insert(nodes.end(), q.begin(), q.end());
                    }
                }
            }

            return 0;
        }

        /// 'Method 2'
        ///
        /// returns the distance from 'p' to the fractal
        /// (after at most 'N' iterations)
        ///
        double Distance(const Point& p) const
        {
            static std::deque<Node> nodes;

            /// distance from [p] to a line segment (p1, p2)
            auto L = [](const Point& p, const Point& p1, const Point& p2)
            {
                auto dp = p2 - p1;

                double t = ((p - p1) & dp)/dp.Sqr2();

                Point q = t < 0 ? p1 :
                          t > 1 ? p2 :
                                  p1 + dp*t;

                return (p - q).Norm2();
            };

            /// initialization

            nodes.clear();

            for(unsigned k = 0; k < data.size()-1; k++)
            {
                nodes.emplace_back(data[k], data[k+1]);
            }

            double D = INFINITY;

            /// search

            while(nodes.size())
            {
                //auto node = nodes.back();
                //nodes.pop_back();

                auto node = nodes.front();
                nodes.pop_front();

                if(node.level > N)
                    break;

                auto R2 = GetSegmentBorder(node.p1, node.p2, a);
                auto p0 = (node.p1 + node.p2)/2;

                /// distance between 'p' and the node center
                double d2 = (p - p0).Sqr2();

                if(d2 < R2)    /// point is inside the border
                {
                    if(node.level == N)
                    {
                        //D2.emplace_back(d2);

                        double d = L(p, node.p1, node.p2);

                        if(!std::isfinite(D) || d < D)
                            D = d;
                    }

                    auto q = node.split(a);

                    ///nodes.reserve(nodes.size() + q.size());
                    nodes.insert(nodes.end(), q.begin(), q.end());
                }
            }

            return D;
        }

        /// 'Method 3'
        ///
        /// uses a pre-computed list of points for 'N' subdivisions
        ///
        double pDistance(const Point& p) const
        {
            /// distance from [p] to a line segment (p1, p2)
            auto L = [](const Point& p, const Point& p1, const Point& p2)
            {
                auto dp = p2 - p1;

                double t = ((p - p1) & dp)/dp.Sqr2();

                Point q = t < 0 ? p1 :
                          t > 1 ? p2 :
                                  p1 + dp*t;

                return (p - q).Norm2();
            };

            ///--------------------------------------------

            double D = INFINITY;

            ///level \in [0, N]
            ///n[level] = (n[0] - 1)*4^level + 1

            unsigned n = ((data.size() - 1) << (2*N)) + 1;

            unsigned idx = 0;
            while(idx < n-1)
            {
                ///Generating a list of nodes from lower levels
                ///that contain the current one and check 'p'
                ///for being within their borders starting from level 0
                ///for quick node discarding

                /// current node:   GetNode(idx,     N)
                /// one level down: GetNode(idx/4,   N-1)
                /// general:        GetNode(idx/4^k, N-k)

                /// checking  if 'p' is within node border (at all levels)
                ///
                for(int k = N; k >= 0; k--)
                {
                    unsigned q0 = idx >> (2*k);     /// base index of the level (N-k)
                    auto node = GetNode(q0, N-k);

                    const auto& p1 = node.first;
                    const auto& p2 = node.second;

                    auto R2 = GetSegmentBorder(p1, p2, a);
                    auto p0 = (p1 + p2)/2;

                    /// distance between 'p' and the node center
                    double d2 = (p - p0).Sqr2();

                    if(d2 >= R2)    /// point is outside of the border
                    {
                        /// discard several nodes
                        /// depending on the current level

                        idx = std::max(idx, q0 + (1 << (2*k)) - 1);

                        break;
                    }

                    if(k == 0)
                    {
                        /// at this points 'p' is inside the borders
                        /// of the node at all levels

                        double d = L(p, p1, p2);

                        if(!std::isfinite(D) || d < D)
                            D = d;
                    }
                }

                idx++;
            }

            return D;
        }

        void Draw(const TChart* chart) const override
        {
            /*auto s2r = [chart](auto x, auto y)
            {
               return chart->ScreenToReal(Point(x, y));
            };*/

            auto r2s = [chart](const Point& p)
            {
               return chart->RealToScreen(p);
            };

            //double M = std::max(chart->K.x, -chart->K.y); /// magnitude

            ///----------------------------------------------------------------

            //Rect Screen = Rect(Point(0, 0), AppSize, false);
            //auto R = GetFractalBorder(data, a);

            /// fractal search area (screen coords)
            //iRect sR = Screen & Rect(r2s(R.A), r2s(R.B));

            ///----------------------------------------------------------------

            //DrawRect(sR, color);

            ///----------------------------------------------------------------

            /**double d;

            for(int x = sR.A.x; x <= sR.B.x; x++)
            {
                for(int y = sR.A.y; y <= sR.B.y; y++)
                {
                    /// Method 1------------------------------

                    //if(TestPixel(chart, x, y))
                        //PutPixel(x, y, color);

                    /// Method 2------------------------------

                    //d = M*Distance(s2r(x, y));

                    //if(d >= 0 && d <= 1)
                        //PutPixel(x, y, (1 - d)*color);

                    /// Method 3------------------------------

                    d = M*pDistance(s2r(x, y));

                    if(d >= 0 && d <= 1)
                        PutPixel(x, y, (1 - d)*color);
                }
            }*/

            /// Method 4------------------------------

            for(unsigned k = 0; k < fdata.size()-1; k++)
            {
                DrawLine(r2s(fdata[k]), r2s(fdata[k+1]), color);
            }
        }
    };

    ///--------------------------------------------------------------------------

    /// draw fractal borders (for each guideline segment)
    ///
    struct FBGraph: public Drawable
    {
        std::vector<Point> data;
        double a, ma;
        Rect R;

        pixel color;

        FBGraph(const std::vector<Point>& data, double a = 1.0, const pixel& color = pixel(200, 0, 0, 0))
               : data(data), a(a), color(color)
        {
            auto q = a*sqrt(2)/4;

            ma = q < 1 ? q/(1 - q) : 0.0;
            ///ma = q < 1 ? 0.25*a/(1 - q) : 0.0;      /// smaller

            //ma = a/4 + ((a/4)*sqrt(2))*a/4 + ((a/4*sqrt(2))*a/4)*sqrt(2)*a/4;

            R = GetFractalBorder(data, a);

            ///auto ca = 0.5 + 0.25*a/(1 - a*sqrt(2)/4);      /// smaller
        }

        void Draw(const TChart* chart) const override
        {
            for(unsigned k = 0; k < data.size()-1; k++)
            {
                const auto& p1 = chart->RealToScreen(data[k]);
                const auto& p2 = chart->RealToScreen(data[k+1]);

                pixel side_color = 0.7*color;

                ///hexagon
                if(0)
                {
                    auto q = a*sqrt(2)/4;
                    double ma = q < 1 ? q/(1 - q) : 0.0;

                    auto l = ma*(p2 - p1);
                    auto n = l.ortho();

                    /// parallel sides
                    DrawDDALine(p1+n, p2+n, side_color);
                    DrawDDALine(p1-n, p2-n, side_color);

                    /// caps
                    //DrawDDALine(p1+n, p1-n, color);
                    //DrawDDALine(p2+n, p2-n, color);

                    /// 'edges'
                    DrawDDALine(p1+n, p1-l, side_color);
                    DrawDDALine(p1-n, p1-l, side_color);

                    DrawDDALine(p2+n, p2+l, side_color);
                    DrawDDALine(p2-n, p2+l, side_color);
                }

                ///circle
                if(1)
                {
                    /*auto q = a*sqrt(2)/4;
                    double ma = q < 1 ? 0.25*a/(1 - q) : 0.0;      /// smaller

                    DrawDDACircle((p1+p2)/2, (ma + 0.5)*(p1 - p2).Norm2(), color);*/

                    auto R2 = GetSegmentBorder(p1, p2, a);

                    DrawDDACircle((p1+p2)/2, sqrt(R2), color);
                }

                ///------------------------------------------

                /// fractal border
                if(0)
                {
                    DrawRect
                    (
                        chart->RealToScreen(R.A),
                        chart->RealToScreen(R.B),
                        color
                    );
                }
            }
        }
    };
}

///-------------------------------------------

#endif // RFRACTAL_H_INCLUDED
