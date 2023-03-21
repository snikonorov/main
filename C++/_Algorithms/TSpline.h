#define _CRT_SECURE_NO_WARNINGS

#ifndef TSPLINE_H_INCLUDED
#define TSPLINE_H_INCLUDED

#include "SGraphics.h"
#include <complex>

struct Interval
{
    struct Range
    {
        double x0, x1;

        double length() const
        {
            return x1 - x0;
        }

        bool contains(double x) const
        {
            return x >= x0 && x <= x1;
        }
    };

    Range range;

    Interval(const Range& range = Range{0.0, 0.0}): range(range)
    {
    }

    const double& x0() const
    {
        return range.x0;
    }

    const double& x1() const
    {
        return range.x1;
    }

    bool contains(double x) const
    {
        return range.contains(x);
    }

    bool operator<(double x) const
    {
        return range.x0 < x;
    }

    bool operator>(double x) const
    {
        return range.x1 > x;
    }
};

struct Cubic: public Interval
{
    ///double d, c, b, a;
    double a, b, c, d;

	Cubic() = default;

	Cubic(double d, const Range& range = Range()): Interval{range}, a(0), b(0), c(0), d(d)
	{
	}

	Cubic(double a, double b, double c, double d, const Range& range = Range()): Interval{range}, a(a), b(b), c(c), d(d)
	{
	}

    double operator()(double x) const
    {
        const double dx = x - x0();
        return d + dx*(c + dx*(b + dx*a));
    }

    std::vector<double> zeroes() const
    {
        std::vector<double> res;

        if(a != 0)
        {
            double q0 = x0() - b/(3*a);

            double r = b*b - 3*a*c;
            double w = b*(2*b*b - 9*a*c) + 27*a*a*d;

            const double c0 = std::pow(2.0, 1./3);
            const std::complex<double> cr(1, sqrt(3));
            const std::complex<double> cl(1, -sqrt(3));

            std::complex<double> q = std::pow(std::sqrt(std::complex<double>(w*w - 4*r*r*r)) - w, 1./3);

            std::complex<double> z;

            z = q0 + (c0*r/q + q/c0)/(3*a);

            if(std::abs(z.imag()) < 1e-12 && contains(z.real()))
                res.push_back(z.real());

            z = q0 - (cr*r/(c0*q) + cl*q/2.0)/(3*a*c0);

            if(std::abs(z.imag()) < 1e-12 && contains(z.real()))
                res.push_back(z.real());

            z = q0 - (cl*r/(c0*q) + cr*q/2.0)/(3*a*c0);

            if(std::abs(z.imag()) < 1e-12 && contains(z.real()))
                res.push_back(z.real());
        }
        else
        {
            double z;

            if(b != 0)
            {
                double q = c*c - 4*b*d;

                if (q >= 0)
                {
                    q = sqrt(q);

                    z = x0() - (c + q)/(2*b);

                    if(contains(z))
                        res.push_back(z);

                    z = x0() - (c - q)/(2*b);

                    if(contains(z))
                        res.push_back(z);
                }
            }
            else
            {
                if(c != 0)
                {
                    z = x0() - d/c;

                    if(contains(z))
                        res.push_back(z);
                }
            }
        }

        return res;
    }

    /// {min, max}
    std::pair<Point, Point> bounds() const
    {
        static std::vector<Point> v(4);

        auto& f = *this;

        ///res.first = min(f(x0()), f(x1()), f(m1), f(m2));
        ///res.second = min(f(x0()), f(x1()), f(m1), f(m2));

        v[0] = {x0(), d};
        v[1] = {x1(), f(x1())};
        v[2] = {x0(), d};
        v[3] = {x0(), d};

        if (a != 0)
        {
            double m1 = x0(), m2 = x0();

            double q = b*b - 3*a*c;

            if (!(q < 0))
            {
                q = sqrt(q);

                m1 = x0() - (b + q)/(3*a);
                m2 = x0() - (b - q)/(3*a);

                if (contains(m1))
                    v[2] = {m1, f(m1)};

                if (contains(m2))
                    v[3] = {m2, f(m2)};
            }
        }
        else if (b != 0)
        {
            double m = x0() - c/(2*b);

            if (contains(m))
                v[2] = {m, f(m)};
        }

        std::sort(v.begin(), v.begin()+4, [](const Point& a, const Point& b){ return a.y < b.y; });

        return {v[0], v[3]};
    }

    /// {{Point(x, f(x)), f''(x)}, ...}
    std::vector<std::pair<Point, double>> extrema() const
    {
        std::vector<std::pair<Point, double>> res;

        double m;
        auto& f = *this;

        if (a != 0)
        {
            double q = b*b - 3*a*c;

            if (!(q < 0))
            {
                if (q == 0)
                {
                    /// saddle point
                    m = x0() - b/(3*a);

                    if (contains(m))
                        res.emplace_back(Point(m, f(m)), 0.0);
                }
                if (q != 0)
                {
                    q = sqrt(q);

                    m = x0() - (b + q)/(3*a);

                    if (contains(m))
                        res.emplace_back(Point(m, f(m)), -2*q);

                    m = x0() - (b - q)/(3*a);

                    if (contains(m))
                        res.emplace_back(Point(m, f(m)), +2*q);
                }
            }
        }
        else if (b != 0)
        {
            m = x0() - c/(2*b);

            if (contains(m))
                res.emplace_back(Point(m, f(m)), 2*b - 3*a*c/b);
        }

        return res;
    }

    std::pair<Cubic, Cubic> split(double x) const
    {
        std::pair<Cubic, Cubic> res;

        auto& q = *this;

        if(x < q.range.x0)
        {
            res.second = q;
        }
        else if(x > q.range.x1)
        {
            res.first = q;
        }
        else
        {
            res.first = q;
            res.first.range.x1 = x;

            double dx = x - q.range.x0;

            res.second.a = q.a;
            res.second.b = q.b + 3*q.a*dx;
            res.second.c = q.c + dx*(2*q.b + 3*q.a*dx);
            res.second.d = q.d + dx*(q.c + dx*(q.b + q.a*dx));

            res.second.range.x0 = x;
            res.second.range.x1 = q.range.x1;
        }

        return res;
    }

    const double& y0() const
    {
        return d;
    }

    const double y1() const
    {
        return (*this)(x1());
    }

	void operator+=(double c0)
	{
		d += c0;
	}

	void operator-=(double c0)
	{
		d -= c0;
	}

	void operator*=(double a0)
	{
		a *= a0;
		b *= a0;
		c *= a0;
		d *= a0;
	}

	void operator/=(double a0)
	{
		a /= a0;
		b /= a0;
		c /= a0;
		d /= a0;
	}

	/*Cubic operator-(double c0) const
	{
		return Cubic(a, b, c, d - c0, range);
	}*/
};

template<class TPoly = Cubic>
struct TSpline : public std::vector<TPoly>
{
    enum TSplineImportMode
    {
        Points = 1,
        Cubics = 2
    };

    enum TSplineStyle
    {
        Finite = 0,
        ConstBorders = 1,
        Periodic = 2,
        Extrapolate = 3
    };

    enum TSplineType
    {
        Linear = 1,
        Cubic = 3,
		Akima = 4
    };

    TSplineStyle style;

	template<class TPoint>
    double Create
    (
        const std::vector<TPoint>& l,
        TSplineStyle s = Finite,
        TSplineType t = Cubic,
        std::pair<double, double> h = {nan(""), nan("")}     /// 2-nd derivative on the boundaries
    )
    {
        style = s;

        auto& res = *this;
        double Amplitude = 0;

        if(l.size() > 2)
        {
            const auto& N = l.size();
            double __min = l[N-1].y, __max = l[N-1].y;

            res.resize(N-1);

            ///-----------------------------------------------------------

			if (t == TSplineType::Linear)
            {
                for (unsigned k = 0; k < N-1; k++)
                {
                    res[k].range = {l[k].x, l[k+1].x};

                    res[k].a = 0;
                    res[k].b = 0;
                    res[k].c = (l[k+1].y - l[k].y)/(l[k+1].x - l[k].x);
                    res[k].d = l[k].y;

                    if(l[k].y < __min) __min = l[k].y;
                    if(l[k].y > __max) __max = l[k].y;
                }

                Amplitude = __max - __min;
            }
			else if(t == TSplineType::Akima && l.size() > 3)
			{
				for (unsigned k = 0; k < N-1; k++)
                {
                    res[k].range = {l[k].x, l[k+1].x};

					double dx = l[k+1].x - l[k].x;
					double dy = l[k+1].y - l[k].y;

					double Dy = dy/dx;

					double Dy_0 = k == 0 ? Dy
										 : (l[k+1].y - l[k-1].y)/(l[k+1].x - l[k-1].x);

					double Dy_1 = k == N-2 ? Dy
										   : (l[k+2].y - l[k].y)/(l[k+2].x - l[k].x);

                    res[k].a = (Dy_1 + Dy_0 - 2*Dy)/(dx*dx);
                    res[k].b = (3*Dy - Dy_1 - 2*Dy_0)/dx;

                    res[k].c = Dy_0;
                    res[k].d = l[k].y;

                    if(l[k].y < __min) __min = l[k].y;
                    if(l[k].y > __max) __max = l[k].y;
                }

				Amplitude = __max - __min;
			}
            else if (t == TSplineType::Cubic)
            {
                /// boundary second derivative values estimation-----------------------------

                double h0 = 0, h1 = 0;

                if(l.size() > 3)
                {
                    if(std::isfinite(h.first))
                    {
                        h0 = h.first;
                    }
                    else
                    {
                        const double& x0 = l[0].x;
                        const double& x1 = l[1].x;
                        const double& x2 = l[2].x;
                        const double& x3 = l[3].x;

                        const double& y0 = l[0].y;
                        const double& y1 = l[1].y;
                        const double& y2 = l[2].y;
                        const double& y3 = l[3].y;

                        h0 = 2*(((3*x0 - x1 - x2 - x3)*y0)/((x0 - x1)*(x0 - x2)*(x0 - x3)) +
                                ((-2*x0 + x2 + x3)*y1)/((x0 - x1)*(x1 - x2)*(x1 - x3)) +
                                ((-2*x0 + x1 + x3)*y2)/((x0 - x2)*(-x1 + x2)*(x2 - x3)) +
                                ((-2*x0 + x1 + x2)*y3)/((x0 - x3)*(-x1 + x3)*(-x2 + x3)));
                    }

                    if(std::isfinite(h.second))
                    {
                        h1 = h.second;
                    }
                    else
                    {
                        const double& x0 = l[N-4].x;
                        const double& x1 = l[N-3].x;
                        const double& x2 = l[N-2].x;
                        const double& x3 = l[N-1].x;

                        const double& y0 = l[N-4].y;
                        const double& y1 = l[N-3].y;
                        const double& y2 = l[N-2].y;
                        const double& y3 = l[N-1].y;

                        h1 = 2*(-(((x1 + x2 - 2*x3)*y0)/((x0 - x1)*(x0 - x2)*(x0 - x3))) +
                                  ((x0 + x2 - 2*x3)*y1)/((x0 - x1)*(x1 - x2)*(x1 - x3)) +
                                  ((x0 + x1 - 2*x3)*y2)/((x0 - x2)*(-x1 + x2)*(x2 - x3)) +
                                  ((x0 + x1 + x2 - 3*x3)*y3)/((x0 - x3)*(-x1 + x3)*(-x2 + x3)));
                    }
                }

                ///--------------------------------------------------------------------------

                res[0].b = h0/2;
                ///res[N-1].b -> h1/2

                /**
                C[0] C[0] ...
                A[1] A[1] B[1] ...
                      ...
                ...  A[n-1] C[n-1]

                n = N-2
                X[0..n-1] <- res[1..n].b
                */

                const unsigned n = N-2;

				auto v = std::vector<double>(n-1);

                const unsigned i = 0;

                double A = l[i+1].x - l[i].x;
                double C = 2*(l[i+2].x - l[i].x);
                double F = 3*( (l[i+2].y - l[i+1].y)/(l[i+2].x - l[i+1].x) -
                               (l[i+1].y - l[i].y)/(l[i+1].x - l[i].x) );

                double q = C;
                res[i+1].b = (F - A*res[0].b)/q;

                for(unsigned i = 1; i < n; i++)
                {
                    double A = l[i+1].x - l[i].x;
                    double C = 2*(l[i+2].x - l[i].x);
                    double F = 3*( (l[i+2].y - l[i+1].y)/(l[i+2].x - l[i+1].x) -
                                   (l[i+1].y - l[i].y)/(l[i+1].x - l[i].x) )
                               - (i == n-1)*(l[i+2].x - l[i+1].x)*h1/2;

                    v[i-1] = A/q;
                    q = C - A*v[i-1];
                    res[i+1].b = (F - A*res[i].b)/q;
                }

                for(unsigned j = n-1; j > 0; j--)
                    res[j].b -= v[j-1]*res[j+1].b;

                ///-----------------------------------------------------------

                for(size_t k = 0; k < N - 1; k++)
                {
                    double dB = (k == N-2) ? (h1/2 - res[k].b)
                                           : (res[k+1].b - res[k].b);

                    double dx = l[k+1].x - l[k].x;

                    res[k].a = dB/(3*dx);
                    res[k].c = (l[k+1].y - l[k].y)/dx - (dB/3 + res[k].b)*dx;
                    res[k].d = l[k].y;

                    res[k].range = {l[k].x, l[k+1].x};

                    if(l[k].y < __min) __min = l[k].y;
                    if(l[k].y > __max) __max = l[k].y;
                }

                Amplitude = __max - __min;
            }
            else
                puts("Unknown spline type");
        }
        else if(l.size() == 2)
        {
            Amplitude = fabs(l[0].y - l[1].y);

            res.resize(1);

            res[0].range.x0 = l[0].x;
            res[0].range.x1 = l[1].x;

            res[0].a = res[0].b = 0;
            res[0].c = (l[1].y - l[0].y)/(l[1].x - l[0].x);
            res[0].d = l[0].y;
        }

        return Amplitude;
    }

    TSpline(): style(Finite)
    {
    }

	TSpline(double c, double x0, double x1): style(Finite)
    {
		auto& res = *this;

		res.resize(1);

		res[0] = c;
		res[0].range.x0 = x0;
		res[0].range.x1 = x1;
    }

	template<class TPoint, class = decltype(std::declval<TPoint>().x)>
	TSpline(const std::vector<TPoint>& l, TSplineStyle s = Finite, TSplineType t = Cubic, std::pair<double, double> h = {nan(""), nan("")})
    {
        Create(l, s, t, h);
    }

	template<class TPoint, class = decltype(std::declval<TPoint>().x)>
	TSpline(const std::vector<TPoint>& l, TSplineType t, std::pair<double, double> h = {nan(""), nan("")})
    {
        Create(l, Finite, t, h);
    }

    TSpline(const std::vector<double>& l, TSplineStyle s = Finite, TSplineType t = Cubic, std::pair<double, double> h = {nan(""), nan("")})
    {
        std::vector<Point> p(l.size());

        for(unsigned k = 0; k < p.size(); k++)
        {
            p[k].x = k;
            p[k].y = l[k];
        }

        Create(p, s, t,h);
    }

    TSpline(const std::initializer_list<Point>& l, TSplineStyle s = Finite, TSplineType t = Cubic, std::pair<double, double> h = {nan(""), nan("")}): style(s)
    {
        Create(std::vector<Point>(l), s, t, h);
    }

    TSpline(unsigned n, TSplineStyle s = Finite): std::vector<TPoly>(n), style(s)
    {
    }

	TSpline(const std::vector<TPoly>& s): style(Finite), std::vector<TPoly>(s)
    {
    }

	template<class Q>
	TSpline(const TSpline<Q>& s): style((TSpline<TPoly>::TSplineStyle)s.style), std::vector<TPoly>(s.size())
    {
		for(unsigned k = 0; k < s.size(); k++)
		{
			(*this)[k] = s[k];
		}
    }

    using std::vector<TPoly>::size;
    using std::vector<TPoly>::front;
    using std::vector<TPoly>::back;
    using std::vector<TPoly>::at;
    using std::vector<TPoly>::clear;
    using std::vector<TPoly>::insert;

    int find(double x) const
    {
        if(size() && std::isfinite(x))
        {
            if(x < front().x0())
                return style == Extrapolate ? 0 : -1;
            else if(x > back().x1())
                return style == Extrapolate ? size()-1 : size();

            //------------------------------------------------------------------

            /// binary
            int a = 0, b = size();

            while(a <= b)
            {
                auto c = a + (b - a)/2;

                if(at(c).contains(x))
                    return c;
                else if(at(c) > x)
                    b = c - 1;
                else
                    a = c + 1;
            }

            //------------------------------------------------------------------
        }

        return size();
    }

    operator const std::vector<Point>() const
    {
        std::vector<Point> res;

        if(size())
        {
            res.resize(size() + 1);

            for(size_t k = 0; k < size(); k++)
            {
                res[k] = {at(k).x0(), at(k).y0()};
            }

            res[size()] = {back().x1(), back().y1()};
        }

        return std::move(res);
    }

    std::vector<Point> insert(const Point& p)
    {
        unsigned idx = p.x < x0() ? 0 : p.x > x1() ? size() + 1 : find(p.x) + 1;

        std::vector<Point> v = *this;
        v.insert(v.begin() + idx, p);

        clear();
        Create(v);

        return v;
    }

    void insert(double x);

    const double& x0() const
    {
        return front().x0();
    }

    const double& x1() const
    {
        return back().x1();
    }

    const double& y0() const
    {
        return front().y0();
    }

    const double y1() const
    {
        return back().y1();
    }

    double operator()(double x) const
    {
        int idx = find(x);

        //return (idx >= 0 && idx < size()) ? at(idx)(x) : 0;

        if(idx >= 0 && (idx < 0 || unsigned(idx) < size()))     /// ~ idx \in [0..size()-1]
        {
            return at(idx)(x);
        }
        else
        {
            switch(style)
            {
                case Finite:
                    return 0;

                case ConstBorders:
                    return idx < 0 ? y0() : y1();

                //case Extrapolate:
                    //return x < x0() ? front()(x) : back()(x);

                case Periodic:
                {
                    const double T = x1() - x0();

                    x -= x > x1() ? +T*int(1 + (x - x1())/T)
                                  : -T*int(1 + (x0() - x)/T);

                    return at(find(x))(x);
                }

                default:
                    break;
            }
        }

		return nan("");
    }

	/*std::vector<Point> GetSample(unsigned N) const
	{
		std::vector<Point> res(N);

		for (unsigned k = 0; k < N; k++)
		{
			const double x = x0() + (x1() - x0())*k/(N - 1);
			res[k] = Point(x, (*this)(x));
		}

		return res;
	}*/

	void operator+=(double c0)
	{
		for(auto& c: *this)
		{
			c += c0;
		}
	}

	void operator-=(double c0)
	{
		for(auto& c: *this)
		{
			c -= c0;
		}
	}

	void operator*=(double a0)
	{
		for(auto& c: *this)
		{
			c *= a0;
		}
	}

	void operator/=(double a0)
	{
		for(auto& c: *this)
		{
			c /= a0;
		}
	}

	TSpline operator-(double c0) const
	{
		TSpline<TPoly> res = *this;

		for(auto& c: res)
		{
			c -= c0;
		}

		return res;
	}

	TSpline operator*(double a0) const
	{
		TSpline<TPoly> res = *this;

		for(auto& c: res)
		{
			c *= a0;
		}

		return res;
	}

	TSpline operator/(double a0) const
	{
		TSpline<TPoly> res = *this;

		for(auto& c: res)
		{
			c /= a0;
		}

		return res;
	}

	void operator=(double c)
    {
		auto& res = *this;

		res[0] = c;
		res[0].range.x0 = x0();
		res[0].range.x1 = x1();

		res.resize(1);
    }
};

using spline = TSpline<Cubic>;

template<>
void TSpline<Cubic>::insert(double x)
{
    if(x < front().x0() || x > back().x1())
    {
        unsigned idx = x < x0() ? 0 : x > x1() ? size() + 1 : find(x) + 1;

        std::vector<Point> v = *this;
        v.insert(v.begin() + idx, Point(x, (*this)(x)));

        clear();
        Create(v);
    }
    else
    {
        unsigned idx = find(x) + 1;

        auto& q = this->at(idx - 1);
        auto res = q.split(x);

        q = res.first;

        this->insert
        (
            this->begin() + idx,
            res.second
        );
    }
}

///--------------------------------------------------------------------------

struct CSpline
{
    spline x;
    spline y;

    template<class TContainer>
    CSpline(const TContainer& data, spline::TSplineType t = spline::Cubic, bool normalize = true)
    {
        Create(data, t, normalize);
    }

    template<class TContainer>
    void Create(TContainer data, spline::TSplineType t = spline::Cubic, bool normalize = true)
    {
        if(data.size())
        {
            /// creating timeline---------------------------

            std::vector<double> T(data.size());

            double L = 0;
            T[0] = 0;

            for(unsigned k = 1; k < data.size(); k++)
            {
                double l = (data[k] - data[k-1]).Norm2();
                T[k] = T[k-1] + l;

                if(normalize)
                    L += l;
            }

            /// separating x and y components---------------

            std::vector<Point> X(data.size());

            for(unsigned k = 0; k < data.size(); k++)
            {
                double t = T[k];

                if(normalize)
                    t /= L;

                X[k].x = t;
                X[k].y = data[k].x;

                data[k].x = t;
            }

            /// spline creation-----------------------------

            x.Create(X, spline::Finite, t);
            y.Create(data, spline::Finite, t);
        }
    }

    Point operator()(double t) const
    {
        return {x(t), y(t)};
    }

    bool empty() const
    {
        return x.empty() || y.empty();
    }

    void clear()
    {
        x.clear();
        y.clear();
    }
};

#endif // TSPLINE_H_INCLUDED
