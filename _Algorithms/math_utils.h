#ifndef MATH_UTILS_H_INCLUDED
#define MATH_UTILS_H_INCLUDED

#include <math.h>
#include "functional_utils.h"

#include "TMatrix.h"

inline double Rand(double a = -1, double b = 1)
{
	return a + (b - a)*rand()/RAND_MAX;
}

template<class T>
inline T sign(const T& x)
{
    return (x > 0) - (x < 0);
}

/// straight line approximation (not necessarily a function);
/// returns 2 points that represent the line
template<class TContainer>
std::pair<Point, Point> line_approx(const TContainer& data)
{
    auto p0 = average(data);

    double C = 0;
    double S = 0;

    for(const auto& p: data)
    {
        double x = p.x - p0.x;
        double y = p.y - p0.y;

        C += (y*y - x*x)/2;
        S += x*y;
    }

    Point l;

    double h = sqrt(C*C + S*S);

    if(h > 1e-15)
    {
        double ca = sqrt((h - C)/(2*h));
        double sa = (S < 0 ? -1 : 1)*sqrt(1 - ca*ca);

        l = Point(ca, sa);
    }

    ///-------------------------------------------------------------------

    double u1 = (data[0] - p0) & l;
    double u2 = u1;

    for(unsigned k = 1; k < data.size(); k++)
    {
        const auto& p = data[k];

        double u = (p - p0) & l;

        if(u > u1)
            u1 = u;

        if(u < u2)
            u2 = u;
    }

    return std::make_pair(p0 + u2*l, p0 + u1*l);
}

/// general approximation
///
///  f: f(n, x) - n-th basis function
///
template<class TContainer, class TFunction>
decltype(auto) approximation(const TContainer& data, const TFunction& f, unsigned N)
{
    auto M = matrix(N, N, [&data, &f](auto row, auto col)
									  {
										  double res = 0;

										  for(const auto& p: data)
										  {
											  res += f(row, p.x)*f(col, p.x);
										  }

										  return res;
									  });

    auto b = matrix(N, 1, [&data, &f](auto row, auto col)
									{
										double res = 0;

										for(const auto& p: data)
										{
											res += p.y*f(row, p.x);
										}

										return res;
									});

    auto a = (!M) * b;

    return a._data;
}

/// general polynomial approximation
template<class TContainer>
decltype(auto) approximation(const TContainer& data, unsigned P = 1)
{
    auto M = matrix(P+1, P+1, [&data](auto row, auto col)
                              {
                                  double res = 0;

                                  for(const auto& p: data)
                                  {
                                      res += pow(p.x, row)*pow(p.x, col);
                                  }

                                  return res;
                              });

    auto b = matrix(P+1, 1, [&data](auto row, auto col)
                            {
                                double res = 0;

                                for(const auto& p: data)
                                {
                                    res += p.y*pow(p.x, row);
                                }

                                return res;
                            });

    auto a = (!M) * b;

    return a._data;
}

decltype(auto) approximate(const std::vector<Point>& data, unsigned P = 1)
{
    auto M = matrix(P+1, P+1, [&data](auto row, auto col)
                              {
                                  double res = 0;

                                  for(const auto& p: data)
                                  {
                                      res += pow(p.x, row)*pow(p.x, col);
                                  }

                                  return res;
                              });

    auto b = matrix(P+1, 1, [&data](auto row, auto col)
                            {
                                double res = 0;

                                for(const auto& p: data)
                                {
                                    res += p.y*pow(p.x, row);
                                }

                                return res;
                            });

    auto a = (!M) * b;
    ///auto a = GrLSolve(M, b);

    return [a](double x)
    {
        double res = 0;

        for(unsigned k = 0; k < a._data.size(); k++)
        {
            res += a._data[k]*pow(x, k);
        }

        return res;
    };
}

/// f0 ~ low-pass cut-off frequency (in Hz)
template<class TPoint>
std::vector<Point> smoothen(const std::vector<TPoint>& data, double f0, bool _linear_boundaries = false)
{
	const double& x0 = data.front().x;
    const double& x1 = data.back().x;

	double T = x1 - x0;
    double dt = T/(data.size() - 1);

    std::vector<Point> r;
    std::vector<Point> res;

    //auto f = spline(data);

    double dx = 1/(4*f0);

    auto x = range(x0 + dx/2, x1 - dx/2, round(T*4*f0) - 1);
	///auto x = range(x0 + dx, x1 - dx, round(T*4*f0) - 2);

    dx = (x.back() - x.front())/(x.size() - 1);

    res.resize(x.size());

    for(unsigned k = 1; k < x.size() - 1; k++)
    {
		unsigned idx0 = round((vbound(x[k] - dx, x0, x1) - x0)/dt);
        unsigned idx1 = round((vbound(x[k] + dx, x0, x1) - x0)/dt);

        r.assign(data.begin() + idx0, data.begin() + idx1);
		//println(r);

        res[k].x = x[k];
        res[k].y = approximate(r, 2)(x[k]);
    }

    r.assign(data.begin() + round((x[0] - x0)/dt), data.begin() + round((x[1] - x0)/dt));

    {
		auto f = approximate(r, 1);

		res[0].x = x0;
		res[0].y = f(x0);

		if(_linear_boundaries)
		{
			double _x = x[0];
			res.insert(res.begin() + 1, Point(_x, f(_x)));

			///_x = (x0 + x[0])/2;
			///res.insert(res.begin() + 1, Point(_x, f(_x)));
		}
	}

    r.assign(data.begin() + round((x[x.size() - 2] - x0)/dt), data.begin() + round((x[x.size() - 1] - x0)/dt));

    {
		auto f = approximate(r, 1);

		res.back().x = x1;
		res.back().y = f(x1);

		if(_linear_boundaries)
		{
			double _x = x.back();
			res.insert(res.end() - 1, Point(_x, f(_x)));

			///_x = (x1 + x.back())/2;
			///res.insert(res.end() - 1, Point(_x, f(_x)));
		}
	}

    return res;
}

/// approximating [data] by line segments
///
inline decltype(auto) linearize(const std::vector<Point>& data, const double r_max = 0.7)
{
    /// distance from [p] to a line (p0, dp)
	auto r = [](const auto& p0, const auto& dp, const Point& p)
    {
        double t = ((p - p0) & dp)/dp.Sqr2();
        return (p - p0 - dp*t).Norm2();
    };

    /// projection of [p] on the [line]
    auto projection = [](const Point& p, const auto& line)
    {
        const Point& p0 = line.first;
        const Point dp = line.second - line.first;

        double t = ((p - p0) & dp)/dp.Sqr2();
        return p0 + dp*t;
    };

    std::vector<Point> points;

    decltype(line_approx(data)) _P;

    unsigned idx = 0, k = 1;
    while(k < data.size())
    {
        auto s = slice(data, idx, k - idx + 1);
        auto P = line_approx(s);

        auto m = max(map([&r, &P](const Point& p)
                         {
                             return r(P.first, P.second - P.first, p);
                         },
                         s));

        if(m >= r_max) /// || (k == data.size() - 1 && !lines.size()))
        {
            if(idx == 0)
            {
                points.push_back(projection(data[idx], _P));
            }

            points.push_back(projection(data[k-1], _P));

            idx = k - 1;
        }
        else
            k++;

        _P = P;
    }

    points.push_back(projection(data.back(), _P));

    return points;
}

std::vector<Point> derivative(const std::vector<Point>& v)
{
    std::vector<Point> res = v;

    if (v.size() > 1)
    {
        for (unsigned k = 0; k < res.size(); k++)
        {
            res[k].y = k == 0              ? (v[k + 1].y - v[k].y) / (v[k + 1].x - v[k].x) :
                       k == res.size() - 1 ? (v[k].y - v[k - 1].y) / (v[k].x - v[k - 1].x) :
                                             (v[k + 1].y - v[k - 1].y) / (v[k + 1].x - v[k - 1].x);

        }
    }

    return res;
}

std::vector<Point> gauss_filter(const std::vector<Point>& v, double s = 10, unsigned short ns = 5)
{
    std::vector<Point> res = v;

    for (unsigned k = 0; k < res.size(); k++)
    {
        double w = 1;

        {
            int p = k - 1;
            double dx;

            while (p >= 0 && (dx = fabs(v[k].x - v[p].x)) <= ns * s)
            {
                double a = exp(-dx * dx / (2 * s*s));

                res[k].y += v[p].y*a;
                w += a;

                p--;
            }
        }

        {
            unsigned p = k + 1;
            double dx;

            while (p < res.size() && (dx = fabs(v[p].x - v[k].x)) <= ns * s)
            {
                double a = exp(-dx * dx / (2 * s*s));

                res[k].y += v[p].y*a;
                w += a;

                p++;
            }
        }

        res[k].y /= w;
    }

    return res;
}

/// spline functions--------------------------------------------
///
#ifdef TSPLINE_H_INCLUDED

/// general---------------------------------------

std::vector<double> zeroes(const spline& s)
{
    std::vector<double> res;

    for(const auto& c: s)
    {
        auto z = c.zeroes();

        for(const auto& x0: z)
        {
            res.push_back(x0);
        }
    }

    return res;
}

std::vector<double> extrema(const spline& s, double saddle_point_error = 0.0)
{
    std::vector<double> res;

    for (unsigned k = 0; k < s.size(); k++)
    {
        auto q = s[k].extrema();

        for (const auto& p: q)
        {
            if (std::abs(p.second) > saddle_point_error)  /// not a saddle point
                res.emplace_back(p.first.x);
        }
    }

    return res;
}

/**std::vector<Point> extrema_points(const spline& s, double saddle_point_error = 0.0)
{
    std::vector<Point> res;

    for (unsigned k = 0; k < s.size(); k++)
    {
        auto q = s[k].extrema();

        for (const auto& p: q)
        {
            if (std::abs(p.second) > saddle_point_error)  /// not a saddle point
                res.emplace_back(p.first);
        }
    }

    return res;
}*/

/// {minima, maxima}
std::pair<std::vector<double>, std::vector<double>> extrema_pair(const spline& s, double saddle_point_error = 0.0)
{
    std::pair<std::vector<double>, std::vector<double>> res;

    for (unsigned k = 0; k < s.size(); k++)
    {
        auto q = s[k].extrema();

        for (const auto& p: q)
        {
            if (std::abs(p.second) > saddle_point_error)  /// not a saddle point
			{
				if(p.second > 0)
					res.first.emplace_back(p.first.x);
				else
					res.second.emplace_back(p.first.x);
			}
        }
    }

    return res;
}

/// approximate gauss filter
///
std::vector<Point> spline_gauss_filter(const std::vector<Point>& v, double s = 10, unsigned short ns = 5)
{
    static spline G(zipmap([](double x){ return exp(-x*x/2); }, range(-9, 9, 101)));

	std::vector<Point> res = v;

    for (unsigned k = 0; k < res.size(); k++)
    {
        double w = 1;

        {
            int p = k - 1;
            double dx;

            while (p >= 0 && (dx = fabs(v[k].x - v[p].x)) <= ns * s)
            {
                double a = G(dx/s);		///exp(-dx * dx / (2 * s*s));

                res[k].y += v[p].y*a;
                w += a;

                p--;
            }
        }

        {
            unsigned p = k + 1;
            double dx;

            while (p < res.size() && (dx = fabs(v[p].x - v[k].x)) <= ns * s)
            {
                double a = G(dx/s);		///exp(-dx * dx / (2 * s*s));

                res[k].y += v[p].y*a;
                w += a;

                p++;
            }
        }

        res[k].y /= w;
    }

    return res;
}

///-----------------------------------------------

struct NotImplemented
{
};

struct InvalidConversion
{
};

/// {a, b, c,...}: ((a*(x-x0) + b)*(x-x0) + c)...
///
struct Polynom: public Interval, public std::vector<double>
{
    Polynom()
    {
    }

    Polynom(const Cubic& c): Interval{c.range}, std::vector<double>(4)
    {
        at(0) = c.a;
        at(1) = c.b;
        at(2) = c.c;
        at(3) = c.d;
    }

	/*Polynom(double c): std::vector<double>(1)
    {
        at(0) = c;
    }*/

    Polynom(const std::initializer_list<double>& l, const Range& range = Range()): Interval{range}, std::vector<double>(l)
    {
    }

    Polynom(const std::vector<double>& v, const Range& range = Range()): Interval{range}, std::vector<double>(v)
    {
    }

	Polynom(unsigned n, const Range& range = Range()): Interval{range}, std::vector<double>(n)
    {
    }

	void operator+=(double c0)
	{
		back() += c0;
	}

	void operator-=(double c0)
	{
		back() -= c0;
	}

    double operator()(double x) const
    {
        double res = 0;
        double dx = x - x0();

        for(unsigned q = 0; q < size(); q++)
        {
            res = res*dx + at(q);
        }

        return res;
    }

    const double& y0() const
    {
        return back();
    }

    const double y1() const
    {
        return (*this)(x1());
    }
};

/// derivative------------------------------------

inline Cubic D(const Cubic& c)
{
	return {0, 3*c.a, 2*c.b, c.c, c.range};
}

double D(const spline& s, double x)
{
    auto idx = s.find(x);
    return idx >= 0 && (idx < 0 || unsigned(idx) < s.size()) ? D(s[idx])(x) : 0;
}

inline Cubic D2(const Cubic& c)
{
	return {0, 0, 6*c.a, 2*c.b, c.range};
}

double D2(const spline& s, double x)
{
    auto idx = s.find(x);
    return idx >= 0 && (idx < 0 || unsigned(idx) < s.size()) ? D2(s[idx])(x) : 0;
}

inline Cubic Derivative(const Cubic& c)
{
	return {0, 3*c.a, 2*c.b, c.c, c.range};
}

inline Polynom Derivative(const Polynom& c)
{
	const unsigned p = c.size();

	Polynom res(std::max(0u, p - 1), c.range);

	for(unsigned k = 0; k < res.size(); k++)
	{
		res[k] = (p - k - 1)*c[k];
	}

	return res;
}

template<class Q>
TSpline<Q> Derivative(TSpline<Q> s)
{
    for(auto c = s.begin(); c != s.end(); c++)
    {
        *c = Derivative(*c);
    }

    return std::move(s);
}

/// integral--------------------------------------

inline Polynom Integral(const Polynom& p)
{
    Polynom res(p.size() + 1, p.range);

	res.back() = 0.0;

	for(unsigned k = 0; k < p.size(); k++)
	{
		res[k] = p[k]/(p.size() - k);
	}

	return res;
}

template<class T>
TSpline<Polynom> Integral(const TSpline<T>& s, bool force_finite = false)
{
	if (force_finite || s.style == TSpline<T>::Finite)
	{
	    TSpline<Polynom> res(s.size());

	    double d = 0;

	    for(unsigned k = 0; k < res.size(); k++)
        {
            res[k] = Integral(s[k]);

            res[k][4] = d;
            d = res[k](res[k].range.x1);
        }

	    return res;
	}
	else
		throw NotImplemented();
}

TSpline<Polynom> Integral(const spline& s, double x0, double x1)
{
	if (s.style == spline::Extrapolate)
	{
	    //TSpline<Polynom> res(s.size());

	    unsigned p = 0;

        while(s[p].range.x1 <= x0 && p < s.size()-1)
        {
            p++;
        }

        TSpline<Polynom> res(s.size() - p, TSpline<Polynom>::Extrapolate);
        double d = 0;

	    for(unsigned k = 0; k < res.size(); k++)
        {
            if (s[k+p].range.x0 >= x1)
            {
                res.resize(k+p+1);
                break;
            }

            res[k] = Integral(s[k+p]);

            if(k == 0)
            {
                res[k].range.x0 = x0;
            }
            else if (k == res.size()-1)
            {
                res[k].range.x1 = x1;
            }

            res[k][4] = d;
            d = res[k](res[k].range.x1);
        }

	    return res;
	}
	else
        throw NotImplemented();
}

///-----------------------------------------------

template<class Q>
inline double mean(const TSpline<Q>& s)
{
	auto q = Integral(s, true);

	return (q(s.x1()) - q(s.x0()))/(s.x1() - s.x0());
}

/// mean of the derivative of `s`
template<class Q>
inline double derivative_mean(const TSpline<Q>& s)
{
	return (s(s.x1()) - s(s.x0()))/(s.x1() - s.x0());
}

/// L2-norm squared
inline double norm2(const Cubic& C, double x1, double x2)
{
    const double& x0 = C.range.x0;
    const double& a = C.a;
    const double& b = C.b;
    const double& c = C.c;
    const double& d = C.d;

    return      x1*(-d*d + x1*(-c*d + x1*(-c*c/3 - 2*b*d/3 +
                x1*(-b*c/2 - a*d/2 + x1*(-b*b/5 - 2*a*c/5 +
                x1*(-a*b/3 - a*a*x1/7)))))) +
                x2*(d*d + x2*(c*d + x2*(c*c/3 + 2*b*d/3 +
                x2*(b*c/2 + a*d/2 + x2*(b*b/5 + 2*a*c/5 +
                x2*(a*b/3 + a*a*x2/7)))))) +
                x0*(x1*(2*c*d + x1*(c*c + 2*b*d + x1*(2*b*c + 2*a*d +
                x1*(b*b + 2*a*c + x1*(2*a*b + a*a*x1))))) +
                x2*(-2*c*d + x2*(-c*c - 2*b*d + x2*(-2*b*c - 2*a*d +
                x2*(-b*b - 2*a*c + x2*(-2*a*b - a*a*x2))))) +
                x0*(x1*(-c*c - 2*b*d + x1*(-3*b*c - 3*a*d + x1*(-2*b*b - 4*a*c +
                x1*(-5*a*b - 3*a*a*x1)))) +
                x2*(c*c + 2*b*d + x2*(3*b*c + 3*a*d + x2*(2*b*b + 4*a*c +
                x2*(5*a*b + 3*a*a*x2)))) +
                x0*(x1*(2*b*c + 2*a*d + x1*(2*b*b + 4*a*c + x1*(20*a*b/3 + 5*a*a*x1))) +
                x2*(-2*b*c - 2*a*d + x2*(-2*b*b - 4*a*c + x2*(-20*a*b/3 - 5*a*a*x2))) +
                x0*(x1*(-b*b - 2*a*c + x1*(-5*a*b - 5*a*a*x1)) +
                x0*(x1*(2*a*b + 3*a*a*x1) + x2*(-2*a*b - 3*a*a*x2) +
                x0*(-a*a*x1 + a*a*x2)) + x2*(b*b + 2*a*c +
                x2*(5*a*b + 5*a*a*x2))))));
}

/// L2-norm
double norm(const spline& s, double x0, double x1)
{
    if (s.style == spline::Extrapolate)
	{
        unsigned k = 0;

        while(s[k].range.x1 <= x0 && k < s.size()-1)
        {
            k++;
        }

        double res = 0;

        bool _f = true;

	    for(; k < s.size(); k++)
        {
            if (s[k].range.x0 >= x1)
                break;

            res += norm2(s[k], _f ? x0 : s[k].range.x0,
                               k == s.size()-1 ? x1 : s[k].range.x1);

            _f = false;
        }

	    return sqrt(res);
	}
	else
        throw NotImplemented();
}

/**
struct TPeak
{
	Point p;			/// peak tip
	double left;		/// f''(left) == 0
	double right;		/// f''(right) == 0

	TPeak(const Point& p, double left, double right): p(p), left(left), right(right)
	{
	}

	operator Point() const
	{
		return p;
	}
};

std::vector<TPeak> peaks(const spline& s)
{
	const double NaN = nan("");
	std::vector<double> z;

	std::vector<TPeak> res;

	for (unsigned k = 0; k < s.size(); k++)
    {
        auto q = s[k].extrema();

		bool _max = false;

        for (const auto& p: q)
        {
            if(p.second < 0)		/// maximum
			{
				_max = true;
				res.emplace_back(p.first, NaN, NaN);

				break;				/// a cubic can have only one maximum
			}
        }

		if(_max)
		{
			int idx = k;

			while(idx >= 0)
			{
				z = D2(s[idx]).zeroes();

				if(z.size())
				{
					if(z[0] < res.back().p.x)
					{
						res.back().left = z[0];
						break;
					}
				}

				idx--;
			}

			idx = k;

			while(idx < s.size())
			{
				z = D2(s[idx]).zeroes();

				if(z.size())
				{
					if(z[0] > res.back().p.x)
					{
						res.back().right = z[0];
						break;
					}
				}

				idx++;
			}
		}
    }

	return res;
}
*/

#endif

#endif
