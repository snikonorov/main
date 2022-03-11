#ifndef FUNCTIONAL_UTILS_H_INCLUDED
#define FUNCTIONAL_UTILS_H_INCLUDED

#include <vector>
#include <algorithm>

/// General-----------------------------------------------------

/// if exists `idx`: !(value < data[idx]) && !(data[idx] < value)
///
///     => return `idx`
///		(~ data[idx] == value)
///
/// else:
///
///     => returns `idx`: value < data[idx] && (data[idx - 1] < value || idx == 0)
///     (~ `idx` is the index of insertion of `value`)
///
template<class TContainer, class T>
unsigned bin_search(const TContainer& data, const T& value)
{
	unsigned a = 0, b = data.size(), c = 0;

	while(a < b)
	{
		c = a + (b - a)/2;

		if(data[c] < value)
			a = c + 1;
		else if(value < data[c])
			b = c;
		else
			break;
	}

	if(a == b)
	{
		///data.insert(data.begin() + a, value);

		return a;
	}
	else
	{
		return c;
	}
}

template<class TContainer, class T, class TFunction>
unsigned bin_search(const TContainer& data, const T& value, const TFunction& f)
{
	unsigned a = 0, b = data.size(), c = 0;
	decltype(f(data[c])) e;

	while(a < b)
	{
		c = a + (b - a)/2;
		e = f(data[c]);

		if(e < value)
			a = c + 1;
		else if(value < e)
			b = c;
		else
			break;
	}

	if(a == b)
	{
		///data.insert(data.begin() + a, value);

		return a;
	}
	else
	{
		return c;
	}
}

///-------------------------------------------------------------

template<class TContainer>
struct $slice
{
	TContainer& v;
	unsigned idx;
	unsigned n;

	$slice(TContainer& v, unsigned idx, unsigned n): v(v), idx(idx), n(n)
	{
		n = std::min(n, unsigned(v.end() - v.begin()) - idx);
	}

	decltype(auto) begin()
	{
		return v.begin() + idx;
	}

	decltype(auto) end()
	{
		return v.begin() + idx  + n;
	}

	decltype(auto) begin() const
	{
		return v.begin() + idx;
	}

	decltype(auto) end() const
	{
		return v.begin() + idx  + n;
	}

	decltype(auto) operator[](unsigned k)
	{
		return v[idx + k];
	}

	decltype(auto) operator[](unsigned k) const
	{
		return v[idx + k];
	}

	decltype(auto) size() const
	{
		return n;
	}

	decltype(auto) front() const
	{
		return v[idx];
	}

	decltype(auto) back() const
	{
		return v[idx + n - 1];
	}

	/**operator TContainer() const
	{
		TContainer res(n);

		for(unsigned k = 0; k < n; k++)
		{
			res[k] = v[idx + k];
		}

		return res;
		///return TContainer(begin(), end());
	}*/
};

template<class TContainer>
decltype(auto) slice(TContainer& v, unsigned idx, unsigned n)
{
	return $slice<TContainer>(v, idx, n);
}

template<class TContainer>
inline void assign(TContainer& v, const $slice<TContainer>& s)
{
	TContainer res(s.n);

	for(unsigned k = 0; k < s.n; k++)
	{
		res[k] = s.v[s.idx + k];
	}

	v = res;

	if(&v == &s.v)
	{
		decltype(auto) q = const_cast<$slice<TContainer>*>(&s);

		q->idx = 0;
		q->n = v.size();
	}
}

///-------------------------------------------------------------

template<class T1, class T2, class T3>
decltype(auto) vbound(T1 v, T2 lower, T3 upper)
{
    return v < lower ? lower : v > upper ? upper : v;
}

decltype(auto) range(double a, double b, unsigned N)
{
    if(N < 2)
		N = 2;

	std::vector<double> x(N);

    for (unsigned k = 0; k < N; k++)
    {
        x[k] = a + (b - a)*k/(N-1);
    }

    return x;
}

unsigned unprime(unsigned N)
{
    double a = log2(N);

    if((1u << unsigned(ceil(a))) == N)
        return N;

    unsigned m = 1 << unsigned(ceil(0.5*a));

    return m*std::max(2u, unsigned(ceil(1.*N/m)));
}

inline decltype(auto) arange(double a, double b, double dx, bool avoid_prime = false)
{
	unsigned N = round((b - a)/dx) + 1;

	if(avoid_prime)
		N = unprime(N);

	dx = (b - a)/(N - 1);

	return std::make_pair(range(a, b, N), dx);
}

decltype(auto) range(int a, int b)
{
    char sgn = b > a ? 1 : -1;
    std::vector<int> x(sgn*(b - a) + 1);

    for (unsigned k = 0; k < x.size(); k++)
    {
        x[k] = a + k*sgn;
    }

    return x;
}

decltype(auto) range(unsigned N)
{
    return N == 0 ? std::vector<int>() : range(0, N-1);
}

template
<
	class F,
	class X,
	class Q = decltype(std::declval<F>()(std::declval<X>()[0])),
	class = typename std::enable_if<!std::is_same<Q, void>::value>::type
>
decltype(auto) map(const F& f, const X& x)
{
    std::vector<decltype(f(x[0]))> res(x.size());

    for (unsigned k = 0; k < x.size(); k++)
    {
        res[k] = f(x[k]);
    }

    return res;
}

template
<
	class F,
	class X,
	class Q = decltype(std::declval<F>()(std::declval<X>()[0])),
	class = typename std::enable_if<std::is_same<Q, void>::value>::type
>
void map(const F& f, const X& x)
{
    for (unsigned k = 0; k < x.size(); k++)
    {
        f(x[k]);
    }
}

template<class T>
decltype(auto) sum(const std::vector<T>& d)
{
    T res = 0;

    for (unsigned k = 0; k < d.size(); k++)
    {
        res += d[k];
    }

    return res;
}

template<class TContainer>
decltype(auto) average(const TContainer& d)
{
    if (d.size())
    {
        auto res = d[0];

        for (unsigned k = 1; k < d.size(); k++)
        {
            res += d[k];
        }

        return res/d.size();
    }
    else
        return typename std::remove_const<typename std::remove_reference<decltype(d[0])>::type>::type();
}

template
<
	class TContainer,
	class T = typename std::remove_const<typename std::remove_reference<decltype(std::declval<TContainer>()[0])>::type>::type
>
decltype(auto) std_sqr(const TContainer& d)
{
    if (d.size() > 1)
    {
        T res = T();
		auto d0 = average(d);

        for (unsigned k = 0; k < d.size(); k++)
        {
            res += (d[k] - d0)*(d[k] - d0);
        }

        return res/(d.size() - 1);
    }
    else
        return T();
}

template<class T>
decltype(auto) norm(const std::vector<T>& d)
{
    if (d.size())
    {
        T res = std::abs(d[0]);

        for (unsigned k = 1; k < d.size(); k++)
        {
			if (std::abs(d[k]) > res)
                res = std::abs(d[k]);
        }

        return res;
    }
    else
        return T(0);
}

template<class TContainer>
decltype(auto) max(const TContainer& d)
{
    if (d.size())
    {
        auto res = d[0];

        for (unsigned k = 1; k < d.size(); k++)
        {
            if (d[k] > res)
                res = d[k];
        }

        return res;
    }
    else
        return typename std::remove_reference<decltype(d[0])>::type();
}

template<class T>
decltype(auto) min(const std::vector<T>& d)
{
    if (d.size())
    {
        auto res = d[0];

        for (unsigned k = 1; k < d.size(); k++)
        {
            if (d[k] < res)
                res = d[k];
        }

        return res;
    }
    else
        return T();
}

template<class TContainer, class TFunction>
decltype(auto) filter(const TContainer& d, const TFunction& predicate)
{
	TContainer res(d.size());

	unsigned idx = 0;
	for(const auto& a: d)
	{
		if(predicate(a))
		{
			res[idx++] = a;
		}
	}

	res.resize(idx);

	return res;
}

template
<
	class TContainer,
	class T = decltype(+std::declval<TContainer>()[0]),
	class = decltype(std::declval<T>() > std::declval<T>())
>
decltype(auto) bounds(const TContainer& d)
{
    if (d.size())
    {
        auto mn = d[0];
		auto mx = d[0];

        for (unsigned k = 1; k < d.size(); k++)
        {
            if (d[k] > mx)
                mx = d[k];

			if (d[k] < mn)
                mn = d[k];
        }

        return std::make_pair(mn, mx);
    }
    else
        return std::make_pair(T(), T());
}

template
<
	class TContainer,
	class Q = decltype(std::declval<TContainer>()[0]),
	class T = decltype(std::declval<Q>().x),
	class T1 = decltype(std::declval<Q>().y)
>
decltype(auto) bounds(const TContainer& d)
{
    if (d.size())
    {
        auto mn = d[0].y;
		auto mx = d[0].y;

        for (unsigned k = 1; k < d.size(); k++)
        {
            if (d[k].y > mx)
                mx = d[k].y;

			if (d[k].y < mn)
                mn = d[k].y;
        }

        return std::make_pair(mn, mx);
    }
    else
        return std::make_pair(T(), T());
}

template<class TContainer>
inline double max_amplitude(const TContainer& d)
{
	auto q = bounds(d);

	return q.second - q.first;
}

template<class T>
decltype(auto) median(std::vector<T> d)
{
    if (d.size())
    {
        std::sort(d.begin(), d.end());

        return d.size() % 2 == 1 ? d[d.size()/2] : (d[d.size()/2] + d[d.size()/2 - 1])/2;
    }
    else
        return T();
}

template<class T>
std::vector<T> diff(const std::vector<T>& d)
{
    if (d.size() > 1)
    {
        std::vector<T> res(d.size()-1);

        for(unsigned k = 0; k < d.size()-1; k++)
        {
            res[k] = d[k+1] - d[k];
        }

        return res;
    }
    else
        return {};
}

template<class T>
std::vector<T> exclude_range(const std::vector<T>& d, double lower_threshold, double upper_threshold)
{
    std::vector<T> res(d.size());
    unsigned n = 0;

    for(unsigned k = 0; k < d.size(); k++)
    {
        if(d[k] > lower_threshold || d[k] < upper_threshold)
            res[n++] = d[k];
    }

    res.resize(n);

    return res;
}

template<class T>
std::vector<T> middle(const std::vector<T>& d)
{
    if (d.size() > 1)
    {
        std::vector<T> res(d.size()-1);

        for(unsigned k = 0; k < d.size()-1; k++)
        {
            res[k] = (d[k+1] + d[k])/2;
        }

        return res;
    }
    else
        return {};
}

/// Point functions---------------------------------------------
///
#ifdef VECTOR3D_H_INCLUDED

decltype(auto) slice(std::vector<Point>& v, double x1, double x2)
{
	unsigned idx1 = bin_search(v, x1, [](const Point& p){ return p.x; });
	unsigned idx2 = bin_search(v, x2, [](const Point& p){ return p.x; });

	return $slice<std::vector<Point>>(v, idx1, idx2 - idx1 + 1);
}

template<>
std::vector<Point> diff(const std::vector<Point>& d)
{
    if (d.size() > 1)
    {
        std::vector<Point> res(d.size()-1);

        for(unsigned k = 0; k < d.size()-1; k++)
        {
            res[k].y = d[k+1].y - d[k].y;
            res[k].x = (d[k+1].x + d[k].x)/2;
        }

        return res;
    }
    else
        return {};
}

template<>
std::vector<Point> exclude_range(const std::vector<Point>& d, double lower_threshold, double upper_threshold)
{
    std::vector<Point> res(d.size());
    unsigned n = 0;

    for(unsigned k = 0; k < d.size(); k++)
    {
        if(d[k].y > lower_threshold || d[k].y < upper_threshold)
            res[n++] = d[k];
    }

    res.resize(n);

    return res;
}

template<class TContainer1, class TContainer2>
std::vector<Point> zip(const TContainer1& a, const TContainer2& b)
{
    if (a.size() == b.size())
    {
        std::vector<Point> res(a.size());

        for (unsigned k = 0; k < a.size(); k++)
        {
            res[k].x = a[k];
            res[k].y = b[k];
        }

        return res;
    }
    else
        return {};
}

template<class F, class X>
decltype(auto) zipmap(const F& f, const X& x)
{
    std::vector<Point> res(x.size());

    for (unsigned k = 0; k < res.size(); k++)
    {
        res[k].x = x[k];
		res[k].y = f(x[k]);
    }

    return res;
}

decltype(auto) average_norm(const std::vector<Point>& d)
{
    if (d.size())
    {
        double res = 0;

        for (unsigned k = 0; k < d.size(); k++)
        {
            res += std::abs(d[k].y);
        }

        return res/d.size();
    }
    else
        return 0.0;
}

decltype(auto) average_norm(const std::vector<double>& d)
{
    if (d.size())
    {
        double res = 0;

        for (unsigned k = 0; k < d.size(); k++)
        {
            res += std::abs(d[k]);
        }

        return res/d.size();
    }
    else
        return 0.0;
}

decltype(auto) average_amplitude(const std::vector<Point>& d)
{
    if (d.size())
    {
        unsigned n = 2;
		double res = std::abs(d[0].y) + std::abs(d.back().y);

        for (unsigned k = 1; k < d.size()-1; k++)
        {
            if((d[k].y - d[k-1].y)*(d[k].y - d[k+1].y) >= 0)
			{
				res += std::abs(d[k].y);
				n++;
			}
        }

        return res/n;
    }
    else
        return 0.0;
}

decltype(auto) bounds(const std::vector<Point>& d)
{
    if (d.size())
    {
        double mn = d[0].y;
		double mx = d[0].y;

        for (unsigned k = 1; k < d.size(); k++)
        {
            if (d[k].y > mx)
                mx = d[k].y;

			if (d[k].y < mn)
                mn = d[k].y;
        }

        return std::make_pair(mn, mx);
    }
    else
        return std::make_pair(double(), double());
}

template<class TContainer>
decltype(auto) bounding_box(const TContainer& d)
{
    if (d.size())
    {
        auto mn = d[0];
		auto mx = d[0];

        for (unsigned k = 1; k < d.size(); k++)
        {
            if (d[k].x > mx.x)
                mx.x = d[k].x;

            if (d[k].y > mx.y)
                mx.y = d[k].y;

			if (d[k].x < mn.x)
                mn.x = d[k].x;

			if (d[k].y < mn.y)
                mn.y = d[k].y;
        }

        return std::make_pair(mn, mx);
    }
    else
        return std::make_pair(decltype(1*d[0])(), decltype(1*d[0])());
}

#endif

/// spline functions--------------------------------------------
///
#ifdef TSPLINE_H_INCLUDED

template<class TPoly>
TSpline<TPoly> slice(const TSpline<TPoly>& s, double x0, double x1)
{
	if(x0 > x1)
		std::swap(x0, x1);

	TSpline<TPoly> res;

	res.style = s.style;

	if(s.size() && x0 < s.x1())
	{
		int idx0 = x0 < s.x0() ? 0 : s.find(x0);
		int idx1 = x1 < s.x0() ? 1 : x1 > s.x1() ? s.size() : s.find(x1) + 1;

		res.assign(s.begin() + idx0, s.begin() + idx1);

		res.front() = res.front().split(x0).second;
		res.back() = res.back().split(x1).first;
	}

	return res;
}

inline decltype(auto) bounds(const spline& s)
{
	double mn = 0, mx = 0;

	for(unsigned k = 0; k < s.size(); k++)
	{
		auto b = s[k].bounds();

		if(k == 0)
		{
			mn = b.first.y;
			mx = b.second.y;
		}
		else
		{
			if(b.first.y < mn)
				mn = b.first.y;

			if(b.second.y > mx)
				mx = b.second.y;
		}
	}

	return std::make_pair(mn, mx);
}

#endif

#endif // FUNCTIONAL_UTILS_H_INCLUDED
