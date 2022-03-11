#ifndef TRECT_H_INCLUDED
#define TRECT_H_INCLUDED

#include "Vector3d.h"

namespace _local
{
    template<class T1, class T2>
    auto min(const T1& a, const T2& b) -> decltype(a + b)
    {
        typedef decltype(a + b) T;
        return a<b ? (T)a : (T)b;
    }

    template<class T1, class T2>
    auto max(const T1& a, const T2& b) -> decltype(a + b)
    {
        typedef decltype(a + b) T;
        return a > b ? (T)a : (T)b;
    }
}

template<class T1, class T2, class T3>
decltype(auto) max3(const T1& a, const T2& b, const T3& c)
{
	return a>b ? (a>c ? a : c) : (b>c ? b : c);
}

template<class T1, class T2, class T3>
decltype(auto) min3(const T1& a, const T2& b, const T3& c)
{
	return a < b ? (a < c ? a : c) : (b < c ? b : c);
}

///---------------------------------------------------------------

template<class Type>
struct TRect
{
    typedef Vector3d<Type> TPoint;

    TPoint A, B;

    TRect()
    {
    }

    template<class T1, class T2, class T3, class T4>
    TRect(const T1& x0, const T2& y0, const T3& width, const T4& height)
    {
        A = TPoint(x0, y0);
        B = TPoint(x0 + width, y0 + height);
    }

    template<class T1, class T2>
    TRect(const T1& width, const T2& height)
    {
        A = TPoint(0, 0);
        B = TPoint(width, height);
    }

    template<class Q>
    TRect(const Vector3d<Q>& a)
    {
        A = a;
        B = a;
    }

    template<class Q1, class Q2>
    TRect(const Vector3d<Q1>& a, const Vector3d<Q2>& b, bool normalize = true)
    {
        if(normalize)
        {
            A.x = _local::min(a.x, b.x);
            A.y = _local::min(a.y, b.y);

            B.x = _local::max(a.x, b.x);
            B.y = _local::max(a.y, b.y);
        }
        else
        {
            A = a;
            B = b;
        }
    }

    template<class Q1, class Q2, class Q3>
    TRect(const Vector3d<Q1>& a, const Vector3d<Q2>& b, const Vector3d<Q3>& c)
    {
        A.x = min3(a.x, b.x, c.x);
        A.y = min3(a.y, b.y, c.y);

        B.x = max3(a.x, b.x, c.x);
        B.y = max3(a.y, b.y, c.y);
    }

    template<class Q>
    TRect& operator=(const TRect<Q>& r)
    {
        A = r.A;
        B = r.B;

        return *this;
    }

	template<class Q>
    void operator=(const Vector3d<Q>& v)
    {
        A = v;
        B = v;
    }

    template<class Q>
    TRect operator+(const Vector3d<Q>& v) const
    {
        return TRect(A + v, B + v);
    }

    template<class Q>
    TRect& operator+=(const Vector3d<Q>& v)
    {
        A += v;
        B += v;

        return *this;
    }

    template<class Q>
    TRect operator-(const Vector3d<Q>& v) const
    {
        return TRect(A - v, B - v);
    }

    template<class Q>
    TRect& operator-=(const Vector3d<Q>& v)
    {
        A -= v;
        B -= v;

        return *this;
    }

    TRect& operator*=(double a)
    {
        A *= a;
        B *= a;

        return *this;
    }

    TRect& Normalize()
    {
        TPoint _A = A;

        A.x = std::min(A.x, B.x);
        A.y = std::min(A.y, B.y);

        B.x = std::max(_A.x, B.x);
        B.y = std::max(_A.y, B.y);

        return *this;
    }

    template<class Q>
    bool Intersects(const TRect<Q>& r, bool normalize = false) const
    {
        return normalize ?
               (std::max(A.x, B.x) >= std::min(r.A.x, r.B.x) &&
                std::max(A.y, B.y) >= std::min(r.A.y, r.B.y) &&
                std::min(A.x, B.x) <= std::max(r.A.x, r.B.x) &&
                std::min(A.y, B.y) <= std::max(r.A.y, r.B.y))

               : (B.x >= r.A.x && B.y >= r.A.y && A.x <= r.B.x && A.y <= r.B.y);
    }

    ///for normalized Rect only!
    template<class Q>
    void Merge(const Vector3d<Q>& p)
    {
        if(p.x < A.x) A.x = p.x;
        else if(p.x > B.x) B.x = p.x;

        if(p.y < A.y) A.y = p.y;
        else if(p.y > B.y) B.y = p.y;
    }

    template<class Q>
    void Merge(const TRect<Q>& r)
    {
        Merge(r.A);
        Merge(r.B);
    }

	/// clipping (both rects are assumed to be normalized)
	template<class Q>
    TRect operator&(const TRect<Q>& r) const
    {
        TPoint A;
		TPoint B;

		A.x = _local::max(this->A.x, r.A.x);
		A.y = _local::max(this->A.y, r.A.y);

		B.x = _local::min(this->B.x, r.B.x);
		B.y = _local::min(this->B.y, r.B.y);

		return (A.x > B.x || A.y > B.y) ? TRect() : TRect(A, B, false);
    }

    TPoint Dim() const
    {
        return B-A;
    }

	template<class Q>
	void Dim(const Q& d)
    {
        B = A + d;
    }

    TPoint Center() const
    {
        return (A+B)/2;
    }

    Type Width() const
    {
        //return (B - A).x;
		return B.x - A.x;
    }

    Type Height() const
    {
        //return (B - A).y;
		return B.y - A.y;
    }

    template<class Q>
    void Width(const Q& w)
    {
        B.x = A.x + w;
    }

    template<class Q>
    void Height(const Q& h)
    {
        B.y = A.y + h;
    }

    template<class Q>
    bool contains(const Vector3d<Q>& p) const
    {
        return p.x >= A.x && p.y >= A.y && p.x <= B.x && p.y <= B.y;
    }

    template<class T1, class T2>
    bool contains(const T1& x, const T2& y) const
    {
        return x >= A.x && y >= A.y && x <= B.x && y <= B.y;
    }

    template<class Q>
    operator TRect<Q>() const
    {
        return TRect<Q>(A, B, false);
    }
};

typedef TRect<double> Rect;
typedef TRect<int> iRect;

#endif // TRECT_H_INCLUDED
