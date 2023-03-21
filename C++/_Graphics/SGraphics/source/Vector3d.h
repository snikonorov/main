#ifndef VECTOR3D_H_INCLUDED
#define VECTOR3D_H_INCLUDED

template<typename T>
struct Vector3d
{
    T x,y,z;

    Vector3d(): x(0), y(0), z(0)
    {
    }

    explicit Vector3d(const T& val): x(val), y(val), z(val)
    {
    }

    Vector3d(const T& X, const T& Y): x(X), y(Y), z(0)
    {
    }

	#ifdef _GLIBCXX_COMPLEX

	template<class Q>
	Vector3d(const std::complex<Q>& a): x(a.real()), y(a.imag()), z(0)
    {
    }

	#endif

    Vector3d(const T& X, const T& Y, const T& Z): x(X), y(Y), z(Z)
    {
    }

    Vector3d(const Vector3d& v): x(v.x), y(v.y), z(v.z)
    {
    }

    template<class Q>
    Vector3d& operator=(const Vector3d<Q>& v)
    {
        x = v.x;
        y = v.y;
        z = v.z;

        return *this;
    }

    /*Vector3d& operator=(const T& a)
    {
        x = a;
        y = a;
        z = a;

        return *this;
    }*/

    Vector3d& Normalize()
    {
        static double r;

        r=(x*x+y*y+z*z);

        if(r!=0 && r!=1)
        {
            r=sqrt(r);

            x/=r;
            y/=r;
            z/=r;
        }

        return *this;
    }

	/*T MaxAbs() const
	{
		return max3(std::abs(x), std::abs(y), std::abs(z));
	}*/

	Vector3d normal() const
	{
		static double r;

		r = (x*x + y * y + z * z);

		if (r != 0 && r != 1)
		{
			r = sqrt(r);

			return {x/r, y/r, z/r};
		}

		return {x, y, z};
	}

	Vector3d ortho() const
	{
		return {y, -x, 0};
	}

    T Sqr() const
    {
        return x*x + y*y + z*z;
    }

    T Sqr2() const
    {
        return x*x + y*y;
    }

    decltype(auto) Norm() const
    {
        return sqrt(x*x + y*y + z*z);
    }

    decltype(auto) Norm2() const
    {
        return sqrt(x*x + y*y);
    }

    template<class Q>
    auto operator+(const Vector3d<Q>& a) const -> Vector3d<decltype(x + a.x)>
    {
        return Vector3d<decltype(x + a.x)>(x+a.x, y+a.y, z+a.z);
    }

    /*template<class Q>
    auto operator+(const Q& a) const -> Vector3d<decltype((T)a + x)>
    {
        return Vector3d<decltype((T)a + x)>(x+a, y+a, z+a);
    }*/

    template<class Q>
    void operator+=(const Vector3d<Q>& a)
    {
        x += a.x;
        y += a.y;
        z += a.z;
    }

    /*template<class Q>
    void operator+=(const Q& a)
    {
        x += a;
        y += a;
        z += a;
    }*/

    template<class Q>
    auto operator-(const Vector3d<Q>& a) const -> Vector3d<decltype(x - a.x)>
    {
        return Vector3d<decltype(x - a.x)>(x-a.x, y-a.y, z-a.z);
    }

    /*template<class Q>
    auto operator-(const Q& a) const -> Vector3d<decltype(x - (T)a)>
    {
        return Vector3d<decltype(x - (T)a)>(x-a, y-a, z-a);
    }*/

    template<class Q>
    void operator-=(const Vector3d<Q>& a)
    {
        x -= a.x;
        y -= a.y;
        z -= a.z;
    }

    /*template<class Q>
    void operator-=(const Q& a)
    {
        x -= a;
        y -= a;
        z -= a;
    }*/

    Vector3d operator-() const
    {
        return Vector3d(-x, -y, -z);
    }

    template<class Q>
    auto operator*(const Vector3d<Q>& a) const -> Vector3d<decltype(a.x*x)>
    {
        return Vector3d<decltype(a.x * x)>(x*a.x, y*a.y, z*a.z);
    }

    template<class Q>
    auto operator*(const Q& k) const -> Vector3d<decltype((T)k * x)>
    {
        return Vector3d<decltype((T)k * x)>(x*k, y*k, z*k);
    }

    template<class Q>
    void operator*=(const Q& k)
    {
        x *= k;
        y *= k;
        z *= k;
    }

    template<class Q>
    void operator*=(const Vector3d<Q>& k)
    {
        x *= k.x;
        y *= k.y;
        z *= k.z;
    }

    template<class Q>
    auto operator&(const Vector3d<Q>& v) const -> decltype(v.x*x)  ///scalar
    {
        return x*v.x + y*v.y + z*v.z;
    }

    template<class Q>
    auto operator|(const Vector3d<Q>& v) const -> decltype(v.y*x)  ///z-component of a vector product
    {
        return x*v.y - y*v.x;
    }

    template<class Q>
    auto operator^(const Vector3d<Q>& v) const -> Vector3d<decltype(v.y*x)>  ///vector
    {
        return Vector3d<decltype(v.y*x)>(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
    }

    template<class Q>
    auto operator/(const Q& k) const -> Vector3d<decltype(x/k)>
    {
        return Vector3d<decltype(x/k)>(x/k, y/k, z/k);
    }

    template<class Q>
    auto operator/(const Vector3d<Q>& v) const -> Vector3d<decltype(x/v.x)>
    {
        return Vector3d<decltype(x/v.x)>(x/v.x, y/v.y, z/v.z);
    }

    template<class Q>
    void operator/=(const Q& k)
    {
        x /= k;
        y /= k;
        z /= k;
    }

    ///----------------------------------------------------------------------

    template<class Q>
    bool operator==(const Vector3d<Q>& v)
    {
        return x == v.x && y == v.y && z == v.z;
    }

    template<class Q>
    bool operator!=(const Vector3d<Q>& v)
    {
        return x != v.x || y != v.y || z != v.z;
    }

    template<class Q>
    bool operator<(const Vector3d<Q>& v)
    {
        return x < v.x;
    }

    ///----------------------------------------------------------------------

    template<class Q>
    operator Vector3d<Q>() const
    {
        return Vector3d<Q>((Q) x, (Q) y, (Q) z);
    }

    #ifdef _GLIBCXX_IOSTREAM

    friend std::ostream& operator<<(std::ostream& out, const Vector3d<T>& v)
    {
        return out << '(' << v.x << ", " << v.y << ", " << v.z << ")";
    }

    #endif
};

///----------------------------------------------------------------------------------------------------

template<class Q, class T>
auto operator+(const Q& a, const Vector3d<T>& v) -> Vector3d<decltype((T)a + v.x)>
{
    return Vector3d<decltype((T)a + v.x)>(a + v.x, a + v.y, a + v.z);
}

template<class Q, class T>
auto operator-(const Q& a, const Vector3d<T>& v) -> Vector3d<decltype((T)a - v.x)>
{
    return Vector3d<decltype((T)a - v.x)>(a - v.x, a - v.y, a - v.z);
}

template<class T, class Q>
auto operator*(const Q& k, const Vector3d<T>& a) -> Vector3d<decltype((T)k * a.x)>
{
    return Vector3d<decltype((T)k * a.x)>(a.x*k, a.y*k, a.z*k);
}

///----------------------------------------------------------------------------------------------------

/*template<class T, class R, class A>
Vector3d<R> operator^(R (*f)(A), const Vector3d<T>& v)
{
    return Vector3d<R>(f(v.x), f(v.y), f(v.z));
}

///more general
template<class T, class Q>
auto operator^(Q f, const Vector3d<T>& v) -> Vector3d<decltype(f(v.x))>
{
    return Vector3d<decltype(f(v.x))>(f(v.x), f(v.y), f(v.z));
}*/

///----------------------------------------------------------------------------------------------------

typedef Vector3d<double> Point;
typedef Vector3d<int> Coord;

#endif // VECTOR3D_H_INCLUDED
