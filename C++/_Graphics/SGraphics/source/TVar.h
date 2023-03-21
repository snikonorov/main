#ifndef TVAR_H_INCLUDED
#define TVAR_H_INCLUDED

#include <functional>

template<template<class> class Obj, class T>
T SubType(const Obj<T>& v);

#define $SubType(obj) decltype(SubType(obj))

///-------------------------------------------------------------

#define INSTANTIATE_ASSIGNMENT \
    template<class Q> \
    void operator+=(const Q& v) \
    { \
        auto& var = *this; \
        var = ($SubType(var)) var +  ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    void operator-=(const Q& v) \
    { \
        auto& var = *this; \
        var = ($SubType(var)) var - ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    void operator*=(const Q& v) \
    { \
        auto& var = *this; \
        var = (($SubType(var)) var) * (($SubType(var)) v); \
    } \
    \
    template<class Q> \
    void operator/=(const Q& v) \
    { \
        auto& var = *this; \
        var = (($SubType(var)) var) / (($SubType(var)) v); \
    }

///-------------------------------------------------------------

/*#define INSTANTIATE_ARIPHMETICS(Obj) \
    template<class Q> \
    Obj operator+(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var + ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    Obj operator-(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var - ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    Obj operator*(const Q& v) const \
    { \
        auto& var = *this; \
        return (($SubType(var)) var) * (($SubType(var)) v); \
    } \
    \
    template<class Q> \
    Obj operator/(const Q& v) const \
    { \
        auto& var = *this; \
        return (($SubType(var)) var) / (($SubType(var)) v); \
    }\
    \
    template<class Q> \
    friend Obj operator+(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) var + ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    friend Obj operator-(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) var - ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    friend Obj operator*(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return (($SubType(var)) var) * (($SubType(var)) v); \
    } \
    \
    template<class Q> \
    friend Obj operator/(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return (($SubType(var)) var) / (($SubType(var)) v); \
    }*/

#define INSTANTIATE_ARIPHMETICS(Obj) \
    template<class Q> \
    Q operator+(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var + ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    Q operator-(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var - ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    Q operator*(const Q& v) const \
    { \
        auto& var = *this; \
        return (($SubType(var)) var) * (($SubType(var)) v); \
    } \
    \
    template<class Q> \
    Q operator/(const Q& v) const \
    { \
        auto& var = *this; \
        return (($SubType(var)) var) / (($SubType(var)) v); \
    }\
    \
    template<class Q> \
    friend Q operator+(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) var + ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    friend Q operator-(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) v - ($SubType(var)) var; \
    } \
    \
    template<class Q> \
    friend Q operator*(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return (($SubType(var)) v) * (($SubType(var)) var); \
    } \
    \
    template<class Q> \
    friend Q operator/(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return (($SubType(var)) v) / (($SubType(var)) var); \
    }

///-------------------------------------------------------------

#define INSTANTIATE_COMPARATORS(Obj) \
    template<class Q> \
    bool operator==(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var == ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    bool operator!=(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var != ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    bool operator>(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var > ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    bool operator<(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var < ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    bool operator>=(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var >= ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    bool operator<=(const Q& v) const \
    { \
        auto& var = *this; \
        return ($SubType(var)) var <= ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    friend bool operator==(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) var == ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    friend bool operator!=(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) var != ($SubType(var)) v; \
    } \
    \
    template<class Q> \
    friend bool operator>(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) v > ($SubType(var)) var; \
    } \
    \
    template<class Q> \
    friend bool operator<(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) v < ($SubType(var)) var; \
    } \
    \
    template<class Q> \
    friend bool operator>=(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) v >= ($SubType(var)) var; \
    } \
    \
    template<class Q> \
    friend bool operator<=(const Q& v, const Obj& o) \
    { \
        auto& var = o; \
        return ($SubType(var)) v <= ($SubType(var)) var; \
    }

///-------------------------------------------------------------

template<class T>
struct IVar             /// basic variable interface
{
    //void (*SetValue)(IVar<T>*, const T&);
    //const T (*GetValue)(const IVar<T>*);

    std::function<void(IVar<T>*, const T&)> SetValue;
    std::function<const T(const IVar<T>*)> GetValue;

    ///------------------------------------------------------

    //template<class Q>
    void operator=(const T& v)
    {
        SetValue(this, v);
    }

	//template<class Q>
	void operator=(T& v)
	{
		SetValue(this, v);
	}

	//template<class Q>
	operator T() const
	{
		return GetValue(this);
	}

	virtual T& GetRef() = 0;

    virtual ~IVar()
    {
    }
};

template<class T>
class TVar : public IVar<T>
{
	protected:

    T val;

    public:

    static void DefSetValue(IVar<T>* obj, const T& v)
    {
        ((TVar<T>*)obj)->val = v;
    }

    static const T DefGetValue(const IVar<T>* obj)
    {
        return ((TVar<T>*)obj)->val;
    }

    ///----------------------------------------------------------

    TVar()
    {
		this->SetValue = DefSetValue;
		this->GetValue = DefGetValue;
    }

    template<class Q>
    explicit TVar(const Q& v)
    {
		this->SetValue = DefSetValue;
		this->GetValue = DefGetValue;

        (*this) = v;
    }

	/*template<class Q>
	explicit TVar(Q& v): val(v)
	{
		Reset();
	}*/

	using IVar<T>::operator=;
	using IVar<T>::operator T;

    INSTANTIATE_ASSIGNMENT
    INSTANTIATE_ARIPHMETICS(TVar<T>)
    INSTANTIATE_COMPARATORS(TVar<T>)

	bool operator!()
	{
		auto& var = *this;
		return !(T) var;
	}

	virtual T& GetRef() override
	{
		return val;
	}
};

#undef INSTANTIATE_ASSIGNMENT
#undef INSTANTIATE_ARIPHMETICS
#undef INSTANTIATE_COMPARATORS

#endif // TVAR_H_INCLUDED
