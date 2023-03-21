#ifndef TO_STRING_H_INCLUDED
#define TO_STRING_H_INCLUDED

#include <string>

template<class T>
inline decltype(auto) to_string(const T& a)
{
    return std::to_string(a);
}

inline decltype(auto) to_string(double a)
{
    static char s[32];

	sprintf(s, "%.15g", a);

	return std::string(s);
}

#ifdef _GLIBCXX_UTILITY

template<class T1, class T2>
inline decltype(auto) to_string(const std::pair<T1, T2>& p)
{
	std::string res = "(";

	res += to_string(p.first);
	res += ", ";
	res += p.second;

	res += ")";

	return res;
}

#endif

template<class TContainer, class = decltype(std::declval<TContainer>().begin(), std::declval<TContainer>().end())>
inline decltype(auto) to_string(const TContainer& v)
{
    std::string res = "{";

    bool _flag = false;

    for(const auto& q: v)
    {
        if (_flag)
            res += ", ";
        else
            _flag = true;

        res += to_string(q);
    }

    res += "}";

	return res;
}

#endif // TO_STRING_H_INCLUDED
