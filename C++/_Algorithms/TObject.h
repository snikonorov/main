#ifndef TOBJECT_H_INCLUDED
#define TOBJECT_H_INCLUDED

#include <utility>
#include <string.h>

#include "to_string.h"

#include <vector>
#include <map>

template<class T>
class InvalidCast
{
};

struct ScalarTypeIndexing
{
};

enum class Type
{
    Void = 0,
    Bool,
    Char,
    Integer,
    Number,
    Array           /// ~ array/map/tuple
};

inline std::string to_string(const Type& type)
{
    switch(type)
    {
        case Type::Void:
            return "void";

        case Type::Bool:
            return "bool";

        case Type::Char:
            return "char";

        case Type::Integer:
            return "integer";

        case Type::Number:
            return "number";

        case Type::Array:
            return "array";

        default:
            return "";
    }
}

struct TObject;
inline std::string to_string(const TObject&);

struct TObject
{
    Type type;

    union
    {
        bool _bool;                                     /// Type::Bool
        unsigned _char;                                 /// Type::Char
        long long _integer;                             /// Type::Integer
        double _number;                                 /// Type::Number
    };

    /// [elements are sorted]
    std::vector<std::pair<TObject, TObject>> data;      /// Type::Array

    TObject(): type(Type::Void)
    {
    }

    TObject(bool a): type(Type::Bool), _bool(a)
    {
    }

    TObject(char a): type(Type::Char), _char(a)
    {
    }

    TObject(int a): type(Type::Integer), _integer(a)
    {
    }

    TObject(long a): type(Type::Integer), _integer(a)
    {
    }

    TObject(long long a): type(Type::Integer), _integer(a)
    {
    }

    TObject(unsigned a): type(Type::Integer), _integer(a)
    {
    }

    TObject(float a): type(Type::Number), _number(a)
    {
    }

    TObject(double a): type(Type::Number), _number(a)
    {
    }

    TObject(const std::string& s): type(Type::Array), data(s.size())
    {
        for(unsigned k = 0; k < data.size(); k++)
        {
            data[k].first = k;
            data[k].second = s[k];
        }
    }

    TObject(const char* s): type(Type::Array), data(strlen(s))
    {
        for(unsigned k = 0; k < data.size(); k++)
        {
            data[k].first = k;
            data[k].second = s[k];
        }
    }

    template<class T>
    TObject(const std::vector<T>& v): type(Type::Array), data(v.size())
    {
        for(unsigned k = 0; k < data.size(); k++)
        {
            data[k].first = k;
            data[k].second = v[k];
        }
    }

    template<class TKey, class TValue>
    TObject(const std::map<TKey, TValue>& m): type(Type::Array), data(m.size())
    {
        unsigned k = 0;
        for(const auto& a: m)
        {
            data[k].first = a.first;
            data[k].second = a.second;

            k++;
        }
    }

    template<class T>
    TObject(const std::initializer_list<T>& v): type(Type::Array), data(v.size())
    {
        unsigned k = 0;
        for(const auto& a: v)
        {
            data[k].first = k;
            data[k].second = a;

            k++;
        }
    }

    bool operator<(const TObject& obj) const
    {
        if(type == obj.type)
        {
            switch(type)
            {
                case Type::Void:
                    return false;
                    break;

                case Type::Bool:
                    return _bool < obj._bool;
                    break;

                case Type::Char:
                    return _char < obj._char;
                    break;

                case Type::Integer:
                    return _integer < obj._integer;
                    break;

                case Type::Number:
                    return _number < obj._number;
                    break;

                default:
                    if(data.size() == obj.data.size())
                    {
                        for(unsigned k = 0; k < data.size(); k++)
                        {
                            if(data[k].second < obj.data[k].second)
                                return true;
                        }

                        return false;
                    }
                    else
                    {
                        return data.size() < obj.data.size();
                    }
                    break;
            }
        }
        else
        {
            return type < obj.type;
        }
    }

    bool operator==(const TObject& obj) const
    {
        if(type == obj.type)
        {
            switch(type)
            {
                case Type::Void:
                    return true;
                    break;

                case Type::Bool:
                    return _bool == obj._bool;
                    break;

                case Type::Char:
                    return _char == obj._char;
                    break;

                case Type::Integer:
                    return _integer == obj._integer;
                    break;

                case Type::Number:
                    return _number == obj._number;
                    break;

                default:
                    if(data.size() == obj.data.size())
                    {
                        for(unsigned k = 0; k < data.size(); k++)
                        {
                            if(!(data[k].second == obj.data[k].second))
                                return false;
                        }

                        return true;
                    }
                    else
                    {
                        return false;
                    }
                    break;
            }
        }

        return false;
    }

    template<unsigned N>
    TObject& operator[](const char (&key)[N])
    {
        return (*this)[TObject(key)];
    }

    TObject& operator[](const TObject& key)
    {
        /**struct proxy
        {
            TObject* obj;
            TObject key;
            unsigned idx;
            bool _insert;

            proxy(): obj(0)
            {
            }

            proxy(TObject* obj, const TObject& key, unsigned idx, bool _insert = false)
                  : obj(obj), key(key), idx(idx), _insert(_insert)
            {
            }

            void operator=(const proxy& P)
            {
                *this = (TObject)P;
            }

            ///decltype(auto) operator=(proxy&&) = delete;

            operator TObject() const
            {
                if(!obj || _insert)
                    return TObject();               /// null
                else
                    return obj->data[idx].second;
            }

            void operator=(const TObject& a)
            {
                if(obj)
                {
                    if(obj->type == Type::Void)
                        obj->type = Type::Array;

                    if(&a == obj)
                    {
                        TObject a0 = a;

                        if(_insert)
                        {
                            obj->data.insert(obj->data.begin() + idx, std::make_pair(key, TObject()));
                            _insert = false;
                        }

                        obj->data[idx].second = a0;
                    }
                    else
                    {
                        if(_insert)
                        {
                            obj->data.insert(obj->data.begin() + idx, std::make_pair(key, TObject()));
                            _insert = false;
                        }

                        obj->data[idx].second = a;
                    }
                }
            }
        };*/

        ///----------------------------------------------------------------

        if(type == Type::Void)
        {
            type = Type::Array;

            data.insert(data.begin(), std::make_pair(key, TObject()));
            return data[0].second;

            ///return proxy(this, key, 0, true);
        }
        else if(type == Type::Array)
        {
            /// binary search
            unsigned a = 0, b = data.size(), c = 0;

            while(a < b)
            {
                c = a + (b - a)/2;

                if(data[c].first < key)
                    a = c + 1;
                else if(key < data[c].first)
                    b = c;
                else
                    break;
            }

            if(a == b)
            {
                data.insert(data.begin() + a, std::make_pair(key, TObject()));
                return data[a].second;

                ///return proxy(this, key, a, true);
            }
            else
            {
                return data[c].second;
                ///return proxy(this, key, c, false);
            }
        }

        throw ScalarTypeIndexing();
    }

	template<unsigned N>
    const TObject& operator[](const char (&key)[N]) const
    {
        return (*this)[TObject(key)];
    }

	const TObject& operator[](const TObject& key) const
    {
        static const TObject null;

        ///----------------------------------------------------------------

        if(type == Type::Void)
        {
            return null;
        }
        else if(type == Type::Array)
        {
            /// binary search
            unsigned a = 0, b = data.size(), c = 0;

            while(a < b)
            {
                c = a + (b - a)/2;

                if(data[c].first < key)
                    a = c + 1;
                else if(key < data[c].first)
                    b = c;
                else
                    break;
            }

            if(a == b)
            {
                return null;
            }
            else
            {
                return data[c].second;
            }
        }

        throw ScalarTypeIndexing();
    }

    operator bool() const
    {
        switch(type)
        {
            case Type::Bool:
                return _bool;

            case Type::Char:
                return _char;

            case Type::Integer:
                return _integer;

            case Type::Number:
                return _number;

            case Type::Array:
                return !data.empty();

            default:
                break;
        }

        throw InvalidCast<bool>();
    }

    operator char() const
    {
        switch(type)
        {
            case Type::Char:
                return _char;

            case Type::Array:
                if(data.size() == 1)
                    return data[0].second;

                break;

            default:
                break;
        }

        throw InvalidCast<char>();
    }

    operator int() const
    {
        switch(type)
        {
            case Type::Char:
                return _char;

            case Type::Integer:
                return _integer;

            case Type::Number:
                return _number;

            case Type::Array:
                if(data.size() == 1)
                    return data[0].second;

                break;

            default:
                break;
        }

        throw InvalidCast<int>();
    }

    operator long long() const
    {
        switch(type)
        {
            case Type::Char:
                return _char;

            case Type::Integer:
                return _integer;

            case Type::Number:
                return _number;

            case Type::Array:
                if(data.size() == 1)
                    return data[0].second;

                break;

            default:
                break;
        }

        throw InvalidCast<long long>();
    }

    operator unsigned() const
    {
        switch(type)
        {
            case Type::Char:
                return _char;

            case Type::Integer:
                return _integer;

            case Type::Number:
                return _number;

            case Type::Array:
                if(data.size() == 1)
                    return data[0].second;

            default:
                break;
        }

        throw InvalidCast<unsigned>();
    }

    operator float() const
    {
        switch(type)
        {
            case Type::Char:
                return _char;

            case Type::Integer:
                return _integer;

            case Type::Number:
                return _number;

            case Type::Array:
                if(data.size() == 1)
                    return data[0].second;

                break;

            default:
                break;
        }

        throw InvalidCast<float>();
    }

    operator double() const
    {
        switch(type)
        {
            case Type::Char:
                return _char;

            case Type::Integer:
                return _integer;

            case Type::Number:
                return _number;

            case Type::Array:
                if(data.size() == 1)
                    return data[0].second;

                break;

            default:
                break;
        }

        throw InvalidCast<double>();
    }

    operator std::string() const
    {
        auto to_printable_ASCII = [](unsigned ch) -> std::string
        {
            if(ch >= 32 && ch <= 126 && ch != '\"')
            {
                return std::string(1, (char)ch);
            }
            else
            {
                if(ch == '\n')
                {
                    return "\\n";
                }
                else if(ch == '\t')
                {
                    return "\\t";
                }
                else if(ch == '\a')
                {
                    return "\\a";
                }
                else if(ch == '\v')
                {
                    return "\\v";
                }
                else if(ch == '\r')
                {
                    return "\\r";
                }
                else if(ch == 0)
                {
                    return "\\0";
                }
                else if(ch == '\b')
                {
                    return "\\b";
                }
                else if(ch == 27)       /// not standard for C
                {
                    return "\\e";
                }
                else if(ch == '\f')
                {
                    return "\\f";
                }
                else if(ch == '\\')
                {
                    return "\\\\";
                }
                else if(ch == '\"')
                {
                    return "\\\"";
                }
                else                        /// using the '\udd...d;' format
                {
                    return "\\u" + to_string(ch) + ';';
                }
            }
        };

        ///----------------------------------------------

        std::string s;

        switch(type)
        {
            case Type::Void:
                s = "null";
                break;

            case Type::Bool:
                s = _bool ? "true" : "false";
                break;

            case Type::Char:
                s += "'";
                s += to_printable_ASCII(_char);
                s += "'";
                break;

            case Type::Integer:
                s = to_string(_integer);
                break;

            case Type::Number:
                s = to_string(_number);
                break;

            case Type::Array:
                {
                    bool _print_keys = false;
                    char _brace = data.size() ? '"' : '{';

                    try
                    {
                        for(unsigned k = 0; k < data.size(); k++)
                        {
                            if((unsigned)data[k].first != k)
                            {
                                _print_keys = true;
                                _brace = '{';
                                break;
                            }
                            else
                            {
                                if(data[k].second.type != Type::Char)
                                {
                                    _brace = '{';
                                }
                            }
                        }
                    }
                    catch(...)
                    {
                        _print_keys = true;
                        _brace = '{';
                    }

                    s += _brace;

                    bool _flag = false;

                    for(unsigned k = 0; k < data.size(); k++)
                    {
                        if(_brace == '{')
                        {
                            if(_flag)
                            {
                                s += ", ";
                            }
                            else
                                _flag = true;
                        }

                        if(_print_keys)
                        {
                            s += to_string(data[k].first);
                            s += ": ";
                        }

                        if(_brace == '"')
                        {
                            s += to_printable_ASCII(data[k].second);
                        }
                        else
                            s += to_string(data[k].second);
                    }

                    if(_brace == '{')
                        s += '}';
                    else
                        s += _brace;
                }
                break;
        }

        return s;
    }

    template<class T>
    operator std::vector<T>() const
    {
        std::vector<T> res;

        switch(type)
        {
            /*case Type::Bool:
                res.emplace_back(_bool);
                break;

            case Type::Char:
                res.emplace_back(_char);
                break;

            case Type::Integer:
                res.emplace_back(_integer);
                break;

            case Type::Number:
                res.emplace_back(_number);
                break;*/

            case Type::Bool:
            case Type::Char:
            case Type::Integer:
            case Type::Number:
                res.emplace_back((T) *this);
                break;

            case Type::Array:
                res.resize(data.size());
                for(unsigned k = 0; k < res.size(); k++)
                {
                    res[k] = (T) data[k].second;
                }

                break;

            default:
                break;
        }

        return res;
    }

    template<class TKey, class TValue>
    operator std::map<TKey, TValue>() const
    {
        std::map<TKey, TValue> res;

        switch(type)
        {
            case Type::Bool:
                res[0] = _bool;
                break;

            case Type::Char:
                res[0] = _char;
                break;

            case Type::Integer:
                res[0] = _integer;
                break;

            case Type::Number:
                res[0] = _number;
                break;

            case Type::Array:
                for(unsigned k = 0; k < data.size(); k++)
                {
                    res[data[k].first] = data[k].second;
                }

                break;

            default:
                break;
        }

        return res;
    }

    decltype(auto) keys() const
    {
        if(type == Type::Array)
        {
            struct TKeys
            {
                const TObject* obj;

                TKeys(const TObject* obj): obj(obj)
                {
                }

                unsigned size() const
                {
                    return obj->data.size();
                }

                const TObject& operator[](unsigned idx) const
                {
                    return obj->data[idx].first;
                }
            };

            return TKeys(this);
        }

        throw;
    }

    decltype(auto) values() const
    {
        if(type == Type::Array)
        {
            struct TValues
            {
                const TObject* obj;

                TValues(const TObject* obj): obj(obj)
                {
                }

                unsigned size() const
                {
                    return obj->data.size();
                }

                const TObject& operator[](unsigned idx) const
                {
                    return obj->data[idx].second;
                }
            };

            return TValues(this);
        }

        throw;
    }

    unsigned size() const
    {
        return data.size();
    }
};

const TObject null;

///--------------------------------------------------------

inline std::string to_string(const TObject& obj)
{
    return obj;
}

template<class T>
struct Array
{
    std::vector<T> data;

    unsigned size() const
    {
        return data.size();
    }

    void add(const T& v)
    {
        data.push_back(v);
    }

    template<class... Args>
    void emplace(Args&&... args)
    {
        data.emplace_back(std::forward<Args>(args)...);
    }

    const T& operator[](unsigned idx) const
    {
        static T t;

        return idx >= data.size() ? t : data[idx];
    }

    T& operator[](unsigned idx)
    {
        static T t;

        return idx >= data.size() ? t : data[idx];
    }
};

struct InvalidObject
{
};

inline TObject from_string(const std::string& s)
{
    struct String
    {
        const std::string& data;

        String(const std::string& data): data(data)
        {
        }

        char operator[](unsigned k) const
        {
            return k >= data.size() ? 0 : data[k];
        }

        unsigned size() const
        {
            return data.size();
        }

        decltype(auto) substr(unsigned idx, unsigned n) const
        {
            return data.substr(idx, n);
        }
    };

    auto data = String(s);

    ///-------------------------------------------------------------

    auto is_space = [](char ch)
    {
        return ch == ' ' || ch == '\t' || ch == 10 || ch == 13;
    };

    auto is_alpha = [](char ch)
    {
        return ch >= 'a' && ch <= 'z';
    };

    auto is_numeric = [](char ch)
    {
        return (ch >= '0' && ch <= '9') || ch == '+' || ch == '-' || ch == '.' || ch == 'e' || ch == 'E';
    };

    auto InvalidObjectError = [&null]() -> TObject
    {
        ///throw InvalidObject();

        return null;
    };

    struct Token
    {
        enum Type: short
        {
            None = 0,           /// <invalid>
            BeginArray,         /// {
            EndArray,           /// }
            Null,
            True,
            False,
            Char,
            Integer,
            Number,
            String,             /// special case as " cannot be nested
            Colon,              /// :
            Separator           /// ,
        };

        Type type;
        unsigned idx0, idx1;

        Token(): type(None)
        {
        }

        Token(Type type, unsigned idx0, unsigned idx1):
              type(type), idx0(idx0), idx1(idx1)
        {
        }

        const char* GetType() const
        {
            switch(type)
            {
                case None:
                    return "None";

                case BeginArray:
                    return "BeginArray";

                case EndArray:
                    return "EndArray";

                case Null:
                    return "Null";

                case True:
                    return "True";
                    break;

                case False:
                    return "False";

                case Char:
                    return "Char";

                case Integer:
                    return "Integer";

                case Number:
                    return "Number";

                case String:
                    return "String";

                case Colon:
                    return "Colon";

                case Separator:
                    return "Separator";

                default:
                    return "";
            }
        }

        void Print() const
        {
            auto str = GetType();
            short n = strlen("BeginArray") - strlen(str) + 1;

            print('[', str, ']');

            for(; n > 0; n--)
            {
                print(' ');
            }
        }
    };

    auto print_tokens = [&data](const auto& tokens)
    {
        int depth = 0;

        for(unsigned k = 0; k < tokens.size(); k++)
        {
            if(tokens[k].type == Token::BeginArray)
            {
                depth++;
            }
            else if(tokens[k].type == Token::EndArray)
            {
                depth--;
            }

            if(depth > 0)
            {
                for(unsigned k = 1; k < depth; k++)
                {
                    print("   ");
                }
            }

            tokens[k].Print();

            unsigned idx0 = tokens[k].idx0;
            unsigned idx1 = tokens[k].idx1;

            if(idx0 <= idx1)
                print("    ", data.substr(idx0, idx1 + 1 - idx0));

            println();
        }
    };

    auto tokenize = [&is_space, &is_alpha, &is_numeric, &InvalidObjectError](const String& s)
    {
        Array<Token> tokens;

        unsigned k = 0;
        while(k < s.size())
        {
            if(s[k] == '{')
            {
                tokens.emplace(Token::BeginArray, k, k);
            }
            else if(s[k] == '}')
            {
                tokens.emplace(Token::EndArray, k, k);
            }
            else if(s[k] == '"' || s[k] == '\'')
            {
                /// locate string or char literal

                const auto& ch = s[k];
                unsigned p = k+1;

                while(s[p] && s[p] != ch)
                {
                    if(s[p] == '\\')
                        p++;

                    p++;
                }

                if(s[p] != ch)
                {
                    println("Closing ", ch, " is missing");

                    InvalidObjectError();
                    return Array<Token>();
                }

                if(ch == '"')
                {
                    tokens.emplace(Token::String, k+1, p-1);
                }
                else
                {
                    tokens.emplace(Token::Char, k+1, p-1);
                }

                k = p;
            }
            else if(s[k] == ':')
            {
                tokens.emplace(Token::Colon, k, k);
            }
            else if(s[k] == ',')
            {
                tokens.emplace(Token::Separator, k, k);
            }
            else if(is_alpha(s[k]))
            {
                /// void | bool

                unsigned p = k;

                while(is_alpha(s[p]))
                {
                    p++;
                }

                auto w = s.substr(k, p-k);

                if(w == "null")
                {
                    tokens.emplace(Token::Null, k, p);
                }
                else if(w == "true")
                {
                    tokens.emplace(Token::True, k, p);
                }
                else if(w == "false")
                {
                    tokens.emplace(Token::False, k, p);
                }
                else
                {
                    println("Unknown identifier: ", w);

                    InvalidObjectError();
                    return Array<Token>();
                }

                k = p-1;
            }
            else if(is_numeric(s[k]))
            {
                unsigned p = k;
                bool _int = true;

                while(is_numeric(s[p]))
                {
                    if(_int && (s[p] == '.' || s[p] == 'e' || s[p] == 'E'))
                    {
                        _int = false;
                    }

                    p++;
                }

                if(_int)
                {
                    tokens.emplace(Token::Integer, k, p-1);
                }
                else
                {
                    tokens.emplace(Token::Number, k, p-1);
                }

                k = p-1;
            }
            else if(is_space(s[k]))
            {
                while(is_space(s[k]))
                {
                    k++;
                }

                k--;
            }
            else
            {
                println("Invalid character encountered: ", s[k]);

                InvalidObjectError();
                return Array<Token>();
            }

            k++;
        }

        return tokens;
    };

    auto parse_token = [&data, &InvalidObjectError](const Token& token) -> TObject
    {
        auto parse_char = [&data, &InvalidObjectError](unsigned& idx)
        {
            auto is_digit = [](auto ch)
            {
                return ch >= '0' && ch <= '9';
            };

            unsigned ch = 0;

            if(data[idx] == '\\')
            {
                idx++;

                if(data[idx] == 'n')
                {
                    ch = '\n';
                }
                else if(data[idx] == 't')
                {
                    ch = '\t';
                }
                else if(data[idx] == 'a')
                {
                    ch = '\a';
                }
                else if(data[idx] == 'v')
                {
                    ch = '\v';
                }
                else if(data[idx] == 'r')
                {
                    ch = '\r';
                }
                else if(data[idx] == '0')
                {
                    ch = 0;
                }
                else if(data[idx] == 'b')
                {
                    ch = '\b';
                }
                else if(data[idx] == 'e')       /// not standard for C
                {
                    ch = 27;
                }
                else if(data[idx] == 'f')
                {
                    ch = '\f';
                }
                else if(data[idx] == '\\')
                {
                    ch = '\\';
                }
                else if(data[idx] == '?')
                {
                    ch = '?';
                }
                else if(data[idx] == '\'')
                {
                    ch = '\'';
                }
                else if(data[idx] == '\"')
                {
                    ch = '\"';
                }
                else if(data[idx] == 'u')       /// '\udd...d;' where `d` is a decimal digit
                {
                    unsigned idx0 = idx-1;

                    ch = 0;
                    idx++;

                    while(is_digit(data[idx]))
                    {
                        ch = 10*ch + (data[idx] - '0');
                        idx++;
                    }

                    if(data[idx] != ';')
                    {
                        println("Invalid character literal: ", data.substr(idx0, idx - idx0 + 1));
                        InvalidObjectError();
                    }
                }
                else                            /// use verbatim
                {
                    ch = data[idx];
                }
            }
            else
            {
                ch = data[idx];
            }

            return ch;
        };

        ///---------------------------------------------

        TObject res;

        switch(token.type)
        {
            case Token::Null:
                ///res = null;
                break;

            case Token::True:
                res = true;
                break;

            case Token::False:
                res = false;
                break;

            case Token::Char:
                {
                    /// parse char

                    const unsigned& idx0 = token.idx0;
                    const unsigned& idx1 = token.idx1;

                    unsigned idx = idx0;

                    res = parse_char(idx);
                    res.type = Type::Char;

                    if(idx != idx1)
                    {
                        println("Invalid character literal: ", data.substr(idx0, idx1 + 1 - idx0));
                        return InvalidObjectError();
                    }
                }
                break;

            case Token::Integer:
                {
                    /// parse integer

                    unsigned idx0 = token.idx0;
                    unsigned idx1 = token.idx1;

                    res = std::stoll(data.substr(idx0, idx1 + 1 - idx0));
                }
                break;

            case Token::Number:
                {
                    /// parse a floating-point literal

                    unsigned idx0 = token.idx0;
                    unsigned idx1 = token.idx1;

                    res = std::stod(data.substr(idx0, idx1 + 1 - idx0));
                }
                break;

            case Token::String:
                {
                    /// parse string

                    std::string s;

                    unsigned idx0 = token.idx0;
                    unsigned idx1 = token.idx1;

                    unsigned k = idx0;
                    while(k <= idx1)
                    {
                        s += parse_char(k);
                        k++;
                    }

                    res = s;
                }
                break;

            default:
                println("A structure token [", token.GetType(), "] cannot be represented by an object");
                return InvalidObjectError();

                break;
        }

        return res;
    };

    ///--------------------------------------------

    TObject res;

    auto tokens = tokenize(s);

    /// {current, idx, value_assigned} | stack
    ///
    std::vector<std::tuple<TObject*, int, bool>> obj;

    TObject* current = &res;
    int idx = -1;
    bool value_assigned = false;

    TObject* override_current = nullptr;
    int depth = 0;

    unsigned k = 0;
    while(k < tokens.size())
    {
        const auto& type = tokens[k].type;

        if(type == Token::None || type == Token::Colon)
        {
            println("Invalid token encountered");
            return InvalidObjectError();
        }
        else if(type == Token::BeginArray)
        {
            obj.emplace_back(current, idx, value_assigned);
            depth++;

            if(idx >= 0)
            {
                if(override_current)
                {
                    current = override_current;
                    override_current = nullptr;
                }
                else
                {
                    if(value_assigned)
                    {
                        println("Missing separator");
                        return InvalidObjectError();
                    }

                    current = &((*current)[idx]);
                }
            }

            idx = 0;
            value_assigned = false;
        }
        else if(type == Token::EndArray)
        {
            if(obj.size())
            {
                current         = std::get<0>(obj.back());
                idx             = std::get<1>(obj.back());
                value_assigned  = std::get<2>(obj.back());

                obj.pop_back();
            }
            else
            {
                println("Unbalanced [EndArray] token encountered");
                return InvalidObjectError();
            }

            depth--;

            ///-----------------------------------

            if(!obj.size())     /// root object reached -> end parsing
            {
                if(k < tokens.size()-1)
                {
                    println("Extra data encountered, skipping");
                }

                break;
            }
        }
        else if(type == Token::Separator)
        {
            if(idx < 0)
            {
                println("Leading comma encountered");
                return InvalidObjectError();
            }
            else
                idx++;

            value_assigned = false;
        }
        else
        {
            if(idx == -1)       /// assign to the root object directly
            {
                if(tokens[k+1].type == Token::Colon)
                {
                    println("Colon outside the array encountered");
                    return InvalidObjectError();
                }

                if(value_assigned)
                {
                    println("Missing separator");
                    return InvalidObjectError();
                }

                *current = parse_token(tokens[k]);
                value_assigned = true;

                ///-----------------------------------

                if(k < tokens.size()-1)
                {
                    println("Extra data encountered, skipping");
                }

                break;
            }
            else
            {
                if(tokens[k+1].type == Token::Colon)
                {
                    auto key = parse_token(tokens[k]);
                    override_current = &((*current)[key]);

                    if(key.type == Type::Integer)
                    {
                        idx = (unsigned)key;
                    }

                    k++;
                }
                else
                {
                    if(value_assigned)
                    {
                        println("Missing separator");
                        return InvalidObjectError();
                    }

                    if(override_current)
                    {
                        *override_current = parse_token(tokens[k]);
                        override_current = nullptr;
                    }
                    else
                        (*current)[idx] = parse_token(tokens[k]);

                    value_assigned = true;
                }
            }
        }

        k++;
    }

    if(depth != 0)
    {
        if(depth == 1)
            println("Missing closing [EndToken] was recovered");
        else
            println(depth, " missing closing [EndToken]-s were recovered");
    }

    return res;
}

#endif // TOBJECT_H_INCLUDED
