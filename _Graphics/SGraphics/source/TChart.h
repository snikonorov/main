#ifndef GTABLE_H_INCLUDED
#define GTABLE_H_INCLUDED

#include <algorithm>

#include "Font.h"
///
///    <- #include "Graphics.h"

#include <tuple>

#define $tuple std::make_tuple

namespace
{
	template<class T>
    T declval();

    template<unsigned N, class T, class... Args>
    auto get_default(const std::tuple<Args...>& t, const T&) -> decltype(declval<typename std::enable_if<N < sizeof...(Args)>::type>(), declval<T>())
    {
        return std::get<N>(t);
    }

    template<unsigned N, class Tuple, class T>
    T get_default(const Tuple&, const T& def_value)
    {
        return def_value;
    }
}

class TChart;

struct Drawable
{
	virtual void Draw(const TChart* chart) const = 0;

    virtual ~Drawable()
    {
    }
};

struct Graph : public std::vector<Point>, public Drawable
{
    enum
    {
        None = 0,
		Points = 1,
        Lines = 2,
		InfiniteLines = 4,
        Bresenham = 8,
        Wu = 16,
        DDA = 32
    };

    class GraphStyle : public TStyle
    {
        public:

        unsigned code;
        bool Visible;

        GraphStyle(const pixel& color = pixel(140, 20, 10, 255), unsigned code = Lines + DDA, bool Visible = true)
                : TStyle(color), code(code), Visible(Visible)
        {
        }

        GraphStyle(unsigned code, const pixel& color = pixel(140, 20, 10, 255), bool Visible = true)
                : TStyle(color), code(code), Visible(Visible)
        {
        }

        operator const unsigned&() const
        {
            return code;
        }

        operator unsigned&()
        {
            return code;
        }
    };

    mutable GraphStyle style;

    Graph(const GraphStyle& style = GraphStyle{})
          : style(style)
    {
    }

    Graph(const std::initializer_list<Point>& l, const GraphStyle& style = GraphStyle{})
          : std::vector<Point>::vector(l), style(style)
    {
    }

    Graph(const std::vector<Point>& v, const GraphStyle& style = GraphStyle{})
          : std::vector<Point>::vector(v), style(style)
    {
    }

    template<class T, class = decltype((double)std::declval<T>())>
    Graph(const std::vector<T>& v, const GraphStyle& style = GraphStyle{})
          : style(style)
    {
        this->resize(v.size());
		for (unsigned k = 0; k < v.size(); k++)
		{
			this->at(k) = Point(k, v[k]);
		}

		//assign(v.begin(), v.end());
    }

	template<class T, class = void, class = decltype(Point(std::declval<T>()))>
    Graph(const std::vector<T>& v, const GraphStyle& style = GraphStyle{})
          : style(style)
    {
        this->resize(v.size());
		for (unsigned k = 0; k < v.size(); k++)
		{
			this->at(k) = Point(v[k]);
		}

		//assign(v.begin(), v.end());
    }

	using std::vector<Point>::assign;

	template<class TContainer>
	void assign(const TContainer& v)
	{
		this->resize(v.size());
		for (unsigned k = 0; k < v.size(); k++)
		{
			this->at(k) = v[k];
		}
	}

	void sort()
	{
		std::sort(begin(), end());
	}

    void Draw(const TChart* chart) const override;

	~Graph()
	{
	}
};

class TChart : public TComponent
{
    protected:

    static const int _l = 100;            ///scale unit length
    static const unsigned char tr = 40;   ///grid transparency

    double _hx, _hy;

    ///----------------------------------------------------------------------------------------

    public:
	
	static sFont* _font_default;

    enum
    {
        None = 0,
        Movable = 2,
        Grid = 4,
		Coords = 8,
		Scalable = 16
    };

    //---------------------------------------------------------------

    //protected:
    Point K;       ///scale

    //public:
    Point Offset;
    Point Cr;      ///axis cross point

    //Coord Sz;      ///(Width, Height)

    std::vector<Drawable*> data;

    sFont* font;
    unsigned style;

    Coord Pos, Size;

	TChart();

    TChart(sFont* font, const Coord& Pos = Coord{0, 0}, const Coord& Size = Coord{-1, -1}, unsigned style = Movable + Grid + Coords + Scalable,
		   const Point& k = Point{1, 1},
		   const Point& o = Point{0, 0},
		   const Point& cr = Point{0, 0});

    ~TChart();

	size_t size() const
	{
		return data.size();
	}

    void push_back(Drawable* g)
    {
        data.push_back(g);
    }

    /*void push_back(const Graph& g)
    {
        data.push_back(new Graph(g));
    }*/

    Drawable* &back()
    {
        return data.back();
    }

    /*template<class... Args>
    void emplace_back(Args&&... args)
    {
        data.push_back(new Graph(std::forward<Args>(args)...));
    }*/

    Drawable* &operator[](int k)
    {
        if(k+1u > data.size())
            data.resize(k+1);

        return data[k >= 0 ? k % data.size() : (k-1) % data.size()];
    }

    Drawable* &at(int k)
    {
        return (*this)[k];
    }

    void Scale(const Point& a);

    void Scale(double a)
    {
        Scale(Point(a));
    }

    void pScale(const Point& a, const Point& p0);

	void pScale(double a, const Point& p0)
	{
		pScale(Point(a), p0);
	}

	Point ScreenToReal(const Point& p) const
	{
		const int& Width = Size.x < 0 ? AppSize.x : Size.x;
		const int& Height = Size.y < 0 ? AppSize.y : Size.y;

		return {(p.x - Width/2 - Offset.x)/K.x + Cr.x, (Height/2 - p.y + Offset.y)/K.y + Cr.y};

		///const Point Sz(Size.x < 0 ? AppSize.x : Size.x, Size.y < 0 ? AppSize.y : Size.y);
		///return Cr + (p - Sz/2 - Offset)/Point(K.x, -K.y);
	}

	Point RealToScreen(const Point& p) const
	{
		const Point Sz(Size.x < 0 ? AppSize.x : Size.x, Size.y < 0 ? AppSize.y : Size.y);

		return (p - Cr)*Point(K.x, -K.y) + Sz/2 + Offset;
	}

    void ProcessEvent(sEvent& event) override;

	/// render chart as an SVG
	///
	std::string SVG() const;
};

#endif // GTABLE_H_INCLUDED
