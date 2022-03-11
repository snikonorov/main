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

	template<typename Functor, class... Args>
	void assign(const std::tuple<Functor, Args...>& F)
	{
		clear();

		///h0 = (b-a)/2^n <= 2^-p

		///const char p = 30;
		///double h0 = (b-a)/pow(2, (int)log2(b-a) + p);

		double h, y0, y, yn;

		const Functor& f = std::get<0>(F);
		double a = std::get<1>(F);
		double b = std::get<2>(F);

		if (a > b)
		{
			y = a;
			a = b;
			b = y;
		}

		const double l0 = b-a;

		const double e = get_default<3>(F, 1e-6);

		//const double h0 = get_default<4>(F, l0*e);
		const double h0 = get_default<4>(F, 1e-6);

		const double hmax = get_default<5>(F, l0/20);

		h = h0;

		const double P = e;
		emplace_back(a, f(a));

		while(a < b)
		{
			y0 = f(a);
			y = f(a + h);
			yn = f(a + 2*h);

			if(!std::isfinite(y0) || !std::isfinite(y))       ///miscontinuity
			{
				bool _cnt = true;
				double a0 = a;

				while(!std::isfinite(y) && h<64*h0 && a0+h < b)
				{
					if(std::isfinite(f(a0+h+h0/2)))
					{
						_cnt = false;

						a += h0/2;
						y = f(a+h);

						break;
					}

					h *= 2;
					y = f(a0+h);
				}

				if(h >= 64*h0)
				{
					h /= 2;

					while(!std::isfinite(y) && a0+h < b)
					{
						h += 32*h0;
						y = f(a0+h);
					}

					if(std::isfinite(y))
					{
						h = 64*h0;
					}
					else break;
				}

				if(_cnt)
				{
					if(!std::isfinite(y)) a0 = b;
					else a0 = a0+h;

					h /= 2;

					while(1)
					{
						do
						{
							h /= 2;

							a0 -= h;

							y = f(a0);
						}
						while(std::isfinite(y) && h>h0);

						if(!std::isfinite(f(a0+h-h0)))
						{
							a = a0;
							y = f(a+h);
							break;
						}

						if(h>h0) a0 += h;
						else break;
					}
				}
			}
			else
			{
				double _y = (y+y0)/2;

				if(fabs( f(a+h/2) - _y ) > P*fabs(_y))
				{
					while(fabs( f(a+h/2) - _y ) > P*fabs(_y) && h > h0)
					{
						h /= 2;
						y = f(a+h);

						_y = (y+y0)/2;
					}
				}
				else
				{
					const double g0 = 1.1;
					yn = f(a + g0*h);

					//while(fabs( y-(yn+y0)/2 ) < P*fabs(y) && h < hmax)// && a+2*h <= b)	//&& h<32*h0
					while(fabs( y - (y0*(g0 - 1) + yn)/g0 ) <= P*fabs(y) && h < hmax)
					{
						h *= g0;
						//h += h0;

						y = yn;
						yn = f(a+g0*h);

						if (a+g0*h > b)
						{
							emplace_back(a+h, y);
							emplace_back(b, f(b));
							return;
						}
					}
				}
			}

			emplace_back(a+h, y);

			a += h;
			if(a >= b)
			{
				emplace_back(b, f(b));
				break;
			}

			h = h0;
		}
	}

    ///template<typename Functor>
    ///Graph(const Functor& f, double a, double b, double e = 1e-3, double h0 = 1e-6)
    template<typename Functor, class... Args>
    Graph(const std::tuple<Functor, Args...>& F, const GraphStyle& style = GraphStyle{}): style(style)
    {
        assign(F);
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
