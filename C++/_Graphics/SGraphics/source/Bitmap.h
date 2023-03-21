#ifndef BITMAP_H_INCLUDED
#define BITMAP_H_INCLUDED

#include <string.h>

#include <vector>
#include <cmath>

template<class T>
struct TPixel
{
    static const unsigned VAL_MAX = 255;

    T a,r,g,b;

    TPixel()
    {
        a = r = g = b = 0;
    }

    explicit TPixel(T A): a(A), r(0), g(0), b(0)
    {
        //a = A;
        //r = g = b = 0;
    }

    TPixel(T R, T G, T B): a(VAL_MAX), r(R), g(G), b(B)
    {
        //a = VAL_MAX;
    }

    TPixel(T A, T R, T G, T B): a(A), r(R), g(G), b(B)
    {
    }

    template<class Q>
    TPixel(const TPixel<Q>& p): a(p.a), r(p.r), g(p.g), b(p.b)
    {
    }

    template<class Q>
    TPixel& operator=(TPixel<Q>& p)
    {
        a = p.a;
        r = p.r;
        g = p.g;
        b = p.b;

        return *this;
    }

    /// #VS Error
	/// void operator=(const TPixel& p) = default;

    double gray() const
    {
        return a*(0.2126*r + 0.7152*g + 0.0722*b)/VAL_MAX/VAL_MAX;
    }

    double Grey() const
    {
        return (0.2126*r + 0.7152*g + 0.0722*b)/VAL_MAX;
    }

    template<class Q>
    void operator+=(const TPixel<Q>& p)
    {
        if(p.a == VAL_MAX || a == 0) *this = p;
        else
        {
            a += p.a - a*p.a/VAL_MAX;
            r -= (r - p.r)*p.a/VAL_MAX;
            g -= (g - p.g)*p.a/VAL_MAX;
            b -= (b - p.b)*p.a/VAL_MAX;
        }
    }

    template<class Q>
    auto operator+(const TPixel<Q>& p) const -> TPixel<decltype(r + p.r)>
    {
        return TPixel<decltype(r + p.r)>(a + p.a - a*p.a/VAL_MAX,
                                         r - (r - p.r)*p.a/VAL_MAX,
                                         g - (g - p.g)*p.a/VAL_MAX,
                                         b - (b - p.b)*p.a/VAL_MAX);
    }

    TPixel operator*(double k) const
    {
        return TPixel((T) round(a*k), r,g,b);
    }

	template<class Q>
	friend TPixel operator*(double k, const TPixel<Q>& p)
	{
		return TPixel((T)round(p.a*k), p.r, p.g, p.b);
	}

    void operator*=(double k)
    {
        a = round(a*k);
    }

    TPixel operator^(double k) const
    {
        return TPixel(a, round(r*k), round(g*k), round(b*k));
    }

    void operator^=(double k)
    {
        r = round(r*k);
        g = round(g*k);
        b = round(b*k);
    }

    TPixel operator/(double k) const
    {
        return TPixel((T) round(a/k), r,g,b);
    }

    void operator/=(double k)
    {
        a = (T) round(a/k);
    }

    bool operator==(const TPixel& p) const
    {
        return (a==p.a && r==p.r && g==p.g && b==p.b);
    }

    bool operator!=(const TPixel& p) const
    {
        return (a!=p.a || r!=p.r || g!=p.g || b!=p.b);
    }

    TPixel operator!() const
    {
        return TPixel(a, VAL_MAX - r, VAL_MAX - g, VAL_MAX - b);
    }

	TPixel bgr() const
    {
        return TPixel(a, b, g, r);
    }
};

typedef TPixel<unsigned char> pixel;

///--------------------------------------------------------------------------------

#include "TRect.h"
///
///    <- #include "Vector3d.h"

///--------------------------------------------------------------------------------

struct mImage
{
    //int X0, Y0;          ///for window crop
    iRect ClipRect;

    std::vector<std::vector<pixel>> matrix;
    int Width, Height;

    mImage(): ClipRect(Coord(0, 0)), Width(0), Height(0)
    {
        //Width = Height = 0;
    }

    mImage(int width, int height);

    ///---------------------------------------------------------------------------

    ~mImage()
    {
        Width = Height = 0;
        ClipRect = Coord(0, 0);
    }

    bool Valid() const
    {
        return Width && Height;
    }

    pixel& operator()(unsigned x, unsigned y)        /// !!! check-free method !!!
    {
        return matrix[x][y];
    }

    pixel GetPixel(int x, int y)
    {
        return (ClipRect.contains(x, y)) ? (*this)(x, y) : pixel();
    }

    template<class Q>
    void PutPixel(int x, int y, const TPixel<Q>& color)
    {
        if(color.a && ClipRect.contains(x, y))  ///x >= 0 && y >= 0 && x < Width && y < Height)
        {
            matrix[x][y] += color;
        }
    }

    template<class Q>
    void PutWuPixel(double x, double y, const TPixel<Q>& c)
    {
        double a = x - floor(x);
        double b = y - floor(y);

        PutPixel((int)x, (int)y, c*(1-a)*(1-b));
        PutPixel((int)x, (int)ceil(y), c*(1-a)*b);
        PutPixel((int)ceil(x), (int)y, c*a*(1-b));
        PutPixel((int)ceil(x), (int)ceil(y), c*a*b);
    }

    template<class Img>
    void Create(const Img& img)
    {
        if(img.Valid())
        {
            resize(img.Width, img.Height);

            for(int x = 0; x < img.Width; x++)
            {
                for(int y = 0; y < img.Height; y++)
                {
                    (*this)(x, y) = img(x, y);
                }
            }
        }

        Width = img.Width;
        Height = img.Height;

        ClipRect.Width(Width - 1);
        ClipRect.Height(Height - 1);
    }

    template<class Img>
    void DrawImage(Img& img, int X, int Y)
    {
        for(int x=0; x<img.Width; x++)
        {
            for(int y=0; y<img.Height; y++)
            {
                PutPixel(x+X, y+Y, img(x, y));
            }
        }
    }

    void Clear()
    {
        for (unsigned x = 0; x < matrix.size(); x++)
        {
            for (unsigned y = 0; y < matrix[0].size(); y++)
            {
                matrix[x][y] = pixel(0);
            }
        }
    }

    ///-----------------------------------------------------------------------------------------

    void resize(int width, int height);

    void Fill(const pixel& color);              /// fill the entire image without clipping!!!
    void FillRect(int x, int y, int x1, int y1, const pixel& color);

    void DrawLine(int x,int y,int x1,int y1,const pixel& color,unsigned char type=0);

    //void DrawSLine(int x,int y,int x1,int y1,const pixel& color,unsigned char style);

	void DrawSmoothLine(double x0, double y0, double x1, double y1, const pixel& color, unsigned char type = 0);

    void DrawDDALine(double x0,double y0,double x1,double y1,const pixel& color,unsigned char type=0);

    void DrawWuLine(double x0, double y0, double x1, double y1,pixel color,unsigned char type=0);

    void DrawRect(int x,int y,int x1,int y1,const pixel& color);

    void DrawCircle(int X,int Y,int R,const pixel& color);

    void DrawDDACircle(double X,double Y,double R,const pixel& color);

    void FillCircle(int X,int Y,int R,const pixel& color);
};

///---------------------------------------------------------------------------------

struct dImage
{
    //int X0, Y0;          ///for window crop
    iRect ClipRect;

    ///pixel BK_COLOR;      <- SWindow::wStyle::color
    unsigned cnst, Cnst;
    unsigned char depth;	///in bytes!

    ///unsigned char* data;
    std::vector<unsigned char> data;
    int Width, Height;

    dImage(): ClipRect(Coord(0, 0)), cnst(0), Cnst(0), Width(0), Height(0)
    {
        //Width = Height = 0;
        //data = 0;
    }

    void Create(int width, int height);

    dImage(int width, int height)
    {
        Create(width, height);
    }

    dImage(const dImage&) = delete;

    ~dImage()
    {
        Cnst = cnst = 0;
        Width = Height = 0;
    }

    bool Valid() const
    {
        return Width && Height;
    }

    TPixel<unsigned char&> operator()(int x, int y)       /// check-free method !!!
    {
        static unsigned char a, r, g, b;
        a = 255;

        /// if(depth == 3) !!!
        return Valid() ? TPixel<unsigned char&>(a, data[cnst*y + depth*x+2], data[cnst*y + depth*x+1], data[cnst*y + depth*x])
                       : TPixel<unsigned char&>(a, r, g, b);
    }

    pixel GetPixel(int x, int y)
    {
        /// if(depth == 3) !!!
        /**return (Valid() && x>=0 && y>=0 && x<Width && y<Height) ? pixel(data[cnst*y + depth*x+2], data[cnst*y + depth*x+1], data[cnst*y + depth*x])
                                                                : pixel();*/

        return (Valid() && ClipRect.contains(x, y)) ? pixel(data[cnst*y + depth*x+2], data[cnst*y + depth*x+1], data[cnst*y + depth*x])
                                                    : pixel();
    }

    template<class Q>
    void PutPixel(int x, int y, const TPixel<Q>& color)
    {
        if(Valid() && color.a && ClipRect.contains(x, y))
        {
            if(color.a == 255)
            {
                data[cnst*y + depth*x] = color.b;
                data[cnst*y + depth*x+1] = color.g;
                data[cnst*y + depth*x+2] = color.r;
            }
            else
            {
                data[cnst*y + depth*x] -= (data[cnst*y + depth*x] - color.b)*color.a/255;
                data[cnst*y + depth*x+1] -= (data[cnst*y + depth*x+1] - color.g)*color.a/255;
                data[cnst*y + depth*x+2] -= (data[cnst*y + depth*x+2] - color.r)*color.a/255;
            }
        }
    }

    template<class Q>
    void PutWuPixel(double x, double y, const TPixel<Q>& c)
    {
        double _a = x - floor(x);
        double _b = y - floor(y);

        PutPixel((int) x, (int) y, c*(1 - _a)*(1 - _b));
        PutPixel((int) x, (int) ceil(y), c*(1 - _a)*_b);
        PutPixel((int) ceil(x), (int) y, c*_a*(1 - _b));
        PutPixel((int) ceil(x), (int) ceil(y), c*_a*_b);
    }

    void Create(mImage& img)
    {
        if(img.Valid())
        {
            resize(img.Width, img.Height);

            for(int x = 0; x < img.Width; x++)
            {
                for(int y = 0; y < img.Height; y++)
                {
                    (*this)(x, y) = img(x, y);
                }
            }
        }

        Width = img.Width;
        Height = img.Height;

        ClipRect.Width(Width - 1);
        ClipRect.Height(Height - 1);
    }

    void Create(dImage& img)
    {
        if(img.Valid())
        {
            resize(img.Width, img.Height);

            memcpy(data.data(), img.data.data(), img.Cnst);

            Cnst = img.Cnst;
            cnst = img.cnst;
        }

        Width = img.Width;
        Height = img.Height;

        ClipRect.Width(Width - 1);
        ClipRect.Height(Height - 1);
    }

    template<class Img>
    void DrawImage(Img& img, int X, int Y)
    {
        if(Valid())
            for(int x=0; x<img.Width; x++)
            {
                for(int y=0; y<img.Height; y++)
                {
					PutPixel(x + X, y + Y, img(x, y));
                }
            }
    }

    //void Clear(const pixel& color=pixel(0), int X0=-1, int Y0=-1, int X1=-1, int Y1=-1);

    void Fill(const pixel& color)       /// fill the entire image without clipping!!!
    {
        unsigned n = 0;

        if(Valid())
        {
            if(color.a == 255)
            {
                if(color.r == color.g && color.r == color.b) memset(data.data(), color.r, Cnst);
                else
                    while(n < Cnst)
                    {
                        data[n] = color.b;
                        data[n+1] = color.g;
                        data[n+2] = color.r;

                        n += depth;
                    }
            }
            else
                while(n < Cnst)
                {
                    data[n] -= (data[n] - color.b)*color.a/255;
                    data[n+1] -= (data[n+1] - color.g)*color.a/255;
                    data[n+2] -= (data[n+2] - color.r)*color.a/255;

                    n += depth;
                }
        }
    }

    void resize(int width, int height);

    //void Fill_alpha(const pixel& color);

    void FillRect(int x,int y,int x1,int y1,const pixel& color);

    void DrawLine(int x,int y,int x1,int y1, const pixel& color,unsigned char type=0);

    void DrawSLine(int x,int y,int x1,int y1,const pixel& color,unsigned char style);

    void DrawDDALine(double x0,double y0,double x1,double y1,const pixel& color,unsigned char type=0);

    void DrawWuLine(double x0, double y0, double x1, double y1,pixel color,unsigned char type=0);

    void DrawRect(int x,int y,int x1,int y1,const pixel& color);

    void DrawCircle(int X,int Y,int R,const pixel& color);

    void DrawDDACircle(double X,double Y,double R,const pixel& color);

    void FillCircle(int X,int Y,int R,const pixel& color);
};

typedef mImage Image;

#endif // BITMAP_H_INCLUDED
