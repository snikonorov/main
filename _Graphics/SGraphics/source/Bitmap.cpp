#include "Bitmap.h"

///# Refine clipping algorithm!
///# "Scale bug" occurs when using this function!
template<class T, class Img>
inline bool _Image_ClipLine(Img* img, T &x0, T &y0, T &x1, T &y1, T& dx, T& dy, const pixel& color)
{
    if(!img->Valid()) return false;

    if(std::abs(dx) < 1 && std::abs(dy) < 1)
    {
        img->PutWuPixel(x0, y0, color);
        //img->PutWuPixel((x0+x1)/2, (y0+y1)/2, color);
        return false;
    }

    if(x0 < -1e-12 || x0 >= img->Width || y0 < -1e-12 || y0 >= img->Height)
    {
		if(x1 < -1e-12 || x1 >= img->Width || y1 < -1e-12 || y1 >= img->Height)
        {
            if((x0<0 && x1<0) || (x0>=img->Width && x1>=img->Width) ||
               (y0<0 && y1<0) || (y0>=img->Height && y1>=img->Height)) return false;
            else
            {
                if(std::abs(dx)<1e-3)
                {
                    y0 = 0;
                    y1 = img->Height-1;
                }
                else if(std::abs(dy)<1e-3)
                {
                    x0 = 0;
                    x1 = img->Width - 1;
                }
                else
                {
                    double y = y1 - (y1-y0)*(x1 - 0)/(x1 - x0);

                    if(y >= 0 && y < img->Height)
                    {
                        y0 = (T)y;
                        x0 = 0;
                    }
                    else
                    {
                        if(y <= 0)
                        {
                            x0 = (T) round(x1 + (0 - y1)*(x1-x0)/(y1-y0));
                            y0 = 0;
                        }
                        else
                        {
                            x0 = (T) round(x1 + (img->Height - 1 - y1)*(x1-x0)/(y1-y0));
                            y0 = img->Height-1;
                        }

                        if(x0 >= img->Width) x0 = img->Width-1;
                        else if(x0 < 0) x0 = 0;
                    }

                    y = y1 - (y1-y0)*(x1 - img->Width + 1)/(x1 - x0);

                    if(y >= 0 && y < img->Height)
                    {
                        y1 = (T) y;
                        x1 = img->Width-1;
                    }
                    else
                    {
                        if(y <= 0)
                        {
                            x1 = (T) round(x1 + (0-y1)*(x1-x0)/(y1-y0));
                            y1 = 0;
                        }
                        else
                        {
                            x1 = (T) round(x1 + (img->Height - 1 - y1)*(x1-x0)/(y1-y0));
                            y1 = img->Height-1;
                        }

                        if(x1 >= img->Width) x1 = img->Width-1;
                        else if(x1 < 0) x1 = 0;
                    }

                    dy = y1-y0;
                    dx = x1-x0;

                    if(std::abs(dx) < 0.5 && std::abs(dy) < 0.5)
                    {
                        img->PutWuPixel((x0+x1)/2, (y0+y1)/2, color);
                        //return;
                    }
                }
            }
        }
        else
        {
            std::swap(x0, x1);
            std::swap(y0, y1);

            dx *= -1;
            dy *= -1;
        }
    }

    return true;
}

///------------------------------------------------------------------------------------

mImage::mImage(int width, int height): ClipRect(width - 1, height - 1), Width(width), Height(height)
{
    //X0 = Y0 = 0;
    //Width = width;
    //Height = height;

    matrix.resize(Width);
    for(int k=0; k<Width; k++) matrix[k].resize(Height);    ///Extend(Height);

    //transparent = true;
}

void mImage::resize(int width, int height)
{
    if(Width < width)
        matrix.resize(width);

    if(Height < height || Width < width)
    {
        for(int x = 0; x < width; x++)
        {
            matrix[x].resize(height);
        }

        Height = height;
    }

	if(Width < width)
		Width = width;

    ClipRect.Width(Width - 1);
    ClipRect.Height(Height - 1);
}

void mImage::Fill(const pixel& color)
{
    for(int y = 0; y<Height; y++)
    {
        for(int x = 0; x<Width; x++)
        {
            matrix[x][y] += color;
        }
    }
}

void mImage::FillRect(int x, int y, int x1, int y1, const pixel& color)
{
    if(!Valid()) return;

    const auto& X0 = ClipRect.A.x;
    const auto& Y0 = ClipRect.A.y;
    const auto& Width = ClipRect.Width();
    const auto& Height = ClipRect.Height();

    if(x < X0) x = X0;
    else if(x >= Width) x = Width-1;

    if(x1 < X0) x1 = X0;
    else if(x1 >= Width) x1 = Width-1;

    if(y < Y0) y = Y0;
    else if(y >= Height) y = Height-1;

    if(y1 < Y0) y1 = Y0;
    else if(y1 >= Height) y1 = Height-1;

    if(x > x1) std::swap(x, x1);
    if(y > y1) std::swap(y, y1);

    int _y = y;

    do
    {
        do
        {
            matrix[x][y] += color;
            y++;
        }
        while(y != y1);

        y = _y;
        x++;
    }
    while(x != x1);
}

void mImage::DrawLine(int x,int y,int x1,int y1, const pixel& color, unsigned char type)
{
    int pdx=0, pdy=0, sx,sy,err,dec,inc;

    inc = std::abs(x1-x);
    dec = std::abs(y1-y);

    if(!(inc|dec))
    {
        if(type==1) PutPixel(x,y,color);
        return;
    }
    else
    {
        if(!_Image_ClipLine(this, x, y, x1, y1, inc, dec, color)) return;
    }
    ///else if( (x<X0 && x1<X0) || (x>Width && x1>Width) || (y<Y0 && y1<Y0) || (y>Height && y1>Height)) return;

    sx=(x1>x) ? 1 : (x1<x) ? -1 : 0;
    sy=(y1>y) ? 1 : (y1<y) ? -1 : 0;

    if (inc > dec) pdx = sx;
    else
    {
        std::swap(inc, dec);
        pdy = sy;
    }

    err = inc>>1;

    if (type != 1)
		PutPixel(x,y,color);

    for(int t=0; t<inc; t++)
    {
        err-=dec;

        if(err<0)
        {
            err+=inc;

            x+=sx;
            y+=sy;
        }
        else
        {
            x+=pdx;
            y+=pdy;
        }

        PutPixel(x,y,color);
    }
}

void mImage::DrawSmoothLine(double x0, double y0, double x1, double y1, const pixel& color, unsigned char type)
{
	double dx = (x1 - x0);
	double dy = (y1 - y0);

	if (!_Image_ClipLine(this, x0, y0, x1, y1, dx, dy, color)) return;

	Point p1(x0, y0);
	Point p2(x1, y1);

	Point dp = p2 - p1;

	if (std::abs(dp.x) > std::abs(dp.y))
	{
		int x1, x2;

		if (dp.x > 0)
		{
			x1 = floor(p1.x);// + 1;
			x2 = ceil(p2.x) - (type == 1);
		}
		else
		{
			x1 = floor(p2.x) + (type == 1);
			x2 = ceil(p1.x);
		}

		const double k = dp.y/dp.x;

		for (; x1 <= x2; x1++)
		{
			double y = p1.y + k*(x1 - p1.x);
			int y1 = round(y);

			Point P0[] = { Point(x1, y1), Point(x1, y1 - 1), Point(x1, y1 + 1) };

			for (const auto& p0: P0)
			{
				double t = (dp & (p0 - p1)) / (dp & dp);

				if (t < 0)
					t = 0;
				else if (t > 1)
					t = 1;

				double d = (p0 - p1 - dp*t).Norm2();
				if (d < 1)
					PutPixel(p0.x, p0.y, color*(1 - d));
			}
		}
	}
	else
	{
		int y1, y2;

		if (dp.y > 0)
		{
			y1 = floor(p1.y);// + 1;
			y2 = ceil(p2.y) - (type == 1);
		}
		else
		{
			y1 = floor(p2.y) + (type == 1);
			y2 = ceil(p1.y);
		}

		const double k = dp.x/dp.y;

		for (; y1 <= y2; y1++)
		{
			double x = p1.x + k*(y1 - p1.y);
			int x1 = round(x);

			Point P0[] = { Point(x1, y1), Point(x1 - 1, y1), Point(x1 + 1, y1) };

			for (const auto& p0 : P0)
			{
				double t = (dp & (p0 - p1)) / (dp & dp);

				if (t < 0)
					t = 0;
				else if (t > 1)
					t = 1;

				double d = (p0 - p1 - dp*t).Norm2();
				if (d < 1)
					PutPixel(p0.x, p0.y, color*(1 - d));
			}
		}
	}
}

void mImage::DrawDDALine(double x0,double y0,double x1,double y1,const pixel& color,unsigned char type)
{
    double dx = (x1-x0);
    double dy = (y1-y0);

    if(!_Image_ClipLine(this, x0, y0, x1, y1, dx, dy, color)) return;

    double L = (std::abs(dx) > std::abs(dy)) ? std::abs(dx)+1 : std::abs(dy)+1;

    dx = dx/L;
    dy = dy/L;

    L=round(L);

    double x,y;

    for(int k=0; k<=L; k++)
    {
        x=x0+k*dx;
        y=y0+k*dy;

        if((type==1 && k==(int)L) || x<0 || x>Width || y<0 || y>Height) break;

        PutWuPixel(x, y, color);
    }
}

void mImage::DrawWuLine(double x0, double y0, double x1, double y1, pixel color,unsigned char type)
{
    double dy=(y1-y0), dx=(x1-x0);

    if(!_Image_ClipLine(this, x0, y0, x1, y1, dx, dy, color)) return;

    int xi=(x1>x0) ? 1 : (x1<x0) ? -1 : 0;
    int yi=(y1>y0) ? 1 : (y1<y0) ? -1 : 0;

    double x=round(x0),y=round(y0);

    x1=round(x1);
    y1=round(y1);

    unsigned char A = color.a;

    if(fabs(dy)>fabs(dx))
    {
        while(yi*(y-y1)<=0)
        {
            if(yi*(y+yi-y1)>0)
            {
                x=x1;//round(x);
                y=y1;//round(y);
            }
            else x=dx*(y-y0)/dy+x0;

            color.a = (unsigned char) (A*(x - (int)x));
            if(color.a) PutPixel((int)x + 1, (int) y, color);

            color.a = (unsigned char) (A*(1 - x + (int)x));
            if(color.a) PutPixel((int)x, (int)y, color);

            y += yi;

            if(type == 1 && yi*(y+yi-y1)>0.6) break;

            //if(x<X0 || x>Width || y<Y0 || y>Height) break;
        }
    }
    else
    {
        while(xi*(x-x1)<=0)
        {
            if(xi*(x+xi-x1)>0)
            {
                x=x1;//round(x);
                y=y1;//round(y);
            }
            else y=dy*(x-x0)/dx+y0;

            color.a = (unsigned char) (A*(y - (int)y));
            if(color.a) PutPixel((int)x, (int)y+1, color);

            color.a = (unsigned char) (A*(1 - y + (int)y));
            if(color.a) PutPixel((int)x, (int)y, color);

            x+=xi;

            if(type==1 && xi*(x+xi-x1)>0.6) break;

            if(x<0 || x>Width || y<0 || y>Height) break;
        }
    }
}

void mImage::DrawRect(int x,int y,int x1,int y1,const pixel& color)
{
    if(x>x1)
    {
        x = x+x1;
        x1 = x-x1;
        x = x-x1;
    }

    if(y>y1)
    {
        y = y+y1;
        y1 = y-y1;
        y = y-y1;
    }

    DrawLine(x,y, x1,y,color);
    DrawLine(x1,y+1, x1,y1-1,color);
    DrawLine(x1,y1, x,y1,color);
    DrawLine(x,y1-1, x,y+1,color);
}

/*void Image::drawCircle(int X,int Y,int R,const pixel& color)
{
    int x=R,y=0,d=1-R;

    while(x >= y)
    {
        putPixel(x+X,y+Y,color);
        putPixel(-x+X,y+Y,color);

        if(x!=y)
        {
            putPixel(y+X,x+Y,color);
            putPixel(y+X,-x+Y,color);
        }

        if(y!=0)
        {
            putPixel(x+X,-y+Y,color);
            putPixel(-x+X,-y+Y,color);

            if(x!=y)
            {
                putPixel(-y+X,x+Y,color);
                putPixel(-y+X,-x+Y,color);
            }
        }

        y++;

        if(d<0) d += 2*y+1;
        else
        {
            x--;
            d += 2*(y-x+1);
        }
    }
}*/

void mImage::DrawCircle(int X,int Y,int R,const pixel& color)
{
    int x=0,y=R,d=3-2*R;

    while(x <= y)
    {
        PutPixel(x+X,y+Y,color);
        PutPixel(x+X,-y+Y,color);
        PutPixel(-x+X,-y+Y,color);
        PutPixel(-x+X,y+Y,color);
        PutPixel(y+X,x+Y,color);
        PutPixel(y+X,-x+Y,color);
        PutPixel(-y+X,-x+Y,color);
        PutPixel(-y+X,x+Y,color);

        if(d<0) d=d+4*x+6;
        else
        {
            d=d+4*(x-y)+10;
            y--;
        }

        x++;
    }
}

void mImage::DrawDDACircle(double X,double Y,double R,const pixel& color)
{
    const double M_SQRT2 = 1.414213562373095;

    double x = R*M_SQRT2/2;
    double y = x;

    PutWuPixel(X+x,Y+y, color);
    PutWuPixel(X+x,Y-y, color);
    PutWuPixel(X-x,Y+y, color);
    PutWuPixel(X-x,Y-y, color);

    x--;
    y = sqrt(R*R - x*x);

    while(x>=0)
    {
        PutWuPixel(X+x,Y+y, color);
        PutWuPixel(X-x,Y-y, color);

        PutWuPixel(X+y,Y+x, color);
        PutWuPixel(X-y,Y-x, color);

        PutWuPixel(X+x,Y-y, color);
        PutWuPixel(X-x,Y+y, color);

        PutWuPixel(X+y,Y-x, color);
        PutWuPixel(X-y,Y+x, color);

        x--;
        y = sqrt(R*R - x*x);
    }
}

void mImage::FillCircle(int X,int Y,int R,const pixel& color)
{
    int x=R,y=0,d=1-R;

    while(x >= y)
    {
        DrawLine(x+X,y+Y,-x+X,y+Y,color);

        if(y!=0)
        {
            DrawLine(x+X,-y+Y,-x+X,-y+Y,color);
        }

        y++;

        if(d<0) d += 2*y+1;
        else
        {
            x--;
            d += 2*(y-x+1);
        }
    }

    int _y=y;

    x=R,y=0,d=1-R;

    while(x >= y)
    {
        if(x!=y)
        {
            DrawLine(y+X,x+Y,y+X,Y+_y,color,1);
            if(y!=0) DrawLine(-y+X,x+Y,-y+X,Y+_y,color,1);

            if(x!=0)
            {
                DrawLine(y+X,-x+Y,y+X,Y-_y,color,1);
                if(y!=0) DrawLine(-y+X,-x+Y,-y+X,Y-_y,color,1);
            }
        }

        y++;

        if(d<0) d += 2*y+1;
        else
        {
            x--;
            d += 2*(y-x+1);
        }
    }
}

///=================================================================================

void dImage::Create(int width, int height)
{
    //X0 = Y0 = 0;

    ClipRect = Rect(width - 1, height - 1);

    Width = width;
    Height = height;

    depth = 3;

    //cnst = depth*(Width + Width%4);
    cnst = depth*Width;

    Cnst = cnst*Height;

    data.resize(Cnst);
}

void dImage::resize(int width, int height)
{
    if(Width < width || Height < height)
    {
        //cnst = depth*(width + width%4);

        cnst = depth*width;
        Cnst = cnst*height;

        data.resize(Cnst);
    }

    Width = width;
    Height = height;

    ClipRect.Width(Width - 1);
    ClipRect.Height(Height - 1);
}

void dImage::FillRect(int x, int y, int x1, int y1, const pixel& color)
{
    if(!Valid()) return;

    const auto& X0 = ClipRect.A.x;
    const auto& Y0 = ClipRect.A.y;
    const auto& Width = ClipRect.Width();
    const auto& Height = ClipRect.Height();

    if(x < X0) x = X0;
    else if(x >= Width) x = Width-1;

    if(x1 < X0) x1 = X0;
    else if(x1 >= Width) x1 = Width-1;

    if(y < Y0) y = Y0;
    else if(y >= Height) y = Height-1;

    if(y1 < Y0) y1 = Y0;
    else if(y1 >= Height) y1 = Height-1;

    if(x > x1) std::swap(x, x1);
    if(y > y1) std::swap(y, y1);

    int _y = y;

    if(color.a == 255)
    {
        if(color.r == color.g && color.r == color.b)
        {
            do
            {
                memset(data.data() + cnst*y + depth*x, color.r, (x1 - x + 1)*depth);
                y++;
            }
            while(y != y1);
        }
        else
        {
            do
            {
                do
                {
                    data[cnst*y + depth*x] = color.b;
                    data[cnst*y + depth*x+1] = color.g;
                    data[cnst*y + depth*x+2] = color.r;

                    y++;
                }
                while(y < y1);

                x++;
                y = _y;
            }
            while(x < x1);
        }
    }
    else
    {
        do
        {
            do
            {
                data[cnst*y + depth*x] -= (data[cnst*y + depth*x] - color.b)*color.a/255;
                data[cnst*y + depth*x+1] -= (data[cnst*y + depth*x+1] - color.g)*color.a/255;
                data[cnst*y + depth*x+2] -= (data[cnst*y + depth*x+2] - color.r)*color.a/255;

                y++;
            }
            while(y < y1);

            x++;
            y = _y;
        }
        while(x < x1);
    }
}

void dImage::DrawLine(int x,int y,int x1,int y1, const pixel& color, unsigned char type)
{
    int pdx=0, pdy=0, sx,sy,err,dec,inc;

    inc = std::abs(x1-x);
    dec = std::abs(y1-y);

    if(!(inc|dec))
    {
        if(type==1) PutPixel(x,y,color);
        return;
    }
    else
    {
        if(!_Image_ClipLine(this, x, y, x1, y1, inc, dec, color)) return;
    }
    ///else if( (x<X0 && x1<X0) || (x>Width && x1>Width) || (y<Y0 && y1<Y0) || (y>Height && y1>Height)) return;

    sx=(x1>x) ? 1 : (x1<x) ? -1 : 0;
    sy=(y1>y) ? 1 : (y1<y) ? -1 : 0;

    if (inc > dec) pdx = sx;
    else
    {
        std::swap(inc, dec);
        pdy = sy;
    }

    err = inc>>1;

    PutPixel(x,y,color);

    for(int t=0; t<inc; t++)
    {
        err-=dec;

        if(err<0)
        {
            err+=inc;

            x+=sx;
            y+=sy;
        }
        else
        {
            x+=pdx;
            y+=pdy;
        }

        PutPixel(x,y,color);
    }
}

/*void dImage::DrawDDALine(double x0,double y0,double x1,double y1,const pixel& color,unsigned char type)
{
    double dx = (x1-x0);
    double dy = (y1-y0);

    _Image_ClipLine(this, x0, y0, x1, y1, dx, dy, color);

    double L = (std::abs(dx) > std::abs(dy)) ? std::abs(dx)+1 : std::abs(dy)+1;

    dx = dx/L;
    dy = dy/L;

    L=round(L);

    double x,y;

    for(int k=0; k<=L; k++)
    {
        x=x0+k*dx;
        y=y0+k*dy;

        if((type==1 && k==(int)L) || x<0 || x>Width || y<0 || y>Height) break;

        PutWuPixel(x, y, color);
    }
}*/

void dImage::DrawDDALine(double x0,double y0,double x1,double y1,const pixel& color,unsigned char type)
{
    double dx = (x1-x0);
    double dy = (y1-y0);

    if(!_Image_ClipLine(this, x0, y0, x1, y1, dx, dy, color)) return;

    unsigned L = (std::abs(dx) > std::abs(dy)) ? unsigned(std::abs(dx) + 1) : unsigned(std::abs(dy) + 1);

    dx = dx/L;
    dy = dy/L;

    double x,y;

    for(unsigned k=0; k<=L; k++)
    {
        x=x0+k*dx;
        y=y0+k*dy;

        ///if((type == 1 && k == L) || x<0 || x>Width || y<0 || y>Height) break;
		if(type == 1 && k == L) break;

        PutWuPixel(x, y, color);
    }
}

void dImage::DrawWuLine(double x0, double y0, double x1, double y1, pixel color,unsigned char type)
{
    double dy=(y1-y0), dx=(x1-x0);

    if(!_Image_ClipLine(this, x0, y0, x1, y1, dx, dy, color)) return;

    int xi=(x1>x0) ? 1 : (x1<x0) ? -1 : 0;
    int yi=(y1>y0) ? 1 : (y1<y0) ? -1 : 0;

    double x=round(x0),y=round(y0);

    x1=round(x1);
    y1=round(y1);

    unsigned char A = color.a;

    if(fabs(dy)>fabs(dx))
    {
        while(yi*(y-y1)<=0)
        {
            if(yi*(y+yi-y1)>0)
            {
                x=x1;//round(x);
                y=y1;//round(y);
            }
            else x=dx*(y-y0)/dy+x0;

            color.a = (unsigned char) (A*(x - (int)x));
            if(color.a) PutPixel((int)x+1, (int)y, color);

            color.a = (unsigned char) (A*(1 - x + (int)x));
            if(color.a) PutPixel((int)x, (int)y, color);

            y += yi;

            if(yi == 0 || (type==1 && yi*(y+yi-y1)>0.6)) break;

            //if(x<X0 || x>Width || y<Y0 || y>Height) break;
        }
    }
    else
    {
        while(xi*(x-x1)<=0)
        {
            if(xi*(x+xi-x1)>0)
            {
                x=x1;//round(x);
                y=y1;//round(y);
            }
            else y=dy*(x-x0)/dx+y0;

            color.a = (unsigned char) (A*(y - (int)y));
            if(color.a) PutPixel((int) x, (int)y+1, color);

            color.a = (unsigned char) (A*(1 - y + (int)y));
            if(color.a) PutPixel((int)x, (int)y, color);

            x += xi;

            if(xi == 0 || (type==1 && xi*(x+xi-x1)>0.6)) break;

            if(x<0 || x>Width || y<0 || y>Height) break;
        }
    }
}

void dImage::DrawRect(int x,int y,int x1,int y1,const pixel& color)
{
    if(x>x1)
    {
        x = x+x1;
        x1 = x-x1;
        x = x-x1;
    }

    if(y>y1)
    {
        y = y+y1;
        y1 = y-y1;
        y = y-y1;
    }

    DrawLine(x,y, x1,y,color);
    DrawLine(x1,y+1, x1,y1-1,color);
    DrawLine(x1,y1, x,y1,color);
    DrawLine(x,y1-1, x,y+1,color);
}

/**void Image::drawCircle(int X,int Y,int R,const pixel& color)
{
    int x=R,y=0,d=1-R;

    while(x >= y)
    {
        putPixel(x+X,y+Y,color);
        putPixel(-x+X,y+Y,color);

        if(x!=y)
        {
            putPixel(y+X,x+Y,color);
            putPixel(y+X,-x+Y,color);
        }

        if(y!=0)
        {
            putPixel(x+X,-y+Y,color);
            putPixel(-x+X,-y+Y,color);

            if(x!=y)
            {
                putPixel(-y+X,x+Y,color);
                putPixel(-y+X,-x+Y,color);
            }
        }

        y++;

        if(d<0) d += 2*y+1;
        else
        {
            x--;
            d += 2*(y-x+1);
        }
    }
}*/

void dImage::DrawCircle(int X,int Y,int R,const pixel& color)
{
    int x=0,y=R,d=3-2*R;

    while(x <= y)
    {
        PutPixel(x+X,y+Y,color);
        PutPixel(x+X,-y+Y,color);
        PutPixel(-x+X,-y+Y,color);
        PutPixel(-x+X,y+Y,color);
        PutPixel(y+X,x+Y,color);
        PutPixel(y+X,-x+Y,color);
        PutPixel(-y+X,-x+Y,color);
        PutPixel(-y+X,x+Y,color);

        if(d<0) d=d+4*x+6;
        else
        {
            d=d+4*(x-y)+10;
            y--;
        }

        x++;
    }
}

void dImage::DrawDDACircle(double X,double Y,double R,const pixel& color)
{
    const double M_SQRT2 = 1.414213562373095;

    double x = R*M_SQRT2/2;
    double y = x;

    PutWuPixel(X+x,Y+y, color);
    PutWuPixel(X+x,Y-y, color);
    PutWuPixel(X-x,Y+y, color);
    PutWuPixel(X-x,Y-y, color);

    x--;
    y = sqrt(R*R - x*x);

    while(x>=0)
    {
        PutWuPixel(X+x,Y+y, color);
        PutWuPixel(X-x,Y-y, color);

        PutWuPixel(X+y,Y+x, color);
        PutWuPixel(X-y,Y-x, color);

        PutWuPixel(X+x,Y-y, color);
        PutWuPixel(X-x,Y+y, color);

        PutWuPixel(X+y,Y-x, color);
        PutWuPixel(X-y,Y+x, color);

        x--;
        y = sqrt(R*R - x*x);
    }
}

void dImage::FillCircle(int X,int Y,int R,const pixel& color)
{
    int x=R,y=0,d=1-R;

    while(x >= y)
    {
        DrawLine(x+X,y+Y,-x+X,y+Y,color);

        if(y!=0)
        {
            DrawLine(x+X,-y+Y,-x+X,-y+Y,color);
        }

        y++;

        if(d<0) d += 2*y+1;
        else
        {
            x--;
            d += 2*(y-x+1);
        }
    }

    int _y=y;

    x=R,y=0,d=1-R;

    while(x >= y)
    {
        if(x!=y)
        {
            DrawLine(y+X,x+Y,y+X,Y+_y,color,1);
            if(y!=0) DrawLine(-y+X,x+Y,-y+X,Y+_y,color,1);

            if(x!=0)
            {
                DrawLine(y+X,-x+Y,y+X,Y-_y,color,1);
                if(y!=0) DrawLine(-y+X,-x+Y,-y+X,Y-_y,color,1);
            }
        }

        y++;

        if(d<0) d += 2*y+1;
        else
        {
            x--;
            d += 2*(y-x+1);
        }
    }
}
