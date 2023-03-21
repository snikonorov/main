#define _CRT_SECURE_NO_WARNINGS

#include "TChart.h"
#include "to_string.h"

///--------------------------------------------------------------------------------------

//sFont* TChart::_font_default = 0;
sFont* TChart::_font_default = new sFont(12);

///--------------------------------------------------------------------------------------

inline double AxisScale(double x, int l)
{
    static int _nx[3]= {2,5,10};
    int _i=0;
    double res=1;

    double k=1;

    while(x<l)
    {
        x*=10;
        k*=10;
    }

    while(x>l)
    {
        if(_i>0)
        {
            x*=_nx[_i-1];
            res/=_nx[_i-1];
        }

        x/=_nx[_i];
        res*=_nx[_i++];

        if(_i==3) _i=0;
    }

    return 1.*res/k;
}

///--------------------------------------------------------------------------------------

void Graph::Draw(const TChart* chart) const
{
	if(!size() || !style.Visible) return;

    const auto& gC = style.color;

    const auto& data = this->data();

	auto r2s = [chart](const auto& p)
	{
		return chart->RealToScreen(p);
	};

	if(style.code & Graph::InfiniteLines)
	{
		for(unsigned p = 0; p < size(); p += 2)
		{
			///DrawInfiniteLine(O + K*data[p], O + K*data[p+1], gC);
			DrawInfiniteLine(r2s(data[p]), r2s(data[p+1]), gC);
		}

		return;
	}

    unsigned p = 0;
    while(p < size() && (std::isnan(data[p].x) || std::isnan(data[p].y)))
		p++;

    ///Point p0 = O + K*data[p];
	Point p0 = r2s(data[p]);

	while(p0.x < 0 && p+1 < size())
	{
		///p0 = O + K*data[++p];
		p0 = r2s(data[++p]);
	}

	if(p)
	{
		///p0 = O + K*data[--p];
		p0 = r2s(data[--p]);
	}

    double& x0 = p0.x;
    double& y0 = p0.y;

    double _x = x0;
    double _y = y0;

    for(; p < size(); p++)
    {
        if(std::isnan(data[p].x) || std::isnan(data[p].y))
		{
			_x = nan("");
			continue;
		}

		if(p+1 < size())
		{
			///Point p_ = O + K*data[p+1];
			Point p_ = r2s(data[p+1]);

			double& x_ = p_.x;
			double& y_ = p_.y;

			if(!std::isnan(_x) && (style.code & Graph::Lines))
			{
				if(std::abs(x_ - _x) >= 1.25 || std::abs(y_ - _y) >= 1.25)
				{
					if(style.code & Graph::Bresenham) DrawLine((int)_x, (int)_y, (int)x_, (int)y_, gC);
					else if(style.code & Graph::Wu) jDrawWuLine(_x, _y, x_, y_,gC);
					else jDrawDDALine(_x, _y, x_, y_,gC);

					_x = x_;
					_y = y_;
				}
			}
			else
			{
				_x = x_;
				_y = y_;
			}
		}

		if(style.code & Graph::Points)
		{
			//x0=round(x0);
			//y0=round(y0);

			const int r = 2;

			FillCircle((int) round(x0), (int) round(y0), r, style.color/2);

			if(style.code & Graph::Bresenham)
				DrawCircle((int) round(x0), (int) round(y0), r, style.color);
			else
				DrawDDACircle(x0, y0, r, style.color);
		}

		x0 = _x;
		y0 = _y;

		///if(x0 > AppSize.x)
			///break;
    }
}

///--------------------------------------------------------------------------------------

TChart::TChart(): TChart(TChart::_font_default)
{
}

TChart::TChart(sFont* font, const Coord& Pos, const Coord& Size, unsigned style, const Point& k, const Point& offset, const Point& cr)
               : K(k), Offset(offset), Cr(cr), font(font), style(style), Pos(Pos), Size(Size)
{
    _hx = 1/AxisScale(K.x, _l);
    _hy = 1/AxisScale(K.y, _l);
}

TChart::~TChart()
{
    for(auto q = data.begin(); q != data.end(); q++)
    {
        delete (*q);
    }

    font = 0;
}

void TChart::Scale(const Point& a)
{
    K *= a;
    Offset *= a;

    _hx = 1/AxisScale(K.x, _l);
    _hy = 1/AxisScale(K.y, _l);
}

void TChart::pScale(const Point& a, const Point& p0)
{
    K *= a;

	const int& Width = Size.x < 0 ? AppSize.x : Size.x;
	const int& Height = Size.y < 0 ? AppSize.y : Size.y;
	const auto& Sz = Point(Width/2, Height/2);

    Offset = Offset*a + (1 - a)*(p0 - Sz);

    _hx = 1/AxisScale(K.x, _l);
    _hy = 1/AxisScale(K.y, _l);
}

void TChart::ProcessEvent(sEvent& event)
{
    static Coord _T;

    const int& X0 = Pos.x;
    const int& Y0 = Pos.y;
    const int& Width = Size.x < 0 ? AppSize.x : Size.x;
    const int& Height = Size.y < 0 ? AppSize.y : Size.y;

	static char _s[128];

    ///---------------------------------------------------------------

    if(!style) return;

    using namespace sEventCodes;

    if(event == Mouse_LeftClick)
    {
        if((style & Movable) && MouseX >= X0 && MouseX <= X0 + Width && MouseY >= Y0 && MouseY <= Y0 + Height)
        {
            _T = Offset - Coord(MouseX, MouseY);
        }
    }
    else if(event == Mouse_Move)
    {
		if(((isLeftPressed && (style & Movable)) || (isRightPressed && (style & Coords))) && MouseX >= X0 && MouseX <= X0 + Width && MouseY >= Y0 && MouseY <= Y0 + Height)
		{
			if(isLeftPressed)
				Offset = _T + Coord(MouseX, MouseY);

			if(!MainWindow->style.DynamicRedraw)
				InvalidateRect(MainWindow->_Internals.hwnd, 0, false);
		}
    }
	else if((event == Mouse_RightClick || event == Mouse_RightUp) && (style & Coords))
	{
		if (!MainWindow->style.DynamicRedraw)
			InvalidateRect(MainWindow->_Internals.hwnd, 0, false);
	}
    else if(event == Mouse_Wheel && (style & Scalable))
    {
        pScale(Point(sKey.Key != 'Y' ? pow(1.1, sMouse.Wheel) : 1,
                     sKey.Key != 'X' ? pow(1.1, sMouse.Wheel) : 1),
               Point(MouseX, MouseY));

		if (!MainWindow->style.DynamicRedraw)
			InvalidateRect(MainWindow->_Internals.hwnd, 0, false);
    }
    else if(event == Window_Paint)
    {
        ///Crop------------------------------------------

        SetClipRect({Pos.x, Pos.y, Width, Height});

        ///----------------------------------------------

        ///const auto Sz = Coord(Width/2, Height/2);
        //const auto Sz = AppSize/2;

        ///Plots.Grid---------------------------------------

        if(style & Grid)
        {
			double cnt, _v, _a;

			DrawLine(X0, (int) round(Offset.y + Height/2.), X0 + Width, (int) round(Offset.y + Height/2.), pixel(100,0,0,0));
            DrawLine((int) round(Offset.x + Width/2.), 0, (int) round(Offset.x + Width/2.), Y0 + Height, pixel(100,0,0,0));

            cnt = Width/2 + Offset.x;

            if(_hx>0)
            {
                _a = round((X0 - cnt)/(K.x*_hx));

                cnt += _a*K.x*_hx;
                _v = _a*_hx;

                _s[0]='0';
                _s[1]=0;

				if(font) Text_out(font, (int)round(Offset.x + Width/2. - (*font)['0'].Width - 5), (int)round(Height/2. + Offset.y + 5), _s);

                while(cnt <= X0 + Width)
                {
                    if(abs(cnt - Offset.x - Width/2)<2)
                    {
                        _s[0]='0';
                        _s[1]=0;
                        _v=0;

                        cnt = Offset.x + Width/2;
                    }
                    else
                    {
                        if(_hx<1e-4 && (int)_v!=0) sprintf(_s,".%d",(int)_v);
                        else sprintf(_s,"%g",(double)_v);

                        DrawLine((int)round(cnt), 0, (int)round(cnt), Y0 + Height, pixel(tr,0,20,50));
                    }

                    if(font && (*_s!='0' || *(_s+1)))
                    {
						const auto& dash = (*font)['-'];

						if(_s[0]=='.')
                        {
							Text_out(font, int(cnt - (Text_Width(font,_s+1) + (_v<0)*dash.Width)/2), int(Height/2 + Offset.y + 5), _s+1);

                            sprintf(_s,"%c%.8g",_v>0 ? '+' : '-',(double)fabs(_v - (int)_v));
                            Text_out(font, int(cnt - (Text_Width(font,_s) + (_v<0)*dash.Width)/2), int(Height/2 + Offset.y + 5 + font->Height+2) ,_s);
                        }
                        else Text_out(font, int(cnt - (Text_Width(font,_s) + (_v<0)*dash.Width)/2), int(Height/2 + Offset.y + 5), _s);
                    }

                    DrawLine((int)round(cnt), int(Offset.y + Height/2 - 4), (int)round(cnt), int(Offset.y + Height/2 + 4), pixel(100,0,0,0));

                    cnt+=K.x*_hx;
                    _v+=_hx;
                }
            }

            cnt = Height/2+Offset.y;

            if(_hy>0)
            {
                _a = round((Y0 - cnt)/(K.y*_hy));

                cnt += _a*K.y*_hy;
                _v = _a*_hy;

                while(cnt <= Y0 + Height)
                {
                    if(abs(cnt - Offset.y - Height/2)<2)
                    {
                        _s[0]='0';
                        _s[1]=0;
                        _v=0;

                        cnt = Offset.y + Height/2;
                    }
                    else
                    {
                        if(_hy<1e-4 && (int)_v!=0) sprintf(_s,"%d%c%.8g",-(int)_v,-_v>0 ? '+' : '-',(double)fabs(_v - (int)_v));
                        else sprintf(_s,"%g",-(double)_v);

                        DrawLine(0, (int)cnt, X0 + Width, (int)cnt, pixel(tr,0,20,50));
                    }

                    if(font && (*_s!='0' || *(_s+1))) Text_out(font, int(Width/2 + Offset.x - Text_Width(font,_s) - 5), (int)round(cnt - font->Height/2. + 1), _s);

                    DrawLine(int(Offset.x + Width/2 - 4), (int)round(cnt), int(Offset.x + Width/2 + 4), (int)round(cnt), pixel(100,0,0,0));

                    cnt+=K.y*_hy;
                    _v+=_hy;
                }
            }
        }

        ///-----------------------------------------------------------------

		if (style & Coords)
		{
			if (isRightPressed && font)
			{
				///Drawing point-coordinates if necessary...

				sprintf(_s, "%g %g", (double)(MouseX-Width/2-Offset.x)/K.x+Cr.x, (double)(Height/2-MouseY+Offset.y)/K.y+Cr.y);
				Text_out(font, MouseX, MouseY - font->Height, _s);
			}
		}

        ///Plots.Drawing----------------------------------------------------

        const auto& N = data.size();

        for(unsigned k=0; k<N; k++)
        {
            data[k]->Draw(this);
        }

        ///Uncrop--------------------------------------------------------------------------------------

        ResetClipRect();

        ///--------------------------------------------------------------------------------------------

        ///DrawRect(Pos, Pos + Coord{Width, Height}, pixel(50));
    }
}

std::string TChart::SVG() const
{
	auto tab = [](unsigned short n)
	{
		return std::string(n, '\t');
	};

	/// render Point to std::string
	///
	auto p2s = [](const Point& p)
	{
		return std::string("(") + to_string(p.x) + ", " + to_string(p.y) + ")";
	};

	std::string res = "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\"?>\n";

	res += "<svg id=\"svg\" width=\"100%\" height=\"100%\" xmlns=\"http://www.w3.org/2000/svg\">\n";

	/// ! `Size` may be partially const !

	res += tab(1) + "<g id=\"chart_size_offset\">\n";
	res += tab(2) + "<g id=\"chart_offset\" transform=\"translate" + p2s(Offset) + "\">\n";
	res += tab(3) + "<g id=\"chart_scale\" transform=\"scale" + p2s(Point(K.x, -K.y)) + "\">\n";

	for(unsigned k = 0; k < data.size(); k++)
	{
		//auto g = data[k]->SVG();
		//return g.to_string();
	}

	res += tab(3) + "</g>\n";
	res += tab(2) + "</g>\n";
	res += tab(1) + "</g>\n";

	res += "</svg>\n";

	return res;
}
