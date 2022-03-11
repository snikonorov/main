#ifndef $GRAPHICS
#define $GRAPHICS

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
    #define _WL_WINDOWS

    #include <windows.h>
#else
    #define _WL_LINUX

    #include <X11/Xlib.h>
    #include <X11/Xutil.h>
    #include <X11/keysym.h>
#endif

#include "Bitmap.h"

#include <functional>
#include "TVar.h"

#include <string>
#include <algorithm>

///---------------------------------------------------------------------

namespace sEventCodes
{
    enum sEventCode
    {
        Null = 0,

        Mouse_Move,
        Mouse_LeftClick,
        Mouse_RightClick,
        Mouse_DoubleClick,
        Mouse_Wheel,
        Mouse_LeftUp,
        Mouse_RightUp,

        Keyboard_Down,
        Keyboard_Char,
        Keyboard_Up,

        Window_Create,
        Window_Show,
        Window_Paint,
        Window_Activate,
        Window_SetFocus,
        Window_KillFocus,
        Window_Move,
        Window_Size,
        Window_Close,

        Unknown
    };
}

namespace sEventNames
{
    extern const char* sEventName[];
}

struct sEvent
{
    sEventCodes::sEventCode code;

    #ifdef _WL_WINDOWS

    MSG NativeEvent;

    #else

    #endif

    sEvent(sEventCodes::sEventCode Code = sEventCodes::Null): code(Code)
                                                            #ifdef _WL_WINDOWS
                                                            , NativeEvent{0, 0, 0, 0, 0, {0, 0}}
                                                            #endif
    {
    }

    operator const sEventCodes::sEventCode&() const
    {
        return code;
    }

    sEvent& operator=(const sEventCodes::sEventCode& Code)
    {
        code = Code;
        return *this;
    }

    bool operator==(const sEventCodes::sEventCode& Code) const
    {
        return code == Code;
    }

    bool operator!=(const sEventCodes::sEventCode& Code) const
    {
        return code != Code;
    }
};

///---------------------------------------------------------------------

struct TStyle
{
    pixel color;

    TStyle() = default;

    TStyle(const pixel& c): color(c)
    {
    }
};

struct wStyle : TStyle
{
    enum wStyleCode
    {
        NORMAL = 0,
        POPUP = 1
    };

    TVar<wStyleCode> code;
    TVar<bool> DynamicRedraw;

    wStyle(const pixel& c = pixel(255, 255, 255), bool dyn_redraw = true, wStyleCode code = NORMAL);

    #ifdef _WL_WINDOWS

    DWORD NativeStyle;

    operator DWORD() const
    {
        return NativeStyle;
    }

    #else

    #endif

	bool operator==(const wStyleCode& Code) const
	{
		return code == Code;
	}

	bool operator!=(const wStyleCode& Code) const
	{
		return code != Code;
	}
};

///---------------------------------------------------------------------

struct TComponent
{
    virtual void ProcessEvent(sEvent& event) = 0;

    virtual ~TComponent()
    {
    }
};

///---------------------------------------------------------------------

class sWindow
{
	protected:

    sWindow();

    sWindow(const sWindow& w) = default;
    ~sWindow();

	public:

    struct _internal
    {
        unsigned char* tdata;
        unsigned char line;
        double X, Y;

        #ifdef _WL_WINDOWS

        HWND hwnd;

        HDC dc;
        PAINTSTRUCT ps;
        BITMAPV4HEADER bi;

        #endif

        _internal();
        ~_internal();
    };

    struct MouseInput
    {
        int X, Y;
        char Wheel;
        bool LeftPressed;
        bool RightPressed;
        bool DoubleClicked;

        MouseInput();

		Coord Pos() const
		{
			return {X, Y};
		}
    };

    struct KeyInput
    {
        unsigned char Char;
        unsigned short Key;
        unsigned short Layout;

        bool Shift, Control;

        KeyInput();
    };

    _internal _Internals;
    KeyInput Key;
    MouseInput Mouse;

    dImage BackBuffer;

    Vector3d<TVar<int>> Size;
    Vector3d<TVar<int>> Pos;

    wStyle style;

    std::function<void(sEvent& event)> ProcessEvent;

    #ifdef _WL_WINDOWS
		///const wchar_t* Caption;
		TVar<std::wstring> Caption;
	#else
		const char* Caption;
		//Display* display;
		//Window window;
	#endif

	std::vector<TComponent*> Components;

	///--------------------------------------------------------------------

	void Add(TComponent* cmp)
	{
		Components.push_back(cmp);
	}

    void operator+=(TComponent* cmp)
	{
        Add(cmp);
	}

	template<class TLambda>
	auto Add(const TLambda& f) -> decltype(f(std::declval<sEvent&>()), void())
	{
		struct Cmp : TComponent
		{
			TLambda f;

			Cmp(const TLambda& f): f(f)
			{
			}

			void ProcessEvent(sEvent& event) override
			{
				f(event);
			}
		};

		Components.push_back(new Cmp(f));
	}

	template<class T>
	void operator+=(T&& cmp)
	{
        Add(std::forward<T>(cmp));
	}

	int Show();

	static sWindow* GetInstance()
	{
		static sWindow* w = new sWindow();
		return w;
	}
};

///---------------------------------------------------------------------------------------------------------------

extern sWindow* MainWindow;

extern decltype(MainWindow->Size)& AppSize;
extern decltype(MainWindow->Pos)& AppPos;

extern const bool& isLeftPressed;
extern const bool& isRightPressed;
extern const bool& isDoubleClicked;
extern const int& MouseX;
extern const int& MouseY;
extern const char& MouseWheel;

extern const sWindow::MouseInput& sMouse;
extern const sWindow::KeyInput& sKey;
extern sWindow::_internal& _sInternals;

///---------------------------------------------------------------------------------------------------------------

inline void SetClipRect(iRect r)
{
    if(r.A.x < 0) r.A.x = 0;
    if(r.A.y < 0) r.A.y = 0;

    if(r.B.x >= MainWindow->BackBuffer.Width)
        r.B.x = MainWindow->BackBuffer.Width - 1;

    if(r.B.y >= MainWindow->BackBuffer.Height)
        r.B.y = MainWindow->BackBuffer.Height - 1;

    MainWindow->BackBuffer.ClipRect = r;
}

inline void ResetClipRect()
{
    MainWindow->BackBuffer.ClipRect = iRect(MainWindow->BackBuffer.Width-1, MainWindow->BackBuffer.Height-1);
}

inline void Clear()
{
    MainWindow->BackBuffer.Fill(MainWindow->style.color);
}

inline void Refresh()
{
    InvalidateRect(MainWindow->_Internals.hwnd, NULL, false);
}

inline void jDrawWuLine(double x, double y, double x1, double y1, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.DrawWuLine(x,y,x1,y1,color,1);
}

inline void jDrawDDALine(double x, double y, double x1, double y1, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.DrawDDALine(x,y,x1,y1,color,1);
}

inline void SaveScreen()
{
    //if(MainWindow->BackBuffer.data)
    //{
        if(!_sInternals.tdata)
            _sInternals.tdata = new unsigned char[MainWindow->BackBuffer.Cnst];

        memcpy(_sInternals.tdata, MainWindow->BackBuffer.data.data(), MainWindow->BackBuffer.Cnst);
    //}
}

inline void PaintSaved()
{
    if(_sInternals.tdata)   ///MainWindow->BackBuffer.data &&
    {
        memcpy(MainWindow->BackBuffer.data.data(), _sInternals.tdata,MainWindow->BackBuffer.Cnst);
    }
}

inline void PutPixel(int x, int y, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.PutPixel(x, y, color);
}

inline void PutWuPixel(double x, double y, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.PutWuPixel(x, y, color);
}

inline pixel GetPixel(int x, int y)
{
    return MainWindow->BackBuffer.GetPixel(x,y);
}

inline void DrawLine(int x, int y, int x1, int y1, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.DrawLine(x, y, x1, y1, color);
}

inline void DrawWuLine(double x, double y, double x1, double y1, const pixel& color = pixel(0,0,0), unsigned char type = 0)
{
    MainWindow->BackBuffer.DrawWuLine(x, y, x1, y1, color, type);
}

inline void DrawDDALine(double x, double y, double x1, double y1, const pixel& color = pixel(0,0,0), unsigned char type = 0)
{
    MainWindow->BackBuffer.DrawDDALine(x, y, x1, y1, color, type);
}

inline void DrawCircle(int x, int y, int r, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.DrawCircle(x, y, r, color);
}

inline void DrawDDACircle(double x, double y, double r, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.DrawDDACircle(x, y, r, color);
}

inline void FillCircle(int x, int y, int r, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.FillCircle(x, y, r, color);
}

inline void MoveTo(double x, double y)
{
    _sInternals.X = x;
    _sInternals.Y = y;
}

inline void Line(double x, double y, const pixel& color = pixel(0,0,0))
{
    switch(_sInternals.line)
    {
        case 0:
            DrawLine((int) round(_sInternals.X), (int) round(_sInternals.Y), (int) round(x), (int) round(y), color);
            break;

        case 1:
            jDrawWuLine(_sInternals.X, _sInternals.Y, x, y, color);
            break;

        default:
            jDrawDDALine(_sInternals.X, _sInternals.Y, x, y, color);
    }
}

inline void Line(double x0, double y0, double x1, double y1, const pixel& color = pixel(0,0,0))
{
    switch(_sInternals.line)
    {
        case 0:
            DrawLine((int) round(x0), (int) round(y0), (int) round(x1), (int) round(y1), color);
            break;

        case 1:
            jDrawWuLine(x0, y0, x1, y1, color);
            break;

        default:
            jDrawDDALine(x0, y0, x1, y1, color);
    }
}

inline void RLine(double dx, double dy, const pixel& color = pixel(0,0,0))
{
    Line(_sInternals.X + dx, _sInternals.Y + dy, color);
}

inline void LineTo(double x, double y, const pixel& color = pixel(0,0,0))
{
    Line(x, y, color);
    MoveTo(x, y);
}

inline void RLineTo(double dx, double dy, const pixel& color = pixel(0,0,0))
{
    Line(_sInternals.X + dx, _sInternals.Y + dy, color);
    MoveTo(_sInternals.X + dx, _sInternals.Y + dy);
}

inline void SetLineType(unsigned char t)
{
    _sInternals.line = t;
}

inline void DrawRect(int x, int y, int x1, int y1, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.DrawRect(x, y, x1, y1, color);
}

inline void FillRect(int x, int y, int x1, int y1, const pixel& color = pixel(0,0,0))
{
    MainWindow->BackBuffer.FillRect(x, y, x1, y1, color);
}

inline void Fill(const pixel& color)
{
    MainWindow->BackBuffer.Fill(color);
}

///inline void FillAlpha(const pixel& color)
///{
///    MainWindow->BackBuffe.FillAlpha(color);
///}

template<class Img>
inline void DrawImage(Img& img, int x, int y)
{
    MainWindow->BackBuffer.DrawImage(img, x, y);
}

template<class Img, class T>
inline void DrawImage(Img& img, const Vector3d<T>& p)
{
    MainWindow->BackBuffer.DrawImage(img, p.x, p.y);
}

///----------------------------------------------------------------------------------------------

template<class T>
inline void PutPixel(const Vector3d<T>& p, const pixel& color = pixel(0,0,0))
{
    PutPixel(p.x, p.y, color);
}

template<class T>
inline void PutWuPixel(const Vector3d<T>& p, const pixel& color = pixel(0,0,0))
{
    PutWuPixel(p.x, p.y, color);
}

template<class T>
inline pixel GetPixel(const Vector3d<T>& p)
{
    return GetPixel(p.x, p.y);
}

template<class T1, class T2>
inline void DrawLine(const Vector3d<T1>& p1, const Vector3d<T2>& p2, const pixel& color = pixel(0,0,0))
{
    DrawLine(p1.x, p1.y, p2.x, p2.y, color);
}

template<class T1, class T2>
inline void DrawWuLine(const Vector3d<T1>& p1, const Vector3d<T2>& p2, const pixel& color = pixel(0,0,0), unsigned char type = 0)
{
    DrawWuLine(p1.x, p1.y, p2.x, p2.y, color, type);
}

template<class T1, class T2>
inline void DrawDDALine(const Vector3d<T1>& p1, const Vector3d<T2>& p2, const pixel& color = pixel(0,0,0), unsigned char type = 0)
{
    DrawDDALine(p1.x, p1.y, p2.x, p2.y, color, type);
}

template<class T1, class T2>
inline void DrawInfiniteLine(const Vector3d<T1>& p1, const Vector3d<T2>& p2, const pixel& color = pixel(0,0,0))
{
    ///const auto& W = MainWindow->BackBuffer.Width;
    ///const auto& H = MainWindow->BackBuffer.Height;

    const auto& W = MainWindow->Size.x;
    const auto& H = MainWindow->Size.y;

    if(p1.x == p2.x)
	{
		DrawDDALine(p1.x, 0, p1.x, H-1, color);
	}
	else if(p1.y == p2.y)
	{
		DrawDDALine(0, p1.y, W-1, p1.y, color);
	}
	else
	{
		double u, q;
		Point dp = p2 - p1;

		static std::vector<double> t(4);
		t.clear();

		u = -p1.y/dp.y;
		q = (p1.x + u*dp.x)/(W - 1);

		if(q >= 0 && q <= 1)
        {
            t.push_back(u);
        }

		u = (H-1 - p1.y)/dp.y;
		q = (p1.x + u*dp.x)/(W - 1);

		if(q >= 0 && q <= 1)
        {
            t.push_back(u);
        }

		u = -p1.x/dp.x;
		q = (p1.y + u*dp.y)/(H - 1);

		if(q >= 0 && q <= 1)
        {
            t.push_back(u);
        }

		u = (W-1 - p1.x)/dp.x;
		q = (p1.y + u*dp.y)/(H - 1);

		if(q >= 0 && q <= 1)
        {
            t.push_back(u);
        }

		if(t.size() > 2)
		{
			std::sort(t.begin(), t.end());
			t[1] = t.back();
			t.resize(2);
		}

        DrawDDALine(p1 + dp*t[0], p1 + dp*t[1], color);
	}
}

template<class T1, class T2>
inline void DrawArrow(const Vector3d<T1>& p1, const Vector3d<T2>& p2, const pixel& color = pixel(0,0,0))
{
	const double M_PI = 3.141592653589793;

	const double a = 8;
    const double alpha = 20*M_PI/180;

    const double sa = sin(alpha);
    const double ca = cos(alpha);

    auto l = Point(p2 - p1).normal();
    auto n = l.ortho();

    auto p3 = p2 - a*l*ca;
    auto sn = a*n*sa;

    DrawDDALine(p1, p2, color);

    DrawDDALine(p2, p3 + sn, color);
    DrawDDALine(p2, p3 - sn, color);
}

template<class T>
inline void DrawCircle(const Vector3d<T>& p, int r, const pixel& color = pixel(0,0,0))
{
    DrawCircle(p.x, p.y, r, color);
}

template<class T>
inline void DrawDDACircle(const Vector3d<T>& p, double r, const pixel& color = pixel(0,0,0))
{
    DrawDDACircle(p.x, p.y, r, color);
}

template<class T>
inline void FillCircle(const Vector3d<T>& p, int r, const pixel& color = pixel(0,0,0))
{
    FillCircle(p.x, p.y, r, color);
}

template<class T>
inline void MoveTo(const Vector3d<T>& p)
{
    _sInternals.X = p.x;
    _sInternals.Y = p.y;
}

template<class T>
inline void Line(const Vector3d<T>& p, const pixel& color = pixel(0,0,0))
{
    Line(p.x, p.y, color);
}

template<class T1, class T2>
inline void Line(const Vector3d<T1>& p0, const Vector3d<T2>& p1, const pixel& color = pixel(0,0,0))
{
    Line(p0.x, p0.y, p1.x, p1.y, color);
}

template<class T>
inline void RLine(const Vector3d<T>& dp, const pixel& color = pixel(0,0,0))
{
    Line(_sInternals.X + dp.x, _sInternals.Y + dp.y, color);
}

template<class T>
inline void LineTo(const Vector3d<T>& p, const pixel& color = pixel(0,0,0))
{
    Line(p, color);
    MoveTo(p);
}

template<class T>
inline void RLineTo(const Vector3d<T>& dp, const pixel& color = pixel(0,0,0))
{
    Line(_sInternals.X + dp.x, _sInternals.Y + dp.y, color);
    MoveTo(_sInternals.X + dp.x, _sInternals.Y + dp.y);
}

template<class T1, class T2>
inline void DrawRect(const Vector3d<T1>& p1, const Vector3d<T2>& p2, const pixel& color = pixel(0,0,0))
{
    DrawRect(p1.x, p1.y, p2.x, p2.y, color);
}

template<class T>
inline void DrawRect(const TRect<T>& r, const pixel& color = pixel(0,0,0))
{
    DrawRect(r.A.x, r.A.y, r.B.x, r.B.y, color);
}

template<class T1, class T2>
inline void FillRect(const Vector3d<T1>& p1, const Vector3d<T2>& p2, const pixel& color = pixel(0,0,0))
{
    FillRect(p1.x, p1.y, p2.x, p2.y, color);
}

template<class T>
inline void FillRect(const TRect<T>& r, const pixel& color = pixel(0,0,0))
{
    FillRect(r.A.x, r.A.y, r.B.x, r.B.y, color);
}

template<class TContainer>
inline void DrawPoints(const TContainer& P, const pixel& color)
{
	const int r = 2;

	for(const auto& p: P)
	{
		FillCircle(p, r, color/2);
		DrawDDACircle(p, r, color);
	}
}

template<class TContainer>
inline void DrawPoints(const TContainer& P, const double r = 2, const pixel& color = pixel(0, 0, 0))
{
	for(const auto& p: P)
	{
		FillCircle(p, r, color/2);
		DrawDDACircle(p, r, color);
	}
}

template<class TContainer>
inline void DrawPolyline(const TContainer& data, const pixel& color = pixel(0, 0, 0))
{
    if(data.size())
    {
        ///SetLineType(2);
        MoveTo(data[0]);

        for(unsigned k = 1; k < data.size(); k++)
        {
            LineTo(data[k], color);
        }
    }
}

///-----------------------------------------------------------------------------------------------

iRect GetDesktopRect();

decltype(MainWindow) SCreateWindow(int width, int height,
                                    #ifdef _WL_WINDOWS
                                        const wchar_t* caption = L"SGraphicsApp",
                                    #else
                                        const char* caption = "SGraphicsApp",
                                    #endif
									const wStyle& style = wStyle{pixel(255, 255, 255), true, wStyle::NORMAL}
                                   );

inline decltype(MainWindow) SCreateWindow(const Coord& Size,
											#ifdef _WL_WINDOWS
												const wchar_t* caption = L"SGraphicsApp",
											#else
												const char* caption = "SGraphicsApp",
											#endif
											const wStyle& style = wStyle{pixel(255, 255, 255), true, wStyle::NORMAL}
										 )
{
    return SCreateWindow(Size.x, Size.y, caption, style);
}

inline decltype(MainWindow) SCreateWindow(
											#ifdef _WL_WINDOWS
                                                const wchar_t* caption=L"SGraphicsApp",
                                            #else
                                                const char* caption = "SGraphicsApp",
                                            #endif
											const wStyle& style = wStyle{pixel(255, 255, 255), true, wStyle::POPUP}
                                          )
{
	iRect rct = GetDesktopRect();
	return SCreateWindow(rct.Width(), rct.Height(), caption, style);
}

///Windows
#ifdef _WL_WINDOWS

inline void quit(int code = 0)
{
    PostQuitMessage(code);
}

///Linux
#else

//int WindowProcedure(XEvent* message);

#endif

#endif
