#include "Graphics.h"

///------------------------------------------------------------------------------------------

sWindow* MainWindow = sWindow::GetInstance();

decltype(MainWindow->Size)& AppSize = MainWindow->Size;
decltype(MainWindow->Pos)& AppPos = MainWindow->Pos;

const bool& isLeftPressed = MainWindow->Mouse.LeftPressed;
const bool& isRightPressed = MainWindow->Mouse.RightPressed;
const bool& isDoubleClicked = MainWindow->Mouse.DoubleClicked;
const int& MouseX = MainWindow->Mouse.X;
const int& MouseY = MainWindow->Mouse.Y;
const char& MouseWheel = MainWindow->Mouse.Wheel;

const sWindow::MouseInput& sMouse = MainWindow->Mouse;
const sWindow::KeyInput& sKey = MainWindow->Key;
sWindow::_internal& _sInternals = MainWindow->_Internals;

///------------------------------------------------------------------------------------------

///setters
namespace
{
    template<class T>
    T& declval();

    ///----------------------------------------------------------------------------

    struct WHSetter
    {
        typedef $SubType(declval<sWindow>().Size.x) T;
        sWindow* w;

        WHSetter(sWindow* W): w(W)
        {
        }

        void operator()(IVar<T>* obj, const T& v) const
        {
            if((T)(*obj) != v)
            {
                //TVar<T>::DefSetValue(obj, v);
                obj->GetRef() = v;

				SetWindowPos(w->_Internals.hwnd, 0, 0, 0, (T) w->Size.x, (T) w->Size.y, SWP_NOMOVE | SWP_NOSENDCHANGING | SWP_NOZORDER);
            }
        }
    };

    ///----------------------------------------------------------------------------

    struct PosSetter
    {
        typedef $SubType(declval<sWindow>().Pos.x) T;
        sWindow* w;

        PosSetter(sWindow* W): w(W)
        {
        }

        void operator()(IVar<T>* obj, const T& v) const
        {
            if((T)(*obj) != v)
            {
                TVar<T>::DefSetValue(obj, v);
                SetWindowPos(w->_Internals.hwnd, 0, (T) w->Pos.x, (T) w->Pos.y, 0, 0, SWP_NOSIZE | SWP_NOSENDCHANGING | SWP_NOZORDER);
            }
        }
    };

    ///----------------------------------------------------------------------------

    struct StyleSetter
    {
        typedef bool T;
        sWindow* w;

        StyleSetter(sWindow* W): w(W)
        {
        }

        void operator()(IVar<T>* obj, const T& v) const
        {
            TVar<T>::DefSetValue(obj, v);
            InvalidateRect(w->_Internals.hwnd, NULL, false);
        }
    };

    ///----------------------------------------------------------------------------

    struct StyleCodeSetter
    {
        typedef wStyle::wStyleCode T;
        sWindow* w;

        StyleCodeSetter(sWindow* W): w(W)
        {
        }

        void operator()(IVar<T>* obj, const T& v) const
        {
            TVar<T>::DefSetValue(obj, v);

            SetWindowLongPtr(w->_Internals.hwnd, GWL_STYLE, v == wStyle::NORMAL ? WS_OVERLAPPEDWINDOW : WS_POPUP);
            SetWindowPos(w->_Internals.hwnd, 0, 0,0,0,0, SWP_NOMOVE | SWP_NOSIZE | SWP_SHOWWINDOW);
        }
    };

    ///----------------------------------------------------------------------------

    struct CaptionSetter
    {
        typedef std::wstring T;
        sWindow* w;

        CaptionSetter(sWindow* W): w(W)
        {
        }

        void operator()(IVar<T>* obj, const T& v) const
        {
            TVar<T>::DefSetValue(obj, v);
            SetWindowTextW(w->_Internals.hwnd, v.data());
        }
    };

    #undef $SubType
}

///-------------------------------------------------------------------------------

wStyle::wStyle(const pixel& c, bool dyn_redraw, wStyleCode code): TStyle(c)
{
    this->code = code;
    DynamicRedraw = dyn_redraw;

    DynamicRedraw.SetValue = StyleSetter(MainWindow);
    this->code.SetValue = StyleCodeSetter(MainWindow);
}

///------------------------------------------------------------

sWindow::sWindow()///: Size{0, 0}, Pos{0, 0}
{
    //Size = {0, 0};
    //Pos = {0, 0};

    ProcessEvent = 0;

    Size.x.SetValue = WHSetter(this);
    Size.y.SetValue = WHSetter(this);

    Pos.x.SetValue = PosSetter(this);
    Pos.y.SetValue = PosSetter(this);

    #ifdef _WL_WINDOWS
        Caption.SetValue = CaptionSetter(this);
    #else
        //Caption = 0;
        //display = 0;
        //window = 0;
    #endif
}

sWindow::~sWindow()
{
    for(auto q = Components.begin(); q != Components.end(); q++)
    {
        if (*q)
		{
			delete (*q);
			*q = 0;
		}
    }
}

sWindow::_internal::_internal(): tdata(0), line(0), X(0), Y(0)
{
    #ifdef _WL_WINDOWS

	hwnd = 0;
	dc = 0;

    ZeroMemory (&bi, sizeof(BITMAPV4HEADER));
    bi.bV4Size = sizeof(BITMAPINFOHEADER);
    bi.bV4BitCount = 24*sizeof(unsigned char);
    bi.bV4Planes = 1;
    bi.bV4V4Compression = BI_RGB;

    #endif
}

sWindow::_internal::~_internal()
{
    #ifdef _WL_WINDOWS

    if(tdata) delete [] tdata;
    tdata = 0;

    #endif
}

sWindow::MouseInput::MouseInput(): X(0), Y(0), Wheel(0), LeftPressed(false), RightPressed(false), DoubleClicked(false)
{
}

sWindow::KeyInput::KeyInput(): Char(0), Key(0), Layout(0), Shift(false), Control(false)
{
}

int sWindow::Show()
{
    #ifdef _WL_WINDOWS

    ShowWindow(_Internals.hwnd, SW_SHOWNORMAL);

    MSG messages;
    while (GetMessage (&messages, NULL, 0, 0))
    {
        TranslateMessage(&messages);
        DispatchMessage(&messages);
    }

    return (int) messages.wParam;

    #else

    #endif
}

///------------------------------------------------------------------------------------------

namespace sEventNames
{
    const char* sEventName[] = {
        "<Null>",

        "Mouse_Move",
        "Mouse_LeftClick",
        "Mouse_RightClick",
        "Mouse_DoubleClick",
        "Mouse_Wheel",
        "Mouse_LeftUp",
        "Mouse_RightUp",

        "Keyboard_Down",
        "Keyboard_Char",
        "Keyboard_Up",

        "Window_Create",
        "Window_Show",
        "Window_Paint",
        "Window_Activate",
        "Window_SetFocus",
        "Window_KillFocus",
        "Window_Move",
        "Window_Size",
        "Window_Close",

        "<Unknown>"
    };
}

///------------------------------------------------------------------------------------------

/// Event receiver
namespace
{
    #ifdef _WL_WINDOWS

    BYTE allKeys[256];

    ///--------------------------------------------------------------------------------------

    inline void SFlush()
    {
        if(_sInternals.dc)
        {
            SetDIBitsToDevice(_sInternals.dc, 0,0, MainWindow->BackBuffer.Width, MainWindow->BackBuffer.Height,
                                              0,0, 0, MainWindow->BackBuffer.Height,
                                              MainWindow->BackBuffer.data.data(),
                                              (BITMAPINFO*)&_sInternals.bi,
                                              DIB_RGB_COLORS);

            ReleaseDC(_sInternals.hwnd, _sInternals.dc);
            _sInternals.dc = 0;

            if(!MainWindow->style.DynamicRedraw) EndPaint(_sInternals.hwnd, &_sInternals.ps);
        }
    }

    inline sEvent TranslateNativeMessage(HWND& hwnd, UINT& message, WPARAM& wParam, LPARAM& lParam)
    {
        using namespace sEventCodes;

        sEvent event;

        event.NativeEvent.hwnd = hwnd;
        event.NativeEvent.message = message;
        event.NativeEvent.lParam = lParam;
        event.NativeEvent.wParam = wParam;

        WINDOWPOS* q;

        switch(message)
        {
            case WM_MOUSEMOVE:
                event = Mouse_Move;

                MainWindow->Mouse.X = LOWORD(lParam);
                MainWindow->Mouse.Y = HIWORD(lParam);

                event.NativeEvent.pt = {MainWindow->Mouse.X, MainWindow->Mouse.Y};
                break;

            case WM_MOUSEWHEEL:
                event = Mouse_Wheel;

                MainWindow->Mouse.Wheel = char(GET_WHEEL_DELTA_WPARAM(wParam) / WHEEL_DELTA);
                break;

            case WM_LBUTTONDOWN:
                event = Mouse_LeftClick;

                MainWindow->Mouse.LeftPressed = true;
                break;

            case WM_RBUTTONDOWN:
                event = Mouse_RightClick;

                MainWindow->Mouse.RightPressed = true;
                break;

            case WM_LBUTTONDBLCLK:
                event = Mouse_DoubleClick;

                MainWindow->Mouse.DoubleClicked = true;
                break;

            case WM_LBUTTONUP:
                event = Mouse_LeftUp;

                MainWindow->Mouse.LeftPressed = false;
                MainWindow->Mouse.DoubleClicked = false;
                break;

            case WM_RBUTTONUP:
                event = Mouse_RightUp;

                MainWindow->Mouse.RightPressed = false;
                break;

            case WM_DESTROY:
                event = Window_Close;
                break;

            case WM_KEYDOWN:
                event = Keyboard_Down;

                GetKeyboardState(allKeys);
                MainWindow->Key.Layout = (short) LOWORD(GetKeyboardLayout(0));
                MainWindow->Key.Shift = (allKeys[VK_SHIFT] & 0x80);
                MainWindow->Key.Control = (allKeys[VK_CONTROL] & 0x80);
                MainWindow->Key.Key = (unsigned short) wParam;
                break;

            case WM_CHAR:
                event = Keyboard_Char;

                MainWindow->Key.Char = (unsigned char) wParam;
                break;

            case WM_KEYUP:
                event = Keyboard_Up;

                MainWindow->Key.Char = 0;
				MainWindow->Key.Key = 0;
				MainWindow->Key.Shift = MainWindow->Key.Control = false;
                break;

            //case WM_CREATE:
                //event = Window_Create;
                //break;

            case WM_SHOWWINDOW:
                event = lParam ? Window_Show : Window_Create;
                break;

            case WM_PAINT:
                event = Window_Paint;

                _sInternals.dc = MainWindow->style.DynamicRedraw ? GetDC(hwnd)
                                                                 : BeginPaint(hwnd, &_sInternals.ps);
                break;

            case WM_WINDOWPOSCHANGED:
                event = Null;

                q = (WINDOWPOS*) lParam;

				if (MainWindow->Pos != Coord{q->x, q->y})
				{
					event = Window_Move;

					MainWindow->Pos.x = q->x;
					MainWindow->Pos.y = q->y;
				}

                if(MainWindow->Size != Coord{q->cx, q->cy})
                {
                    event = Window_Size;

                    MainWindow->Size.x = q->cx;
                    MainWindow->Size.y = q->cy;
                }

                if(!MainWindow->style.DynamicRedraw) InvalidateRect(hwnd, NULL, false);
                break;

            case WM_ACTIVATE:
                event = Window_Activate;
                break;

            case WM_SETFOCUS:
                event = Window_SetFocus;
                break;

            case WM_KILLFOCUS:
                event = Window_KillFocus;
                break;

            case WM_QUIT:
                //event = Null;
                break;

            default:
                event = sEventCodes::Unknown;
        }

        return event;
    }

    inline LRESULT CALLBACK SDefWindowProc (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
    {
        sEvent event;

        ///if(MainWindow->Components.size())
        event = TranslateNativeMessage(hwnd, message, wParam, lParam);

        if(!MainWindow->ProcessEvent)
        {
            bool _ret = false;

            switch(message)
            {
                case WM_PAINT:
                    if(!_sInternals.dc)
                    {
                        _sInternals.dc = MainWindow->style.DynamicRedraw ? GetDC(hwnd)
                                                                         : BeginPaint(hwnd, &_sInternals.ps);
                    }

					Clear();
                    break;

                case WM_KEYDOWN:
                    if(wParam == 27) PostQuitMessage(0);
                    break;

                case WM_DESTROY:
                    PostQuitMessage(0);
                    break;

                //case WM_CREATE:
                    //puts("!");
                    //break;

                default:
                    _ret = true;
                    ///return DefWindowProc(hwnd, message, wParam, lParam);
            }

            for(auto q = MainWindow->Components.begin(); q != MainWindow->Components.end(); q++)
            {
				if(*q)
                    (*q)->ProcessEvent(event);
				else
                    puts("WTF?! [Graphics.cpp: SDefWindowProc()]");
            }

            if(event == sEventCodes::Window_Paint)
                SFlush();

            if(_ret) return DefWindowProc(hwnd, message, wParam, lParam);
        }
        else
        {
            if(event == sEventCodes::Null)
                event = TranslateNativeMessage(hwnd, message, wParam, lParam);

            ///------------------------------------------------------------------------------------

            MainWindow->ProcessEvent(event);

            for(auto q = MainWindow->Components.begin(); q != MainWindow->Components.end(); q++)
            {
                if(*q)
                    (*q)->ProcessEvent(event);
				else
                    puts("WTF?! [Graphics.cpp: SDefWindowProc()]");
            }

            if(event == sEventCodes::Window_Paint)
                SFlush();

            return (event == sEventCodes::Unknown) ? DefWindowProc(hwnd, message, wParam, lParam) : 0;
        }

        return 0;
    }

    #else

    #endif
}

iRect GetDesktopRect()
{
    iRect res(Coord(0, 0));

    #ifdef _WL_WINDOWS
        RECT dskt;
        GetWindowRect(GetDesktopWindow(), &dskt);

        res.Width(dskt.right - dskt.left);
        res.Height(dskt.bottom - dskt.top);
    #else
        XWindowAttributes attr;
        Display* D = XOpenDisplay(0);
        XGetWindowAttributes(D, DefaultRootWindow(D), &attr);

        res.Width(attr.width);
        res.Height(attr.height);
    #endif

    return res;
}

decltype(MainWindow) SCreateWindow(int width, int height,
                                   #ifdef _WL_WINDOWS
								       const wchar_t* caption,
                                   #else
                                       const char* caption,
                                   #endif
								   const wStyle& style
                                  )
{
	MainWindow->style = style;
	MainWindow->Size = Coord(width, height);
	MainWindow->Caption = caption;

	///------------------------------------------------------------------------

	///Windows
    #ifdef _WL_WINDOWS

        WNDCLASSEX wincl;

        wincl.hInstance = NULL;
        wincl.lpszClassName = TEXT("ñWindowsApp");
        wincl.lpfnWndProc = SDefWindowProc;
        wincl.style = CS_DBLCLKS;
        wincl.cbSize = sizeof (WNDCLASSEX);

        wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
        wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
        wincl.hCursor = LoadCursor (NULL, IDC_ARROW);
        wincl.lpszMenuName = NULL;
        wincl.cbClsExtra = 0;
        wincl.cbWndExtra = 0;
        wincl.hbrBackground = NULL;

        if (!RegisterClassEx (&wincl)) return 0;

        iRect R = GetDesktopRect();
        MainWindow->BackBuffer.Create((R.Width() % 4) ? R.Width() + 4 - R.Width()%4 : R.Width(),
                                       R.Height());

        _sInternals.bi.bV4Width = MainWindow->BackBuffer.Width;
        _sInternals.bi.bV4Height = -MainWindow->BackBuffer.Height;

        Clear();

        //_sInternals.tdata = new unsigned char[MainWindow->BackBuffer.Cnst];
        //printf("%d\n", MainWindow->BackBuffer.Cnst/1024/1024);

        DWORD WS = WS_OVERLAPPEDWINDOW;
        if(MainWindow->style.code == wStyle::POPUP)
        {
            WS = WS_POPUP;
            //MainWindow->Size = R.Dim();
        }

        MainWindow->style.NativeStyle = WS;

        _sInternals.hwnd = CreateWindowEx (
                   0,                    //Extended possibilites for variation
                   wincl.lpszClassName,         //Classname
                   NULL,       //Title Text
                   WS,
                   CW_USEDEFAULT,       //Windows decides the position
                   CW_USEDEFAULT,       //where the window ends up on the screen
                   MainWindow->Size.x,                 //The programs width
                   MainWindow->Size.y,                 //and height in pixels
                   HWND_DESKTOP,        //The window is a child-window to desktop
                   NULL,                //No menu
                   NULL,       //Program Instance handler
                   NULL                 //No Window Creation data
                   );

        SetWindowTextW(_sInternals.hwnd, MainWindow->Caption.GetRef().data());

    ///Linux
    #else

    int run()
    {
        MainWindow->display = XOpenDisplay(0);
        if(!MainWindow->display) return 1;

        int scr = DefaultScreen(MainWindow->display);

        int blackColor = BlackPixel(MainWindow->display, scr);

        MainWindow->window = XCreateSimpleWindow(MainWindow->display,
                                                 RootWindow(MainWindow->display, scr),
                                                 0, 0, MainWindow->Width, MainWindow->Height,
                                                 0, blackColor, blackColor);
        XSetWindowBackgroundPixmap(MainWindow->display, MainWindow->window, None);
        XStoreName(MainWindow->display, MainWindow->window, MainWindow->Caption);

        XSelectInput(MainWindow->display, MainWindow->window,
                     StructureNotifyMask | ExposureMask | PointerMotionMask |
                     ButtonPressMask | KeyPressMask | ButtonReleaseMask | KeyReleaseMask);

        XMapWindow(MainWindow->display, MainWindow->window);

        XWindowAttributes attr;
        XGetWindowAttributes(MainWindow->display, DefaultRootWindow(MainWindow->display), &attr);

        MainWindow->BackBuffer->Width=attr.width;
        MainWindow->BackBuffer->Height=attr.height;

        MainWindow->BackBuffer->transparent=false;
        MainWindow->BackBuffer->depth = 4;

        MainWindow->BackBuffer->cnst = MainWindow->BackBuffer->depth*attr.width;
        MainWindow->BackBuffer->Cnst = MainWindow->BackBuffer->cnst*attr.height;
        MainWindow->BackBuffer->data=new unsigned char[MainWindow->BackBuffer->Cnst];
        MainWindow->BackBuffer->clear();

        MainWindow->Containers._tdata=new unsigned char[MainWindow->BackBuffer->Cnst];

        sMouse.LeftPressed=sMouse.RightPressed=false;

        setbuf(stdout, NULL);

        XEvent e;

        while(1)
        {
            int p = XPending(MainWindow->display);

            if(p) XNextEvent(MainWindow->display, &e);
            else e.type = 0;

            if(!WindowProcedure(&e)) break;
        }

        XCloseDisplay(MainWindow->display);
    }

    #endif
	///------------------------------------------------------------------------

	return MainWindow;
}
