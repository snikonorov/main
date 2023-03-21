#define _CRT_SECURE_NO_WARNINGS

#include "Font.h"

const Char& sFont::operator[](unsigned key_code)
{
	auto q = data.find(key_code);

	if (q != data.end())
	{
		return q->second;
	}
	else
	{
		Char& chr = data[key_code];

		#ifdef _WL_WINDOWS

		HDC hdc = GetDC(MainWindow->_Internals.hwnd);
		HDC buffDC = CreateCompatibleDC(hdc);

		HFONT hFnt = CreateFontA(Height, 0, 0, 0, FW_REGULAR, false, false, false,
												ANSI_CHARSET, OUT_OUTLINE_PRECIS, CLIP_DEFAULT_PRECIS, PROOF_QUALITY,
												DEFAULT_PITCH | FF_DONTCARE, _family);

		SelectObject(buffDC, hFnt);
		SetTextColor(buffDC, RGB(0, 0, 0));
		SetBkMode(buffDC, TRANSPARENT);

		char txt[5] = { 0, 0, 0, 0, 0 };
		//*((int*)txt) = key_code;

		reinterpret_cast<unsigned&>(txt) = key_code;

		SIZE sz;
		GetTextExtentPoint32W(buffDC, (LPCWSTR) txt, 1, &sz);

		RECT rc = {0, 0, sz.cx, sz.cy + 2};

		auto& Width = rc.right;
		auto& Height = rc.bottom;

		auto pWidth = 3*Width;
		pWidth += (pWidth%4 ? (4 - pWidth%4) : 0);

		_data.resize(pWidth*Height);

		BITMAPINFO info;
		ZeroMemory(&info.bmiHeader, sizeof(info.bmiHeader));

		info.bmiHeader.biSize = sizeof(info.bmiHeader);
		info.bmiHeader.biHeight = -Height;
		info.bmiHeader.biPlanes = 1;
		info.bmiHeader.biBitCount = 24;
		info.bmiHeader.biCompression = BI_RGB;

		HBITMAP buffBitmap = CreateCompatibleBitmap(hdc, Width, Height);
		SelectObject(buffDC, buffBitmap);

		int savedDC = SaveDC(buffDC);

		HBRUSH hBrush = CreateSolidBrush(RGB(255, 255, 255));

		///--------------------------------------------------------

		GetTextExtentPoint32W(buffDC, (LPCWSTR) txt, 1, &sz);
		rc.right = sz.cx;

		FillRect(buffDC, &rc, hBrush);
		DrawTextW(buffDC, (LPCWSTR) txt, 1, &rc, DT_NOCLIP);

		info.bmiHeader.biWidth = Width;
		chr.Width = (unsigned short) Width;

		pWidth = 3*Width;
		pWidth += (pWidth%4 ? (4 - pWidth%4) : 0);

		GetDIBits(buffDC, buffBitmap, 0, Height, _data.data(), &info, DIB_RGB_COLORS);

		chr.matrix.resize(Width);
		for(int x = 0; x < Width; x++)
			chr.matrix[x].resize(Height);

		for(int y = 0; y < Height; y++)
		{
			for(int x = 0; x < Width; x++)
			{
				auto cnst = pWidth*y + 3*x;
				pixel p = {_data[cnst + 2], _data[cnst + 1], _data[cnst]};
				chr.matrix[x][y] = (unsigned char) round(255*(p.Grey()));
			}
		}

		RestoreDC(buffDC, savedDC);

		DeleteObject(hBrush);
		DeleteObject(buffBitmap);
		DeleteDC(buffDC);

		DeleteObject(hFnt);

		_data.clear();

		ReleaseDC(MainWindow->_Internals.hwnd, hdc);

		#endif

		return chr;
	}
}
