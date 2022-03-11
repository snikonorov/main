#ifndef FONT_H_INCLUDED
#define FONT_H_INCLUDED

#include "Graphics.h"

#include <map>

struct Char
{
	std::vector<std::vector<unsigned char>> matrix;
	unsigned short Width;
};

class sFont
{
	private:

	std::vector<unsigned char> _data;

	#ifdef _WL_WINDOWS

	const char* _family;

	#endif

	public:

    std::map<unsigned, Char> data;

	///unsigned short Height;
	TVar<unsigned short> Height;

	sFont(unsigned short height, const char* family = "Lucida console"): _family(family), Height(height)
	{
	    Height.SetValue = [this](auto* obj, const auto& val)
	    {
	        obj->GetRef() = val;
	        this->data.clear();
        };
	}

	sFont(const sFont&) = delete;       /// #ToFix: gotta deal with it

	/*void SetHeight(unsigned short height)
	{
		data.clear();
		Height = height;
	}

	unsigned short getHeight()
	{
		return Height;
	}*/

	const Char& operator[](unsigned key_code);

	const Char& chr(unsigned key_code)
	{
		return (*this)[key_code];
	}

	const char* GetFontFamily() const
	{
		return _family;
	}
};

///--------------------------------------------------------------------------------------------------

template<class TChar>
inline int Text_out(sFont* font, int X, int Y, const TChar* ch, const pixel& color = pixel(0, 0, 0))
{
	int x0 = X;

	if(font)
	{
		int h = 0;

		pixel c = color;

		while (*ch)
		{
			const auto& p = *ch;

			if (p == 13)
			{
				ch++;
				continue;
			}
			else if (p == 10)
			{
				h += font->Height + 5;
				x0 = X;
				ch++;
				continue;
			}

			const auto& chr = (*font)[p];

			for (unsigned short x = 0; x < chr.Width; x++)
			{
				for (unsigned short y = 0; y < font->Height; y++)
				{
					if (chr.matrix[x][y] != 255)
					{
						c.a = (255 - chr.matrix[x][y])*color.a/255;
						PutPixel(x + x0, y + Y + h, c);
					}
				}
			}

			x0 += chr.Width;

			ch++;
		}
	}

	return x0;
}

///unsigned Text_Width(sFont* font, const char* s, char terminate = 0, int wt = -1);

template<class TChar>
unsigned Text_Width(sFont* font, const TChar* s, TChar terminate = 0, int wt = -1)
{
    unsigned x0 = 0;

    if(font)
	{
		const auto* _s = s;

		while(*s != terminate && *s != 0 && !(wt > 0 && s-_s >= wt))
		{
			x0 += (*font)[*s].Width;
			s++;
		}
	}

    return x0;
}

#endif // FONT_H_INCLUDED
