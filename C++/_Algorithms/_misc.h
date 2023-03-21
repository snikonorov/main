#ifndef _MISC_H_INCLUDED
#define _MISC_H_INCLUDED

#include <stdio.h>
#include <string>

void execute(const std::wstring& Cmd, bool halt = false, DWORD flags = CREATE_NO_WINDOW)
{
    static std::wstring cmd;

    STARTUPINFOW si;
    PROCESS_INFORMATION pi;

    ZeroMemory(&si, sizeof(si));
    ZeroMemory(&pi, sizeof(pi));

    cmd = Cmd;

    CreateProcessW(NULL, (wchar_t*) cmd.data(), NULL, NULL, FALSE, flags, NULL, NULL, &si, &pi);

    if(halt)
        WaitForSingleObject(pi.hProcess, INFINITE);

    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);
}

///--------------------------------------------------------------------

/// returns `true` if successful
///
bool str_to_file(const std::string& s, const char* path)
{
	FILE* f = fopen(path, "w");
	
	if(f)
	{
		fprintf(f, "%s", s.data());
		fclose(f);
		
		return true;
	}
	
	return false;
}

std::string str_from_file(const char* path)
{
	std::string str;

	FILE* f = fopen(path, "r");

	if(f)
	{
		fseek(f, 0, SEEK_END);
		unsigned N = ftell(f);
		fseek(f, 0, SEEK_SET);

		str.resize(N);
		fread((void*)str.data(), sizeof(char), str.size(), f);

		fclose(f);		
	}

	return str;
}

bool str_to_file(const std::string& s, const std::string& path)
{
	return str_to_file(s, path.data());
}

std::string str_from_file(const std::string& path)
{
	return str_from_file(path.data());
}

///--------------------------------------------------------------------

template<class T>
std::string Dec2Hex(T a)
{
    static char hex[] = {'0', '1', '2', '3', '4', '5', '6', '7',
                         '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};

    std::string s;
    s.resize(2*sizeof(a));

    for(int k = s.size()-1; k >= 0; k--)
    {
        s[k] = hex[a % 16];
        a /= 16;
    }

    return s;
}

#endif // _MISC_H_INCLUDED
