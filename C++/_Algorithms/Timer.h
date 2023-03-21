#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED

#include <stdio.h>
#include <time.h>
#include <windows.h>

struct Timer
{
    clock_t T;

    Timer()
    {
        T = clock();
    }

    ~Timer()
    {
        T = clock() - T;

        const auto WHITE = FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED;

        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), WHITE | FOREGROUND_INTENSITY);

        printf("%g seconds [%d ticks] passed\n", 1.*T/CLOCKS_PER_SEC, (int) T);

        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), WHITE);
    }
};

#endif // TIMER_H_INCLUDED
