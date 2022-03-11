@echo off

:: compiler and `make` setup

set COMPILER_PATH=C:\Program Files\mingw-w64\x86_64-8.1.0-win32-seh-rt_v6-rev0\mingw64\bin\

set MAKE=mingw32-make.exe

set PATH=%COMPILER_PATH%;%PATH%

::-------------------------------------------------------------------------------------

:: build using `make`

"%COMPILER_PATH%\%MAKE%"