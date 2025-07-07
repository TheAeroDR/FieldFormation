@echo off
echo Compiling CUDA double precision field computation...
echo.

REM Set up Visual Studio environment (try multiple possible locations)
echo Setting up Visual Studio environment...

REM Try VS 2022 Community first
if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" (
    echo Found VS 2022 Community
    call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
    goto compile
)

REM Try VS 2022 Professional
if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" (
    echo Found VS 2022 Professional
    call "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
    goto compile
)

REM Try VS 2019 Community
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat" (
    echo Found VS 2019 Community
    call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
    goto compile
)

REM Try VS Build Tools
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
    echo Found VS 2019 Build Tools
    call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
    goto compile
)

REM Try common MSVC locations for older versions
if exist "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" (
    echo Found VS 2015
    call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x64 > nul 2>&1
    goto compile
)

REM Check if cl.exe is already in PATH
where cl.exe > nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo Found cl.exe in PATH
    goto compile
)

REM If no Visual Studio found, show error
echo ERROR: Could not find Visual Studio environment setup script
echo.
echo Please install Visual Studio with C++ tools or run from:
echo "Developer Command Prompt for VS" or "x64 Native Tools Command Prompt"
echo.
echo You can also manually run:
echo call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
echo.
goto end

:compile
echo.
echo Compiling with nvcc...

REM Compile CUDA code with double precision (no --use_fast_math for better precision)
nvcc -O3 -arch=sm_86 -std=c++17 fields_cuda.cu -o fields_cuda.exe

if %ERRORLEVEL% EQU 0 (
    echo.
    echo ===== COMPILATION SUCCESSFUL =====
    echo Executable: fields_cuda.exe
    echo Double precision enabled for accurate small magnetic field values
    echo --use_fast_math removed to preserve precision
    echo.
    echo To run: fields_cuda.exe
    echo.
) else (
    echo.
    echo ===== COMPILATION FAILED =====
    echo Error code: %ERRORLEVEL%
    echo.
    echo Troubleshooting:
    echo 1. Make sure Visual Studio is installed with C++ tools
    echo 2. Make sure CUDA Toolkit is installed
    echo 3. Try running from "Developer Command Prompt for VS"
    echo 4. Check if cl.exe is in PATH: where cl
    echo.
)

pause

:end
