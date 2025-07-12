@echo off
echo Compiling CUDA Sensitivity Study...
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
pause
exit /b 1

:compile
echo Environment setup complete.
nvcc -O3 -arch=sm_89 -std=c++17 sensittivity_cuda.cu -o sensitivity_cuda.exe

if errorlevel 1 (
    echo.
    echo Compilation FAILED!
    echo Please check:
    echo   1. CUDA toolkit is installed
    echo   2. GPU compute capability matches -arch flag
    echo   3. All source files are present
    pause
    exit /b 1
) else (
    echo.
    echo Compilation SUCCESSFUL!
    echo Executable: sensitivity_cuda.exe
    echo.
    echo To run sensitivity study: run_sensitivity_study.bat
    echo To run individual test: sensitivity_cuda.exe [test_name]
    echo.
    pause
)
