@echo off
setlocal enabledelayedexpansion

echo Starting Comprehensive Sensitivity Study
echo ==========================================

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
pause
exit /b 1

:compile
echo Environment setup complete.

:: Compile the CUDA program
echo Compiling sensittivity_cuda.cu...
nvcc -O3 -arch=sm_89 -std=c++17 sensittivity_cuda.cu -o sensitivity_cuda.exe
if errorlevel 1 (
    echo Compilation failed!
    pause
    exit /b 1
)
echo Compilation successful.
echo.

:: Create results directory
if not exist "results" mkdir results

:: List of all sensitivity tests
set tests=baseline rad_profile_inner_10 rad_profile_outer_10 theta_plus_10 theta_minus_10 height_plus_10 height_minus_10 charge_plus_10 charge_minus_10 electron_density_plus_10 electron_density_minus_10 eyewall_velocity_plus_10 eyewall_velocity_minus_10 dust_loading_plus_10 dust_loading_minus_10 dust_density_plus_10 dust_density_minus_10 avg_particle_size_plus_10 avg_particle_size_minus_10 translational_velocity_plus_10 translational_velocity_minus_10

:: Counter for progress tracking
set count=0
for %%t in (%tests%) do set /a count+=1
set total=%count%
set current=0

:: Run each sensitivity test
for %%t in (%tests%) do (
    set /a current+=1
    echo.
    echo [!current!/%total%] Running sensitivity test: %%t
    echo ============================================
    
    :: Run the simulation
    sensitivity_cuda.exe %%t
    
    if errorlevel 1 (
        echo Error running test %%t
        echo Continuing with next test...
    ) else (
        echo Test %%t completed successfully.
        
        :: Move output files to results directory
        if exist "sensitivity*_timeseries_%%t.csv" (
            move "sensitivity*_timeseries_%%t.csv" "results\"
            echo Moved timeseries file to results directory.
        )
        
        if exist "sensitivity_metrics_%%t.txt" (
            move "sensitivity_metrics_%%t.txt" "results\"
            echo Moved metrics file to results directory.
        )
    )
)

echo.
echo ==========================================
echo Sensitivity Study Complete!
echo ==========================================
echo All results are stored in the 'results' directory.
echo.

:: Create summary file
echo Creating summary of all tests...
echo Sensitivity Study Summary > results\sensitivity_study_summary.txt
echo Generated on: %date% %time% >> results\sensitivity_study_summary.txt
echo. >> results\sensitivity_study_summary.txt

echo Tests performed: >> results\sensitivity_study_summary.txt
for %%t in (%tests%) do (
    echo   - %%t >> results\sensitivity_study_summary.txt
)

echo. >> results\sensitivity_study_summary.txt
echo Total tests: %total% >> results\sensitivity_study_summary.txt

echo Summary created at: results\sensitivity_study_summary.txt
echo.
pause
