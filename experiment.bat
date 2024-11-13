@echo off

:: Compile the C++ program with g++
g++ -O3 -ffast-math -o experiment.exe main.cpp

:: Check if the compilation was successful
if %errorlevel% neq 0 (
    echo Compilation failed.
    exit /b %errorlevel%
)

:: Run the program with seed values from 1 to 50
for /L %%i in (1,1,50) do (
    echo Running experiment with seed %%i
    experiment.exe ABBBBBBABBBAB -seed %%i -target -5.6104 -nfesLmt 1000000 -runtimeLmt 60.0 -Np 300
)