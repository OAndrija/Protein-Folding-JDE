@echo off

g++ -O3 -ffast-math -o experiment.exe main.cpp
if %errorlevel% neq 0 (
    echo Compilation failed.
    exit /b %errorlevel%
)

echo seed,E,runtime,nfes,speed > results.csv

for /L %%i in (1,1,50) do (
    echo Running experiment with seed %%i
    experiment.exe ABBBBBBABBBAB -seed %%i -target -5.6104 -nfesLmt 1000000 -runtimeLmt 60.0 -Np 300
)
