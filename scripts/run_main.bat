@echo off

g++ -O3 -ffast-math -o main.exe src\main.cpp
if %errorlevel% neq 0 (
    echo Compilation failed.
    exit /b %errorlevel%
)

echo seed,E,runtime,nfes,speed > data\results.csv

for /L %%i in (1,1,50) do (
    echo Running experiment with seed %%i
    experiment.exe ABBBBBBABBBAB -seed %%i -target -5.6104 -nfesLmt 1000000 -runtimeLmt 60.0 -Np 300
)

echo Running data analysis...
python data\analyze_runs.py
if %errorlevel% neq 0 (
    echo Data analysis script failed.
    exit /b %errorlevel%
)

echo All tasks completed successfully.
