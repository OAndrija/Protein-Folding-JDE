@echo off

g++ -O3 -ffast-math -o main.exe src\main.cpp
if %errorlevel% neq 0 (
    echo Compilation failed.
    exit /b %errorlevel%
)

echo run,seed,E,runtime,nfes,speed> data\results.csv

main.exe ABBBBBBABBBAB -seed 1 -target -5.6104 -nfesLmt 1000000 -runtimeLmt 60.0 -Np 300 -expRuns 50 -expThreads 16 >> data\results.csv


echo Running data analysis...
python scripts\analyze_runs.py
if %errorlevel% neq 0 (
    echo Data analysis script failed.
    exit /b %errorlevel%
)

echo All tasks completed successfully.
