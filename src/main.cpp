#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <chrono>
#include <random>
#include <fstream>

struct Solution {
    std::vector<double> angles;
    double F;
    double Cr;
};

struct Point {
    double x, y, z;
};

struct RunTask {
    unsigned int seed;
    unsigned int runNumber;
};

// Shared queue and synchronization primitives
std::queue<RunTask> taskQueue;
std::mutex queueMutex;
std::condition_variable queueCV;
bool stop = false;

double randomDouble(std::mt19937& generator, std::uniform_real_distribution<>& distribution) {
    return distribution(generator);
}

double cCoefficient(char si, char sj) {
    if (si == 'A' && sj == 'A') return 1.0;
    if (si == 'B' && sj == 'B') return 0.5;
    return -0.5;
}

double calculateDistance(const Point& p1, const Point& p2) {
    return sqrt((p2.x - p1.x) * (p2.x - p1.x) +  (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z));
}

//  Xj,i,G = rand[0,1] * (XjU - XjL) + XjL ; i=1, ...,Np;  j=1, ..., D;  Fi,1=0.5;  Cri,1=0.9
std::vector<Solution> initializePopulation(const unsigned int& Np, const unsigned int& D, std::mt19937& generator, std::uniform_real_distribution<>& distribution) {
    std::vector<Solution> population(Np);

    for (unsigned int i = 0; i < Np; ++i) {
        Solution solution;

        for (unsigned int j = 0; j < D; ++j) {
            double angle = randomDouble(generator, distribution) * 2 * M_PI - M_PI;
            solution.angles.push_back(angle);
        }

        solution.F = 0.5;
        solution.Cr = 0.9;

        population[i] = solution;
    }

    return population;
}

std::vector<Point> calculatePositions(const std::vector<double>& angles, unsigned int L) {
    std::vector<Point> positions(L);

    positions[0] = {0, 0, 0};
    positions[1] = {0, 1, 0};
    positions[2] = {cos(angles[0]), 1 + sin(angles[0]), 0};

    for (unsigned int i = 3; i < L; ++i) {
        double theta = angles[i - 2];
        double beta = angles[L - 2 + i - 3];

        double cos_theta = cos(theta);
        double sin_theta = sin(theta);
        double cos_beta = cos(beta);
        double sin_beta = sin(beta);

        positions[i].x = positions[i - 1].x + cos_theta * cos_beta;
        positions[i].y = positions[i - 1].y + sin_theta * cos_beta;
        positions[i].z = positions[i - 1].z + sin_beta;
    }

    return positions;
}

double calculateEnergy(const Solution& solution, const std::string& S) {
    double E1 = 0.0, E2 = 0.0;
    unsigned int L = S.length();

    // E1 = 1/4 * SUM (1 - cos(theta_i))
    for (unsigned int i = 0; i < L - 2; ++i) {
        E1 += (1 - std::cos(solution.angles[i]));
    }
    E1 /= 4.0;

    // We calculate the positions of angles in a Solution
    std::vector<Point> positions = calculatePositions(solution.angles, L);

    // E2 = 4 * SUM ( SUM (d(pi,pj)^-12 - c_ij * d(pi,pj)^-6) )
    for (unsigned int i = 0; i < L - 2; ++i) {
        for (unsigned int j = i + 2; j < L; ++j) {
            double d_ij = calculateDistance(positions[i], positions[j]);
            double d_ij_6 = 1.0 / (d_ij * d_ij * d_ij * d_ij * d_ij * d_ij);
            double d_ij_12 = d_ij_6 * d_ij_6;
            double c_ij = cCoefficient(S[i], S[j]);

            E2 += (d_ij_12 - c_ij * d_ij_6);
        }
    }
    E2 *= 4;

    return E1 + E2;
}

Solution jdeAlgorithm(std::vector<Solution>& population, std::mt19937& generator,
                      std::uniform_real_distribution<>& distribution, const std::string& S,
                      const unsigned int &Np, const unsigned int &D,
                      double target, double epsilon, unsigned int nfesLmt, double runtimeLmt,
                      unsigned int &nfes)
{
    auto startTime = std::chrono::high_resolution_clock::now();
    Solution bestSolution = population[0];
    double bestEnergy = calculateEnergy(bestSolution, S);

    while (true) {
        for (unsigned int i = 0; i < Np; ++i) {
            // MUTACIJA: if(rand[0, 1] < 0.1) Fi,G = 0.1 + rand[0, 1] * 0.9 else  Fi,G = Fi,G

            if (randomDouble(generator, distribution) < 0.1) {
                population[i].F = 0.1 + randomDouble(generator, distribution) * 0.9;
            }

            // Vj,i,G=Xj,r1,G+Fi,G*(Xj,r2,G - Xj,r3,G), i=1, ..., Np; j=1, ..., D; r= rand{1, ..., Np} r1 != r2 != r3 != i
            unsigned int r1, r2, r3;
            do { r1 = std::uniform_int_distribution<>(0, Np - 1)(generator); } while (r1 == i);
            do { r2 = std::uniform_int_distribution<>(0, Np - 1)(generator); } while (r2 == i || r2 == r1);
            do { r3 = std::uniform_int_distribution<>(0, Np - 1)(generator); } while (r3 == i || r3 == r1 || r3 == r2);

            Solution V;
            V.angles.resize(D);
            for (unsigned int j = 0; j < D; ++j) {
                V.angles[j] = population[r1].angles[j] + population[i].F * (population[r2].angles[j] - population[r3].angles[j]);
            }

            // KRIZANJE: if(rand[0, 1] < 0.1) Cri,G = rand[0, 1] else Cri,G = Cri,G
            if (randomDouble(generator, distribution) < 0.1) {
                population[i].Cr = randomDouble(generator, distribution);
            }

            Solution U;
            U.angles.resize(D);
            unsigned int jrand = std::uniform_int_distribution<>(0, D - 1)(generator);

            //if(rand[0,1] < Cri,G  ||  ji,rand == j) Uj,i,G= Vj,i,G else Uj,i,G = Xj,i,G
            for (unsigned int j = 0; j < D; ++j) {
                if (randomDouble(generator, distribution) < population[i].Cr || jrand == j) {
                    U.angles[j] = V.angles[j];
                } else {
                    U.angles[j] = population[i].angles[j];
                }

                // POPRAVLJANJE
                if (U.angles[j] < -M_PI) U.angles[j] = -M_PI;
                if (U.angles[j] > M_PI) U.angles[j] = M_PI;
            }

            double energy_U = calculateEnergy(U, S);
            double energy_Xi = calculateEnergy(population[i], S);
            nfes++;

            // SELEKCIJA
            if (energy_U < energy_Xi) {
                population[i] = U;
                if (energy_U < bestEnergy) {
                    bestEnergy = energy_U;
                    bestSolution = U;
                }
            }

            auto currentTime = std::chrono::high_resolution_clock::now();
            double elapsedTime = std::chrono::duration<double>(currentTime - startTime).count();

            if (bestEnergy <= target + epsilon || elapsedTime >= runtimeLmt || nfes >= nfesLmt) {
                return bestSolution;
            }
        }
    }
}

void parseArguments(int argc, char* argv[], std::string &S, unsigned int &seed, double &target, unsigned int &nfesLmt,
    double &runtimeLmt, unsigned int &Np, unsigned int &D, unsigned int &expRuns, unsigned int &expThreads, unsigned int &algThreads) {
    if (argc != 16) {
        throw std::invalid_argument("Invalid number of arguments. Expected number: 15.");
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (i == 1) {
            S = argv[i];
            D = 2 * S.length() - 5;
        } else if (arg == "-seed" && i + 1 < argc) {
            seed = std::strtoul(argv[++i], nullptr, 10);
        } else if (arg == "-target" && i + 1 < argc) {
            target = std::strtod(argv[++i], nullptr);
        } else if (arg == "-nfesLmt" && i + 1 < argc) {
            nfesLmt = std::strtoul(argv[++i], nullptr, 10);
        } else if (arg == "-runtimeLmt" && i + 1 < argc) {
            runtimeLmt = std::strtod(argv[++i], nullptr);
        } else if (arg == "-Np" && i + 1 < argc) {
            Np = std::strtoul(argv[++i], nullptr, 10);
        } else if (arg == "-expRuns" && i + 1 < argc) {
            expRuns = std::strtoul(argv[++i], nullptr, 10);
        } else if (arg == "-expThreads" && i + 1 < argc) {
            expThreads = std::strtoul(argv[++i], nullptr, 10);
        } else if (arg == "-algThreads" && i + 1 < argc) {
            algThreads = std::strtoul(argv[++i], nullptr, 10);
        } else {
            throw std::invalid_argument("Unknown argument: " + arg);
        }
    }
}

void printPopulation(const std::vector<Solution>& population) {
    std::cout << "Population size: " << population.size() << "\n";
    for (unsigned int i = 0; i < population.size(); ++i) {
        std::cout << "Solution " << i + 1 << ":\n";
        std::cout << "Angles: ";
        for (double angle : population[i].angles) {
            std::cout << angle << " ";
        }
        std::cout << "\n";
        std::cout << "F: " << population[i].F << ", Cr: " << population[i].Cr << "\n\n";
    }
}

// Producer function: generates tasks and pushes them to the queue
void producer(unsigned int expRuns, unsigned int seed) {
    for (unsigned int i = 0; i < expRuns; ++i) {
        {
            std::lock_guard<std::mutex> lock(queueMutex);
            taskQueue.push({seed + i, i + 1}); // Seed and run number
        }
        queueCV.notify_one(); // Notify a consumer thread
    }

    // Notify all threads to stop after the last task
    {
        std::lock_guard<std::mutex> lock(queueMutex);
        stop = true;
    }
    queueCV.notify_all();
}

// Consumer function: processes tasks from the queue
void consumer(const std::string& S, unsigned int Np, unsigned int D, double target,
              unsigned int nfesLmt, double runtimeLmt, double epsilon) {
    while (true) {
        RunTask task;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            queueCV.wait(lock, [] { return !taskQueue.empty() || stop; });

            // Exit if there are no tasks and stop is true
            if (stop && taskQueue.empty()) {
                return;
            }

            task = taskQueue.front();
            taskQueue.pop();
        }

        // Initialize the generator for this run
        std::mt19937 generator(task.seed);
        std::uniform_real_distribution<> distribution(0.0, 1.0);
        std::vector<Solution> population = initializePopulation(Np, D, generator, distribution);

        unsigned int nfes = 0;

        auto startTime = std::chrono::high_resolution_clock::now();
        Solution bestSolution = jdeAlgorithm(population, generator, distribution, S, Np, D, target, epsilon, nfesLmt, runtimeLmt, nfes);
        auto endTime = std::chrono::high_resolution_clock::now();

        double bestEnergy = calculateEnergy(bestSolution, S);
        double runtime = std::chrono::duration<double>(endTime - startTime).count();
        double speed = nfes / runtime;

        // Print results for this run
        std::cout << task.runNumber << ",";
        std::cout << task.seed << ",";
        std::cout << bestEnergy << ",";
        std::cout << runtime << ",";
        std::cout << nfes << ",";
        std::cout << speed << "\n";
        // std::cout << "Solution Angles: ";
        // for (double angle : bestSolution.angles) {
        //     std::cout << angle << " ";
        // }
        // std::cout << "\n\n";
    }
}

int main(int argc, char* argv[]) {
    try {
        std::string S;
        unsigned int seed = 0, nfesLmt = 0, Np = 0, D = 0, expRuns = 0, expThreads = 0, algThreads = 0;
        double target = 0.0, runtimeLmt = 0.0;

        parseArguments(argc, argv, S, seed, target, nfesLmt, runtimeLmt, Np, D, expRuns, expThreads, algThreads);

        double epsilon = 1e-6;

        auto startTime = std::chrono::high_resolution_clock::now();
        // Start producer thread
        std::thread producerThread(producer, expRuns, seed);

        // Start consumer threads
        std::vector<std::thread> consumerThreads;
        for (unsigned int i = 0; i < expThreads; ++i) {
            consumerThreads.emplace_back(consumer, S, Np, D, target, nfesLmt, runtimeLmt, epsilon);
        }

        // Join threads
        producerThread.join();
        for (auto& thread : consumerThreads) {
            thread.join();
        }
        auto endTime = std::chrono::high_resolution_clock::now();

        // Calculate elapsed time
        std::chrono::duration<double> elapsed = endTime - startTime;
        double executionTime = elapsed.count();

        // Save to CSV
        std::ofstream csvFile("data/timing_results.csv", std::ios::app);
        if (csvFile.is_open()) {
            csvFile << expThreads << "," << executionTime << "\n";
            csvFile.close();
        } else {
            std::cerr << "Error: Unable to open file for writing.\n";
        }

    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
