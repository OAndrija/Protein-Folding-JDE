#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include<chrono>

struct Solution {
    std::vector<double> angles;
    double F;
    double Cr;
};

struct Point {
    double x, y, z;
};

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

        positions[i].x = positions[i - 1].x + cos(theta) * cos(beta);
        positions[i].y = positions[i - 1].y + sin(theta) * cos(beta);
        positions[i].z = positions[i - 1].z + sin(beta);
    }

    return positions;
}

double calculateEnergy(const Solution& solution, const std::string& S) {
    double E1 = 0.0, E2 = 0.0;
    unsigned int L = S.length();

    for (unsigned int i = 0; i < L - 2; ++i) {
        E1 += (1 - std::cos(solution.angles[i]));
    }
    E1 /= 4.0;

    std::vector<Point> positions = calculatePositions(solution.angles, L);

    for (unsigned int i = 0; i < L - 2; ++i) {
        for (unsigned int j = i + 2; j < L; ++j) {
            double d_ij = calculateDistance(positions[i], positions[j]);
            double d_ij_6 = pow(d_ij, -6);
            double d_ij_12 = pow(d_ij, -12);
            double c_ij = cCoefficient(S[i], S[j]);

            E2 += 4 * (d_ij_12 - c_ij * d_ij_6);
        }
    }

    return E1 + E2;
}

Solution jdeAlgorithm(std::vector<Solution>& population, std::mt19937& generator,
                      std::uniform_real_distribution<>& distribution, const std::string& S,
                      const unsigned int &Np, const unsigned int &D,
                      double target, double epsilon, unsigned int nfesLmt, double runtimeLmt)
{
    auto startTime = std::chrono::high_resolution_clock::now();
    unsigned int nfes = 0;
    Solution bestSolution = population[0];
    double bestEnergy = calculateEnergy(bestSolution, S);

    while (true) {
        for (unsigned int i = 0; i < Np; ++i) {
            if (randomDouble(generator, distribution) < 0.1) {
                population[i].F = 0.1 + randomDouble(generator, distribution) * 0.9;
            }

            if (randomDouble(generator, distribution) < 0.1) {
                population[i].Cr = randomDouble(generator, distribution);
            }

            unsigned int r1, r2, r3;
            do { r1 = std::uniform_int_distribution<>(0, Np - 1)(generator); } while (r1 == i);
            do { r2 = std::uniform_int_distribution<>(0, Np - 1)(generator); } while (r2 == i || r2 == r1);
            do { r3 = std::uniform_int_distribution<>(0, Np - 1)(generator); } while (r3 == i || r3 == r1 || r3 == r2);

            Solution V;
            V.angles.resize(D);
            for (unsigned int j = 0; j < D; ++j) {
                V.angles[j] = population[r1].angles[j] + population[i].F * (population[r2].angles[j] - population[r3].angles[j]);
            }

            Solution U;
            U.angles.resize(D);
            unsigned int jrand = std::uniform_int_distribution<>(0, D - 1)(generator);
            for (unsigned int j = 0; j < D; ++j) {
                if (randomDouble(generator, distribution) < population[i].Cr || jrand == j) {
                    U.angles[j] = V.angles[j];
                } else {
                    U.angles[j] = population[i].angles[j];
                }

                if (U.angles[j] < -M_PI) U.angles[j] = -M_PI;
                if (U.angles[j] > M_PI) U.angles[j] = M_PI;
            }

            double energy_U = calculateEnergy(U, S);
            double energy_Xi = calculateEnergy(population[i], S);
            nfes += 2;

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

void parseArguments(int argc, char* argv[], std::string &S, unsigned int &seed, double &target, unsigned int &nfesLmt, unsigned int &Np, unsigned int &D) {
    if (argc != 11) {
        throw std::invalid_argument("Invalid number of arguments. Expected number: 10.");
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "S" && i + 1 < argc) {
            S = argv[++i];
            D = 2 * S.length() - 5;
        } else if (arg == "-seed" && i + 1 < argc) {
            seed = std::strtoul(argv[++i], nullptr, 10);
        } else if (arg == "-target" && i + 1 < argc) {
            target = std::strtod(argv[++i], nullptr);
        } else if (arg == "-nfesLmt" && i + 1 < argc) {
            nfesLmt = std::strtoul(argv[++i], nullptr, 10);
        } else if (arg == "-Np" && i + 1 < argc) {
            Np = std::strtoul(argv[++i], nullptr, 10);
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

int main(int argc, char* argv[]) {
    try {
        std::string S;
        unsigned int seed = 0, nfesLmt = 0, Np = 0, D = 0;
        double target = 0.0;

        parseArguments(argc, argv, S, seed, target, nfesLmt, Np, D);

        std::mt19937 generator(seed);
        std::uniform_real_distribution<> distribution(0.0, 1.0);
        std::vector<Solution> population = initializePopulation(Np, D, generator, distribution);

        double epsilon = 1e-6;
        double runtimeLmt = 60.0;

        auto startTime = std::chrono::high_resolution_clock::now();
        Solution bestSolution = jdeAlgorithm(population, generator, distribution, S, Np, D, target, epsilon, nfesLmt, runtimeLmt);
        auto endTime = std::chrono::high_resolution_clock::now();

        double bestEnergy = calculateEnergy(bestSolution, S);
        unsigned int nfes = nfesLmt;
        double runtime = std::chrono::duration<double>(endTime - startTime).count();
        double speed = nfes / runtime;

        std::cout << "S: " << S << "\n";
        std::cout << "seed: " << seed << "\n";
        std::cout << "nfes: " << nfes << "\n";
        std::cout << "runtime: " << runtime << "\n";
        std::cout << "speed: " << speed << "\n";
        std::cout << "E: " << bestEnergy << "\n";
        std::cout << "solution angles: ";
        for (double angle : bestSolution.angles) {
            std::cout << angle << " ";
        }
        std::cout << "\n";

    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
