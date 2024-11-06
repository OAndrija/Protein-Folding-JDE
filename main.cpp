#include <iostream>
#include <random>
#include <vector>
#include <cmath>

struct Solution {
    std::vector<double> angles;
    double F;
    double Cr;
};

double randomDouble(std::mt19937& generator, std::uniform_real_distribution<>& distribution) {
    return distribution(generator);
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

Solution jdeAlgorithm(std::vector<Solution>& population, std::mt19937& generator, std::uniform_real_distribution<>& distribution, unsigned int &Np, unsigned int &D) {
    for (unsigned int i = 0; i < Np; ++i) {
        if (randomDouble(generator, distribution) < 0.1) {
            population[i].F = 0.1 + randomDouble(generator, distribution) * 0.9;
        }

        if (randomDouble(generator, distribution) < 0.1) {
            population[i].Cr = randomDouble(generator, distribution);
        }

        unsigned int r1, r2, r3;
        do {
            r1 = std::uniform_int_distribution<>(0, Np - 1)(generator);
        } while ( r1 == i);

        do {
            r2 = std::uniform_int_distribution<>(0, Np - 1)(generator);
        } while (r2 == i || r2 == r1);

        do {
            r3 = std::uniform_int_distribution<>(0, Np - 1)(generator);
        } while (r3 == i || r3 == r1 || r3 == r2);

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

            if (U.angles[j] < -M_PI) {
                U.angles[j] = -M_PI;
            } else if (U.angles[j] > M_PI) {
                U.angles[j] = M_PI;
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

        std::cout << "S: " << S << "\n";
        std::cout << "Seed: " << seed << "\n";
        std::cout << "Target: " << target << "\n";
        std::cout << "NFES Limit: " << nfesLmt << "\n";
        std::cout << "Np: " << Np << "\n";
        std::cout << "D: " << D << "\n";

        std::mt19937 generator(seed);
        std::uniform_real_distribution<> distribution(0.0, 1.0);

        std::vector<Solution> population = initializePopulation(Np, D, generator, distribution);

        printPopulation(population);
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
