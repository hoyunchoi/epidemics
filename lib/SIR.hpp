#pragma once

#include <map>
#include <random>  //uniform distribution
#include <set>
#include <string>
#include <vector>

#include "CSV.hpp"
#include "Networks.hpp"
#include "pcg_random.hpp"
#include "stringFormat.hpp"
#include "Epidemics.hpp"

/*
    SIR Model simulation
    S+I -> I+I with rate SI_II
    I -> R with rate I_R
*/
namespace SIR {
//* pre-defined parameters
const std::string rootDirectory = "../data/epidemics/SIR/";
std::map<std::string, int> stateToInt = {{"S", 0}, {"I", 1}, {"R", 2}};
unsigned seedSize;

//* Network parameters
std::string networkType;
unsigned networkSize;
unsigned meanDegree;
std::vector<Node_Epidemic<unsigned>> nodes;

//* SIR rate parameter
double SI_II;
double I_R;

//* Random parameter
int randomEngineSeed;
pcg32 randomEngine;
std::uniform_real_distribution<double> probabilityDistribution(0, 1);

//* File name convention
const std::string fileName() {
    const std::string fileName = "N" + to_stringWithExponent((double)networkSize, 1) + ",M" + std::to_string(meanDegree) + ",SIII" + to_stringWithPrecision(SI_II, 2) + ",IR" + to_stringWithPrecision(I_R, 2);

    return randomEngineSeed == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(randomEngineSeed) + ".txt";
}

//* Set Network
void setNetwork(const Network<unsigned>& t_network) {
    networkType = t_network.type;
    networkSize = t_network.size;
    meanDegree = t_network.meanDegree;

    nodes.clear();
    nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        Node_Epidemic node(index, "S");
        node.neighbors = t_network.adjacency[index];
        nodes.emplace_back(node);
    }
}

//* Set rate parameters
void setRate(const double& t_SI_II, const double& t_I_R, const unsigned& t_seedSize = 1) {
    SI_II = t_SI_II;
    I_R = t_I_R;
    seedSize = t_seedSize;
}

//* Simulate SIR model using 4th-order Runge-Kutta method
namespace RK4 {
//* pre-defined parameters
const int recentR_length = 5;
const double err = 1e-7;
const double deltaT = 1e-2;

//* Parameters used for RK4
unsigned iteration;
std::vector<double> recentR;
double currentTime;
double ratioS;
double ratioI;
double ratioR;

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngine;
    iteration = 0;
    currentTime = 0.0;
    recentR.assign(recentR_length, -1);
    ratioS = 1.0 - (double)seedSize / networkSize;
    ratioI = (double)seedSize / networkSize;
    ratioR = 0.0;
}

//* RK-4
const double dotS(const double& t_S, const double& t_I) { return -SI_II * meanDegree * t_S * t_I; }
const double dotI(const double& t_S, const double& t_I) { return SI_II * meanDegree * t_S * t_I - I_R * t_I; }
const double dotR(const double& t_I) { return I_R * t_I; }

//* update one step
void update() {
    const double ratioS1 = dotS(ratioS, ratioI);
    const double ratioI1 = dotI(ratioS, ratioI);
    const double ratioR1 = dotR(ratioI);

    const double ratioS2 = dotS(ratioS + ratioS1 * deltaT / 2, ratioI + ratioI1 * deltaT / 2);
    const double ratioI2 = dotI(ratioS + ratioS1 * deltaT / 2, ratioI + ratioI1 * deltaT / 2);
    const double ratioR2 = dotR(ratioI + ratioI1 * deltaT / 2);

    const double ratioS3 = dotS(ratioS + ratioS2 * deltaT / 2, ratioI + ratioI2 * deltaT / 2);
    const double ratioI3 = dotI(ratioS + ratioS2 * deltaT / 2, ratioI + ratioI2 * deltaT / 2);
    const double ratioR3 = dotR(ratioI + ratioI2 * deltaT / 2);

    const double ratioS4 = dotS(ratioS + ratioS3 * deltaT, ratioI + ratioI3 * deltaT);
    const double ratioI4 = dotI(ratioS + ratioS3 * deltaT, ratioI + ratioI3 * deltaT);
    const double ratioR4 = dotR(ratioI + ratioI3 * deltaT);

    ratioS += (ratioS1 + 2 * ratioS2 + 2 * ratioS3 + ratioS4) * deltaT / 6;
    ratioI += (ratioI1 + 2 * ratioI2 + 2 * ratioI3 + ratioI4) * deltaT / 6;
    ratioR += (ratioR1 + 2 * ratioR2 + 2 * ratioR3 + ratioR4) * deltaT / 6;

    ++iteration;
    currentTime += deltaT;
    recentR[iteration % recentR_length] = ratioR;
}  //* End of function SIR::RK4::update

//* Check if the state reached equilibrium
bool equilibrium() {
    if (recentR[recentR_length - 1] < 0) {
        return false;
    }
    if (ratioI < 1.0 / networkSize || ratioR > 1.0 - 1.0 / networkSize) {
        return true;
    }
    const double average = std::accumulate(recentR.begin(), recentR.end(), 0.0) / recentR_length;
    for (auto& r : recentR) {
        if (fabs(r - average) > err) {
            return false;
        }
    }
    return true;
}  //* End of function SIR::RK4::equilibrium

double run(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/RK4/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Write the result
    writeFile << currentTime << "," << ratioS << "," << ratioI << "," << ratioR << "\n";
    while (!equilibrium()) {
        update();
        writeFile << currentTime << "," << ratioS << "," << ratioI << "," << ratioR << "\n";
    }
    writeFile.close();
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SIR::RK4::run
}  // namespace RK4

//* Simulate SIR model using Gillespi Algorithm
namespace GA {
//* Parameters used for GA
std::set<unsigned> reactingIndex;
double deltaT, currentTime;
int numS;
int numI;
int numR;

//* Update transition rate of single node
void updateTransitionRate(const unsigned& t_index) {
    const int intState = stateToInt[nodes[t_index].state];
    switch (intState) {
        //* S process
        case 0: {
            unsigned infectiousNeighbor = 0;
            for (const unsigned& neighbor : nodes[t_index].neighbors) {
                if (nodes[neighbor].state == "I") {
                    ++infectiousNeighbor;
                }
            }
            nodes[t_index].transitionRate = SI_II * infectiousNeighbor;
            break;
        }
        //* I process
        case 1: {
            nodes[t_index].transitionRate = I_R;
            break;
        }
        //* R Process
        case 2: {
            nodes[t_index].transitionRate = 0.0;
            break;
        }
    }
}

//* Update transition rate of all reacting nodes
void updateTransitionRate() {
    for (const unsigned& index : reactingIndex) {
        updateTransitionRate(index);
    }
}

//* Initiallize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngine;

    currentTime = 0.0;
    numS = networkSize - seedSize;
    numI = seedSize;
    numR = 0;
    reactingIndex.clear();
    for (unsigned index = 0; index < seedSize; ++index) {
        nodes[index].state = "I";
        reactingIndex.emplace(index);
        for (const unsigned& neighbor : nodes[index].neighbors) {
            reactingIndex.emplace(neighbor);
        }
    }
    updateTransitionRate();
}

//* update one step for every nodes
void syncUpdate() {
    deltaT = 1e-2;

    //* Do reactions according to each transition rate and add I into reacting
    std::set<unsigned> newReactingIndex;
    for (const unsigned& index : reactingIndex) {
        const int intState = stateToInt[nodes[index].state];
        const double transitionProb = 1.0 - std::exp(-1.0 * nodes[index].transitionRate * deltaT);
        switch (intState) {
            //* S process
            case 0: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].state = "I";
                    --numS;
                    ++numI;
                    newReactingIndex.emplace(index);
                }
                break;
            }
            //* I process
            case 1: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].state = "R";
                    --numI;
                    ++numR;
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
        }
    }

    //* Add neighbor S of I node into reacting nodes
    reactingIndex = newReactingIndex;
    for (const unsigned& index : newReactingIndex) {
        for (const unsigned& neighbor : nodes[index].neighbors) {
            if (nodes[neighbor].state == "S") {
                reactingIndex.emplace(neighbor);
            }
        }
    }

    //* Update time and transition rate
    currentTime += deltaT;
    updateTransitionRate();
}  //* End of function syncUpdate

//* update one step for single node
void asyncUpdate() {
    //* Calculate total tarnsition rate and delta time
    double totalTransitionRate = 0.0;
    for (const unsigned& index : reactingIndex) {
        totalTransitionRate += nodes[index].transitionRate;
    }
    deltaT = std::log(1.0 / probabilityDistribution(randomEngine)) / totalTransitionRate;
    totalTransitionRate *= probabilityDistribution(randomEngine);

    //* Choose target node to be reacted
    std::vector<unsigned> shuffledReactingIndex(reactingIndex.begin(), reactingIndex.end());
    std::shuffle(shuffledReactingIndex.begin(), shuffledReactingIndex.end(), randomEngine);
    double cumulativeTransitionRate = 0.0;
    unsigned target = 0;
    for (unsigned i = 0; i < shuffledReactingIndex.size(); ++i) {
        target = shuffledReactingIndex[i];
        cumulativeTransitionRate += nodes[target].transitionRate;
        if (cumulativeTransitionRate > totalTransitionRate) {
            break;
        }
    }

    //* React target node and update transition rate
    const int intState = stateToInt[nodes[target].state];
    switch (intState) {
        //* S process
        case 0: {
            nodes[target].state = "I";
            --numS;
            ++numI;
            updateTransitionRate(target);
            for (const unsigned& neighbor : nodes[target].neighbors) {
                if (nodes[neighbor].state == "S") {
                    reactingIndex.emplace(neighbor);
                    updateTransitionRate(neighbor);
                }
            }
            break;
        }
        //* I process
        case 1: {
            nodes[target].state = "R";
            --numI;
            ++numR;
            reactingIndex.erase(target);
            for (const unsigned& neighbor : nodes[target].neighbors) {
                if (nodes[neighbor].state == "S") {
                    nodes[neighbor].transitionRate -= SI_II;
                    if (!nodes[neighbor].transitionRate) {
                        reactingIndex.erase(neighbor);
                    }
                }
            }
            break;
        }
    }

    //* Update time
    currentTime += deltaT;

}  //* End of function SIR::GA::asyncUpdate

double syncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_sync/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        syncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SIR::GA::syncRun

double asyncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_async/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        asyncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SIR::GA::asyncRun
}  // namespace GA
}  // namespace SIR
