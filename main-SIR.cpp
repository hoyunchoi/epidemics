#include <chrono>
#include <iostream>
#include <random>   //random_device()

#include "../library-Git/Epidemics.hpp"

int main(int argc, char *argv[]){
    const unsigned networkSize = std::stoul(argv[1]);
    const unsigned meanDegree = std::stoul(argv[2]);
    const unsigned linkSize = networkSize*meanDegree/2;
    const double SI_II = std::stod(argv[3]);
    const double I_R = std::stod(argv[4]);
    const int randomSeed = std::stoi(argv[5]);
    auto start=std::chrono::system_clock::now();

    //* Pre-defined variables
    const bool deletion = false;
    const int seedSize = 1;

    //* Generate Random Engine
    const int randomEngineSeed = randomSeed==-1 ? (std::random_device())() : randomSeed;
    pcg32 randomEngine;

    //* declare network parameters
    Network network;
    double orderParameter;
    const double degreeExponent = 2.5;

    using namespace SIR;
    //* RR-RK4
    randomEngine.seed(randomEngineSeed); network = RR::generate(networkSize, meanDegree, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = RK4::run(randomEngineSeed, randomEngine, deletion);

    //* RR-GA_sync
    randomEngine.seed(randomEngineSeed); network = RR::generate(networkSize, meanDegree, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::syncRun(randomEngineSeed, randomEngine, deletion);

    //* RR-GA_async
    randomEngine.seed(randomEngineSeed); network = RR::generate(networkSize, meanDegree, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::asyncRun(randomEngineSeed, randomEngine, deletion);

    //* ER-GA_sync
    randomEngine.seed(randomEngineSeed); network = ER::generate(networkSize, linkSize, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::syncRun(randomEngineSeed, randomEngine, deletion);

    //* ER-GA_async
    randomEngine.seed(randomEngineSeed); network = ER::generate(networkSize, linkSize, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::asyncRun(randomEngineSeed, randomEngine, deletion);

    //* SF-GA_sync
    randomEngine.seed(randomEngineSeed); network = SF::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::syncRun(randomEngineSeed, randomEngine, deletion);

    //* SF-GA_async
    randomEngine.seed(randomEngineSeed); network = SF::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::asyncRun(randomEngineSeed, randomEngine, deletion);

    //* CL-GA_sync
    randomEngine.seed(randomEngineSeed); network = CL::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::syncRun(randomEngineSeed, randomEngine, deletion);

    //* CL-GA_async
    randomEngine.seed(randomEngineSeed); network = CL::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SI_II, I_R, seedSize);
    orderParameter = GA::asyncRun(randomEngineSeed, randomEngine, deletion);

    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf("%0.10f second for SIR model.\n", sec.count());

    return 0;
}