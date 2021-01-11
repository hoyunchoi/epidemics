#include <chrono>
#include <iostream>
#include <random>   //random_device()

#include "GSEIR.hpp"



int main(int argc, char *argv[]){
    const unsigned networkSize = std::stoul(argv[1]);
    const unsigned meanDegree = std::stoul(argv[2]);
    const unsigned linkSize = networkSize*meanDegree/2;
    const double SE_EE = std::stod(argv[3]);
    const double SI_EI = std::stod(argv[3]);
    const double E_I = std::stod(argv[4]);
    const double I_R = std::stod(argv[5]);
    const int randomEngineSeed = std::stoul(argv[6]);
    auto start=std::chrono::system_clock::now();

    //* Pre-defined variables
    const bool deletion = true;
    const int seedSize = 1;

    //* Generate Random Engine
    const int seed = randomEngineSeed==-1 ? (std::random_device())() : randomEngineSeed;
    pcg32 randomEngine;

    //* Network
    Network network;
    double orderParameter;
    const double degreeExponent = 2.5;

    using namespace GSEIR;
    //* RR-RK4
    randomEngine.seed(seed); network = RR::generate(networkSize, meanDegree, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = RK4::run(seed, randomEngine, deletion);

    //* RR-GA_sync
    randomEngine.seed(seed); network = RR::generate(networkSize, meanDegree, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::syncRun(seed, randomEngine, deletion);

    //* RR-GA_async
    randomEngine.seed(seed); network = RR::generate(networkSize, meanDegree, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::asyncRun(seed, randomEngine, deletion);

    //* ER-GA_sync
    randomEngine.seed(seed); network = ER::generate(networkSize, linkSize, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::syncRun(seed, randomEngine, deletion);

    //* ER-GA_async
    randomEngine.seed(seed); network = ER::generate(networkSize, linkSize,  randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::asyncRun(seed, randomEngine, deletion);

    //* SF-GA_sync
    randomEngine.seed(seed); network = SF::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::syncRun(seed, randomEngine, deletion);

    //* SF-GA_async
    randomEngine.seed(seed); network = SF::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::asyncRun(seed, randomEngine, deletion);

    //* CL-GA_sync
    randomEngine.seed(seed); network = CL::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::syncRun(seed, randomEngine, deletion);

    //* CL-GA_async
    randomEngine.seed(seed); network = CL::generate(networkSize, linkSize, degreeExponent, randomEngine);
    setNetwork(network); setRate(SE_EE, SI_EI, E_I, I_R, seedSize);
    orderParameter = GA::asyncRun(seed, randomEngine, deletion);

    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf("%0.10f second\n", sec.count());

    return 0;
}