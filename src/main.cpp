#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <random>
//#include <armadillo>
#include <time.h>
#include "galacticcluster.h"
#include "integrator.h"
#include "vec3.h"
#include "celestialbody.h"

using namespace std;
//using namespace arma;

/* TODO:
 * - cmd-arguments: N, R0?, M0?, dt, eps, numsteps
 * - dr^3 + eps^2*dr, check if correct in situation
 * - remove method dependence in integrator?
 * - fix why file initialization doesnt work first call
 * - should write masses to file somehow
 *
 *
 */

int main(int argc, char *argv[])
{
    int N  = 10; // number of objects
    int numsteps = 1000;
    double dt = 0.01;
    double R0 = 20.0;
    double eps = 0.0;
    char const *method = "Verlet";

    // Initializing filename
    string filename = "../benchmarks/pos_dt"+to_string(dt)+"_steps"+to_string(numsteps)+".xyz";

    // Initializing random device and random number generators
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> uniform_RNG(0.0,R0); // Uniform probability distribution
    std::normal_distribution<double> gaussian_RNG(10,1); // Gaussian probability distribution

    // Initializing system and adding bodies
    GalacticCluster galacticCluster(eps);
    for (int i=0; i<N; i++)
    {
        vec3 pos(uniform_RNG(gen),uniform_RNG(gen),uniform_RNG(gen));
        galacticCluster.createCelestialBody(pos, vec3(0.0,0.0,0.0), gaussian_RNG(gen));
    }

    // Integration loop
    Integrator integrator(dt,method);
    galacticCluster.writeToFile(filename);
    for (int step=0; step<numsteps; step++)
    {
        integrator.integrateOneStep(galacticCluster);

        if (step%100 == 0) galacticCluster.writeToFile(filename);
    }

    return 0;
}
