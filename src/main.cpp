#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <time.h>
#include "galacticcluster.h"
#include "integrator.h"
#include "vec3.h"
#include "celestialbody.h"

using namespace std;

/* TODO:
 * - cmd-arguments: N, R0?, M0?, dt, eps, numsteps             semi DONE
 * - dr^3 + eps^2*dr, check if correct in situation
 * - fix gravitational constant G.
 * - remove method dependence in integrator?                        DONE
 * - fix why file initialization doesnt work first call
 * - should write masses to file somehow                            DONE
 * - Should we have task structure? if so, how would we set it up?
 * - particles ejected
 * - atm initializing in a cube of R0^3 volume, not sphere
 *
 * Tasks:
 *  a) extend code from 3, need to find G in these units! so some math
 *    - add function to galacticcluster to determine G after adding all bodies?
 *    - find t_coll
 *  b) run system for a few t_coll, check for equilibrium time
 *  c) calculate kinetic and potential energy of the system, check conservation
 *    - check ejection from system, and its dependence on N
 *  d) add smoothing factor, try different values for energy conservation
 *    - compare with previous results
 *  e) test of virial theorem
 *  f) density of particles, mostly not c++ related
 *
 * Questions:
 *  - should mu = Mmean just be M0? Or is it important to calculate it more precisely?
 *  - rng inside galacticcluster::createcelestialbody, or just in main?
 *
 * Potential task structure:
 *  b) Runs program with just storing positions (and masses) to file
 *  c) Runs with checks for potential and kinetic energy. Calculates energy removed by ejection
 *  - both b) and c) can be run to check for the effects of varying eps, so don't need unique for that one?
 *  e) virial theorem, do we need to do anything in c++?
 *  f) density, can use b) since we write masses to file there? Yeah, python plotting and running b for various N
 */

vec3 to_cart(double r, double theta, double phi, int limits);

int main(int argc, char *argv[])
{
    if (argc<5)
    {
        cout<<"Usage: "<<argv[0]<<" N"<<" t_coll"<<" dt"<<" eps"<<endl;
        cout << "N sets the number of stars in our galactic cluster" << endl;
        cout << "Number of collision times we run the simulation" << endl; //Can this be formulated better?
        cout << "dt is the time step used" << endl;
        cout << "eps is the value of the smoothing factor in the force calculation. Should be set to a small value." << endl;
        exit(1);
    }
    int N  = atoi(argv[1]);           // number of objects
    double t_coll = atof(argv[2]);    // number of t_coll we run the simulation
    double dt  = atof(argv[3]);       // time step
    double eps = atof(argv[4]);       // smoothing factor to newton's gravitation
    int numsteps = (int)(t_coll/dt);  // number of time steps
    double masscale = 100.0/N;       // mass scaling factor for constant total mass

    // Initializing filename and determining how often we want to write to file
    string filename = "../benchmarks/pos_N"+to_string(N)+"_dt"+to_string(dt)+"_tcoll_"+to_string(t_coll)+"_eps"+to_string(eps)+".xyz";
    string filenameEnergy = "../benchmarks/energy_N"+to_string(N)+"_dt"+to_string(dt)+"_tcoll_"+to_string(t_coll)+"_eps"+to_string(eps)+".xyz";
    int numfilesteps = 100;  // max number of steps we write to file, might change later
    int iter = 1;
    if (numsteps>numfilesteps) iter = numsteps/numfilesteps;

    // Initializing random device and random number generators
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> uniform_RNG(0.0,1.0); // Uniform probability distribution
    std::normal_distribution<double> gaussian_RNG(10*masscale,1*masscale); // Gaussian probability distribution

    // Initializing system and adding bodies
    GalacticCluster galacticCluster(eps);
    for (int i=0; i<N; i++)
    {
        vec3 pos = to_cart(uniform_RNG(gen),uniform_RNG(gen),uniform_RNG(gen), 20);
        //vec3 pos(uniform_RNG(gen),uniform_RNG(gen),uniform_RNG(gen));
        galacticCluster.createCelestialBody(pos, vec3(0.0,0.0,0.0), gaussian_RNG(gen));
    }
    galacticCluster.initializer(20.0,10.0);

    // Integration loop
    Integrator integrator(dt);
    galacticCluster.writeToFilePerParticle(filename);
    galacticCluster.writeToFileEnergyPerTime(filenameEnergy);
    galacticCluster.writeToFilePerParticle(filename);
    galacticCluster.writeToFileEnergyPerTime(filenameEnergy);
    for (int step=0; step<numsteps; step++)
    {
        integrator.integrateOneStep(galacticCluster);

        if (step%(iter*10) == 0)
        {
            cout<<galacticCluster.totalEnergy()<<endl;
        }
        if (step%iter == 0)
        {
            galacticCluster.calculateEnergyPerParticle();
            galacticCluster.writeToFilePerParticle(filename);  // only write to file every iter steps
            galacticCluster.writeToFileEnergyPerTime(filenameEnergy);
        }
    }
    return 0;
}

vec3 to_cart(double r, double theta, double phi, int limits)
{
    // to ensure a sphere in cartesian coordinates

    double theta1 = theta*M_PI;
    double phi1 = phi*2*M_PI;
    double r1 = r*limits;

    double x = r1*sin(theta1)*cos(phi1);
    double y = r1*sin(theta1)*sin(phi1);
    double z = r1*cos(theta1);

    vec3 cartesian(x, y, z);

    return cartesian;
}
