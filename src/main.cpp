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
    
    // Writing the first values to file. Doubled up due to the file opening method not working on first call for unknown reason.
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
