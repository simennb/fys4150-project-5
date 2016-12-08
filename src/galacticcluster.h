#ifndef GALACTICCLUSTER_H
#define GALACTICCLUSTER_H

#include "celestialbody.h"
#include <vector>
#include <string>
#include <fstream>

class GalacticCluster
{
public:
    GalacticCluster(double epsilon);
    CelestialBody &createCelestialBody(vec3 position, vec3 velocity, double mass);
    void calculateForcesAndEnergy();
    void calculateEnergyPerParticle();
    void initializer(double R0, double M0);
    int numberOfBodies() const;

    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    void writeToFile(std::string filename);

    vec3 angularMomentum() const;
    vec3 momentum() const;
    std::vector<CelestialBody> &bodies();
    std::vector<CelestialBody> m_bodies;

    double eps;

private:
    vec3 m_angularMomentum;
    vec3 m_momentum;
    std::vector<double> vec_potentialEnergy;
    std::vector<double> vec_kineticEnergy;
    std::ofstream m_file;
    double m_kineticEnergy;
    double m_potentialEnergy;
    double G;
};

#endif // GALACTICCLUSTER_H
