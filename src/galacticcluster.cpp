#include "galacticcluster.h"
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

GalacticCluster::GalacticCluster(double epsilon) :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
    eps = epsilon;
}

CelestialBody& GalacticCluster::createCelestialBody(vec3 position, vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) );
    return m_bodies.back(); // Return reference to the newest added celestial body
}

void GalacticCluster::initializer(double R0, double M0)
{
    G = M_PI*M_PI*R0*R0*R0/(8.0*numberOfBodies()*M0); // assuming M0 is mean, which is at least satisfied-ish for large N
    vec_potentialEnergy.resize(numberOfBodies(), 0);
    vec_kineticEnergy.resize(numberOfBodies(), 0);
}

void GalacticCluster::calculateForcesAndEnergy()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();
    m_momentum.zeros();

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            double dr = deltaRVector.length();

            // Calculate the force
            vec3 force = -G*body1.mass * body2.mass / (dr*dr*dr + eps*eps*dr) * deltaRVector; //
            body1.force += force;
            body2.force -= force;

            // Potential energy
            m_potentialEnergy -= G * body1.mass * body2.mass / dr;
        }

        //m_momentum += body1.mass*body1.velocity;
        //m_angularMomentum += body1.position.cross(m_momentum);

        m_kineticEnergy   += 0.5*body1.mass*body1.velocity.lengthSquared();
    }
}

void GalacticCluster::calculateEnergyPerParticle()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        double temp_potentialEnergy = 0; // Temporary sum of potential energy
        double temp_kineticEnergy =   0;   // Temporary kinetic energy

        for(int j=0; j<numberOfBodies(); j++)
        {
            if (i != j)
            {
                CelestialBody &body2 = m_bodies[j];
                vec3 deltaRVector = body1.position - body2.position;
                double dr = deltaRVector.length();

                // Potential energy
                temp_potentialEnergy -= G * body1.mass * body2.mass / dr;
            }
         }

    temp_kineticEnergy   = 0.5*body1.mass*body1.velocity.lengthSquared();
    vec_kineticEnergy[i] = temp_kineticEnergy;

    vec_potentialEnergy[i] =  temp_potentialEnergy;
    }
}

int GalacticCluster::numberOfBodies() const
{
    return m_bodies.size();
}

double GalacticCluster::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double GalacticCluster::potentialEnergy() const
{
    return m_potentialEnergy;
}

double GalacticCluster::kineticEnergy() const
{
    return m_kineticEnergy;
}

void GalacticCluster::writeToFile(string filename)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    m_file << numberOfBodies() << endl;
    m_file << "Comment line that needs to be here. Balle." << endl;
    double counter = 0;
    for(CelestialBody &body : m_bodies) {
        m_file<<body.mass<<" "<< body.position.x() <<" "<< body.position.y() <<" "<< body.position.z();
        m_file<< " " << vec_potentialEnergy[counter] << " " << vec_kineticEnergy[counter] << "\n";
        counter += 1;
    }
}

vec3 GalacticCluster::angularMomentum() const
{
    return m_angularMomentum;
}

vec3 GalacticCluster::momentum() const
{
    return m_momentum;
}

std::vector<CelestialBody> &GalacticCluster::bodies()
{
    return m_bodies;
}

