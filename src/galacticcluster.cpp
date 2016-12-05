#include "galacticcluster.h"
#include <iostream>
#include <iomanip>
using namespace std;

GalacticCluster::GalacticCluster(double epsilon) :
    m_kineticEnergy(0),
    m_potentialEnergy(0)
{
    eps = epsilon;
    pi = 3.14159265358979323846;
}

CelestialBody& GalacticCluster::createCelestialBody(vec3 position, vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) );
    return m_bodies.back(); // Return reference to the newest added celestial body
}

void GalacticCluster::gravitationalConstant(int N, double R0, double M0)
{
    G = pi*pi*R0*R0*R0/(8.0*N*M0); // assuming M0 is mean, which is at least satisfied-ish for large N
    cout<<G<<endl;
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

        m_momentum += body1.mass*body1.velocity;
        m_angularMomentum += body1.position.cross(m_momentum);
        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
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
    for(CelestialBody &body : m_bodies) {
        m_file<<body.mass<<" "<< body.position.x() <<" "<< body.position.y() <<" "<< body.position.z() << "\n";
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

