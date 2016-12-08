#include "integrator.h"
#include "galacticcluster.h"
#include <string.h>
#include <cmath>

Integrator::Integrator(double dt) :
    m_dt(dt)
{
}

void Integrator::integrateOneStep(GalacticCluster &system)
{
    system.calculateForcesAndEnergy();

    /* Velocity Verlet */
    for(CelestialBody &body : system.bodies()) {
        body.velocity += (m_dt/2)*(body.force / body.mass); //Calculate velocity at half step
        body.position += body.velocity*m_dt;
    }
    system.calculateForcesAndEnergy(); // updating forces at x(i+1)

    for(CelestialBody &body : system.bodies())
    {
        body.velocity += (m_dt/2)*(body.force / body.mass);
    }

}

