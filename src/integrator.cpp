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
        vec3 accel1 = body.force/body.mass; // acceleration before updating x(i) to x(i+1)
        body.position += body.velocity*m_dt + 0.5 * accel1 * pow(m_dt,2);

        system.calculateForcesAndEnergy(); // updating forces at x(i+1)
        body.velocity += 0.5* (accel1 + body.force / body.mass) * m_dt;
    }
}

