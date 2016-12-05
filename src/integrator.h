#ifndef INTEGRATOR_H
#define INTEGRATOR_H


class Integrator
{
public:
    double m_dt;
    Integrator(double dt);
    void integrateOneStep(class GalacticCluster &system);
};

#endif // INTEGRATOR_H
