#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "specincld.h"

class Integrator
{
 public:
  Integrator(vector<double> GMin, vector<double> stin, double tin, double dtin);

  vector<double> GM;
  vector<double> state;
  double t;
  double official_delta_t;

  int n_equations;
  int n_planets;
  double epsilon;
  int max_n_planets;
  int max_n_equations;
  vector<double> deriv;

  void integrate();
  vector<double> getState();
  double getTime();
  double getDeltaT();
  void setDeltaT(double newdeltat);
  vector<double> equations_of_motion(vector<double> local_state);
};

#endif // INTEGRATOR_H
