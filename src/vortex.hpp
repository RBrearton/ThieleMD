#pragma once
#include "particle.hpp"
#include <array>
#include <cmath>

class Vortex : public Particle
{
protected:
  double forceExponentialfactor = 1;

  // The xIsotropy variable introduces isotropy into the system, defining
  // a preferred direction. This reduces domain formation.
  // 1.0 means no isotropy, other than that there will be isotropy.
  // It is multipled to dx before calculating the distance.
  double xIsotropy = 1.1;

public:
  using Particle::Particle;

  Vortex(std::array<double, 2> pos);
  Vortex(double x, double y);

  void setParameters(int *rows, int *cols, double lambda);

  //virtual void actForce(double dt);
  
  virtual void interact(Particle *particle) override;
  virtual void interactWall() override;
};