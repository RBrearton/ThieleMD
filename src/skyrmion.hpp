#pragma once
#include "vortex.hpp"
#include <array>
#include <cmath>

class Skyrmion : public Vortex
{
private:
    // // Determines the offset of the Skyrmions trajectory
    // double hallAngle = 0; 

    // Sin and Cos of the hallAngle. These are saved after initial parameters
    // are set so that they need not be calculated at run time for each step.
    double sinTh;
    double cosTh;

public:
  using Vortex::Vortex;

  Skyrmion(std::array<double, 2> pos);
  Skyrmion(double x, double y);

  // Copy of setParameters to allow hall angle.
  void setParameters(int *rows, int *cols, double lambda, double hallAngle);

  // Change force to push skyrmion in direction of hall angle.
  virtual void actForce(double dt) override;
};