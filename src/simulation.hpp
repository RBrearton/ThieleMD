#pragma once

#include <atomic>
#include <fstream>
#include <memory>
#include <random>

#include "datamanager.hpp"
#include "grid.hpp"

class Simulation
{
private:
  DataManager dataManager;
  Grid grid;
  std::vector<std::function<void(Particle *)>> constraints;

  double temperature = 0;
  double minimumTemp = 0;
  // Value by which temp is multiplied each timestep.
  // Set to zero for no cooling, <1 for cooling, or >1 for heating.
  double coolingRate = 0.4;
  double drag;
  double numParticles; // the number of particles in the grid
  double shearForceMagnitude = 0;
  // The step number at which the shear force is introduced.
  int shearForceStartTime = 0;
  int shearForceEndTime = 0;

  // Not currently in use.
  std::atomic<int> threadNum;

  std::ofstream smtest;

  // Generates a random number to adjust position by, based on the temperature.
  double thermalNoise();

public:
  Simulation(Grid grid);

  void setSaveCounter(int savecounter);

  void setTemperature(double temperature);
  double getTemperature();
  void setMinimumTemp(double minTemp);
  void setCoolingRate(double cool);
  void coolTemperature(); // Adjusts temperature by coolingRate.

  void initThreads(int numThreads);
  void exitThreads();

  void setShearForce(double shearForce, int shearForceStart, int shearForceEnd);
  void shearForce(Particle *particle);

  void addConstraint(std::function<void(Particle *)> constraint);
  void relax();
  void relax(int numSteps);

  void runTimeStep(int stepNumber);
};