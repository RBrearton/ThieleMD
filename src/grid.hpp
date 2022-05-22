#pragma once

#include <array>
#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <vector>

#include "cell.hpp"
#include "particle.hpp"

class Grid
{
private:
  // checks to see if particles part_1 and part_2 are close enough to interact
  bool is_close(Particle *part_1, Particle *part_2);
  // as above, but with periodic boundary conditions
  bool is_close_periodic(Particle *part_1, Particle *part_2);
  // as above, but in a system with circular symmetry
  bool is_close_circular(Particle *part_1, Particle *part_2);

  //holds which function to use, is_close or is_close_periodic depending on
  // whether the system is periodic.
  std::function<bool(Particle *, Particle *)> is_close_function =
      std::bind(&Grid::is_close_periodic, this, std::placeholders::_1, std::placeholders::_2);

protected:
  int rows;
  int cols;
  std::vector<std::vector<Cell>> cells;

  std::vector<std::array<Particle *, 2>> interactionList;

public:
  Grid(int rows, int cols);
  Grid(std::vector<std::vector<Cell>> cells);

  void setPeriodic(bool isperiodic);

  int getRows();         // get number of rows in the grid
  int getCols();         // get number of cols in the grid
  int getNumParticles(); // calculates the number of particles in the grid
  std::vector<std::vector<Cell>> getCells();

  // adds a particle to the grid, throws an error if particle is out of bounds
  void addParticle(Particle *particle);
  // updates the particle's cell if it left its previous cell.
  void updateParticle(Particle *particle);
  // work out which particles are close enough to interact with each other
  void makeInteractionList();
  std::vector<std::array<Particle *, 2>> getInteractionList();

  // gets all the particles in the grid
  std::vector<Particle *> getParticles();

  // runs through the interaction list and works out the force acting on each
  // particle due to other particles in the simulation
  void calculateForces();
};