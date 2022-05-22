#pragma once

#include "particle.hpp"
#include <array>
#include <vector>

class Cell
{
private:
  std::array<int, 2> coords;         // x, y coordinates of this Cell
  std::vector<Particle *> particles; // pointers to particles in this Cell

public:
  Cell(int x, int y);
  Cell(std::array<int, 2> coordinates);
  Cell(int x, int y, std::vector<Particle *> particles);
  Cell(std::array<int, 2> coordinates, std::vector<Particle *> particles);

  std::array<int, 2> getCoords();
  std::vector<Particle *> getParticles();

  void addParticle(Particle *particle);
  void removeParticle(Particle *particle);
  void removeParticle(int particleIndex);
};