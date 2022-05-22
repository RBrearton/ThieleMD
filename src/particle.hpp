#pragma once
#include <array>
#include <atomic>

class Particle
{
protected:
  int *ROWS;
  int *COLS;

public:
  int id; //index for this particle. Used when saving.

  // atomic variables have no copy constructor -> they must be public
  // std collections require elements to be copyable, so  cant store in arrays
  std::atomic<double> pos_x;       // current x position
  std::atomic<double> pos_y;       // current y position
  std::atomic<double> old_pos_x;   // previous x position, currently unused
  std::atomic<double> old_pos_y;   // previous y position, currently unused
  std::atomic<double> force_x;     // current force along x direction
  std::atomic<double> force_y;     // current force along y direction
  std::atomic<double> old_force_x; // previous force along x direction
  std::atomic<double> old_force_y; // previous force along y direction

  Particle();                          // unused
  Particle(std::array<double, 2> pos); // inits position
  Particle(double x, double y);        // inits position

  virtual void interact(Particle *particle) = 0;
  virtual void interactWall() = 0;

  virtual void actForce(double dt);
};