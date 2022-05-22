#include "particle.hpp"

Particle::Particle() {}
Particle::Particle(std::array<double, 2> pos)
    : pos_x{pos[0]}, pos_y{pos[1]}
{
}
Particle::Particle(double x, double y)
    : pos_x{x}, pos_y{y} {}

void Particle::actForce(double dt)
{
    // using dx = Fdt as the Reichhardts do, we have X_new = X + Fdt
    // TODO: make sure this compiles down to lock-free machine code
    old_pos_x = static_cast<double>(pos_x);
    pos_x = pos_x + dt * force_x;
    old_force_x = static_cast<double>(force_x);
    force_x = 0; // reset the force ready for the next iteration

    old_pos_y = static_cast<double>(pos_y);
    pos_y = pos_y + dt * force_y;
    old_force_y = static_cast<double>(force_y);
    force_y = 0; // reset the force ready for the next iteration
}