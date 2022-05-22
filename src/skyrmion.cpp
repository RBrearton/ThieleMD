#include "skyrmion.hpp"

Skyrmion::Skyrmion(std::array<double, 2> position)
    : Vortex(position)
{
}

Skyrmion::Skyrmion(double x, double y)
    : Vortex(x, y)
{
}

void Skyrmion::setParameters(int *rows, int *cols, double lambda, double hallAngle)
{
    this->ROWS = rows;
    this->COLS = cols;
    this->forceExponentialfactor = 1 / lambda;
    this->sinTh = sin((M_PI/180.0)*hallAngle);
    this->cosTh = cos((M_PI/180.0)*hallAngle);
}

void Skyrmion::actForce(double dt)
{
    // The Skyrmion should travel on a trajectory that makes the hall angle with
    // the force it experiences.


    // using dx = Fdt as the Reichhardts do, we have X_new = X + Fdt
    // TODO: make sure this compiles down to lock-free machine code

    // Rotate the force on the skyrmion depending on its hall angle.
    double rotated_force_x = (force_x * cosTh) - (force_y * sinTh);
    double rotated_force_y = (force_x * sinTh) + (force_y * cosTh);

    old_pos_x = static_cast<double>(pos_x);
    pos_x = pos_x + dt * rotated_force_x; 
    old_force_x = static_cast<double>(force_x);
    force_x = 0; // reset the force ready for the next iteration

    old_pos_y = static_cast<double>(pos_y);
    pos_y = pos_y + dt * rotated_force_y;
    old_force_y = static_cast<double>(force_y);
    force_y = 0; // reset the force ready for the next iteration
}