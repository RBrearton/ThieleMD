#include "vortex.hpp"

Vortex::Vortex(std::array<double, 2> position)
    : Particle(position)
{
}

Vortex::Vortex(double x, double y)
    : Particle(x, y)
{
}

void Vortex::setParameters(int *rows, int *cols, double lambda)
{
    this->ROWS = rows;
    this->COLS = cols;
    this->forceExponentialfactor = 1 / lambda;
}

void Vortex::interact(Particle *particle)
{
    /*
        Calculates the interaction force between two particles, assuming that
        the force looks like exp(x1 - x2). As this models a vortex and not a
        skyrmion, there is no force rotation
        In grid.cpp, we assume a particle is close if the separation is < 1. Here,
        we use that to impose periodice boundary conditions: if the particles
        have a separation >1, then we take the negative of the force.
    */
    std::array<double, 2> force;

    // work out the distance between the two particles
    double dx = this->pos_x - particle->pos_x;
    dx = dx * xIsotropy; // Include interaction isotropy.
    double dx_2 = dx * dx;
    // Account for periodic distances by observing whether separation is >1.
    if (dx_2 > 1.0)
    {
        // Too far away so calculate periodic distance.
        // Here we want dx to flip sign, so that the forces are in the correct directions.
        // Note that the magnitude of max_x should always be bigger than dx.
        // Also, account for xIsotropy.
        dx = (dx > 0) ? dx - (xIsotropy * (*ROWS)) : dx + (xIsotropy * (*ROWS));
        dx_2 = dx * dx;
    }

    double dy = this->pos_y - particle->pos_y;
    double dy_2 = dy * dy;
    // Do same for periodic y.
    if (dy_2 > 1.0)
    {
        dy = (dy > 0) ? dy - (*COLS) : dy + (*COLS);
        dy_2 = dy * dy;
    }

    double dist = std::sqrt(dx_2 + dy_2); // distance between particles

    // magnitude of the force looks like exp(-distance_between_particles)
    double mag_force = forceExponentialfactor * std::exp(-dist * forceExponentialfactor);
    // For an input of lambda, just set forceExponentialfactor to 1/lambda input.

    // from this construct the force vector
    force[0] = mag_force * dx / dist;
    force[1] = mag_force * dy / dist;

    // append this force to this->force and particle-> force
    this->force_x = this->force_x + force[0];
    particle->force_x = particle->force_x - force[0];

    this->force_y = this->force_y + force[1];
    particle->force_y = particle->force_y - force[1];
}

void Vortex::interactWall()
{
    // Add force dependent on displacement from boundary.
    // This force is used in the following timestep, so we actually miss
    // the wall interaction on the first time step with this method...
    // Using 1.0 as max interaction range, just like the rest of the code.

    if (this->pos_x < 1.0)
    {
        double force_x = std::exp(-this->pos_x * forceExponentialfactor);
        this->force_x = this->force_x + force_x;
    }
    else if (this->pos_x > *ROWS - 1.0)
    {
        double dx = *ROWS - this->pos_x;
        double force_x = std::exp(-dx * forceExponentialfactor);
        this->force_x = this->force_x - force_x;
    }

    if (this->pos_y < 1.0)
    {
        double force_y = std::exp(-this->pos_y * forceExponentialfactor);
        this->force_y = this->force_y + force_y;
    }
    else if (this->pos_y > *COLS - 1.0)
    {
        double dy = *COLS - this->pos_y;
        double force_y = std::exp(-dy * forceExponentialfactor);
        this->force_y = this->force_y - force_y;
    }
}