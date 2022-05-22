#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <experimental/filesystem>
#include <fstream>
#include <math.h>

#include "datamanager.hpp"
#include "grid.hpp"
#include "particle.hpp"
#include "simulation.hpp"
#include "vortex.hpp"
#include "skyrmion.hpp"

int ROWS = 6;
int COLS = 6;

class Timer
{
private:
    // type aliases to make accessing nested type easier
    using clock_t = std::chrono::high_resolution_clock;
    using ns_t = std::chrono::high_resolution_clock::duration;
    using second_t = std::chrono::duration<double, std::ratio<1>>;

    std::chrono::time_point<clock_t> m_beg;

public:
    Timer() : m_beg(clock_t::now())
    {
    }

    void reset()
    {
        m_beg = clock_t::now();
    }

    double elapsed() const
    {
        return std::chrono::duration_cast<ns_t>(clock_t::now() - m_beg).count();
    }
};

void stickyBoundary(Particle *particle)
{
    // particle sticks to wall if it should pass it.
    if (particle->pos_x > COLS)
    {
        // Can't set it to COLS, due to how we label the cell sizes.
        // If its position was =COLS then it would be outside the grid.
        particle->pos_x = COLS - 0.0001;
    }
    else if (particle->pos_x < 0)
    {
        particle->pos_x = 0;
    }

    if (particle->pos_y > ROWS)
    {
        particle->pos_y = ROWS - 0.0001;
    }
    else if (particle->pos_y < 0)
    {
        particle->pos_y = 0;
    }
}

void bouncyBoundary(Particle *particle)
{
    // particle bounces from the wall simply, quite crude this way though.
    if (particle->pos_x > COLS)
    {
        particle->pos_x = (2 * COLS) - particle->pos_x;
    }
    else if (particle->pos_x < 0)
    {
        particle->pos_x = -particle->pos_x;
    }

    if (particle->pos_y > ROWS)
    {
        particle->pos_y = (2 * ROWS) - particle->pos_y;
    }
    else if (particle->pos_y < 0)
    {
        particle->pos_y = -particle->pos_y;
    }
}

void periodicBoundary(Particle *particle)
{
    // This boundary condition is recursively called until
    // the given particle is back within the box.
    if(particle->pos_x >= 0 && particle->pos_x < COLS
    && particle->pos_y >= 0 && particle->pos_y < ROWS)
    {
        // Stop the recrusive calling.
        return;
    }

    if (particle->pos_x >= COLS)
    {
        particle->pos_x = particle->pos_x - COLS;
    }
    else if (particle->pos_x < 0)
    {
        particle->pos_x = particle->pos_x + COLS;
    }

    if (particle->pos_y >= ROWS)
    {
        particle->pos_y = particle->pos_y - ROWS;
    }
    else if (particle->pos_y < 0)
    {
        particle->pos_y = particle->pos_y + ROWS;
    }

    // Recursively call this constraint until the particle is 
    // back within the box.
    periodicBoundary(particle);
}

void twistBoundary(Particle *particle)
{
    // particle sticks to wall if it should pass it.
    if (particle->pos_x > COLS)
    {
        // Can't set it to COLS, due to how we label the cell sizes.
        // If its position was =COLS then it would be outside the grid.
        particle->pos_x = COLS - 0.0001;
    }
    else if (particle->pos_x < 0)
    {
        particle->pos_x = 0;
    }

    if (particle->pos_y > ROWS)
    {
        particle->pos_y = ROWS - 0.0001;
    }
    else if (particle->pos_y < 0)
    {
        particle->pos_y = 0;
    }

    particle->interactWall();
}

//Method to add any number of random particles to system.
void addRandomVortices(Grid *grid, int numParticles, double lambda)
{
    /* 
    We create an array of vortex addresses for the number of random particles we want.
    Then we create that many particles at random x and y positions, then add them to 
    the grid.
    */
    srand(time(NULL));
    Vortex *vortices[numParticles];
    for (int i = 0; i < numParticles; ++i)
    {
        //Vortex v(((float)std::rand() / (float) RAND_MAX) * ROWS, ((float)std::rand() / (float) RAND_MAX) * COLS);
        vortices[i] = new Vortex(((float)std::rand() / (float)RAND_MAX) * ROWS,
                                 ((float)std::rand() / (float)RAND_MAX) * COLS);
        vortices[i]->setParameters(&ROWS, &COLS, lambda);
        grid->addParticle(vortices[i]);
    }
}

void addRandomSkyrmions(Grid *grid, int numParticles, double lambda, double hallAngle)
{
    /* 
    We create an array of vortex addresses for the number of random particles we want.
    Then we create that many particles at random x and y positions, then add them to 
    the grid.
    */
    srand(time(NULL));
    Skyrmion *skyrmions[numParticles];
    for (int i = 0; i < numParticles; ++i)
    {
        //Vortex v(((float)std::rand() / (float) RAND_MAX) * ROWS, ((float)std::rand() / (float) RAND_MAX) * COLS);
        skyrmions[i] = new Skyrmion(((float)std::rand() / (float)RAND_MAX) * ROWS,
                                 ((float)std::rand() / (float)RAND_MAX) * COLS);
        skyrmions[i]->setParameters(&ROWS, &COLS, lambda, hallAngle);
        grid->addParticle(skyrmions[i]);

        //TEMP: Square grid for 576 particles
        //         skyrmions[i] = new Skyrmion(((i%24)/4.0)+0.1,
        //                          ((i/24)/4.0)+0.1);
        //                          //std::cout << i%ROWS << " " << i/COLS << std::endl;
        // skyrmions[i]->setParameters(&ROWS, &COLS, lambda, hallAngle);
        // grid->addParticle(skyrmions[i]);
    }
}

int main()
{
    Timer timer; // for debugging
    // putting particles into the grid here

    // make some vortices
    // Vortex v1(3.05, 3.05);

    // add them to the grid
    // grid.addParticle(&v1);

    // take in the inputs from inputs.text file
    std::ifstream input("../inputs.txt");
    if (!input)
    {
        std::cout << "Could not find inputs." << std::endl;
        return 0;
    }
    std::string buffer;
    int simulationLength;
    int numberparticles;
    double temperature;
    double minTemp;
    double coolingRate;
    double lambda;
    int savecounter;
    int periodic;
    double hallAngle;
    double shearForceMag;
    int shearForceStart;
    int shearForceEnd;

    input.ignore(256, ' ');
    input >> ROWS;
    input.ignore(256, ' ');
    input >> COLS;
    input.ignore(256, ' ');
    input >> simulationLength;
    input.ignore(256, ' ');
    input >> numberparticles;
    input.ignore(256, ' ');
    input >> temperature;
    input.ignore(256, ' ');
    input >> minTemp;
    input.ignore(256, ' ');
    input >> coolingRate;
    input.ignore(256, ' ');
    input >> lambda;
    input.ignore(256, ' ');
    input >> savecounter;
    input.ignore(256, ' ');
    input >> periodic;
    input.ignore(256, ' ');
    input >> hallAngle;
    input.ignore(256, ' ');
    input >> shearForceMag;
    input.ignore(256, ' ');
    input >> shearForceStart;
    input.ignore(256, ' ');
    input >> shearForceEnd;
    input.close();
    // need lambda for vortex, savetimer in datamanager?

    Grid grid(ROWS, COLS);

    //addRandomVortices(&grid, numberparticles, lambda);
    addRandomSkyrmions(&grid, numberparticles, lambda, hallAngle);
    
    std::cout << "There are "
              << grid.getNumParticles() << " particles in the system."
              << std::endl;

    // grid.makeInteractionList();

    if (periodic == 1)
    {
        // Grid needs to be set before simulation is made.
        grid.setPeriodic(true);
    }
    else
    {
        grid.setPeriodic(false);
    }

    // testing threading
    Simulation simulation(grid);
    simulation.setSaveCounter(savecounter);
    simulation.setShearForce(shearForceMag, shearForceStart, shearForceEnd);

    if (periodic == 1)
    {
        //constraints need to be added after simulation is made.
        simulation.addConstraint(&periodicBoundary);
    }
    else
    {
        simulation.addConstraint(&twistBoundary);
    }
    // std::cout << grid.getInteractionList().size() << std::endl;

    // assign whether the simulation is periodic, including constraints.

    // simulation.initThreads(2);
    // this is where the simulation goes

    simulation.setTemperature(temperature);
    simulation.setMinimumTemp(minTemp);
    simulation.setCoolingRate(coolingRate);
    simulation.relax(simulationLength);

    // this is where the simulation ends
    //simulation.exitThreads();

    // simulation.lockCout.lock();
    std::cout << "Executed in " << timer.elapsed() / 1e9 << "s" << std::endl;
    // simulation.lockCout.unlock();

    return 0;
}