#include "simulation.hpp"

Simulation::Simulation(Grid grid)
    : grid(grid),
      threadNum{0},
      io_service{},
      work(new boost::asio::io_service::work(io_service)),
      smtest("../data/smtest.dat") {}

void Simulation::setSaveCounter(int savecounter) { this->dataManager.saveCounter = savecounter; }

void Simulation::setTemperature(double temp) { this->temperature = temp; }
double Simulation::getTemperature() { return temperature; }
void Simulation::setMinimumTemp(double minTemp) { this->minimumTemp = minTemp;}
void Simulation::setCoolingRate(double cool) { this->coolingRate = cool; }
void Simulation::coolTemperature()
{
    if (temperature <= minimumTemp)
        return;
    setTemperature(temperature * coolingRate);
}
double Simulation::thermalNoise()
{
    // Setting the temperature to <= zero turns off temperature fluctuations.
    if (temperature <= 0)
        return 0;
    double tempWeighting = 0.001;
    std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> dist(0, getTemperature() * tempWeighting);
    // Create a random number from a normal distribution, where the sd is based on temp.
    return dist(generator);
}

void Simulation::initThreads(int numThreads)
{
    // initialize the threads to run in the threadpool - only call once!
    // maybe this would be better off in the constructor...
    for (int i = 0; i < numThreads; ++i)
    {
        threadpool.create_thread(std::bind(&Simulation::threadWork, this));
    }
}

void Simulation::addConstraint(std::function<void(Particle *)> constraint)
{
    // add a constraint to the system, used to define inaccessible locations
    constraints.push_back(constraint);
}

void Simulation::threadWork()
{
    // run through all the threads in the threadpool and make them do idle work
    lockCout.lock();
    std::cout << "Initialized thread number " << threadNum++ << std::endl;
    lockCout.unlock();

    io_service.run();

    lockCout.lock();
    std::cout << "Exited thread number " << --threadNum << std::endl;
    lockCout.unlock();
}

void Simulation::exitThreads()
{
    // reset the work
    work.reset();

    // join_all should now hang until all actual work has finished processing
    threadpool.join_all();
}


void Simulation::setShearForce(double shearForce, int shearForceStart, int shearForceEnd) 
{  
    this->shearForceMagnitude = shearForce;
    this->shearForceStartTime = shearForceStart;
    this->shearForceEndTime = shearForceEnd;
}

void Simulation::shearForce(Particle *particle)
{
    // Fx = 1 / (y - midpoint)^2. 
    // This leeps the force symmetrical and continuous at boundaries.
    // Also, avoid singularities (though crudely). It is almost impossible for this
    // to happen, however.
    if(particle->pos_y == grid.getRows()*0.5) particle->pos_y = particle->pos_y + 0.000001;
    double shear_force_x = 1 / (particle->pos_y - ((double)grid.getRows() * 0.5));
    shear_force_x *= shearForceMagnitude * shear_force_x;
    particle->force_x = particle->force_x + shear_force_x;

    if(particle->pos_x == grid.getCols()*0.5) particle->pos_x = particle->pos_x + 0.000001;
    double shear_force_y = 1 / (particle->pos_x - ((double)grid.getCols() * 0.5));
    shear_force_y *= shearForceMagnitude * shear_force_y;
    particle->force_y = particle->force_y + shear_force_y;
    
}

void Simulation::runTimeStep(int stepNumber)
{
    /*
    Makes interaction list, calculates forces and then moves the particles 
    according to the force acting on them. Also takes into account temperature
    effects and other external forces, as well as taking simulation geometry 
    into account by iterating over constraints.
    */
    grid.makeInteractionList();
    grid.calculateForces();

    // get all the particles in the simulation
    std::vector<Particle *> particles = grid.getParticles();

    double cols = static_cast<double>(grid.getCols());
    double rows = static_cast<double>(grid.getRows());

    // Previously, thermalNoise was in a separate function.
    // It is much faster to create the distribution here instead.
    double dt = 0.01;
    std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> dist(0, getTemperature() * dt);
    // Create a random number from a normal distribution, where the sd is based on temp.

    for (Particle *particle : particles)
    {
        // Only apply the shear force after a given time step.
        if(stepNumber >= shearForceStartTime && stepNumber < shearForceEndTime)  shearForce(particle);

        // TEMPORARY - setting dt = 0.001.
        particle->actForce(dt);

        // wiggle them due to temperature here, after the force is applied.
        // Apply gaussian noise to the particle's position due to thermal noise.
        particle->pos_x = particle->pos_x + dist(generator);
        particle->pos_y = particle->pos_y + dist(generator);

        // impose all boundary and simulation shape conditions
        for (auto constraint : this->constraints)
        {
            constraint(particle);
        }

        // Need to update each particle's cell if it has moved cell.
        // If it has moved, it needs to be removed from the previous and added to the new cell.
        // Update its cell, since it has moved.
        grid.updateParticle(particle);
    }
}

void Simulation::relax(int numSteps)
{
    // Initialize all forces and boundary conditions on the particles before the first time-step.
    std::vector<Particle *> particles = grid.getParticles();
    for (Particle *particle : particles)
    {
        // impose all boundary and simulation shape conditions
        for (auto constraint : this->constraints)
        {
            constraint(particle);
        }
        grid.updateParticle(particle);
    }

    // runTimeStep numSteps times
    for (int i = 0; i < numSteps; ++i)
    {
        runTimeStep(i);
        // Change temperature due to the cooling parameter.
        coolTemperature();
        // Save at a constant number of timesteps. Change to an input number.
        if (i % dataManager.saveCounter == 0)
        {
            std::cout <<  "Saving to file: Step " << i << ". " << std::endl;
            // Currently saves just the binary data at each given timestep.
            dataManager.saveFrame(i, grid.getParticles());
        }
    }

    // This gives the option on whether or not to convert to images.
    // For now, left out to be able to leave simulations on.

    // std::string convertToImage;
    // std::cout << "Convert Data to Images? [y/N]" << std::endl;
    // if (std::getline(std::cin, convertToImage))
    // {
    //     if (convertToImage.length() == 1 && convertToImage[0] == 'y')
    //     {
    //         // At the end of the simulation, the binary data file is reinterpreted as a series of images.
    //         dataManager.convertToImages(&grid);
    //     }
    // }
    
    dataManager.convertToImages(&grid);

    //End of Simulation
}