#include "grid.hpp"

bool Grid::is_close_periodic(Particle *part_1, Particle *part_2)
{
    // THIS USES "NORMAL" PERIODIC BOUNDARY CONDITIONS!!!!
    // checks to see if particles part_1 and part_2 are close enough to interact

    // check if distance between particles is 1
    // impose periodic boundary conditions!
    double dx = part_1->pos_x - part_2->pos_x;
    if (dx < 0)
        dx = -dx;
    if (dx > static_cast<double>(cols) / 2.0)
        dx = static_cast<double>(cols) - dx;

    // also impose periodicity in the y direction
    double dy = part_1->pos_y - part_2->pos_y;
    if (dy < 0)
        dy = -dy;
    if (dy > static_cast<double>(rows) / 2.0)
        dy = static_cast<double>(rows) - dy;

    double dist_squared = dx * dx + dy * dy;

    // no need to sqrt! If particles are too far, return false
    if (dist_squared > 1.0 || dist_squared == 0)
        return false;

    return true;
}

bool Grid::is_close(Particle *part_1, Particle *part_2)
{
    // THIS USES NO PARTICULAR BOUNDARY CONDITIONS
    double dx = part_1->pos_x - part_2->pos_x;
    double dy = part_1->pos_y - part_2->pos_y;

    double dist_sq = dx * dx + dy * dy;

    // no need to sqrt! If particles are further apart than 1, return false
    if (dist_sq > 1.0 || dist_sq == 0)
        return false;
    return true;
}

void Grid::setPeriodic(bool periodic)
{
    if (periodic)
        is_close_function =
            std::bind(&Grid::is_close_periodic, this, std::placeholders::_1, std::placeholders::_2);
    else
        is_close_function =
            std::bind(&Grid::is_close, this, std::placeholders::_1, std::placeholders::_2);
}

Grid::Grid(int rows, int cols) : rows{rows}, cols{cols}
{
    cells.resize(cols);
    for (int i = 0; i < cols; ++i)
    {
        cells[i].reserve(rows);
        for (int j = 0; j < rows; ++j)
        {
            cells[i].push_back(Cell(i, j));
        }
    }
}

Grid::Grid(std::vector<std::vector<Cell>> cells) : cells{cells}
{
    cols = cells.size();
    rows = cells[0].size();
}

int Grid::getRows() { return rows; }
int Grid::getCols() { return cols; }
int Grid::getNumParticles()
{
    int counter = 0;
    for (auto i : cells)
        for (auto j : i)
        {
            counter += j.getParticles().size();
        }
    return counter;
}

std::vector<std::vector<Cell>> Grid::getCells() { return cells; }

void Grid::addParticle(Particle *particle)
{
    // first work out if the particle is in a legal position
    std::array<double, 2> position{particle->pos_x, particle->pos_y};

    // work out which cell the particle is in
    std::array<int, 2> cellPos = {static_cast<int>(position[0]),
                                  static_cast<int>(position[1])};

    // throw an invalid_argument error if the particle is in the wrong place
    if (cellPos[0] >= cols || cellPos[0] < 0 ||
        cellPos[1] >= rows || cellPos[1] < 0)
    {
        throw std::invalid_argument("Particle passed to addParticle has invalid"
                                    " position");
    }

    // the particle should be in a legit position if execution reaches here
    // Give the particle a unique id, once it has been added to the grid.
    particle->id = getNumParticles();
    cells[cellPos[0]][cellPos[1]].addParticle(particle);
}

void Grid::updateParticle(Particle *particle)
{
    // first work out if the particle is in a legal position
    std::array<double, 2> position{particle->pos_x, particle->pos_y};

    // work out which cell the particle is in now
    std::array<int, 2> cellPos = {static_cast<int>(position[0]),
                                  static_cast<int>(position[1])};
    // work out its previous cell position.
    std::array<int, 2> prevCellPos = {static_cast<int>(particle->old_pos_x),
                                      static_cast<int>(particle->old_pos_y)};
    if (cellPos[0] == prevCellPos[0] && cellPos[1] == prevCellPos[1])
        return;

    // throw an invalid_argument error if the particle is in the wrong place
    if (cellPos[0] >= cols || cellPos[0] < 0 ||
        cellPos[1] >= rows || cellPos[1] < 0)
    {
        throw std::invalid_argument("Particle passed to updateParticle has invalid"
                                    " position");
    }
    // the particle should be in a legitimate position if execution reaches here
    // Remove the particle from its previous cell now and add it to its new one.
    cells[cellPos[0]][cellPos[1]].addParticle(particle);
    cells[prevCellPos[0]][prevCellPos[1]].removeParticle(particle);
}

void Grid::makeInteractionList()
{
    /*
    Iterate over each cell in the grid. If the cell contains a particle, check 
    if any nearest neighbour cells also contain a particle. If they do, add the
    particle pair to the interaction list.
    */

    // for now the implementation is single-threaded
    // TODO: optimization, multithreading, different boundary conditions

    // loop over all the cells
    int num_cols{cells.size()};
    for (int i = 0; i < num_cols; ++i)
    {
        int num_rows{cells[i].size()};
        for (int j = 0; j < num_rows; ++j)
        {
            // further optimizations are possible! currently not "remembering"
            // previously empty cells

            // check if there's a particle in the cell
            std::vector<Particle *> particles(cells[i][j].getParticles());

            // if the cell is empty, continue
            if (particles.empty())
                continue;

            // currently assuming periodic boundary conditions!
            // check neighbouring cells for particles

            // define up, down, left, right cells separately
            int right = i + 1;
            int left = i - 1;
            int up = j + 1;
            int down = j - 1;

            // avoid overflow!
            if (right == num_cols)
                right = 0;
            if (left == -1)
                left == num_cols - 1;
            if (up == num_rows)
                up = 0;
            if (down == -1)
                down = num_rows - 1;

            std::vector<Particle *> nearParticles; // init empty storage vectors
            std::vector<Particle *> temp;

            // check "top right" neighbour
            temp = cells[right][up].getParticles();
            nearParticles.insert(nearParticles.end(), temp.begin(), temp.end());

            // check "right" neighbour
            temp = cells[right][j].getParticles();
            nearParticles.insert(nearParticles.end(), temp.begin(), temp.end());

            // check "bottom centre" neighbour
            temp = cells[i][down].getParticles();
            nearParticles.insert(nearParticles.end(), temp.begin(), temp.end());

            // check "bottom right" neighbour
            temp = cells[right][down].getParticles();
            nearParticles.insert(nearParticles.end(), temp.begin(), temp.end());

            // if nearParticles is empty, dont loop
            if (!nearParticles.empty())
            {
                // loop through all nearParticles and all particles, pairing them
                // and adding them to the interaction list
                for (auto part_1 : particles)
                {
                    for (auto part_2 : nearParticles)
                    {
                        // make sure the particles are close enough to interact
                        // Here we can choose whether the simulation is periodic or not,
                        // by choosing either is_close or is_close_periodic since
                        // all other periodicness depends on this.
                        if (is_close_function(part_1, part_2))
                        {
                            // make the particle pair
                            std::array<Particle *, 2>
                                interact_pair{part_1, part_2};

                            // add the particle pair to the interaction list
                            interactionList.push_back(interact_pair);
                        }
                    }
                }
            }

            // now deal with interactions between particles in the same cell
            // dont bother if there's only one particle in there
            if (particles.size() == 1)
                continue;

            // loop over all particles in the cell avoiding repeats
            for (int i = 0; i < particles.size(); ++i)
            {

                for (int j = i; j < particles.size(); ++j)
                {
                    // make sure they're close enough to actually interact
                    if (is_close_periodic(particles[i], particles[j]))
                    {
                        // make the particle pair
                        std::array<Particle *, 2>
                            interact_pair{particles[i], particles[j]};

                        interactionList.push_back(interact_pair);
                    }
                }
            }
        }
    }
}

std::vector<std::array<Particle *, 2>> Grid::getInteractionList()
{
    return interactionList;
}

void Grid::calculateForces()
{
    // runs through the interaction list and works out the force acting on each
    // particle due to other particles in the simulation

    // make sure that the interaction list isn't empty
    if (interactionList.size() == 0)
        std::cout << "Empty interaction list!" << std::endl;

    // loop over every particle pair in the interaction list
    for (std::array<Particle *, 2> pair : interactionList)
    {
        // only call interact on the first particle, passing the 2nd particle
        // as the argument to interact
        pair[0]->interact(pair[1]); // this affects pair[0]->force_x,
                                    // pair[1]->force_x, etc.
    }

    // the last thing that should be done is to clear the interaction list
    interactionList.clear();
}

std::vector<Particle *> Grid::getParticles()
{
    std::vector<Particle *> particles;
    // iterate through all cells in grid and grab all their particles
    for (std::vector<Cell> cell_row : cells)
    {
        for (Cell cell : cell_row)
        {
            std::vector<Particle *> temp = cell.getParticles();
            particles.insert(particles.end(), temp.begin(), temp.end());
        }
    }

    return particles;
}