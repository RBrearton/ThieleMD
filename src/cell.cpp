#include "cell.hpp"
#include <iostream>
Cell::Cell(int x, int y)
    : coords{std::array<int, 2>{x, y}} {}
Cell::Cell(std::array<int, 2> coords)
    : coords{coords} {}
Cell::Cell(int x, int y, std::vector<Particle *> particles)
    : coords{std::array<int, 2>{x, y}}, particles{particles} {}
Cell::Cell(std::array<int, 2> coords, std::vector<Particle *> particles)
    : coords{coords}, particles{particles} {}

std::array<int, 2> Cell::getCoords() { return coords; }
std::vector<Particle *> Cell::getParticles() { return particles; }

void Cell::addParticle(Particle *particle)
{
    particles.push_back(particle);
}
void Cell::removeParticle(int particleIndex)
{
    particles.erase(particles.begin() + particleIndex);
}
void Cell::removeParticle(Particle *particle)
{
    for (int i = 0; i < particles.size(); ++i)
    {
        if (particles[i]->id == particle->id)
        {
            removeParticle(i);
            return;
        }
    }
}