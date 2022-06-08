#pragma once

#include "grid.hpp"
#include "particle.hpp"
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <math.h>
#include <fftw3.h>

class DataManager
{
private:
    int getFileSize(std::string filePath);

    // put a particle's position into a consistent format to save frame.
    void saveParticle(std::ofstream *stream, Particle *particle);
    int imageResolution = 64; // the side length of each cell in pixels
    int particleRadius = 3;   // the graphical radius (number of pixels) of each particle in the simulation

    std::string filePath = "../bin/data.bin";
    std::string imageDir = "../bin/images/";

public:
    DataManager();

    int saveCounter = 200; // the number of timesteps until a frame is saved

    void fourierToImage(Grid *grid, int byteIndex, int frameIndex);
    // Read a file created by WriteFile(), to give skyrmion positions at each time step.
    void readFile();
    // for converting a single frame to a .ppm file
    void convertFrameToImage(Grid *grid, int byteIndex, int frameIndex);
    // save all data from current simulation as images for each frame.
    void convertToImages(Grid *grid);
    // Save skyrmion positions and save to a binary file.
    void writeFile();
    // Save the positions of all skyrmions at the current timestep.
    void saveFrame(int frameIndex, std::vector<Particle *> particles);
};