#include "datamanager.hpp"

DataManager::DataManager() {}

void DataManager::readFile()
{
    std::ifstream in(filePath, std::ios::binary);
    if (!in)
    {
        std::cout << "Load failed." << std::endl;
        return;
    }

    // get the length of the file
    in.seekg(0, in.end);
    int fileLength = in.tellg(); // length of file
    in.seekg(0, in.beg);

    int numberParticles = 0;
    in.read(reinterpret_cast<char *>(&numberParticles), 4);

    int frameSize = 8 + (20 * numberParticles);

    std::cout << "File Length: " << fileLength / frameSize << std::endl;

    double d = 0;
    int n = 0;
    // read through each number of the file
    // 4 comes from each entry being 4 bytes long
    in.seekg(0, in.beg);
    for (int i = 0; i < fileLength; i += frameSize)
    {
        // read through each frame
        // number of vortices
        in.read(reinterpret_cast<char *>(&n), 4);
        std::cout << "Number Particles: " << n << std::endl;
        // frame index
        in.read(reinterpret_cast<char *>(&n), 4);
        std::cout << "Frame Index: " << n << std::endl;

        // each skyrmion's id, posx, posy
        for (int j = 0; j < numberParticles; ++j)
        {
            // particle ID
            in.read(reinterpret_cast<char *>(&n), 4);
            std::cout << "Particle ID: " << n << std::endl;
            // posx
            in.read(reinterpret_cast<char *>(&d), 8);
            std::cout << "X Position: " << d << std::endl;
            // posy
            in.read(reinterpret_cast<char *>(&d), 8);
            std::cout << "Y Position: " << d << std::endl;
        }
    }
}

void DataManager::convertFrameToImage(Grid *grid, int byteIndex, int frameIndex)
{
    // creates an image with a number of pixels given by the grid COLS and ROWS
    // a white pixel denotes no skyrmion, a black pixel denotes a skyrmion

    // Currently we are writing as a .pbm, but if we want colour we could use a
    // .ppm instead (though the code should be change to allow that).

    // here we append zeroes to the frame labels so that the video encoder does not misorder frames
    std::string imagePath = imageDir + "image";
    if (frameIndex < 10)
        imagePath += "000" + std::to_string(frameIndex);
    else if (frameIndex < 100)
        imagePath += "00" + std::to_string(frameIndex);
    else if (frameIndex < 1000)
        imagePath += "0" + std::to_string(frameIndex);
    else
        imagePath += std::to_string(frameIndex);
    imagePath += ".pbm";

    std::ofstream img(imagePath);
    std::ifstream file(filePath, std::ios::binary);

    if (!img || !file)
        return;

    int width = grid->getRows() * imageResolution;
    int height = grid->getCols() * imageResolution;

    img << "P1" << std::endl;
    img << width << " " << height << std::endl;
    // img << "255" << std::endl;

    file.seekg(byteIndex, file.beg);

    int numberParticles = 0;
    file.read(reinterpret_cast<char *>(&numberParticles), 4);
    file.seekg(4, file.cur); // skip the timeframe int

    // currently a boolean array, noting whether a particle exists at that position
    int imageMatrix[width * height];
    memset(imageMatrix, 0, sizeof imageMatrix);

    double x = 0;
    double y = 0;
    float squareParticleRadius = particleRadius * particleRadius;

    for (int i = 0; i < numberParticles; ++i)
    {
        // fill the boolean matrix which holds particle positions
        file.seekg(4, file.cur); // skip the particle ID

        // x position
        file.read(reinterpret_cast<char *>(&x), 8);
        // y position
        file.read(reinterpret_cast<char *>(&y), 8);

        int x_pos = (int)rint(x * imageResolution);
        int y_pos = (int)rint(y * imageResolution);

        for (int a = -particleRadius; a <= particleRadius; ++a)
        {
            // ensure we don't access points outside of boundaries
            if (x_pos + a < 0 || x_pos + a >= width)
                continue;
            for (int b = -particleRadius; b <= particleRadius; ++b)
            {
                // ensure we don't access points outside of boundaries
                if (y_pos + b < 0 || y_pos + b >= height)
                    continue;
                // ensure we draw a circle.
                if (a * a + b * b <= squareParticleRadius)
                {
                    /*
                    For every pixel in the imageMatrix, note true if within particle radius
                    of a particle. Each true point in the matrix will be drawn as a black pixel
                    */
                    imageMatrix[x_pos + a + (width * (y_pos + b))] = 1;
                }
            }
        }
    }

    for (int j = height - 1; j >= 0; --j)
    {
        for (int i = 0; i < width; ++i)
        {
            img << imageMatrix[i + (width * j)];
        }
    }

    file.close();
    img.close();
}

void DataManager::fourierToImage(Grid *grid, int byteIndex, int frameIndex)
{
    std::ifstream file(filePath, std::ios::binary);

    if (!file)
        return;

    int width = grid->getRows() * imageResolution;
    int height = grid->getCols() * imageResolution;

    file.seekg(byteIndex, file.beg);

    int numberParticles = 0;
    file.read(reinterpret_cast<char *>(&numberParticles), 4);
    file.seekg(4, file.cur); // skip the timeframe int

    fftw_complex *ftMatrix;
    ftMatrix = (fftw_complex *)fftw_malloc(height * width * sizeof(fftw_complex));
    // Initialize as all zeroes.
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            ftMatrix[i + (width * j)][0] = 0;
            ftMatrix[i + (width * j)][1] = 0;
        }
    }

    double x = 0;
    double y = 0;

    for (int i = 0; i < numberParticles; ++i)
    {
        // fill the boolean matrix which holds particle positions
        file.seekg(4, file.cur); // skip the particle ID

        // x position
        file.read(reinterpret_cast<char *>(&x), 8);
        // y position
        file.read(reinterpret_cast<char *>(&y), 8);

        int x_pos = (int)rint(x * imageResolution);
        int y_pos = (int)rint(y * imageResolution);

        if (x_pos >= width || x_pos < 0)
            continue;
        if (y_pos >= height || y_pos < 0)
            continue;

        ftMatrix[x_pos + (width * y_pos)][0] = pow(-1.0, x_pos + y_pos);
    }

    //     double realOut = 0;
    //    double imgOut = 0;
    //    double ampOut[width][height] = {0};
    //    for (int yWave = 0; yWave < height; yWave++) {
    //       for (int xWave = 0; xWave < width; xWave++) {
    //          for (int ySpace = 0; ySpace < height; ySpace++) {
    //             for (int xSpace = 0; xSpace < width; xSpace++) {
    //                realOut = (ftMatrix[ySpace][xSpace] * cos(2 *
    //                   M_PI * ((1.0 * xWave * xSpace / width) + (1.0 * yWave * ySpace /
    //                   height)))) / sqrt(width * height);
    //                imgOut = (ftMatrix[ySpace][xSpace] * sin(2 * M_PI
    //                   * ((1.0 * xWave * xSpace / width) + (1.0 * yWave * ySpace / height)))) /
    //                   sqrt( width * height);
    //                ampOut[yWave][xWave] = sqrt(
    //                   (realOut * realOut) +
    //                   (imgOut * imgOut));
    //             }
    //          }
    //       }
    //       std::cout << "Fourier Line " << yWave << "finished." <<std::endl;
    //    }

    fftw_complex *ftOut;
    ftOut = (fftw_complex *)fftw_malloc(height * width * sizeof(fftw_complex));
    fftw_plan plan;
    plan = fftw_plan_dft_2d(width, height, ftMatrix, ftOut, FFTW_FORWARD, FFTW_ESTIMATE);
    // std::cout << "FT Planned." << std::endl;
    fftw_execute_dft(plan, ftMatrix, ftOut);
    // std::cout << "FT Done." << std::endl;
    fftw_destroy_plan(plan);

    fftw_free(ftMatrix);

    // rfftwnd_create_plan rfftwnd_create_plan(int rank, const int *n,
    //                        fftw_direction dir, int flags);
    std::string imagePath = "../bin/fourierimages/image";
    if (frameIndex < 10)
        imagePath += "000" + std::to_string(frameIndex);
    else if (frameIndex < 100)
        imagePath += "00" + std::to_string(frameIndex);
    else if (frameIndex < 1000)
        imagePath += "0" + std::to_string(frameIndex);
    else
        imagePath += std::to_string(frameIndex);
    imagePath += ".pgm";

    std::ofstream fourier(imagePath);

    if (fourier)
    {
        fourier << "P2" << std::endl;
        fourier << width << " " << height << std::endl;
        fourier << 255 << std::endl; // 255 shades between white and black.
        // double imagePixels[width * height] = {0};
        double maxPixelBrightness = 0;

        for (int j = height - 1; j >= 0; --j)
        {
            for (int i = 0; i < width; ++i)
            {
                // fourier << ampOut[i][j];
                ftOut[i + (height * j)][0] = ((ftOut[i + (height * j)][0] * ftOut[i + (height * j)][0]) +
                                              (ftOut[i + (height * j)][1] * ftOut[i + (height * j)][1])) /
                                             sqrt(height * width);

                // // Push it all into the unit square at the centre of the image.
                // int x_pos = (width / 2) - (10 * (int)remainder((i - (width/2)), (2 * M_PI * imageResolution)));
                // //i - (((i - (int)(width / 2.0))/imageResolution)*imageResolution);
                // if(x_pos < 0 || x_pos >= width) x_pos = 0;
                // int y_pos = (height / 2) - (10 * (int)remainder((j - (height/2)), (2 * M_PI * imageResolution)));
                // if(y_pos < 0 || y_pos >= height) y_pos = 0;
                // ftOut[x_pos + (height * y_pos)][3] += ftOut[i + (height * j)][0];

                ftOut[i + (height * j)][0] = log(1 + ftOut[i + (height * j)][0]);

                // value = log(1 + value);
                // imagePixels[i + (height * j)] = value;
                if (ftOut[i + (height * j)][0] > maxPixelBrightness)
                    maxPixelBrightness = ftOut[i + (height * j)][0];
                // fourier << (int)(value * 255) << " ";
            }
        }

        // for (int j = height - 1; j >= 0; --j)
        // {
        //     for (int i = 0; i < width; ++i)
        //     {
        //         ftOut[i + (height * j)][1] += ftOut[i + (height * j)][0];
        //         //ftOut[i + (height * j)][0] = log(1 + (ftOut[i + (height * j)][0]));

        //         //value = log(1 + value);
        //         //imagePixels[i + (height * j)] = value;
        //         //if(ftOut[i + (height * j)][0] > maxPixelBrightness) maxPixelBrightness = ftOut[i + (height * j)][0];
        //         //fourier << (int)(value * 255) << " ";
        //     }
        // }

        for (int j = height - 1; j >= 0; --j)
        {
            for (int i = 0; i < width; ++i)
            {
                // Low cutoff, for better visualization (removes some noise)
                if (ftOut[i + (height * j)][0] / maxPixelBrightness <= 0.01)
                    ftOut[i + (height * j)][0] = 0;
                // fourier << 0 << " ";
                fourier << (int)(255 * ftOut[i + (height * j)][0] / maxPixelBrightness) << " ";
            }
        }
    }
    fourier.close();
    fftw_free(ftOut);
    // std::cout << "Fourier Image Saved. Step " << frameIndex << std::endl;
}

void DataManager::convertToImages(Grid *grid)
{
    std::ifstream file(filePath, std::ios::binary);
    if (!file)
        return;

    // first, delete all frames in the folder, so that there is no overlap or excess frames
    system("rm -r ../bin/images/*.pbm");
    system("rm -r ../bin/fourierimages/*.pgm");

    file.seekg(0, file.beg);
    int numberParticles = 0;
    file.read(reinterpret_cast<char *>(&numberParticles), 4);
    int frameSize = 8 + numberParticles * 20; // See below
    // get the length of the file
    file.seekg(0, file.end);
    int fileLength = file.tellg(); // length of file
    file.close();

    // fileLength/frameSize gives the number of frames in the file
    // note that currently it assumes a constant number of particles in the system
    for (int i = 0; i < fileLength / frameSize; ++i)
    {
        std::cout << "Converting frame " << i << std::endl;
        convertFrameToImage(grid, i * frameSize, i);
        fourierToImage(grid, i * frameSize, i);
    }

    // // Save fourier snapshot of the last image.
    // int j = fileLength / frameSize;
    //         //fourierToImage(grid, 0, 0);
    //         fourierToImage(grid, fileLength - frameSize, j - 1);

    // THIS ISN'T A GOOD WAY TO DO THIS
    // currently just a simple quick way to save time for now
    system("ffmpeg -pattern_type glob -framerate 5 -i \"../bin/images/*.pbm\" ../bin/realspace.avi -y");
    system("ffmpeg -pattern_type glob -framerate 5 -i \"../bin/fourierimages/*.pgm\" ../bin/reciprocalspace.avi -y");
}

void DataManager::writeFile()
{
}

void DataManager::saveParticle(std::ofstream *stream, Particle *particle)
{
    // Int so 4 bytes long.
    stream->write(reinterpret_cast<char *>(&(particle->id)), 4);
    // Doubles so 8 bytes long.
    stream->write(reinterpret_cast<char *>(&(particle->pos_x)), 8);
    stream->write(reinterpret_cast<char *>(&(particle->pos_y)), 8);
}

void DataManager::saveFrame(int frameIndex, std::vector<Particle *> particles)
{
    // Each skyrmion is 20 bytes. Save method:
    // (int)#vortices + (int)frameTime + ((int)id + (double)x + (double)y)) *#vortices
    // Each frame is therefore, 8 + 20 * #vortices bytes long.

    // If this is a new run then overwrite the file (binary), otherwise append new data (binary | app).
    std::ofstream out(filePath, (frameIndex == 0) ? std::ios::binary : std::ios::binary | std::ios::app);

    if (!out)
    {
        std::cout << "Save failed." << std::endl;
        return;
    }

    int numParticles = particles.size();
    out.write(reinterpret_cast<char *>(&numParticles), 4); // save number of vortices
    out.write(reinterpret_cast<char *>(&frameIndex), 4);   // save frame number

    for (Particle *i : particles)
    {
        // For every particle in the simulation, append its id, pos_x, and pos_y data.
        DataManager::saveParticle(&out, i);
    }

    out.close();
}