// Modified by Nathan Devlin 10/22/2020 for CIS 563 Project 1

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>
#include <unordered_set>

#include "SimulationDriver.h"

int main(int argc, char* argv[])
{
    std::cout << "Starting" << std::endl;

    using T = double;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;

    SimulationDriver<T,dim> driver;

    // set up mass spring system

    T youngs_modulus = 1.0;
    T damping_coeff = 1.0;
    T dt = 0.00001;
    T totalMass = 1.0;

    // node data
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<bool> node_is_fixed;

    // segment data
    std::vector<Eigen::Matrix<int,2,1>> segments;
    std::vector<T> rest_length;

    if (argc < 2)
    {
        std::cout << "Please indicate test case number: 0 (cloth) or 1 (volumetric bunny)" << std::endl;
        exit(0);
    }

    if (strcmp(argv[1], "0") == 0) // cloth case
    {
        // TODO
        /*
            1. Create node data: position, mass, velocity
            2. Fill segments and rest_length, including struct springs, shearing springs and bending springs.
            3. Choose proper youngs_modulus, damping_coeff and dt.
            4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
            5. Generate quad mesh for rendering.
        */

        // Set physical values to change visual properties of cloth
        // A large disparity in values usually necessitates a smaller dt time interval
        youngs_modulus = 4.0;
        damping_coeff = 0.0001;
        dt = 0.0001;
        totalMass = 1.0;

        int clothWidthIndices = 65;
        int clothHeightIndices = 65;

        T clothWidth = 2.0;
        T clothHeight = 2.0;

        T totalPoints = (T)(clothWidthIndices * clothHeightIndices);

        std::vector<Eigen::Matrix<int,4,1> > faces;

        TV position = TV(0.0, 0.0, 2.0);

        // Create node data
        for(int r = 0; r < clothHeightIndices; r++)
        {
            for(int c = 0; c < clothWidthIndices; c++)
            {
                x.emplace_back(position);

                node_is_fixed.emplace_back(false);

                m.emplace_back(totalMass / totalPoints);

                v.emplace_back(0.0, 0.0, 0.0);

                position[0] += clothWidth / ((T)clothWidthIndices - 1.0);
            }
            position[0] = 0.0;
            position[2] -= clothHeight / ((T)clothHeightIndices - 1.0);
        }

        int r = 0;
        int c = 0;

        // Create segment data
        for(r = 0; r < clothHeightIndices - 1; r++)
        {
            for(c = 0; c < clothWidthIndices - 1; c++)
            {
                //Add face for obj
                faces.emplace_back(c + r * clothHeightIndices,
                                   c + r * clothHeightIndices + 1,
                                   c + (r + 1) * clothHeightIndices + 1,
                                   c + (r + 1) * clothHeightIndices);

                //Horizontal Structural Springs
                segments.emplace_back(c + r * clothHeightIndices,
                                      r * clothHeightIndices + c + 1);

                //Vertical Structural Springs
                segments.emplace_back(c + r * clothHeightIndices,
                                      (r + 1) * clothHeightIndices + c);

                // Right Shear Springs
                segments.emplace_back(c + r * clothHeightIndices,
                                      (r + 1) * clothHeightIndices + c + 1);

                if(c != 0)
                {
                    // Left Shear Springs
                    segments.emplace_back(c + r * clothHeightIndices,
                                          (r + 1) * clothHeightIndices + c - 1);
                }
                if(c < clothWidthIndices - 2)
                {
                    // Horizontal Flexion Springs
                    segments.emplace_back(c + r * clothHeightIndices,
                                          r * clothHeightIndices + c + 2);
                }
                if(r < clothHeightIndices - 2)
                {
                    // Vertical Flexion Springs
                    segments.emplace_back(c + r * clothHeightIndices,
                                          (r + 2) * clothHeightIndices + c);
                }
            }
            // Vertical Structural Springs
            segments.emplace_back(c + r * clothHeightIndices,
                                  (r + 1) * clothHeightIndices + c);
            // Left Shear Springs
            segments.emplace_back(c + r * clothHeightIndices,
                                  (r + 1) * clothHeightIndices + c - 1);
            if(r < clothHeightIndices - 2)
            {
                // Vertical Flexion Springs
                segments.emplace_back(c + r * clothHeightIndices,
                                      (r + 2) * clothHeightIndices + c);
            }
        }

        // Last Row
        for(c = 0; c < clothWidthIndices - 2; c++)
        {
            // Horizontal Structural Springs
            segments.emplace_back(c + clothHeightIndices * (clothHeightIndices - 1),
                                  clothHeightIndices * (clothHeightIndices - 1) + c + 1);
            // Horizontal Flexion Springs
            segments.emplace_back(c + r * clothHeightIndices,
                                  r * clothHeightIndices + c + 2);
        }
        // Horizontal Structural Springs
        segments.emplace_back(c + clothHeightIndices * (clothHeightIndices - 1),
                              clothHeightIndices * (clothHeightIndices - 1) + c + 1);

        // Populate rest_length
        for(Eigen::Matrix<int,2,1> seg : segments)
        {
            TV vec = x[seg[0]] - x[seg[1]];
            rest_length.push_back(vec.norm());
        }

        // Add fixed handles
        int fixedPt1 = 0;
        int fixedPt2 = clothWidthIndices - 1;

        node_is_fixed[fixedPt1] = true;
        node_is_fixed[fixedPt2] = true;


        // Write out obj file
        std::string filename = "data/writtenCloth.obj";
        std::ofstream objBuffer;
        objBuffer.open(filename);
        if(!objBuffer.is_open())
        {
            std::cout << filename << " failed to open." << std::endl;
        }
        else
        {
            std::cout << filename << " successfully opened." << std::endl;

            objBuffer << "# Object File writtenCloth.obj\n";
            for (TV X : x)
            {
                objBuffer << "v ";
                for (int i = 0; i < dim; i++)
                {
                    objBuffer << " " << X(i);
                }
                if (dim == 2)
                {
                    objBuffer << " 0";
                }
                objBuffer << "\n";
            }
            objBuffer << "\n";
            for (const Eigen::Matrix<int, 4, 1> &face : faces)
            {
                objBuffer << "f " << (T) (face(0) + 1.0) << " " << (T) (face(1) + 1.0) <<
                          " " << (T) (face(2) + 1.0) << " " << (T) (face(3) + 1.0)
                          << "\n"; // poly segment mesh is 1-indexed
            }
            objBuffer << "\n#End of File\n";
            objBuffer.close();
        }


        driver.helper = [&](T t, T dt)
        {
            // TODO

            TV pt1Vel = TV(0.0, 0.0, 0.0);
            TV pt2Vel = TV(0.0, 0.0, 0.0);

            // Move handles along z until t = 2.3
            if(t < 2.3)
            {
                T displace = -4 * cos(2 * t);
                pt1Vel[2] = displace;
                pt2Vel[2] = displace;
            }
            else
            {
                // Stretch fabric along z afterwards
                T displace = 0.5 * cos(2 * t);

                pt1Vel[0] = -1 * displace;
                pt2Vel[0] = displace;
            }

            // Update positions and velocities

            driver.ms.x[fixedPt1] += pt1Vel * dt;
            driver.ms.x[fixedPt2] += pt2Vel * dt;

            driver.ms.v[fixedPt1] = pt1Vel;
            driver.ms.v[fixedPt2] = pt2Vel;
        };

        driver.test="cloth";
    }


    else if (strcmp(argv[1], "1") == 0) // volumetric bunny case
    {
        // TODO
        /*
            1. Create node data from data/points: The first line indicates the number of points and dimension (which is 3).
            2. Fill segments and rest_length from data/cells: The first line indicates the number of tetrahedra and the number of vertices of each tet (which is 6). Each edge in this tetrahedral mesh will be a segment. Be careful not to create duplicate edges.
            3. Choose proper youngs_modulus, damping_coeff, dt;
            4. Set boundary condition (node_is_fixed) and helper function (to achieve moving boundary condition).
        */

        // Set physical values to change visual properties of bunny
        // A large disparity in values usually necessitates a smaller dt time interval
        youngs_modulus = 10.0;
        damping_coeff =  2.0;
        dt = 0.00005;
        totalMass = 2.0;


        // Open points file

        std::string filename = "data/points";
        std::ifstream pointsFile(filename);

        if(!pointsFile.is_open())
        {
            std::cout << filename << " failed to open." << std::endl;
        }
        else
        {
            std::cout << filename << " successfully opened." << std::endl;

            std::string line = "";

            bool hasRun = false;
            int numPoints = 0;

            std::getline(pointsFile, line);

            do
            {
                std::vector<std::string> splitLine;
                std::string buffer = "";
                std::stringstream sStream(line);

                while (sStream >> buffer)
                {
                    splitLine.push_back(buffer);
                }

                if (!hasRun)
                {
                    numPoints = std::stoi(splitLine[0]);
                    std::cout << "Number of points: " << numPoints << std::endl;

                    if (std::stoi(splitLine[1]) != 3)
                    {
                        std::cout << "Error: Improper points file" << std::endl;
                        break;
                    }

                    hasRun = true;
                    continue;
                }

                T xPos = (T) std::stod(splitLine[0]);
                T yPos = (T) std::stod(splitLine[1]);
                T zPos = (T) std::stod(splitLine[2]);

                x.emplace_back(xPos, yPos, zPos);

                node_is_fixed.emplace_back(false);

                m.emplace_back(totalMass / numPoints);

                v.emplace_back(0.0, 0.0, 0.0);

            } while (std::getline(pointsFile, line));

            pointsFile.close();
        }


        // To facilitate hashing of std::pair objects
        struct pairHash
        {
            std::size_t operator()(const std::pair<int,int> & p) const
            {
                return p.first * 31 + p.second;
            }
        };

        // Hash Set to prevent duplicate springs
        std::unordered_set<std::pair<int, int>, pairHash> springs;


        // Open cells file

        filename = "data/cells";
        std::ifstream cellsFile (filename);

        if(!cellsFile.is_open())
        {
            std::cout << filename << " failed to open." << std::endl;
        }
        else
        {
            std::cout << filename << " successfully opened." << std::endl;

            std::string line = "";
            int tetNumber = 1;
            bool hasRun = false;
            int numTets = 0;

            std::getline(cellsFile, line);

            do
            {
                std::vector<std::string> splitLine;
                std::string buffer = "";
                std::stringstream sStream(line);

                while (sStream >> buffer)
                {
                    splitLine.push_back(buffer);
                }

                if (!hasRun)
                {
                    numTets = std::stoi(splitLine[0]);
                    std::cout << "Number of tets: " << numTets << std::endl;

                    if (std::stoi(splitLine[1]) != 4)
                    {
                        std::cout << "Error: Improper tets file" << std::endl;
                        break;
                    }

                    hasRun = true;
                    continue;
                }

                int pt1 = std::stoi(splitLine[0]);
                int pt2 = std::stoi(splitLine[1]);
                int pt3 = std::stoi(splitLine[2]);
                int pt4 = std::stoi(splitLine[3]);

                // Add Springs
                springs.insert(std::pair<int, int>(pt1, pt2));
                springs.insert(std::pair<int, int>(pt1, pt3));
                springs.insert(std::pair<int, int>(pt1, pt4));
                springs.insert(std::pair<int, int>(pt2, pt3));
                springs.insert(std::pair<int, int>(pt2, pt4));
                springs.insert(std::pair<int, int>(pt3, pt4));

                tetNumber++;

            } while (std::getline(cellsFile, line));

            cellsFile.close();
        }
        
        // Populate rest_length
        for(std::pair<int, int> p : springs)
        {
            Eigen::Matrix<int,2,1> seg = Eigen::Matrix<int,2,1>(p.first, p.second);
            segments.emplace_back(p.first, p.second);
            TV vec = x[seg[0]] - x[seg[1]];
            rest_length.emplace_back(vec.norm());
        }

        // Fix handle points
        int ear1 = 2140;
        int ear2 = 2346;
        int tail = 1036;

        node_is_fixed[ear1] = true;
        node_is_fixed[ear2] = true;
        node_is_fixed[tail] = true;

        driver.helper = [&](T t, T dt)
        {
            // TODO

            // Pull back tail for 1 second then release
            if(t < 1.0)
            {
                T pointXVel = 0.4 * cos(t);

                driver.ms.v[tail][0] = pointXVel;

                driver.ms.x[tail][0] += pointXVel * dt;
            }
            else
            {
                driver.ms.node_is_fixed[tail] = false;
            }

        };
        driver.test="bunny";
    }

    else
    {
        std::cout << "Wrong case number!" << std::endl;
        exit(0);
    }

    // simulate

    driver.dt = dt;
    driver.ms.segments = segments;
    driver.ms.m = m;
    driver.ms.v = v;
    driver.ms.x = x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = node_is_fixed;
    driver.ms.rest_length = rest_length;

    driver.run(120);

    std::cout << "Done" << std::endl;

    return 0;
}

