// Modified by Nathan Devlin for CIS 563 Project 1

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

    // Bunny Working decently well at;
    /*
    T youngs_modulus = 1000;
    T damping_coeff = 0.2;
    T dt = 0.000001;
    T totalMass = 550.0;

     Also at
      T youngs_modulus = 10.0;
        T damping_coeff = 2.0;
        T dt = 0.000001;
        T totalMass = 2.0;
     */

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
    std::vector<Eigen::Matrix<int,2,1> > segments;
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

        youngs_modulus = 15000.0;
        damping_coeff = 10.0;
        dt = 0.0001;
        totalMass = 1000.0;


        int pointNumber = 1;

        int clothWidthIndices = 100;
        int clothHeightIndices = 100;

        T totalPoints = (T)(clothWidthIndices * clothHeightIndices);

        T clothWidth = 2.0;
        T clothHeight = 2.0;

        TV position = TV(0.0, 0.0, 2.0);

        for(int r = 0; r < clothHeightIndices; r++)
        {
            for(int c = 0; c < clothWidthIndices; c++) {

                /*
                std::cout << pointNumber << ": " <<
                          position[0] << "," << position[1] << "," << position[2]
                          << std::endl;
                */

                pointNumber++;


                x.push_back(TV(position));
                node_is_fixed.push_back(false);

                //
                m.push_back(totalMass / totalPoints);

                v.push_back(TV(0.0, 0.0, 0.0));

                position[0] += clothWidth / ((T) clothWidthIndices - 1.0);

            }
            position[0] = 0.0;
            position[2] -= clothHeight / ((T)clothHeightIndices - 1.0);
        }

        int r = 0;
        int c = 0;

        for(r = 0; r < clothHeightIndices - 1; r++)
        {
            for(c = 0; c < clothWidthIndices - 1; c++)
            {
                //Horizontal Structural Springs
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                          r * clothHeightIndices + c + 1));

                //Vertical Structural Springs
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                          (r + 1) * clothHeightIndices + c));

                // Right Shear Springs
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                          (r + 1) * clothHeightIndices + c + 1));

                if(c != 0)
                {
                    // Left Shear Springs
                    segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                              (r + 1) * clothHeightIndices + c - 1));
                }
                if(c < clothWidthIndices - 2)
                {
                    // Horizontal Flexion Springs
                    segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                              r * clothHeightIndices + c + 2));
                }
                if(r < clothHeightIndices - 2)
                {
                    // Vertical Flexion Springs
                    segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                              (r + 2) * clothHeightIndices + c));
                }
            }
            // Vertical Structural Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                      (r + 1) * clothHeightIndices + c));
            // Left Shear Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                      (r + 1) * clothHeightIndices + c - 1));
            if(r < clothHeightIndices - 2)
            {
                // Vertical Flexion Springs
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                          (r + 2) * clothHeightIndices + c));
            }
        }

        for(c = 0; c < clothWidthIndices - 2; c++)
        {
            // Horizontal Structural Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + clothHeightIndices * (clothHeightIndices - 1),
                                                      clothHeightIndices * (clothHeightIndices - 1) + c + 1));
            // Horizontal Flexion Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices,
                                                      r * clothHeightIndices + c + 2));
        }
        // Horizontal Structural Springs
        segments.push_back(Eigen::Matrix<int,2,1>(c + clothHeightIndices * (clothHeightIndices - 1),
                                                  clothHeightIndices * (clothHeightIndices - 1) + c + 1));

        for(Eigen::Matrix<int,2,1> seg : segments)
        {
            TV vec = x[seg[0]] - x[seg[1]];
            rest_length.push_back(vec.norm());
        }

        int fixedPt1 = 0;
        int fixedPt2 = clothWidthIndices - 1;

        node_is_fixed[fixedPt1] = true;
        node_is_fixed[fixedPt2] = true;


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

        /*
        youngs_modulus = 1000.0;
        damping_coeff = 2.0;
        dt = 0.00001;
        totalMass = 2.0;
         */

        youngs_modulus = 10.0;
        damping_coeff = 2.0;
        dt = 0.00001;
        totalMass = 2.0;


        // Open obj
        /*
        std::ifstream bunnyObj ("data/bunny.obj");

        if(bunnyObj.is_open())
        {
            std::cout << "File successfully opened." << std::endl;
        }
        else
        {
            std::cout << "File failed to open." << std::endl;
        }

        std::string line;

        int pointNumber = 1;

        while(std::getline(bunnyObj, line))
        {
            if(line == "" || line[0] == ' ' || line[0] == '#')
            {
                continue;
            }

            std::vector<std::string> splitLine;

            std::stringstream tokenize(line);

            std::string intermediate;

            while(std::getline(tokenize, intermediate, ' '))
            {
                splitLine.push_back(intermediate);
            }

            if(splitLine[0][0] == 'v')
            {
                T xPos = (T)std::stod(splitLine[1]);
                T yPos = (T)std::stod(splitLine[2]);
                T zPos = (T)std::stod(splitLine[3]);

                TV position = TV(xPos, yPos, zPos);

                x.push_back(position);
                node_is_fixed.push_back(false);
                m.push_back(0.1);
                v.push_back(TV(0.0, 0.0, 0.0));

                pointNumber++;
            }
        }
        bunnyObj.close();
        */



        // Open points

        std::string filename = "data/points";
        std::ifstream pointsFile(filename);

        if(pointsFile.is_open())
        {
            std::cout << filename << " successfully opened." << std::endl;
        }
        else
        {
            std::cout << filename << " failed to open." << std::endl;
        }

        std::string line = "";

        int pointNumber = 1;
        bool hasRun = false;
        int numPoints = 0;

        std::getline(pointsFile, line);

        do
        {
            std::vector<std::string> splitLine;
            std::string buffer = "";
            std::stringstream sStream(line);

            while(sStream >> buffer)
            {
                splitLine.push_back(buffer);
            }

            if(!hasRun)
            {
                numPoints = std::stoi(splitLine[0]);
                std::cout << "Number of points: " << numPoints << std::endl;

                if(std::stoi(splitLine[1]) != 3)
                {
                    std::cout << "Error: Improper points file" << std::endl;
                    break;
                }

                hasRun = true;
                continue;
            }

            T xPos = (T)std::stod(splitLine[0]);
            T yPos = (T)std::stod(splitLine[1]);
            T zPos = (T)std::stod(splitLine[2]);

            TV position = TV(xPos, yPos, zPos);

            x.push_back(position);

            node_is_fixed.push_back(false);

            m.push_back(totalMass / numPoints);

            v.push_back(TV(0.0, 0.0, 0.0));

            pointNumber++;

        } while(std::getline(pointsFile, line));

        pointsFile.close();


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

        if(cellsFile.is_open())
        {
            std::cout << filename << " successfully opened." << std::endl;
        }
        else
        {
            std::cout << filename << " failed to open." << std::endl;
        }

        line = "";
        int tetNumber = 1;
        hasRun = false;
        int numTets = 0;

        std::getline(cellsFile, line);

        do
        {
            std::vector<std::string> splitLine;
            std::string buffer = "";
            std::stringstream sStream(line);

            while(sStream >> buffer)
            {
                splitLine.push_back(buffer);
            }

            if(!hasRun)
            {
                numTets = std::stoi(splitLine[0]);
                std::cout << "Number of tets: " << numTets << std::endl;

                if(std::stoi(splitLine[1]) != 4)
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


            /*
            segments.push_back(Eigen::Matrix<int,2,1>(pt1, pt2));
            segments.push_back(Eigen::Matrix<int,2,1>(pt1, pt3));
            segments.push_back(Eigen::Matrix<int,2,1>(pt1, pt4));
            segments.push_back(Eigen::Matrix<int,2,1>(pt2, pt3));
            segments.push_back(Eigen::Matrix<int,2,1>(pt2, pt4));
            segments.push_back(Eigen::Matrix<int,2,1>(pt3, pt4));
            */


            tetNumber++;

        } while(std::getline(cellsFile, line));

        cellsFile.close();


        // Populate rest_length
        for(std::pair<int, int> p : springs)
        {
            Eigen::Matrix<int,2,1> seg = Eigen::Matrix<int,2,1>(p.first, p.second);
            segments.push_back(seg);
            TV vec = x[seg[0]] - x[seg[1]];
            rest_length.push_back(vec.norm());
        }

        /*
        for(Eigen::Matrix<int,2,1> seg : segments)
        {
            TV vec = x[seg[0]] - x[seg[1]];
            rest_length.push_back(vec.norm());
        }
        */


        //Ear1 orig point is 2140;
        //Ear2 orig point is 2346;
        //Tail orig point is 1036;

        // Ear 1
        node_is_fixed[2140] = true;
        /*
        node_is_fixed[962] = true;
        node_is_fixed[2130] = true;
        node_is_fixed[1230] = true;
        node_is_fixed[2144] = true;
        node_is_fixed[950] = true;
         */

        // Ear 2
        node_is_fixed[2346] = true;
        /*
        node_is_fixed[1202] = true;
        node_is_fixed[2773] = true;
        node_is_fixed[1995] = true;
        node_is_fixed[1203] = true;
        node_is_fixed[1200] = true;
        */

        // Tail
        node_is_fixed[1036] = true;
        /*
        node_is_fixed[1067] = true;
        node_is_fixed[1068] = true;
        node_is_fixed[1026] = true;
        node_is_fixed[1014] = true;
        node_is_fixed[1000] = true;
        node_is_fixed[1025] = true;
        node_is_fixed[688] = true;
         */


        driver.helper = [&](T t, T dt)
        {
            // TODO

            if(t < 1.1)
            {
                T pointXVel = 0.5 * cos(1 * t);

                driver.ms.v[1036][0] = pointXVel;
                /*
                driver.ms.v[1067][0] = pointXVel;
                driver.ms.v[1068][0] = pointXVel;
                driver.ms.v[1026][0] = pointXVel;
                driver.ms.v[1014][0] = pointXVel;
                driver.ms.v[1000][0] = pointXVel;
                driver.ms.v[1025][0] = pointXVel;
                driver.ms.v[688][0] = pointXVel;
                 */


                driver.ms.x[1036][0] += pointXVel * dt;
                /*
                driver.ms.x[1067][0] += pointXVel * dt;
                driver.ms.x[1068][0] += pointXVel * dt;
                driver.ms.x[1026][0] += pointXVel * dt;
                driver.ms.x[1014][0] += pointXVel * dt;
                driver.ms.x[1000][0] += pointXVel * dt;
                driver.ms.x[1025][0] += pointXVel * dt;
                driver.ms.x[688][0] += pointXVel * dt;
                 */

            }
            else
            {
                driver.ms.node_is_fixed[1036] = false;
                /*
                driver.ms.node_is_fixed[1067] = false;
                driver.ms.node_is_fixed[1068] = false;
                driver.ms.node_is_fixed[1026] = false;
                driver.ms.node_is_fixed[1014] = false;
                driver.ms.node_is_fixed[1000] = false;
                driver.ms.node_is_fixed[1025] = false;
                driver.ms.node_is_fixed[688] = false;
                 */

            }

        };
        driver.test="bunny";
    }

    else {
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

    // Orignal frame #s: 120
    driver.run(120);

    std::cout << "Done" << std::endl;

    return 0;
}