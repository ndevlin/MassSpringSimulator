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
    T youngs_modulus = 1000;
    T damping_coeff = 0.2;
    // Working well at 0.0001;
    T dt = 0.0001;

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

        int pointNumber = 1;

        int clothWidthIndices = 100;
        int clothHeightIndices = 100;

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
                m.push_back(0.1);
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
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, r * clothHeightIndices + c + 1));

                //Vertical Structural Springs
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, (r + 1) * clothHeightIndices + c));

                // Right Shear Springs
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, (r + 1) * clothHeightIndices + c + 1));

                if(c != 0)
                {
                    // Left Shear Springs
                    segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, (r + 1) * clothHeightIndices + c - 1));
                }
                if(c < clothWidthIndices - 2)
                {
                    // Horizontal Flexion Springs
                    segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, r * clothHeightIndices + c + 2));
                }
                if(r < clothHeightIndices - 2)
                {
                    // Vertical Flexion Springs
                    segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, (r + 2) * clothHeightIndices + c));
                }
            }
            // Vertical Structural Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, (r + 1) * clothHeightIndices + c));
            // Left Shear Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, (r + 1) * clothHeightIndices + c - 1));
            if(r < clothHeightIndices - 2)
            {
                // Vertical Flexion Springs
                segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, (r + 2) * clothHeightIndices + c));
            }
        }

        for(c = 0; c < clothWidthIndices - 2; c++)
        {
            // Horizontal Structural Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + clothHeightIndices * (clothHeightIndices - 1),
                                                      clothHeightIndices * (clothHeightIndices - 1) + c + 1));
            // Horizontal Flexion Springs
            segments.push_back(Eigen::Matrix<int,2,1>(c + r * clothHeightIndices, r * clothHeightIndices + c + 2));
        }
        // Horizontal Structural Springs
        segments.push_back(Eigen::Matrix<int,2,1>(c + clothHeightIndices * (clothHeightIndices - 1),
                                                  clothHeightIndices * (clothHeightIndices - 1) + c + 1));


        for(Eigen::Matrix<int,2,1> seg : segments)
        {
            TV vec = x[seg[0]] - x[seg[1]];
            rest_length.push_back(vec.norm());
        }

        int fixedPt1 = 101;
        int fixedPt2 = 198;

        node_is_fixed[fixedPt1] = true;
        node_is_fixed[fixedPt2] = true;

        //std::cout<< "x[0]: " << x[0][0] << ", " << x[0][1] << ", " <<  x[0][2] << std::endl;
        //std::cout<< "x[9]: " << x[9][0] << ", " << x[9][1] << ", " <<  x[9][2] << std::endl;



        driver.helper = [&](T t, T dt)
        {
            // TODO

            //driver.ms.x[fixedPt1][2] += t * dt * -2.0;
            //driver.ms.x[fixedPt2][2] += t * dt * -2.0;


            if(t < 2.3)
            {
                // Should be driver.ms.v
                v[fixedPt1][2] = 4 * cos(2 * t);
                v[fixedPt2][2] = 4 * cos(2 * t);
            }
            else
            {
                // Should be driver.ms.v
                v[fixedPt1][2] = 0;
                v[fixedPt2][2] = 0;
                v[fixedPt1][0] = -0.8 * cos(2 * t);
                v[fixedPt2][0] = 0.8 * cos(2 * t);
            }



            driver.ms.x[fixedPt1][0] += v[fixedPt1][0] * dt;
            driver.ms.x[fixedPt2][0] += v[fixedPt2][0] * dt;

            driver.ms.x[fixedPt1][2] += v[fixedPt1][2] * dt;
            driver.ms.x[fixedPt2][2] += v[fixedPt2][2] * dt;

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

        std::ifstream pointsFile ("data/points");

        if(pointsFile.is_open())
        {
            std::cout << "File successfully opened." << std::endl;
        }
        else
        {
            std::cout << "File failed to open." << std::endl;
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
            m.push_back(0.1);
            v.push_back(TV(0.0, 0.0, 0.0));

            //std::cout << pointNumber << ": " << position[0] << ", " << position[1] << ", " << position[2] << std::endl;


            pointNumber++;

        } while(std::getline(pointsFile, line));

        pointsFile.close();



        struct pairHash
        {
            std::size_t operator()(const std::pair<int,int> & p) const
            {
                return p.first * 31 + p.second;
            }
        };

        std::unordered_set<std::pair<int, int>, pairHash> springs;


        // Open cells file

        std::ifstream cellsFile ("data/cells");

        if(cellsFile.is_open())
        {
            std::cout << "File successfully opened." << std::endl;
        }
        else
        {
            std::cout << "File failed to open." << std::endl;
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

            std::cout << splitLine[0] << " " << splitLine[1] << " " << splitLine[2] << " " <<splitLine[3] << std::endl;

            int pt1 = std::stoi(splitLine[0]);
            int pt2 = std::stoi(splitLine[1]);
            int pt3 = std::stoi(splitLine[2]);
            int pt4 = std::stoi(splitLine[3]);

            /*
            std::cout << pt1[0] << ", " << pt1[1] << ", " << pt1[2] << std::endl;
            std::cout << pt2[0] << ", " << pt2[1] << ", " << pt2[2] << std::endl;
            std::cout << pt3[0] << ", " << pt3[1] << ", " << pt3[2] << std::endl;
            std::cout << pt4[0] << ", " << pt4[1] << ", " << pt4[2] << std::endl;
            */

            // Add Spring
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


        int ear1 = 2140;
        int ear2 = 2346;
        int tail = 1036;

        node_is_fixed[ear1] = true;
        node_is_fixed[ear2] = true;
        node_is_fixed[tail] = true;



        driver.helper = [&](T t, T dt)
        {
            // TODO

            if(t < 1.1)
            {
                driver.ms.v[tail][0] = 0.5 * cos(1 * t);

                driver.ms.x[tail][0] += driver.ms.v[tail][0] * dt;
            }
            else
            {
                driver.ms.node_is_fixed[tail] = false;
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