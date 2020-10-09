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
    T youngs_modulus = 100000.0;
    T damping_coeff = 0.1;
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

        node_is_fixed[101] = true;
        node_is_fixed[198] = true;

        driver.helper = [&](T t, T dt){
            // TODO


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

        driver.helper = [&](T t, T dt) {
            // TODO
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
