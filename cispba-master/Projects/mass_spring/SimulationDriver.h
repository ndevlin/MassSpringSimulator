// Modified by Nathan Devlin for CIS 563 Project 1

#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MassSpringSystem.h"
#include <functional>


template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MassSpringSystem<T,dim> ms;
    T dt;
    TV gravity;
    std::string test;
    T frameRate;

    std::function<void(T, T)> helper = [&](T, T) {};

    SimulationDriver()
    : dt((T)0.00001) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.80665;
        frameRate = 24.0;
    }

    void run(const int max_frame)
    {
        T accumulate_t = 0;
        mkdir("output/", 0777);
        std::string output_folder = "output/" + test;
        mkdir(output_folder.c_str(), 0777);
        std::string filename = output_folder + "/" + std::to_string(0) + ".poly";

        ms.dumpPoly(filename);

        int N_substeps = (int)((1.0 / frameRate) / dt) + 1;
        std::cout << "Num steps: " << N_substeps << std::endl;

        for(int frame = 1; frame <= max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int count = 1;
            int tenPercent = (int)((T)N_substeps / 10.0);

            for (int step = 1; step <= N_substeps; step++) {
                //std::cout << "Step " << step << std::endl;
                helper(accumulate_t, dt);

                if(step > tenPercent)
                {
                    std::cout << count * 2 << "0% " << std::endl;
                    tenPercent += tenPercent;
                    count++;
                }

                advanceOneStepExplicitIntegration();
                accumulate_t += dt;
            }
            std::cout << std::endl;

            mkdir("output/", 0777);
            std::string output_folder = "output/" + test;
            mkdir(output_folder.c_str(), 0777);
            std::string filename = output_folder + "/" + std::to_string(frame) + ".poly";
            ms.dumpPoly(filename);
            std::cout << std::endl;
        }
    }

    void advanceOneStepExplicitIntegration()
    {
        std::vector<TV> f_spring;
        ms.evaluateSpringForces(f_spring);
        std::vector<TV> f_damping;
        ms.evaluateDampingForces(f_damping);

        // TODO: update position and velocity according to Newton's law.

        for(int i = 0; i < ms.x.size(); i++)
        {
            // int version of node_is_fixed to avoid branching
            T isFree = (T)(int)(!ms.node_is_fixed[i]);

            T dtByMass = dt / ms.m[i];

            // Add gravitational force
            TV addedVelocity = gravity * dt;

            addedVelocity += f_spring[i] * dtByMass;

            addedVelocity += f_damping[i] * dtByMass;

            // Update velocity and position if this point is free

            ms.x[i] += isFree * (addedVelocity + ms.v[i]) * dt;

            ms.v[i] += isFree * addedVelocity;
        }

    }
};


