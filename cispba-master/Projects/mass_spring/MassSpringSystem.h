// Modified by Nathan Devlin for CIS 563 Project 1

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class MassSpringSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    T youngs_modulus;
    T damping_coeff;
    std::vector<bool> node_is_fixed;
    std::vector<T> rest_length;

    MassSpringSystem()
    {}

    void evaluateSpringForces(std::vector<TV >& f)
    {
        // TODO: evaluate spring force

        // Initialize f to number of points
        for(int i = 0; i < x.size(); i++)
        {
            f.push_back(TV(0.0, 0.0, 0.0));
        }

        for(int i = 0; i < segments.size(); i++)
        {
            int point1 = segments[i][0];
            int point2 = segments[i][1];

            TV springVec = x[point1] - x[point2];

            TV springVecDir = springVec.normalized();

            T springVecLength = springVec.norm();

            T forceAmount = -1 * youngs_modulus * (springVecLength / rest_length[i] - 1.0);

            f[point1] += forceAmount * springVecDir;

            f[point2] += -1.0 * forceAmount * springVecDir;
        }

    }

    void evaluateDampingForces(std::vector<TV >& f)
    {
        // TODO: evaluate damping force

        // Initialize f to number of points
        for(int i = 0; i < x.size(); i++)
        {
            f.push_back(TV(0.0, 0.0, 0.0));
        }

        for(int i = 0; i < segments.size(); i++)
        {
            int point1 = segments[i][0];
            int point2 = segments[i][1];

            TV pt1Vel = v[point1];
            TV pt2Vel = v[point2];

            // Dashpot damping
            TV springVecDir = (x[point1] - x[point2]).normalized();

            T vRelative = (pt1Vel - pt2Vel).dot(springVecDir);

            T fScalar = -1 * damping_coeff * vRelative;

            //Dashpot Damping
            f[point1] = fScalar * springVecDir;
            f[point2] = -1.0 * fScalar * springVecDir;

            //Small Drag Damping for Air resistance
            T dragDampingCoef = 0.0001;
            f[point1] += - dragDampingCoef * pt1Vel;
            f[point2] += - dragDampingCoef * pt2Vel;
        }

    }


    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        for (const Eigen::Matrix<int, 2, 1>& seg : segments)
            fs << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
        fs << "END\n";
        fs.close();
    }
};
