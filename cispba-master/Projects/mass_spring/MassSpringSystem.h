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

        for(int i = 0; i < x.size(); i++)
        {
            f.push_back(TV(0.0, 0.0, 0.0));
        }

        for(int i = 0; i < segments.size(); i++)
        {
            TV springVec = x[segments[i][0]] - x[segments[i][1]];

            TV springVecDir = springVec.normalized();

            T springVecLength = springVec.norm();

            T forceAmount = -1 * youngs_modulus * (springVecLength / rest_length[i] - 1.0);

            f[segments[i][0]] += forceAmount * springVecDir;

            f[segments[i][1]] += -1.0 * forceAmount * springVecDir;
        }


    }

    void evaluateDampingForces(std::vector<TV >& f)
    {
        // TODO: evaluate damping force

        for(int i = 0; i < x.size(); i++)
        {
            f.push_back(TV(0.0, 0.0, 0.0));
        }

        for(int i = 0; i < segments.size(); i++)
        {
            TV springVec = x[segments[i][0]] - x[segments[i][1]];

            TV springVecDir = springVec.normalized();

            T vRelative = (v[segments[i][0]] - v[segments[i][1]]).dot(springVecDir);

            T fScalar = -1 * damping_coeff * vRelative;

            f[segments[i][0]] += fScalar * springVecDir;

            f[segments[i][1]] += -1.0 * fScalar * springVecDir;
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
