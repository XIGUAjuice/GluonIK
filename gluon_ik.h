#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <limits>

using namespace std;
using Eigen::MatrixXd;
using Eigen::Matrix4d;
using Eigen::Matrix3d;
using Eigen::Matrix;
using Eigen::Vector3d;
using Eigen::RowVector4d;

namespace GluonLib{
    class Gluon
    {
    public:
        static double rad(double degree);
        static double deg(double radian);
        static Matrix4d poseTotrans(double euler_x, double euler_y, double euler_z, double x, double y, double z);
        static vector<vector<double>> IkSolver(Matrix4d T70);
        static vector<double> findNearestSolution(vector<vector<double>> solutions, vector<double> angles_now);
    };
}