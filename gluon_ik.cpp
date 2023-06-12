#include "gluon_ik.h"

namespace GluonLib{
        const double PI = 3.14159265358979323846;
        double Gluon::rad(double degree)
        {
            return degree * PI / 180;
        }
        double Gluon::deg(double radian)
        {
            return radian * 180 / PI;
        }
        Matrix4d Gluon::poseTotrans(double euler_x, double euler_y, double euler_z, double x, double y, double z)
        {
            /*
            将Unity下的位姿表示为齐次变换矩阵
            Unity的角度表示为ZXY固定角
            Unity中采用左手坐标系，因此坐标需取反
            */
            euler_x = rad(euler_x);
            euler_y = rad(euler_y);
            euler_z = rad(euler_z);
            x = -x;
            y = -y;
            z = -z;

            auto Rx = Matrix3d{
                {1, 0, 0},
                {0, cos(euler_x), -sin(euler_x)},
                {0, sin(euler_x), cos(euler_x)}
            };
            auto Ry = Matrix3d{
                {cos(euler_y), 0, sin(euler_y)},
                {0, 1, 0},
                {-sin(euler_y), 0, cos(euler_y)}
            };
            auto Rz = Matrix3d{
                {cos(euler_z), -sin(euler_z), 0},
                {sin(euler_z), cos(euler_z), 0},
                {0, 0, 1}
            };
            auto R = Ry * Rx * Rz;

            auto pos = Vector3d{x, y, z};
            Matrix<double, 3, 4> temp_T;
            temp_T << R, pos;

            auto row_last = RowVector4d{0, 0, 0, 1};
            Matrix4d T;
            T << temp_T, 
                row_last;

            return T;
        }
        vector<vector<double>> Gluon::IkSolver(Matrix4d T70)
        {   
            double d2 = 116.5;
            double a2 = 203.5;
            double a3 = 173;
            double d5 = 79.2;
            double d7 = 42;

            auto T76_i = Matrix4d{
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, -d7},
                {0, 0, 0, 1}
            };
            auto T60 = T70 * T76_i;

            double r11 = T60(0, 0);
            double r12 = T60(0, 1);
            double r13 = T60(0, 2);
            double px = T60(0, 3);
            double r21 = T60(1, 0);
            double r22 = T60(1, 1);
            double r23 = T60(1, 2);
            double py = T60(1, 3);
            double r31 = T60(2, 0);
            double r32 = T60(2, 1);
            double r33 = T60(2, 2);
            double pz = T60(2, 3);

            vector<vector<double>> solutions;

            // 第一个解(theta1+ theta5+ theta3+)
            auto theta1 = atan2(py, px) - atan2(d2, sqrt(pow(px, 2) + pow(py, 2) - pow(d2, 2)));
            auto theta5 = acos(-r13 * sin(theta1) + r23 * cos(theta1));
            auto theta6 = atan2((-r12 * sin(theta1) + r22 * cos(theta1)) / -sin(theta5), (-r11 * sin(theta1) + r21 * cos(theta1)) / sin(theta5));
            
            auto m = d5 * ((r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6)) + px * cos(theta1) + py * sin(theta1);
            auto n = d5 * (r31 * sin(theta6) + r32 * cos(theta6)) + pz;
            auto p = (r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6);
            auto q = r31 * sin(theta6) + r32 * cos(theta6);

            auto theta3 = acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            auto k1 = a2 + a3 * cos(theta3);
            auto k2 = -a3 * sin(theta3);

            auto theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            auto theta4 = atan2(-p, -q) - theta2 - theta3;

            vector<double> solution{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            // 第二个解(theta1+ theta5+ theta3-)
            theta3 = -acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            k1 = a2 + a3 * cos(theta3);
            k2 = -a3 * sin(theta3);

            theta2 = theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            theta4 = theta4 = atan2(-p, -q) - theta2 - theta3;

            solution = vector<double>{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            // 第三个解(theta1+ theta5- theta3+)
            theta5 = -acos(-r13 * sin(theta1) + r23 * cos(theta1));
            theta6 = theta6 = atan2((-r12 * sin(theta1) + r22 * cos(theta1)) / -sin(theta5), (-r11 * sin(theta1) + r21 * cos(theta1)) / sin(theta5));
            
            m = d5 * ((r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6)) + px * cos(theta1) + py * sin(theta1);
            n = d5 * (r31 * sin(theta6) + r32 * cos(theta6)) + pz;
            p = (r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6);
            q = r31 * sin(theta6) + r32 * cos(theta6);

            theta3 = acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            k1 = a2 + a3 * cos(theta3);
            k2 = -a3 * sin(theta3);

            theta2 = theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            theta4 = theta4 = atan2(-p, -q) - theta2 - theta3;

            solution = vector<double>{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            // 第四个解(theta1+ theta5- theta3-)
            theta3 = -acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            k1 = a2 + a3 * cos(theta3);
            k2 = -a3 * sin(theta3);

            theta2 = theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            theta4 = theta4 = atan2(-p, -q) - theta2 - theta3;

            solution = vector<double>{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            // 第五个解(theta1- theta5+ theta3+)
            theta1 = theta1 = atan2(py, px) - atan2(d2, -sqrt(pow(px, 2) + pow(py, 2) - pow(d2, 2)));
            theta5 = acos(-r13 * sin(theta1) + r23 * cos(theta1));
            theta6 = atan2((-r12 * sin(theta1) + r22 * cos(theta1)) / -sin(theta5), (-r11 * sin(theta1) + r21 * cos(theta1)) / sin(theta5));
            
            m = d5 * ((r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6)) + px * cos(theta1) + py * sin(theta1);
            n = d5 * (r31 * sin(theta6) + r32 * cos(theta6)) + pz;
            p = (r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6);
            q = r31 * sin(theta6) + r32 * cos(theta6);

            theta3 = acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            k1 = a2 + a3 * cos(theta3);
            k2 = -a3 * sin(theta3);

            theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            theta4 = atan2(-p, -q) - theta2 - theta3;

            solution = vector<double>{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            // 第六个解(theta1- theta5+ theta3-)
            theta3 = -acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            k1 = a2 + a3 * cos(theta3);
            k2 = -a3 * sin(theta3);

            theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            theta4 = atan2(-p, -q) - theta2 - theta3;

            solution = vector<double>{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            //第七个解(theta1- theta5- theta3+)
            theta5 = -acos(-r13 * sin(theta1) + r23 * cos(theta1));
            theta6 = atan2((-r12 * sin(theta1) + r22 * cos(theta1)) / -sin(theta5), (-r11 * sin(theta1) + r21 * cos(theta1)) / sin(theta5));
            
            m = d5 * ((r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6)) + px * cos(theta1) + py * sin(theta1);
            n = d5 * (r31 * sin(theta6) + r32 * cos(theta6)) + pz;
            p = (r11 * cos(theta1) + r21 * sin(theta1)) * sin(theta6) + (r12 * cos(theta1) + r22 * sin(theta1)) * cos(theta6);
            q = r31 * sin(theta6) + r32 * cos(theta6);

            theta3 = acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            k1 = a2 + a3 * cos(theta3);
            k2 = -a3 * sin(theta3);

            theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            theta4 = atan2(-p, -q) - theta2 - theta3;

            solution = vector<double>{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            // 第八个解(theta1- theta5- theta3-)
            theta3 = -acos((pow(m, 2) + pow(n, 2) - pow(a2, 2) - pow(a3, 2)) / (2 * a2 * a3));

            k1 = a2 + a3 * cos(theta3);
            k2 = -a3 * sin(theta3);

            theta2 = theta2 = atan2(m * k2 - n * k1, m * k1 + n * k2);
            theta4 = theta4 = atan2(-p, -q) - theta2 - theta3;

            solution = vector<double>{deg(theta1), deg(theta2), deg(theta3), deg(theta4), deg(theta5), deg(theta6)};
            solutions.push_back(solution);

            return solutions;
        }
        vector<double> Gluon::findNearestSolution(vector<vector<double>> solutions, vector<double> angles_now)
        {
            int n = solutions.size();
            int nearest = -1;
            double min_distance = numeric_limits<double>::infinity();

            for(int i = 0; i < n; i++)
            {
                double distance = 0;
                for(int j = 0; j < 6; j++)
                {   
                    if(isnan(solutions[i][j]))
                    {
                        distance = numeric_limits<double>::infinity();
                        break;
                    }
                    distance = distance + pow((solutions[i][j] - angles_now[j]), 2);
                }
                if(distance < min_distance)
                {
                    min_distance = distance;
                    nearest = i;
                }
            }

            if(nearest != -1)
                return solutions[nearest];
            else
                return vector<double>{};
        }
}
// double rad(double degree);
// double deg(double radian);
// Matrix4d poseTotrans(double euler_x, double euler_y, double euler_z, double x, double y, double z);
// vector<vector<double>> IkSolver(Matrix4d T70);
// vector<double> findNearestSolution(vector<vector<double>> solutions, vector<double> angles_now);

// int main()
// {
//     auto T = poseTotrans(-32.242, -72.808, 137.192, 159.9763, 268.8852, 341.8257);
//     auto solutions = IkSolver(T);
//     for(auto solution : solutions)
//     {
//         for(auto theta : solution)
//         {
//             cout << theta << " ";
//         }
//         cout << endl;
//     }

//     auto nearest = findNearestSolution(solutions, vector<double>{30, -32, -29, 90, 29, 29});
//     cout << "neatest: " << endl;
//     for(auto theta : nearest)
//     {
//         cout << theta << " ";
//     }
//     cout << endl;

//     return 0;
// }