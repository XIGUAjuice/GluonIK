%module gluon_ik
%include "std_vector.i"

%{
#include "Eigen/Dense"
#include "gluon_ik.h"
%}

%include "Eigen/Dense"
%include "gluon_ik.h"


namespace std{
    %template(DoubleVector) vector<double>;
    %template(DoubleVector2d) vector<vector<double>>;
}
