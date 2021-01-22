//
//  numerical.hpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 25/11/20.
//

#ifndef numerical_hpp
#define numerical_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include <cmath>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>

typedef std::vector<double> vec;


double square(double);

class Cubic_interpolator
{
public:
    typedef boost::math::cubic_b_spline<double> interp;
    interp spline;
    ~Cubic_interpolator();
    int free_cubic();
    int set_cubic_spline(vec,vec,bool);
    double interpolate(double);
    
};


class Bary_interpolator
{
public:
    typedef boost::math::barycentric_rational<double> non_interp;
    non_interp bary;
  //  ~Bary_interpolator();
 //   Bary_interpolator();
    int set_barycentric(vec,vec,bool);
    double interpolate(double);
    
    
};


class two_linear_interp
{
private:
    bool logx=false;
    bool logy=false;
    double **f_grid;
    vec data_x,data_y;
    int x_dim,y_dim;
    int make_grid();
    int convert_logx();
    int convert_logy();
    double interpolate(double,double);
    int free_grid();
    
public:

    vec *data;
    int input_vectors(vec,vec,vec);
    double interp2d(double,double);
    int set_logscale_x();
    int set_logscale_y();

};


#endif /* numerical_hpp */
