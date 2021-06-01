//
//  minimizer.hpp
//
//
//  Created by Dipyaman Pramanik on 28/03/21.
//

#ifndef minimizer_hpp
#define minimizer_hpp

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <functional>


typedef std::vector<double> vec;

#define SHFT(a,b,c,d) (a) = (b); (b) = (c) ; (c) = (d);

double SIGN(double a, double b);
double FMAX(double x,double y);

class line_minimizer
{
    private:
        double GOLD,GLIMIT,TINY;
        std::function <double(double)> _func;

    public:
        line_minimizer();
        int bracket(double,double,std::function <double(double)>);
        int minimize(double,double, std::function <double(double)>);
        int brent(std::function <double(double)>);

        int verbosity;
        double TOL;
        double x,a,b,c,fa,fb,fc;
        double u,v,w,m,fu,fv,fw,fx,fm;
};

template<class T>
class multi_minimizer
{
    private:
        int dimension;
        double one_d(double);
        int iterate(int);
        double metric(double*,double*);
        int add_vec(double*,double*,double*);
        int subtract_vec(double*,double*,double*);
        double dot_product(double*,double*);
        int scalar_product(double,double*,double*);
        int copyvec(double*,double*);
        int copymatrix(double**,double**);
        int normalize(double*);

        double distance;
        double *ini_vec,*direc;
        double **Identity;

        int bracket(double,double);
        int minimize(double,double);
        int brent();

        int maximum_index(double*);

        double GOLD,GLIMIT,TINY;
        T _func;
        int verbosity;
        double x,a,b,c,fa,fb,fc;
        double u,v,w,m,fu,fv,fw,fx,fm;

    public:
        int init_minimizer(int,vec,T);
        int linmin(double*,double*,double*);
        int powell();
        int set_tolerance(double);
        int free_minimizer();
        int verbosity_m;
       
        double *pos;
        double TOL;
        double** P;
        double** n;
        double* change;

        double curr_val;
};








#endif /* minimizer_hpp  */
