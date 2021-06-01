#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <functional>
#include <vector>

#include "minimizer.hpp"
#include "chisqmin.hpp"

line_minimizer::line_minimizer()
{
    GOLD = (-1+sqrt(5))/2.0;
    GLIMIT = 100;
    TINY = 1.0e-20;
    verbosity = 0;
    TOL = 1e-8;

}


double SIGN(double a, double b)
{
    if(b>0)
    {
        return a;
    }
    else
    {
        return -a;
    }

    return 0;

}

double FMAX(double x,double y)
{
    double m;

    if(x<y)
    {
        m = y;
    }
    else
    {
        m = x;
    }

    return m;
}



int line_minimizer::bracket(double _a,double _b, std::function <double(double)> func)
{
    double ulim,u,r,q,fu,dum;

    a = _a;
    b = _b;

    fa = func(a);
    fb = func(b);

    /* We always want to go from a to b down hill so if f(a)
     * < f(b), we switch a and b;
     */

    if(fb>fa)    
    {
        SHFT(dum,a,b,dum);
        SHFT(dum,fa,fb,dum);

    }

    c = b + GOLD*(b-a);
    fc = func(c);

    while(fb>fc)
    {
        r = (b-a)*(fb-fc);
        q = (b-c)*(fb-fa);
/*
 * u is the extremum of the fitted parabola
 */
        u = b - ((b-c)*q - (b-a)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));

        ulim = b + GLIMIT*(c-b);

        if((b-u)*(u-c)>=0.0)  /* u is between b and c */
        {
            fu = func(u);
            if(fu<fc)    /* minimum lies betwen b and c */
            {
                a = b;
                b = u;
                fa = fb;
                fb = fu;
                return 0;
            }
            else if(fu>fb)  /* minimum lies between a and u */
            {
                c = u;
                fc = fu;
                return 0;

            }

            u = c + GOLD*(c-b);  /* parabolic fit was no use, Use default magnification */
            fu = func(u);

        }
        else if ((c-u)*(u-ulim)>0.0)
        {
            fu = func(u);
            if(fu<fc)
            {
                SHFT(b,c,u,c+GOLD*(c-b));
                SHFT(fb,fc,fu,func(u));
            }
        }
        else if ((u-ulim)*(ulim-c)>= 0.0)
        {
            u = ulim;
            fu = func(u);
        }
        else
        {
            u = c + GOLD*(c-b);
            fu = func(u);
        }

        SHFT(a,b,c,u);
        SHFT(fa,fb,fc,fu);

    }



    

    return 0;
}

int line_minimizer::minimize(double _a,double _b,std::function <double (double)> func)
{
    m = 0;

    if(verbosity!=0)
    {
        std::cout<<"Verbosity is set to True, everything will be printed..."<<"\n";
    }

    bracket(_a,_b,func);

    if(verbosity!=0)
    {
        std::cout<<"The initital bracket is :"<<a<<" "<<c<<std::endl;
    }

    brent(func);


    if(verbosity!=0)
    {
        std::cout<<"The minimum is :"<<m<<std::endl;
    }

    fm = func(m);


    return 0;
}

int line_minimizer::brent(std::function <double(double)> func)
{
    double dum;
    b = c;

    
    if(a>b)
    {
        SHFT(dum,a,b,dum);
    }


    if(verbosity!=0)
    {
        std::cout<<"Initializing brent method..."<<std::endl;
    }
    m = (a+b)/2;

    v = a + (1-GOLD)*(b-a);
    w = v;
    x = v;

    fv = func(v);
    fw = fv;
    fx = fv;

    if(verbosity!=0)
    {
        std::cout<<"Check initial tolerance..."<<std::endl;
        std::cout<<"Tolerance = "<<TOL<<"\n";
    }
    
    if(fabs(b-a)<TOL)
    {

        if(verbosity!=0)
        {
            std::cout<<"found minimum, exiting..."<<std::endl;
            std::cout<<"minimum = "<<m<<"\n";
        }

        return 0;
    }
    
    if(verbosity!=0)
    {
        std::cout<<"Enetering loop..."<<std::endl;
    }

    int k =0;
    int ITMAX = 100;

    double etemp,e,d,p,q;

    e = 0.0;
    
   


    u = x;

    while(fabs(b-a)>=TOL && k<=ITMAX)
    {

        p = (w-x)*(w-x)*(fx-fv)+(v-x)*(v-x)*(fw-fx);
        q = (w-x)*(fx-fv)+(v-x)*(fw-fx);

        etemp = e;
        e = d;
        
        

        if(q!=0 && (a-u)*(u-b)>=0.0 && fabs(p/q)<0.5*fabs(etemp))
        {

 
            d = 0.5*p/q;
            u = x + d;


            if(verbosity!=0)
            {
                std::cout<<"k = "<<k<<" parabolic ... m = "<<m<<std::endl;
            }

        }
        else
        {
            if(x>=m)
            {
                u = GOLD*x + (1-GOLD)*a;
                d = u-x;
                    
            }
            else if(x<m)
            {

                
                u = GOLD*x + (1-GOLD)*b;
                d = u-x;
            }
            
            
            if(verbosity!=0)
            {
                std::cout<<"k = "<<k<<" "<<" golden section ... m = "<<m<<std::endl;
            }



        }

        fu = func(u);

        if(fu<=fx)
        {
            if(u<=x)
            {
                b=x;
            }
            else
            {
                a=x;
            }
            SHFT(v,w,x,u);
            SHFT(fv,fw,fx,fu);
        }
        else
        {
            if(u<=x)
            {
                a=u;
            }
            else
            {
                b=u;
            }
            if(fu<=fw || w==x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if(fu<=fv || v==x||x==w)
            {
                v=u;
                fv=fu;
            }
        }



        k = k+1;
        m = (a+b)/2.0;
  
           if(k>100)
           {
               return 0;
           }

    }


    




    return 0;
}



template<class T>
int multi_minimizer<T>::bracket(double _a,double _b)
{
    double ulim,u,r,q,fu,dum;

    a = _a;
    b = _b;

    fa = one_d(a);
    fb = one_d(b);

    /* We always want to go from a to b down hill so if f(a)
     * < f(b), we switch a and b;
     */

    if(fb>fa)    
    {
        SHFT(dum,a,b,dum);
        SHFT(dum,fa,fb,dum);

    }

    c = b + GOLD*(b-a);
    fc = one_d(c);

    while(fb>fc)
    {
        r = (b-a)*(fb-fc);
        q = (b-c)*(fb-fa);
/*
 * u is the extremum of the fitted parabola
 */
        u = b - ((b-c)*q - (b-a)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));

        ulim = b + GLIMIT*(c-b);

        if((b-u)*(u-c)>=0.0)  /* u is between b and c */
        {
            fu = one_d(u);
            if(fu<fc)    /* minimum lies betwen b and c */
            {
                a = b;
                b = u;
                fa = fb;
                fb = fu;
                return 0;
            }
            else if(fu>fb)  /* minimum lies between a and u */
            {
                c = u;
                fc = fu;
                return 0;

            }

            u = c + GOLD*(c-b);  /* parabolic fit was no use, Use default magnification */
            fu = one_d(u);

        }
        else if ((c-u)*(u-ulim)>0.0)
        {
            fu = one_d(u);
            if(fu<fc)
            {
                SHFT(b,c,u,c+GOLD*(c-b));
                SHFT(fb,fc,fu,one_d(u));
            }
        }
        else if ((u-ulim)*(ulim-c)>= 0.0)
        {
            u = ulim;
            fu = one_d(u);
        }
        else
        {
            u = c + GOLD*(c-b);
            fu = one_d(u);
        }

        SHFT(a,b,c,u);
        SHFT(fa,fb,fc,fu);

    }



    

    return 0;
}

template<class T>
int multi_minimizer<T>::minimize(double _a,double _b)
{
    m = 0;

    if(verbosity!=0)
    {
        std::cout<<"Verbosity is set to True, everything will be printed..."<<"\n";
    }

    bracket(_a,_b);

    if(verbosity!=0)
    {
        std::cout<<"The initital bracket is :"<<a<<" "<<c<<std::endl;
    }

    brent();


    if(verbosity!=0)
    {
        std::cout<<"The minimum is :"<<m<<std::endl;
    }

    fm = one_d(m);


    return 0;
}

template<class T>
int multi_minimizer<T>::brent()
{
    double dum;
    b = c;

    
    if(a>b)
    {
        SHFT(dum,a,b,dum);
    }


    if(verbosity!=0)
    {
        std::cout<<"Initializing brent method..."<<std::endl;
    }
    m = (a+b)/2;

    v = a + (1-GOLD)*(b-a);
    w = v;
    x = v;

    fv = one_d(v);
    fw = fv;
    fx = fv;

    if(verbosity!=0)
    {
        std::cout<<"Check initial tolerance..."<<std::endl;
        std::cout<<"Tolerance = "<<TOL<<"\n";
    }
    
    if(fabs(b-a)<TOL)
    {

        if(verbosity!=0)
        {
            std::cout<<"found minimum, exiting..."<<std::endl;
            std::cout<<"minimum = "<<m<<"\n";
        }

        return 0;
    }
    
    if(verbosity!=0)
    {
        std::cout<<"Enetering loop..."<<std::endl;
    }

    int k =0;
    int ITMAX = 100;

    double etemp,e,d,p,q;

    e = 0.0;
    
   


    u = x;

    while(fabs(b-a)>=TOL && k<=ITMAX)
    {

        p = (w-x)*(w-x)*(fx-fv)+(v-x)*(v-x)*(fw-fx);
        q = (w-x)*(fx-fv)+(v-x)*(fw-fx);

        etemp = e;
        e = d;
        
        

        if(q!=0 && (a-u)*(u-b)>=0.0 && fabs(p/q)<0.5*fabs(etemp))
        {

 
            d = 0.5*p/q;
            u = x + d;


            if(verbosity!=0)
            {
                std::cout<<"k = "<<k<<" parabolic ... m = "<<m<<std::endl;
            }

        }
        else
        {
            if(x>=m)
            {
                u = GOLD*x + (1-GOLD)*a;
                d = u-x;
                    
            }
            else if(x<m)
            {

                
                u = GOLD*x + (1-GOLD)*b;
                d = u-x;
            }
            
            
            if(verbosity!=0)
            {
                std::cout<<"k = "<<k<<" "<<" golden section ... m = "<<m<<std::endl;
            }



        }

        fu = one_d(u);

        if(fu<=fx)
        {
            if(u<=x)
            {
                b=x;
            }
            else
            {
                a=x;
            }
            SHFT(v,w,x,u);
            SHFT(fv,fw,fx,fu);
        }
        else
        {
            if(u<=x)
            {
                a=u;
            }
            else
            {
                b=u;
            }
            if(fu<=fw || w==x)
            {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            }
            else if(fu<=fv || v==x||x==w)
            {
                v=u;
                fv=fu;
            }
        }



        k = k+1;
        m = (a+b)/2.0;
  
           if(k>100)
           {
               return 0;
           }

    }


    




    return 0;
}

template<class T>
int multi_minimizer<T>::init_minimizer(int dim,vec start_value, T obj)
{
    
    GOLD = (-1+sqrt(5))/2.0;
    GLIMIT = 100;
    TINY = 1.0e-20;
    verbosity = 0;
    verbosity_m=0;
    TOL = 1e-8;
    _func = obj;

    dimension = dim;
    

    pos = new double[dimension];

    P = new double*[dimension+1];
    n = new double*[dimension];
    Identity = new double*[dimension];
    change = new double[dimension];

    ini_vec = new double[dimension];
    direc = new double[dimension];

    for(int i=0;i<dimension;i++)
    {
        n[i] = new double[dimension];
        Identity[i] = new double[dimension];
    }

    for(int i=0;i<dimension+1;i++)
    {
        P[i] = new double[dimension];
    }

    for(int i =0;i<dimension+1;i++)
    {
        for(int j=0;j<dimension;j++)
        {

            if(i==0)
            {
                P[i][j] = start_value[j];
            }
            else
            {
                P[i][j] = 0.0;
            }
        }
    }


    for(int i =0;i<dimension;i++)
    {
        for(int j=0;j<dimension;j++)
        {

            if(i==j)
            {
                Identity[i][j] = 1.0;
            }
            else
            {
                Identity[i][j] = 0.0;
            }
        }
    }

    copymatrix(Identity,n);

    return 0;
}

template<class T>
int multi_minimizer<T>::powell()
{
    int k=0;
    int ITMAX = 10000;

    /* take initial guess*/

    /* initialize by minimizing along the basis vectors */

    iterate(k);


    /* Minimizing loop if the tolerance is not met */

    while(distance >=TOL && k<=ITMAX)
    {
        iterate(k);

        if(verbosity_m!=0)
        {
            std::cout<<k<<"\t";
            for(int j=0;j<dimension;j++)
            {
                std::cout<<pos[j]<<"\t";

            } 
            std::cout<<TOL<<"\t"<<curr_val<<"\t"<<distance<<"\t"<<"\n";
        }
//        cout<<k<<"\t"<<pos[0]<<"\t"<<pos[1]<<"\t"<<curr_val<<"\n";
        k=k+1;
    }

    





    return 0;
}

template<class T>
int multi_minimizer<T>::iterate(int k)
{
   double previous; 

    for(int i =0;i<dimension;i++)
    {
        previous = curr_val;

        linmin(P[i],P[i+1],n[i]);

        change[i] = previous-curr_val;
    }

    for(int i =0;i<dimension-1;i++)
    {
        copyvec(n[i+1],n[i]);

    }
    subtract_vec(P[dimension],P[0],n[dimension-1]);

    normalize(n[dimension-1]);


    distance  = metric(P[0],P[dimension]);


   // cout<<distance<<"\n";

    int ind = maximum_index(change);


    for(int i=0;i<dimension;i++)
    {
   //     cout<<n[dimension-1][i]<<"\t"<<n[ind][i]<<"\n";
    }

  //  cout<<ind<<endl;
//    copyvec(n[dimension-1],n[ind]);






    if(k%(2*dimension+1)==0)
    {
        /* after N iterations, reset directions along the
         * basis vectors. */

        copymatrix(Identity,n);
   
    }

    linmin(P[dimension],P[0],n[dimension-1]);


    for(int i=0;i<dimension;i++)
    {
        pos[i] = P[0][i];
    }



    return 0;
}

template<class T>
int multi_minimizer<T>::linmin(double* P1,double* P2,double* _n)
{
   copyvec(P1,ini_vec); 
   copyvec(_n,direc);


  
   TOL = 1e-8;


   minimize(-0.4,0.4);
   double lambda = m;

   double *temp;
   temp = new double[dimension];
   scalar_product(lambda,_n,temp);
   add_vec(P1,temp,P2);

   curr_val =fm;


   
 //  for(int i=0;i<dimension;i++)
 //  {
 //      cout<<P1[i]<<"\t"<<_n[i]<<"\n";
 //  }
 //
 
    
   delete temp;

//   cout<<P2[0]<<"\t"<<P2[1]<<"\t"<<curr_val<<"\n";

    return 0;
}

template<class T>
double multi_minimizer<T>::one_d(double lambda)
{
    double f=0;

    double *lam_n, *temp;

    lam_n = new double[dimension];
    temp = new double[dimension];


    scalar_product(lambda,direc,lam_n);

    add_vec(ini_vec,lam_n,temp);


    vec v;

    for(int i=0;i<dimension;i++)
    {

        v.push_back(temp[i]);
    }

    f = _func.Function(v);


  //  cout<<_func.Function(v)<<"\t"<<"here\n";

    delete lam_n;
    delete temp;

    return f;
}

template<class T>
int multi_minimizer<T>::set_tolerance(double _tol)
{
    TOL = _tol;
    return 0;
}

template<class T>
int multi_minimizer<T>::free_minimizer()
{
    delete [] P;
    delete [] n;
    delete [] Identity;
    delete change;
    delete ini_vec;
    delete direc;

    return 0;
}

template<class T>
double multi_minimizer<T>::metric(double* x,double* y)
{
    double f = 0;

    for(int i=0;i<dimension;i++)
    {
        f = f + pow((x[i]-y[i]),2.0);
    }

    f = sqrt(f);

    return f;
}

template<class T>
int multi_minimizer<T>::copyvec(double* x,double* y)
{
    for(int i=0;i<dimension;i++)
    {
        y[i] = x[i];
    }

    return 0;
}

template<class T>
int multi_minimizer<T>::subtract_vec(double* x,double* y,double *z)
{
 //   z = new double[dimension];

    for(int i=0;i<dimension;i++)
    {
        z[i] = x[i]-y[i];
    }

    return 0;
}

template<class T>
int multi_minimizer<T>::add_vec(double* x,double* y,double *z)
{
 //   z = new double[dimension];

    for(int i=0;i<dimension;i++)
    {
        z[i] = x[i]+y[i];
    }

    return 0;
}

template<class T>
double multi_minimizer<T>::dot_product(double* x,double* y)
{
    double f = 0;

    for(int i=0;i<dimension;i++)
    {
        f = f + x[i]*y[i];
    }


    return f;
}

template<class T>
int multi_minimizer<T>::scalar_product(double s,double* x,double *y)
{
  //  y = new double[dimension];

    for(int i=0;i<dimension;i++)
    {
        y[i] = s*x[i];
    }


    return 0;
}

template<class T>
int multi_minimizer<T>::copymatrix(double** A,double** B)
{
    for(int i=0;i<dimension;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            B[i][j] = A[i][j];

        }
    }

    return 0;
}

template<class T>
int multi_minimizer<T>::normalize(double* V)
{
    double N = dot_product(V,V);

    for(int i=0;i<dimension;i++)
    {
        if(N!=0)
        {
            V[i] = 1/sqrt(N)*V[i];
        }
    }

    return 0;

}

template<class T>
int multi_minimizer<T>::maximum_index(double *x)
{
    int ind=0;

    double max=x[0];

    for(int i=1;i<ind;i++)
    {
        if(x[i]>max)
        {
            max = x[i];
            ind = i;
        }

    }


    return ind;
}

double Rosenbrock(vec x)
{


    double f = pow((10.0-x[0]),2.0)+100.0*pow((x[1]-x[0]*x[0]),2.0);


    return f;
}

double Rosenbrock_multi(vec x)
{


    double f=0.0;

        for(int i=0;i<2;i++)
        {
            double term = 100*pow((pow(x[2*(i+1)-1],2.0)-x[2*(i+1)]),2.0)+pow((x[2*(i+1)-1]-1.0),2.0);

            f = f + term;

        }

    return f;
}


double test_func(double x)
{
    double f = x*x*x/3.0 - x*x/2.0-x-1;



    return f;
}

double spherical(vec x)
{
    double f=0.0;

    for(int i=0;i<30;i++)
    {
        f = f + (x[i]*x[i]);
    }

    return f;

}


template class multi_minimizer<chisq::SK_IV>;
