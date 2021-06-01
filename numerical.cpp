//
//  numerical.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 25/11/20.
//

#include "numerical.hpp"


//interpolator::interpolator()
//{
    
//}

/* ***********************************************************************************************************************/

/* ***********************************************************************************************************************/


int Cubic_interpolator::free_cubic()
{
    delete &spline;
    
    return 0;
}

/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/


Cubic_interpolator::~Cubic_interpolator()
{
    
    
}
/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/

int Cubic_interpolator::set_cubic_spline(vec data_x, vec data_y,bool logarithmic)
{

    
    if(logarithmic==false)
    {
        double step0,step1;

        for(int i=0;i<data_x.size()-3;i++)
        {
            step0 = data_x[i+1]-data_x[i];
            step1 = data_x[i+2]-data_x[i+1];
            
      /*      if(fabs(step0-step1)>1e-16)
            {
                std::cout<<"ERROR!! Cubic spline not valid for non-equispaced data!!\n";
                exit(-1);
                
            }  */
            

        }
        
        
        double step = data_x[1]-data_x[0];

        interp tempo_spline(data_y.data(), data_y.size(),data_x[0],step);
        spline = tempo_spline;
        

    }
    if(logarithmic != false)
    {
        vec log_x;
        for(int i = 0;i<data_x.size();i++)
        {
            log_x.push_back(log(data_x[i]));
            
        }
        
        double step0,step1;

        for(int i=0;i<data_x.size()-3;i++)
        {
            step0 = data_x[i+1]-data_x[i];
            step1 = data_x[i+2]-data_x[i+1];
            
         /*   if(fabs(step0-step1)>1e-16)
            {
                std::cout<<"ERROR!! Cubic spline not valid for non-equispaced data!!\n";
                exit(-1);
                
            }  */
        }
        
        
        double step = log_x[1]-log_x[0];
        interp tempo_spline(data_y.data(),data_y.size(),log_x[0],step);
        spline = tempo_spline;
        
        
        
        

    }

    return 0;

}

double Cubic_interpolator::interpolate(double x)
{

    double p = 0;
    if(std::isnan(x)!=true && std::isinf(x)!=true)
    {   
        p = spline(x);
    }
    else
    {
        std::cout<<"WARNING!!Encountered infinity or NaN so set interpoalte to 0\n";
    }


    return p;

}
/* ***********************************************************************************************************************/
/* ***********************************************************************************************************************/
/* ***********************************************************************************************************************/


/*Bary_interpolator::~Bary_interpolator()
{
    
    
}
*/

int Bary_interpolator::set_barycentric(vec data_x,vec data_y,bool logarithmic)
{
    if(logarithmic==false)
    {

        non_interp bary_inter(data_y.data(), data_x.data(),data_x.size());
        bary = bary_inter;

        
    }
    if(logarithmic!=false)
    {
        vec log_x;
        for(int i = 0;i<data_x.size();i++)
        {
            log_x.push_back(log(data_x[i]));
            
        }
        
        non_interp bary_inter(data_y.data(), data_x.data(),data_x.size());
        bary = bary_inter;
        
    }
    
    
    
    return 0;
}
/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/

double Bary_interpolator::interpolate(double x)
{
    double p = bary(x);

    return p;
}

/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/

int two_linear_interp::make_grid()
{
    double temp_x0 = data[0][0];
    int i=0;
    while(data[0][i]==temp_x0)
    {
        data_y.push_back(data[1][i]);
        i++;
        
    }
    
    
    i=0;
    
    while(i<data[0].size())
    {
        if(data[0][i]!=temp_x0)
        {
            data_x.push_back(data[0][i]);
            temp_x0 = data[0][i];
        }
        
       i++; 
    }
   


    x_dim = int(data_x.size());
    y_dim = int(data_y.size());
  
    
    
    f_grid = new double*[x_dim];
    for(int j=0;j<x_dim;j++)
    {
        f_grid[j] = new double[y_dim];
        
    }

    
    
    for(int j=0;j<x_dim;j++)
    {
        for(int k=0;k<y_dim;k++)
        {
            f_grid[j][k] = data[2][y_dim*j+k];

        }
        
    }
    
//    std::cout<<i<<"\n";
    

    return 0;
}







/* ***********************************************************************************************************************/

double two_linear_interp::interpolate(double x,double y)
{
    int i = 0;
    
    double x_lo=0.0,x_ho=0.0;
    double y_lo=0.0,y_ho=0.0;
    
    while(i< x_dim && data_x[i]<x)
    {
        x_lo = data_x[i];
        x_ho = data_x[i+1];
        i++;
        
    }
    
    
    int j=0;
    
    while(j<y_dim && data_y[j]<=y)
    {
        y_lo = data_y[j];
        y_ho = data_y[j+1];
      //  printf("%d\t%g\t%g\n",j,y,data_y[j]);

        j++;
    }
    

    
 //   printf("%d\t%d\n",i,j);


    
    double f_x1y1;// = f_grid[i][j];
    double f_x1y2;// = f_grid[i][j+1];
    double f_x2y1;// = f_grid[i+1][j];
    double f_x2y2;// = f_grid[i+1][j+1];
    
    double f;
 
    if(i==x_dim)
    {
        x_lo = data_x[i-1];
        x_ho = data_x[i];
        
        f_x1y1 = f_grid[i-1][j];
        f_x1y2 = f_grid[i-1][j+1];
        f_x2y1 = f_grid[i][j];
        f_x2y2 = f_grid[i][j+1];

        
    }
        
    if(j>=y_dim-1)
    {
 //       std::cout<<"here\n";
        y_lo = data_y[y_dim-2];
        y_ho = data_y[y_dim-1];
        
        
        f_x1y1 = f_grid[i][y_dim-2];
        f_x1y2 = f_grid[i][y_dim-1];
        f_x2y1 = f_grid[i+1][y_dim-2];
        f_x2y2 = f_grid[i+1][y_dim-1];
        

        
        double f_x1y=f_x1y1+(f_x1y2-f_x1y1)*(y_ho-y_lo)/(y_ho-y_lo);
        double f_x2y=f_x2y1+(f_x2y2-f_x2y1)*(y_ho-y_lo)/(y_ho-y_lo);
        
        f = f_x1y+(f_x2y-f_x1y)*(x-x_lo)/(x_ho-x_lo);
        
 //       std::cout<<y<<"\t"<<f<<"\n";

        
      //  std::cout<<f_x1y<<"\t"<<f_x2y<<"\t"<<f<<"\n";
 
    }
    else
    {
    
   // std::cout<<f_x1y1<<"\t"<<f_x1y2<<"\t"<<f_x2y1<<"\t"<<f_x2y2<<"\n";
        
        f_x1y1 = f_grid[i][j];
        f_x1y2 = f_grid[i][j+1];
        f_x2y1 = f_grid[i+1][j];
        f_x2y2 = f_grid[i+1][j+1];

    
    double f_x1y=f_x1y1+(f_x1y2-f_x1y1)*(y-y_lo)/(y_ho-y_lo);
    double f_x2y=f_x2y1+(f_x2y2-f_x2y1)*(y-y_lo)/(y_ho-y_lo);
    
    f = f_x1y + (f_x2y-f_x1y)*(x-x_lo)/(x_ho-x_lo);
    
    
  //  if(j==y_dim)
  //  {
   //     std::cout<<x_lo<<"\t"<<x_ho<<"\t"<<y_lo<<"\t"<<y_ho<<"\n";
  //      std::cout<<f_x1y1<<"\t"<<f_x1y2<<"\t"<<f_x2y1<<"\t"<<f_x2y2<<"\n";
    //    std::cout<<f<<"\n";

//std::cout<<y<<"\t"<<f<<"\n";
        
    }
    
//    std::cout<<y<<"\t"<<f<<"\n";

    
    return f;
}


/* ***********************************************************************************************************************/

/* ***********************************************************************************************************************/
double two_linear_interp::interp2d(double x, double y)
{
    if(logx==true)
    {
        convert_logx();
        
    }
    if(logy==true)
    {
        convert_logy();
        
    }
    
    
    
    make_grid();
    double f = interpolate(x, y);
    free_grid();
    return f;
}



/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/


int two_linear_interp::free_grid()
{
    for (int i = 0; i < x_dim; ++i)
    {
        delete [] f_grid[i];
    }
    
    
    delete [] f_grid;

    delete [] data;

    
    return 0;
}


/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/

int two_linear_interp::input_vectors(vec x, vec y, vec z)
{
    data = new vec[3];
    
    for(int i=0;i<x.size();i++)
    {
        data[0].push_back(x[i]);
        data[1].push_back(y[i]);
        data[2].push_back(z[i]);

        
    }
    
    
    return 0;
}

int two_linear_interp::set_logscale_x()
{
    logx = true;
    
    return 0;
    
}

int two_linear_interp::set_logscale_y()
{
    logy = true;
    return 0;
    
}

int two_linear_interp::convert_logx()
{
    for(int i=0;i<data[0].size();i++)
    {
        data[0][i] = log10(data[0][i]);
        
    }
    
    
    
    return 0;
}

int two_linear_interp::convert_logy()
{
    for(int i=0;i<data[1].size();i++)
    {
        data[1][i] = log10(data[1][i]);
        
    }
    
    
    
    return 0;
}


double square(double x)
{
    return x*x;
    
}
