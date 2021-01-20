//
//  visiible_anti.cpp
//  Solar_neutrino
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#include "visiible_anti.hpp"
#include "interactions.hpp"
#include "probability.hpp"
#include "read_files.hpp"
#include "numerical.hpp"
#include <fstream>
#include<cmath>

visible_anti_prob::visiible_anti_prob(std::string type_name)
{
    which_type = type_name;
    
}

visible_anti_prob::~visible_anti_prob()
{
    
}

int visible_anti_prob::interpolate_data()
{
    /* Read the production distributions inside the sun */
    
        production distri;
        distri.read_prod_dist("prod_dist/bs05opflux2.dat");
    
    /* Stores the data from the file into eight different vectors for corresponding
     production mechanism**/
  
        vec *prod_data;
        int distri_n_col = distri.n_col;
        prod_data = new vec[distri_n_col];
    
        for(int i=0;i<distri.n_row;i++)
        {
            for(int j=0;j<distri.n_col;j++)
            {
                prod_data[j].push_back(distri.prod.data[j][i]);

            }
        


        }
    
    distri.clean_memory();
    
    /* Next we define the interpolation functions for various production mechanism **/
 
        Cubic_interpolator *inter_prod;
        inter_prod = new Cubic_interpolator[distri.n_col];
    
        for(int i =1;i<distri.n_col;i++)
        {
            inter_prod[i].set_cubic_spline(prod_data[0],prod_data[i],false);

        }
    
    
    /* Read the solar density profile */
  
        solar_density rho_sol;
        rho_sol.read_density("solar_density/bs2005op.dat");
        
    /* Define the interpoltion function for the solar density*/
    
    Bary_interpolator::non_interp Baryc(rho_sol.density.data[0].data(),rho_sol.density.data[1].data(),rho_sol.density.data[0].size());
    
    
 
    rhosol = new double[prod_data[0].size()];
    prod_data_n_row = int(prod_data[0].size());
    
    ff = new vec[8];


    
    
    for(int i=0;i<prod_data[0].size();i++)
    {
        rhosol[i] = Baryc(prod_data[0][i]);
      
    
        for(int j=0;j<8;j++)
        {
            double x = prod_data[0][i];
            ff[j].push_back(inter_prod[j+1].interpolate(fabs(x)));
            
        }
    
    
    
    }
    
    

    rho_sol.clean_memory();
    
    
    delete[] inter_prod;
    delete[] prod_data;
    
    
    return 0;
}

double visible_anti_prob::Propagation();
{
    double Gamma = osc_params[3];
    double p = (1.0-exp(-Gamma_i*L));
    P_ij = p;
    
    return p;
}

int visible_anti_prob::Calculate_probability(vec E)
{
    std::cout<<"Initializing Probability ...\n";
    
    interpolate_data();
//    read_regen();

    std::cout<<"Calculating Probabilities ...\n";
    
    for(int i=0;i<E.size();i++)
    {
        Energy = E[i];
        E_vec.push_back(E[i]);
        prob_inside_sun();
        regeneration_earth();
        for(int j=0;j<8;j++)
        {
            E_p_day[j].push_back(pday[j]);
            for(int k=0;k<9;k++)
            {
                E_p_night[j][k].push_back(pnight[j][k]);
                
            }
            
        }
    }
    
    if(print_prob == 1)
    {
        std::cout<<"Writing the probability output in files\n";
      
        Print_probability();
        
    }
    
    

    return 0;
}

/* We were writing here, it is not done. Import the Integration class next to do the
integration over alpha. */

double visible_anti_prob::integrand(double E)
{
    double res;
    
    prob_inside_sun();
    
    double c13 = cos(osc_params[1]);

    double c13_4 = c13*c13*c13*c13;
    
    weighted_differential W_rate(which_type);
    
    double res = func
    
    
    
    
    return res;
}

double visible_anti_prob::Calculate_probability(double E)
{
    double prob;
    
    Energy = E;
   
    
    double w = integrand(E);
    
    for(int j=0;j<8;j++)
    {
        
        
    }
    

    
    
    
    return prob;
}

int visible_anti_prob::Probability_curve(double E_l,double E_h,double E_step,bool loga)
{
    
    
    return 0;
}

int visible_anti_prob::prob_inside_sun()
{
    /* Oscillation parameters*/
    
    double c12 = cos(osc_params[0]);
    double s12 = sin(osc_params[0]);

    double c_2th12 = c12*c12-s12*s12;
    double s_2th12 = 2*s12*c12;

//    double s2_th12 = s12*s12;

    double c13 = cos(osc_params[1]);
    double s13 = sin(osc_params[1]);

    double s2_th13 = s13*s13;

    
    double fmed[8];

    for(int i=0;i<8;i++)
    {
    
        fmed[i] = 0.0;
        pmed[i] = 0.0;
    
    }



/* Average over the production point */

    for(int i=0;i<prod_data_n_row-1;i++)
    {
        double temp_rhosol = (1.-s2_th13)*pow(10,rhosol[i]);

        double aa = temp_rhosol-0.25e-6*osc_params[2]/Energy*c_2th12;
        double bb = 0.25e-6*osc_params[2]/Energy*s_2th12;
        
        
        
  //      double aa = temp_rhosol-Energy*c_2th12;
  //      double bb = Energy*s_2th12;
 //       double s_2tm = bb/(sqrt(aa*aa+bb*bb));
        double c_2tm = -aa/(sqrt(aa*aa+bb*bb));
    
        double pee = 0.5*(1+c_2tm*c_2th12);
    //    std::cout<<aa<<"\n";

    
    
        for(int j=0;j<8;j++)
        {


            fmed[j] = fmed[j]+ff[j][i];
            pmed[j] = pmed[j]+ff[j][i]*pee;
        
        }
    
    
    
    }
    
    double pmed1[8];

    for(int i=0;i<8;i++)
    {
        pmed[i] = pmed[i]/fmed[i];
        pday[i] = s2_th13*s2_th13+(1-s2_th13)*(1-s2_th13)*pmed[i];

    }
    
    
    
    
    
    return 0;
}

int visible_anti_prob::free_data();
{
    delete[] ff;
    delete[] rhosol;
    
    for(int i=0;i<4;i++)
    {
        exp[i].clean_memory();
        
    }
    
    delete[] exp;
    
    return 0;
}

int visible_anti_prob::Print_probability();
{
    return 0;
}
