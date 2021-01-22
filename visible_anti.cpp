//
//  visiible_anti.cpp
//  Solar_neutrino
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#include "visible_anti.hpp"
#include "interactions.hpp"
#include "probability.hpp"
#include "read_files.hpp"
#include "numerical.hpp"
#include <fstream>
#include<cmath>

visible_anti_prob::visible_anti_prob(std::string type_name,double e_max,int num_param)
{
    which_type = type_name;
    E_max = e_max;
    
    init_interpolate_flux();
    interpolate_data();
    osc_params = new double[num_param];
}

visible_anti_prob::~visible_anti_prob()
{
    
}

int visible_anti_prob::Wrap_oscparams()
{
    osc_params[0] = atan(sqrt(Tan_Th12));
    osc_params[1] = Th13;
    osc_params[3] = Delta;
    osc_params[4] = Dm21;
    osc_params[6] = Tau;
    
    
    return 0;
    
    
    
}

int visible_anti_prob::init_interpolate_flux()
{
    flux.read_flux("flux/b8spec-2006.dat");
    
    vec flux_x,flux_y;
    
    for(int i =0;i<flux.Flux.data[0].size();i++)
    {
        flux_x.push_back(flux.Flux.data[0][i]);
        flux_y.push_back(flux.Flux.data[1][i]);
        
    }
    
    flux_interpolator.set_cubic_spline(flux_x,flux_y,false);
    
    return 0;
    
    
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

double visible_anti_prob::Propagation(double E)
{
    double Gamma_i = 1/(E*osc_params[6]);
    double p = (1.0-exp(-Gamma_i*L));
    P_ij = p;
    
    return p;
}
/*
int visible_anti_prob::Calculate_decayed_flux(vec E)
{
    std::cout<<"Initializing Decayed Flux Calculator...\n";
    
    interpolate_data();
//    read_regen();

    std::cout<<"Calculating anti-neutrino flux profile due to visible decay ...\n";
    
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
        std::cout<<"Writing the anti-neutrino flux profile output in files\n";
      
        Print_probability();
        
    }
    
    

    return 0;
}

*/

double visible_anti_prob::integrand(double E)
{
    
    prob_inside_sun();
    
    
    
    weighted_differential W_rate(which_type);
    
    /*ONLY BORON IS BEING CONSIDERED FOR SIMPLICITY. LATER WE CAN ADD ALL SOURCES*/
    
    double res = flux_interpolator.interpolate(E)*pday[4]*square(sin(osc_params[0]))*W_rate.weighted_rate(osc_params[3],E,Energy);
    

    
    
    return res;
}

double visible_anti_prob::integrate()
{
    double h = (E_max-Energy);
    
    double integral= 3.0*h/8.0*(integrand(Energy)+3.0*integrand((2.0*Energy+E_max)/3.0)+3.0*integrand((Energy+2.0*E_max)/3.0)+integrand(E_max));
    
    return integral;
}

double visible_anti_prob::Calculate_decayed_flux(double E)
{
    Wrap_oscparams();
    
    Energy = E;
    
    double c13 = cos(osc_params[1]);

    double c13_4 = c13*c13*c13*c13;
    
    double P_1e_earth = square(cos(osc_params[0]));
   
    
    double Flux_dec = c13_4*integrate()*Propagation(E)*P_1e_earth;
    
    return Flux_dec;
}

int visible_anti_prob::Decayed_Flux_curve(double E_l,double E_h,double E_step,bool loga)
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

int visible_anti_prob::free_data()
{
    delete[] ff;
    delete[] rhosol;
    delete osc_params;
  /*
    for(int i=0;i<4;i++)
    {
        exp[i].clean_memory();
        
    }
    
    delete[] exp;
    */
    return 0;
}
/*
int visible_anti_prob::Print_probability();
{
    return 0;
}
*/



