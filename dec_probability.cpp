//
//  visiible_anti.cpp
//  Solar_neutrino
//
//  Created by Dipyaman Pramanik on 19/01/21.
//

#include "dec_probability.hpp"
#include <fstream>
#include <cmath>
#include <iostream>
#include <mutex>
#include <string>



int dec_prob::Init_prob(std::string type_name,int num_param)
{
    which_type = type_name;
    
    regeneration = SOL_NO;

    n_of_params = num_param;
    
    init_interpolate_flux();
    
    calculate_total_flux();
    
    interpolate_data(); 
    osc_params = new double[num_param];
    
    return 0;
}

dec_prob::~dec_prob()
{
    
}


int dec_prob::Wrap_oscparams()
{
    osc_params[0] = atan(sqrt(Tan_Th12));
    osc_params[1] = asin(sqrt(Th13));

    osc_params[3] = Delta;
    osc_params[4] = Dm21;
    osc_params[6] = Tau1;
    osc_params[7] = Tau2;
    
    
    return 0;
    
    
    
}

int dec_prob::init_interpolate_flux()
{
    flux.read_flux("flux/b8spectrum.txt");
    
    vec flux_x,flux_y;
    
    for(int i =0;i<flux.Flux.data[0].size();i++)
    {
        flux_x.push_back(flux.Flux.data[0][i]);
        flux_y.push_back(flux.Flux.data[1][i]);
        
    }
    
    flux_interpolator.set_cubic_spline(flux_x,flux_y,false);
    
    
   	x_max = flux_x[flux_x.size()-2];

    
    return 0;
    
    
}

int dec_prob::interpolate_data()
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

double dec_prob::Propagation(double E)
{
    double Gamma_i;
  
 

    if(life>1e-16)
    {
        Gamma_i = 1/(E*life);
  
    }
    else
    {   
        Gamma_i = 0.0;
    }
   
    
    
    
    double p = (1.0-exp(-Gamma_i*L*4.98e-4));



    P_ij = p;

 //   std::cout<<Gamma_i*L*4.98e-4<<"\t"<<p<<"\n";
    
    return p;
}

double dec_prob::survival(double E)
{
    double Gamma_i;

    if(life>1e-16)
    {
        Gamma_i = 1/(E*life);
    }
    else
    {
        Gamma_i = 0;
    }
   
    double p = exp(-Gamma_i*L*4.98e-4);



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

int dec_prob::calculate_total_flux()
{
	
	
	return 0;
}

double dec_prob::integrand(double E)
{
    pe2[4] = 0;
    pe1[4] = 0;
    prob_inside_sun(E);
    
    double E_beta = Energy; 
    double E_alpha = E;
    
    weighted_differential W_rate(which_type);
    
    /*ONLY BORON IS BEING CONSIDERED FOR SIMPLICITY. LATER WE CAN ADD ALL SOURCES*/
    
    double res;
    
    double w;

  /*  if(particle==1)
    {
        w = E_beta/square(E_alpha);
    }
    else
    {
        w = 1/E_alpha*(1-E_beta/E_alpha);
    }
*/
    if(E<x_max)
    {
        if(_channel==1)
        {
		    res = flux_interpolator.interpolate(E)*pe2[4]*Propagation(E)*W_rate.weighted_rate(osc_params[3],E_alpha,E_beta,particle);
        }
        else if(_channel==2)
        {
		    res = flux_interpolator.interpolate(E)*pe2[4]*Propagation(E)*W_rate.weighted_rate(osc_params[3],E_alpha,E_beta,particle);
        }
        else if(_channel==3)
        {
        
     //       std::cout<<"Energy  =====  "<<E<<"\n";

		    res = flux_interpolator.interpolate(E)*Propagation(E)*W_rate.weighted_rate(osc_params[3],E_alpha,E_beta,particle);
        }
        else if (_channel==4)
        {
            res = flux_interpolator.interpolate(E)*pe2[4]*Propagation(E)*W_rate.weighted_rate(osc_params[3],E_alpha,E_beta,particle);


        }
        else if(_channel==5)
        {

            res = flux_interpolator.interpolate(E)*pe1[4]*Propagation(E)*W_rate.weighted_rate(osc_params[3],E_alpha,E_beta,particle);
        }
       // std::cout<<w<<"\t"<<W_rate.weighted_rate(osc_params[3],E_alpha,E_beta,particle)<<"\n";
	}
    else
	{
		res = 0.0;
	}
    //~ std::cout<<osc_params[3]<<"\n";
    
    
    
    return res;
}

double dec_prob::integrate()
{
	E_max = Energy/square(osc_params[3]);

    if(E_max>x_max)
    {
        E_max = x_max;
    }
	

    	
    double h = (E_max-Energy)/4.0;
    
//    std::cout<<x_max<<"\t"<<E_max<<"\t"<<Energy<<"\n\n";

  //  double integral= 3.0*h/8.0*(integrand(Energy)+3.0*integrand((2.0*Energy+E_max)/3.0)+3.0*integrand((Energy+2.0*E_max)/3.0)+integrand(E_max));


    double integral = (integrand(Energy)+integrand(E_max));

    for(int i=0;i<3;i++)
    {
        double term = 2.0*integrand(Energy+(i+1)*h);
        integral = integral + term;
    //    std::cout<<term<<std::endl;
        
    }


    
    return h/2.0*integral;
}



double dec_prob::integrate_flux()
{
    E_max = 14.8;

    double E_min = 4.0;

    double h = (E_max-E_min)/4.0;

    double integral = (flux_interpolator.interpolate(E_min)+flux_interpolator.interpolate(E_max));

    for( int i =0; i<3; i++)
    {
        double term = 2.0*(flux_interpolator.interpolate(E_min+(i+1)*h));
        integral = integral + term;
    }

    return h/2.0*integral;
 


}




double dec_prob::Calculate_decayed_flux(double E,int fin_flav,int parti,int channel)
{
    Wrap_oscparams();

    _channel = channel;
    

 //   std::cout<<"Energy ============================="<<E<<"\n";

    Energy = E;
    particle = parti;


    prob_inside_sun(Energy);


    double c13 = cos(osc_params[1]);
    double s13 = sin(osc_params[1]);
    double s12 = sin(osc_params[0]);
    double c12 = cos(osc_params[0]);


    double s13_4 = s13*s13*s13*s13;
    double c13_4 = c13*c13*c13*c13;

    double s13_2 = s13*s13;
    double c13_2 = c13*c13;

    double s12_2 = s12*s12;
    double c12_2 = c12*c12;


    double Flux_dec;

    if(channel==1)
    {
        if(fin_flav==0)
        {
            if(parti==1)
            {
                life = osc_params[6];  /* Lifetime for the decay into neutrino */
                

                double P_inv_prime = pe1[4]*square(c12) + pe2[4]*square(s12)*survival(Energy);

                double flux_0 = flux_interpolator.interpolate(Energy);

                double flux_inv_prime = flux_0*P_inv_prime;
                


                double flux_inv = c13_4*flux_inv_prime;

                double flux_vis;

                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_4*integrate()*square(c12);
                }

                
          //         Flux_dec = flux_vis;              

                Flux_dec = flux_inv_prime + flux_vis + s13_4*flux_0;

          //      std::cout<<flux_inv_prime<<"\t"<<flux_vis<<"\t"<<Propagation(E)<<"\t"<<survival(E)<<std::endl;

                
            }
            if(parti==-1)
            {
                
                life = osc_params[7];   /* Lifetime for the deacy into anti-neutrino */
                 
        
                double P_1e_earth = c13_2*square(cos(osc_params[0]));
   
                if(regeneration!=SOL_YES)
                {
                    Flux_dec = c13_4*integrate()*P_1e_earth;
                }


            }
        }
    }
    else if(channel==2)
    {

        if(fin_flav==0)
        {
            if(parti==1)
            {

//                std::cout<<channel<<std::endl;

                life = osc_params[6];

                double P_inv_prime = pe1[4]*s12 + pe2[4]*c12*survival(Energy);

                double flux_0 = flux_interpolator.interpolate(Energy);

                double flux_inv_prime = flux_0*P_inv_prime;

                double flux_inv = c13_4*flux_inv_prime;



                double c13_2_s13_2 = c13*c13*s13*s13;

                double flux_vis;

                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2_s13_2*integrate();
                }

                Flux_dec = flux_inv + flux_vis;

            }
            if(parti==-1)
            {

                life = osc_params[7];


                double c13_2_s13_2 = c13*c13*s13*s13;
    
                
                Flux_dec = c13_2_s13_2*integrate();

            }
        }
    }
    else if(channel==3)
    {
        if(fin_flav==0)
        {
            if(parti==1)
            {
                life = osc_params[6];

                double P_inv_prime = pe1[4]*c12_2*c13_2+pe2[4]*c13_2*s12_2+s13_2*s13_2*survival(Energy);

                double flux_0 = flux_interpolator.interpolate(Energy);

                double flux_inv_prime = flux_0*P_inv_prime;


                double flux_inv = c13_4*flux_inv_prime;


                double c13_2_s13_2 = c13*c13*s13*s13;

                double flux_vis;

                if(regeneration!=SOL_YES)
                {
                    flux_vis = c13_2*s13_2*c12_2*integrate();
                }

                Flux_dec = flux_inv + flux_vis;



            }
 
            if(parti==-1)
            {
                life = osc_params[7];

   
             //   std::cout<<s13_2<<std::endl;

                double P_1e_earth = c13_2*square(cos(osc_params[0]));

                if(regeneration!=SOL_YES)
                {

                    Flux_dec = s13_2*integrate()*P_1e_earth;
                }


            }


        }

    }

    else if(channel==4)
    {
        if(fin_flav==0)
        {
            if(parti==-1)
            {
                life= osc_params[7];


                if(regeneration!=SOL_YES)
                {
                    Flux_dec = c13_2*integrate()*s13_2;

                }

            }

        }

    }
    else if(channel==5)
    {
        if(fin_flav==0)
        {
            if(parti==-1)
            {
                life = osc_params[7];

                if(regeneration!=SOL_YES)
                {
                    Flux_dec = c13_2*integrate()*s13_2;
                }


            }

        }
    }

    else
    {
        std::cout<<"ERROR!! Invalid Channel..."<<std::endl;
        exit(-1);
    }


    return Flux_dec;
}



int dec_prob::Decayed_Flux_curve(double E_l,double E_h,double E_step,bool loga)
{
    
    
    return 0;
}


int dec_prob::prob_inside_sun(double E)
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
    double pe1med[8];
    double pe2med[8];
    double pmed1[8];

    for(int i=0;i<8;i++)
    {
    
        fmed[i] = 0.0;
        pe1med[i] = 0.0;
        pe2med[i] = 0.0;
        pmed1[i] = 0.0;
    }



/* Average over the production point */

    for(int i=0;i<prod_data_n_row-1;i++)
    {
        double temp_rhosol = square(c13)*pow(10,rhosol[i]);

 //       double aa = temp_rhosol-0.25e-6*osc_params[4]/Energy*c_2th12;
 //       double bb = 0.25e-6*osc_params[4]/Energy*s_2th12;
  
        double aa = osc_params[4]*c_2th12-4.0*E*1e6*temp_rhosol;
        double bb = osc_params[4]*s_2th12;
 
        
  //      double aa = temp_rhosol-Energy*c_2th12;
  //      double bb = Energy*s_2th12;
 
  //      double s_2tm = bb/(sqrt(aa*aa+bb*bb));
        double c_2tm = aa/(sqrt(aa*aa+bb*bb));
   
        double s_tm = 0.5*(1-c_2tm);
        double c_tm = 0.5*(1+c_2tm);



        double pee = 0.5*(1+c_2tm*c_2th12);
            

    
        for(int j=0;j<8;j++)
        {


            fmed[j] = fmed[j] + ff[j][i];
            pe1med[j] = pe1med[j] + ff[j][i]*c_tm;
            pe2med[j] = pe2med[j] + ff[j][i]*s_tm;
            pmed1[j] = pmed1[j] + ff[j][i]*pee;
        
        }
    
    
    
    }
    

    for(int i=0;i<8;i++)
    {
        pe1[i] = pe1med[i]/fmed[i];
        pe2[i] = pe2med[i]/fmed[i];


    }
    
    
    
    
    
    return 0;
}


int dec_prob::free_data()
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


