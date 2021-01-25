//
//  probability.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 21/11/20.
//

#include "probability.hpp"
#include "read_files.hpp"
#include "numerical.hpp"
#include <fstream>
#include<cmath>

Probability::~Probability()
{
  //  delete[] ff;
  //  delete rhosol;
    
    
}

/* ***********************************************************************************************************************/

/* This function reads the solar density and production distribution data and interpolates for averaging out
    in a later stage
 The function returns the arrays rhosol and ff. */
/* ***********************************************************************************************************************/

int Probability::interpolate_data()
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
        
    /* Define the interpoltion function for the solar density**/
    
  //  Cubic_interpolator inter_rho_sol;
  //  inter_rho_sol.set_cubic_spline(rho_sol.density.data[0],rho_sol.density.data[1],false);
    
    Bary_interpolator::non_interp Baryc(rho_sol.density.data[0].data(),rho_sol.density.data[1].data(),rho_sol.density.data[0].size());
    
    
 
    rhosol = new double[prod_data[0].size()];
    prod_data_n_row = int(prod_data[0].size());
    
    ff = new vec[8];
    std::ofstream ofil;
    ofil.open("debug_data/test.dat");

    
    
    for(int i=0;i<prod_data[0].size();i++)
    {
        rhosol[i] = Baryc(prod_data[0][i]);
      
    
        for(int j=0;j<8;j++)
        {
            double x = prod_data[0][i];
            ff[j].push_back(inter_prod[j+1].interpolate(fabs(x)));
       //     if(j==0)
            ofil<<j<<"\t"<<inter_prod[j+1].interpolate(fabs(x))<<"\t"<<x<<"\n";


            
        }
    
    
    
    }
    
    

    rho_sol.clean_memory();
    
    
    delete[] inter_prod;
    delete[] prod_data;
    
    
    return 0;
}

/* ***********************************************************************************************************************/

/* This function calculates the the probability for a array fo energy values is a range and return the corresponding probability vector in a data file */

/* ***********************************************************************************************************************/


int Probability::Probability_curve(double E_l, double E_h, double E_step, int flav,bool loga)
{
    vec Ener;
    for(double E=E_l;E<=E_h;E=E+E_step)
    {
        if(loga==false)
        {
            Ener.push_back(E);
        }
        else
        {
       //     std::cout<<pow(10,E)<<"\n";

            Ener.push_back(pow(10,E));
        }
    }
    
    Calculate_probability(Ener,flav);
    
    
    
    return 0;
}

/* ***********************************************************************************************************************/
/* This function calculates the probability for both day and night */

/* ***********************************************************************************************************************/



double Probability::Calculate_probability(double E,int time_in,int flav)
{
    Energy = E;
    fin_flav = flav;
    prob_inside_sun();
    read_regen();
    regeneration_earth();

        
    return 0;
    
}
/* ***********************************************************************************************************************/

/* This function calculates the probability for a vector of array for both day and night */

/* ***********************************************************************************************************************/





int Probability::Calculate_probability(vec E,int flav)
{
    std::cout<<"Initializing Probability ...\n";
    
    fin_flav = flav;
    interpolate_data();
    read_regen();

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


/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/

int Probability::Print_probability()
{
    
   
    
    std::string prod_name[8] = {"pp","pep","hep","7Be","8B","N","O","F"};
    
    std::string exp_name[4] = {"SK","SAGE","SNO","Gallex"};
    
    std::ofstream **ofl;
    
    ofl = new std::ofstream*[8];
    
    for(int i=0;i<8;i++)
    {
        ofl[i] = new std::ofstream[4];
    }
    
    for(int i=0;i<8;i++)
    {
        for(int j=0;j<4;j++)
        {
            std::string out_name = file_path + outfile + "_" + prod_name[i]+ "_" + exp_name[j] + ".dat";
            
            ofl[i][j].open(out_name);
            
            for(int k=0;k<E_vec.size();k++)
            {
                if(j==0)
                {
                  
                    ofl[i][j]<<E_vec[k]<<"  "<<E_p_day[i][k]<<"  ";
                        for(int l=0;l<7;l++)
                        {
                        ofl[i][j]<<E_p_night[i][l][k]<<"  ";

                        
                        }
                    
                    ofl[i][j]<<"\n";
                    
                    
                }
                else if(j>0)
                {
                    int bin_ind = j+6;
                    ofl[i][j]<<E_vec[k]<<"  "<<E_p_day[i][k]<<"  "<<E_p_night[i][bin_ind][k]<<"\n";
                }
                
            }
            
            
            ofl[i][j].close();
            
        }
        
        
    }
    
    std::cout<<"Writing Complete\n";
    
    for(int i =0;i<8;i++)
    {
        delete [] ofl[i];
    }
    
    
    delete [] ofl;
    
    
    return 0;
    
}





/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/




int Probability::prob_inside_sun()
{
    
   
    /* Oscillation parameters */
    
        double c12 = cos(osc_params[0]);
        double s12 = sin(osc_params[0]);
    
        double c_2th12 = c12*c12-s12*s12;
        double s_2th12 = 2*s12*c12;

    //    double s2_th12 = s12*s12;
    
        double c13 = cos(osc_params[1]);
        double s13 = sin(osc_params[1]);
    
        double s2_th13 = s13*s13;
    
    
    
    
        double *fmed;
        fmed = new double[8];
    
        for(int i=0;i<8;i++)
        {
        
            fmed[i] = 0.0;
            pmed[i] = 0.0;
        
        }
    

    
    /* Average over the production point */
    
        for(int i=0;i<prod_data_n_row-1;i++)
        {
            double temp_rhosol = (1.-s2_th13)*pow(10,rhosol[i]);

            double aa = temp_rhosol-0.25e-6*osc_params[4]/Energy*c_2th12;
            double bb = 0.25e-6*osc_params[4]/Energy*s_2th12;
            
            
            
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
            E_p_day.push_back(s2_th13*s2_th13+(1-s2_th13)*(1-s2_th13)*pmed[i]);
            
    
        }
        
 //   std::cout<<"one\n";
    
    return 0;
}


/* ***********************************************************************************************************************/

/* This function calculates the regeneration due to the earth matters */

/* ***********************************************************************************************************************/


int Probability::regeneration_earth()
{
    
    double tg = pow(tan(osc_params[0]),2.0);
    double logdm = log10(0.25e-6*osc_params[4]/Energy);//0.15e-13;//log10(1.27*osc_params[4]/Energy);
    
    
  
    vec *preg;
    vec cords[2];
    preg = new vec[10];
        
   
    
    for(int i=0;i<exp[0].n_row;i++)
    {
        cords[0].push_back(exp[0].regen.data[0][i]);
        cords[1].push_back(exp[0].regen.data[1][i]);
        
        preg[0].push_back(exp[0].regen.data[2][i]);
        preg[1].push_back(exp[0].regen.data[3][i]);
        preg[2].push_back(exp[0].regen.data[4][i]);
        preg[3].push_back(exp[0].regen.data[5][i]);
        preg[4].push_back(exp[0].regen.data[6][i]);
        preg[5].push_back(exp[0].regen.data[7][i]);
        preg[6].push_back(exp[0].regen.data[8][i]);

        preg[7].push_back(exp[1].regen.data[2][i]);
        preg[8].push_back(exp[2].regen.data[2][i]);
        preg[9].push_back(exp[3].regen.data[2][i]);


        
    }
    

    
    
    two_linear_interp *inter;
    
    inter = new two_linear_interp[10];

    double p2e[9];
    
    for(int i =0;i<10;i++)
    {
        inter[i].set_logscale_y();
        inter[i].input_vectors(cords[0],cords[1], preg[i]);
        p2e[i] = inter[i].interp2d(tg,logdm);
    //    std::cout<<tg<<"\t"<<logdm<<"\t"<<p2e[i]<<"\n";

    }
    
    double s2_13 = sin(osc_params[1])*sin(osc_params[1]);
    double s2_12 = sin(osc_params[0])*sin(osc_params[0]);
    double c_212 = cos(2*osc_params[0]);
    
    for(int i=0;i<8;i++)
    {
        for(int j=0;j<10;j++)
        {
            pnight[i][j] = pow(s2_13,2) + pow((1-s2_13),2)*(pmed[i]-s2_12+p2e[j]*(1.0-2.0*pmed[i]))/c_212;
       //     std::cout<<logdm<<"\t"<<pnight[i][j]<<"\n";
            
            
        }
        
    }
    
    

    delete [] inter;
    
    delete[] preg;
    
    return 0;
}


/* ***********************************************************************************************************************/
        
/* This function reads the data files for the regeneration effect */

/* ***********************************************************************************************************************/
int Probability::read_regen()
{
    
    exp = new regeneration[4];
    
    exp[0].read_regeneration("regeneration/pregsk.dat");   /* SK has 6+1 bins   */
    exp[1].read_regeneration("regeneration/pregsg.dat");   /* Sage has 1 bin */
    exp[2].read_regeneration("regeneration/pregsno.dat");   /* SNO has 1 bin */
    exp[3].read_regeneration("regeneration/preggno.dat");    /*                  */

    
    
    
    return 0;
}







/* ***********************************************************************************************************************/
/*This function clears the memory for the temporary variables */

/* ***********************************************************************************************************************/






int Probability::free_data()
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


