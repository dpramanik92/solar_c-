//
//  main.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 20/11/20.
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <string.h>
#include "event.hpp"
#include "dec_probability.hpp"
#include "chisqmin.hpp"
#include "numerical.hpp"
#include "minimizer.hpp"


using namespace std;

typedef vector<double> vec;


double min(double,double);


double prior(double val,double bf,double sigma)
{
    double p = pow(((val-bf)/sigma),2.0);

    return p;
}

int main(int argc, const char * argv[]) {

    
    
    double E_max = 10.0;

    string file_prefix = string(argv[1]);

    string Experiment = string(argv[2]);

    string which_type = string(argv[3]);

    const char* chan = argv[4];

  //  string delta = string(argv[4]);
    
  //  double _delta = atof(argv[4]);

    string file_name = file_prefix +"_" + chan + "_" + Experiment + "_" + which_type + ".dat";


     cout<<"The output is being written to: "<<file_name<<endl;

    dec_prob Prob;
    Prob.Tan_Th12 = tan(33.56*(M_PI/180.0));
    Prob.Th13 = 0.02;
    Prob.Dm21 = 7.5e-5;
    Prob.Delta = 0;
    Prob.Tau1 = 6.6e-3;
    Prob.Tau2 = 6.6e-3;
    Prob.L = 1.0;  // Baseline in A.U. 
    Prob.Init_prob(which_type,8);   
   


     Event_generator _event;
    
     _event.resolution[0] = 0.0;
     _event.resolution[1] = 0.06;
     _event.resolution[2] = 0.0;
    
     _event.e_min = 1.8;
     _event.e_max = E_max;
     _event.n_bins = 100;
    

    
     _event.Set_probability_engine(Prob);
    
     file_reader exp_data,bkg_data;


      

     if(Experiment=="SK-IV")
     {
        string data_file = "exp_data/SuperK-IV_data.dat";
        string bkg_file = "exp_data/SuperK-IV_bkg.dat";

        exp_data.read_file(data_file);
        bkg_data.read_file(bkg_file);

        double Time = 2970.1; /* days */
        Time = Time*86400;
        double fiducial = 22.5; /* ktons */

        _event.efficiency = 0.05;

        _event.Normalization = 15.0*1e-41*1e31;

        _event.Exposure = Time*fiducial;

     }
     else if(Experiment=="Kamland")
     {

        string data_file = "exp_data/KamLAND_data.dat";
        string bkg_file = "exp_data/KamLAND_bkg.dat";

        exp_data.read_file(data_file);
        bkg_data.read_file(bkg_file);

        double Time = 2343; /* days */
        Time = Time*86400;
        double fiducial = 1.0; /* ktons */

        _event.efficiency = 0.85;

        _event.Normalization = 8.5*1e-41*1e32;

        _event.Exposure = Time*fiducial;
     }
     else if(Experiment=="Borexino")
     {
         

        string data_file = "exp_data/Borexino_data.dat";
        string bkg_file = "exp_data/Borexino_bkg.dat";

        exp_data.read_file(data_file);
        bkg_data.read_file(bkg_file);

        double Time = 2485; /* days */
        Time = Time*86400;
        double Np = 1.32e31; /* ktons */

        _event.efficiency = 0.85;

        _event.Normalization = 30.0*1e-41;

        _event.Exposure = Time*Np;
     }
     else
     {
         cerr<<"ERROR!!Invalid Experiment name!!"<<endl;
         exit(-1);
     }



    
     _event.res_stat = SOL_NO;


     chisq::SK_IV chi2;
     chi2.fin_flav = 0;
     chi2.channel = atoi(chan);
     chi2.particle = -1;

     

     vec fit_params;

     fit_params.push_back(0.31);
     fit_params.push_back(0.02);
     fit_params.push_back(0);
     fit_params.push_back(0.9);
     fit_params.push_back(7.5e-5);
     fit_params.push_back(2.5e-3);
     fit_params.push_back(1e-4);
     fit_params.push_back(1e-4);

     ofstream ofl;
     ofl.open(file_name);

     vec start_values;

     for(int i=0;i<16;i++)
     {
         start_values.push_back(0.0);
     }

     chisq::sys_minimizer<chisq::SK_IV> _minimizer;
     

     chi2.Set_sys(0.2,0.2,0.05,0.05);
     chi2.Init(_event,exp_data,bkg_data,SOL_NO);
     chi2.Set_binned_systematics(SOL_NO);

     chi2.statistics(start_values);

     double tau = 1e-6;
     double delt = 0.01;


     cout<<"Calculating chi2s...\n";

 //    cout.precision(4);
   /*  for(double tau=-6;tau<=-1;tau=tau+0.1)
     {
         fit_params[6] = pow(10,tau);

         for(double delt = 0.01;delt<=0.99;delt=delt+0.01)
         {
             fit_params[3] = delt;
*/
             double res = 100000000000.0;

             for(double the13=0.02;the13<0.024;the13=the13+0.0002)
             {
                 fit_params[1] = the13;

                 for(double the12=31.0;the12<=36.0;the12=the12+0.2)
                 {
                     fit_params[0] = tan(the12*(M_PI/180.0));

                     for(double ldm=6.8;ldm<=8.0;ldm=ldm+0.04)
                     {
                         fit_params[4] = ldm*1e-5;

                         double chi = _minimizer.Minimize(chi2,start_values,fit_params) + prior(the13,0.02221,0.0006)+prior(the12,33.44,0.78)+prior(ldm,7.42,0.21);

                         res = min(res,chi);

                      cout<<pow(10,tau)<<"\t"<<delt<<"\t"<<the13<<"\t"<<the12<<"\t"<<ldm<<"\t"<<chi<<endl;

                     }


                 }
                 

             }

             cout<<pow(10,tau)<<"\t"<<delt<<"\t"<<res<<endl;
             ofl<<pow(10,tau)<<"\t"<<delt<<"\t"<<res<<endl;
  //       }
  //   }


/*
     for(double tau = -6;tau<=-1;tau=tau+0.1)
     {
         fit_params[6] = pow(10,tau);

         for(double delt = 0.01;delt<=0.99;delt=delt+0.01)
         {
             fit_params[3] = delt;

             double res = 10000000000.0;

             for(double the13=0.0;the13<0.1;the13=the13+0.001)
             {

                 fit_params[1] = the13;

                 double chi = _minimizer.Minimize(chi2,start_values,fit_params);
  //                              + square(the13-0.02221)/square(0.0006);

//                 res = min(res,chi);


  //           }


             cout<<pow(10,tau)<<"\t"<<delt<<"\t"<<chi<<"\n";
             ofl<<pow(10,tau)<<"\t"<<delt<<"\t"<<chi<<"\n";


         }
     }
*/
    
 /* 

    for(double tau = -6;tau<=-1;tau=tau+0.1)
    {
        fit_params[6] = pow(10,tau);

        double res=10000000000.0;

        for(double the13=0.0;the13<0.1;the13=the13+0.001)
        {
            fit_params[1] = the13;

            double chi = _minimizer.Minimize(chi2,start_values,fit_params);
                            //+ square(the13-0.02221)/square(0.0006);

      //      cout<<pow(10,tau)<<"\t"<<the13<<"\t"<<chi<<endl;

      //      ofl<<pow(10,tau)<<"\t"<<the13<<"\t"<<chi<<endl;
      //
              res = min(chi,res);

        }

        cout<<pow(10,tau)<<"\t"<<res<<"\n";
        ofl<<pow(10,tau)<<"\t"<<res<<"\n";
    }

    ofl.close();
*/
/*    for(double th12 = 0.2;th12<=0.7;th12=th12+0.01)
    {
        fit_params[0] = 0.0;

//        double res = 10000000000000.0;

        for(double tau=-7;tau<=-1;tau=tau+0.1)
          {
            fit_params[6] = pow(10,tau);

            double res = 1000000000000000.0;
        
            for(double dm21=2.0;dm21<=15.0;dm21=dm21+0.1)
            {
                fit_params[4] = dm21*1e-5;

                double chi = chi2.Calc_chi2_nosys(fit_params);

  //              cout<<fit_params[0]<<"\t"<<fit_params[6]<<"\t"<<dm21<<"\t"<<chi<<"\n";

                res = min(chi,res);

            }

            cout<<fit_params[0]<<"\t"<<fit_params[6]<<"\t"<<res<<"\n";

        }

//        ofl<<fit_params[6]<<"\t"<<chi<<endl;

 //       cout<<fit_params[0]<<"\t"<<chi<<endl;

    //    cout<<endl;     
    }
    
    ofl.close();
     
    cout<<"Calculation of chi2 complete...\n";
    


*/
/*
    for(double dm21 = 2.0;dm21<=15.0;dm21=dm21+0.01)
    {
        fit_params[4] = dm21*1e-5;

        for(double tau = -5.0;tau<=0.0;tau=tau+0.1)
        {
            fit_params[6] = pow(10,tau);

            double chi = chi2.Calc_chi2_nosys(fit_params);

            ofl<<dm21<<"\t"<<tau<<"\t"<<chi<<endl;
            cout<<dm21<<"\t"<<tau<<"\t"<<chi<<endl;



        }
    }

*/

/*
     for(double delt = 0;delt<=1.0;delt=delt+0.01 )
     {
         fit_params[3] = delt;

         for(double tau = -5.0;tau<=0;tau=tau+0.1)
         {
             fit_params[6] = pow(10,tau);

             double chi = chi2.Calc_chi2_nosys(fit_params);

             ofl<<delt<<"\t"<<tau<<"\t"<<chi<<"\n";

             cout<<delt<<"\t"<<tau<<"\t"<<chi<<"\n";
         }
     }


*/

    return 0;
}


double min(double x,double y)
{
    double m;

    if(x<=y)
    {
        m = x;
    }
    else if(x>y)
    {
        m = y;
    }
        


    return m;
}
