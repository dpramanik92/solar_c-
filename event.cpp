//
//  event.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 21/11/20.
//

#include "dec_probability.hpp"
#include "event.hpp"
#include "read_files.hpp"
#include "numerical.hpp"
#include <fstream>
#include <cmath>



int Event_generator::Init_evgen(int flav,int parti, int channel)
{

    _flavour = flav;
    _channel = channel;
    _particle = parti;
	
	smearing_matrix=SOL_NO;
	eff_vector=SOL_NO;
	//Man_bins=SOL_NO;
	sys_stat=SOL_YES;
//	res_stat=SOL_YES;
	FAST_GEN = SOL_NO;

	norm = Normalization*Exposure;


    create_bins();

    //~ find_sigma();
    //~ create_smearing_matrix();
    
    init_interpolate_flux();
    init_interpolate_cross();
    


	
	
	
	std::cout<<"The Event generator is initialized...\n";
	
	return 0;
	
}

int Event_generator::Init_fast_generator()
{
    sampling_space = (samp_max-samp_min)/double(n_samplings);
    
    for(int i = 0;i<n_samplings;i++)
    {
        double s=samp_min+(0.5+i)*sampling_space;;
        samplings.push_back(s);
                
    }
    
    create_lookup_matrix();

    return 0;
    
}

int Event_generator::Set_fast_event_generator(int what,double _min,double _max,int _num)
{
    FAST_GEN = what;
    res_stat = what;
    
    
    samp_min = _min;
    samp_max = _max;
    n_samplings = _num;
    
    
    
    
    
    return 0;
}

double Event_generator::resFunc(double E,double E_p)
{
    find_sigma(E);

    double res = 1.0/(sqrt(2.0*M_PI)*Sigma)*exp(-0.5*((E-E_p)/Sigma)*((E-E_p)/Sigma));
    return res;


}

int Event_generator::create_lookup_matrix()
{
    
    smear_mat = new double*[n_bins];
    
    for(int i=0;i<n_bins;i++)
    {
        smear_mat[i] = new double[n_samplings];
        
        
    }
    
    for(int i=0;i<n_bins;i++)
    {
        double h = (bin_f[i]-bin_i[i]);
        double a = bin_i[i];
        double b = bin_f[i];
        
        for(int j=0;j<n_samplings;j++)
        {
            find_sigma(samplings[j]);
            
            
				smear_mat[i][j] = h/3.0*(resFunc(samplings[j],a)+3.0*resFunc(samplings[j],(a+2.0*b)/3.0)+3.0*resFunc(samplings[j],(2.0*a+b)/3.0)+resFunc(samplings[j],b));
        }
        
    }
    
    
    
    return 0;
}

int Event_generator::Set_probability_engine(dec_prob Prob)
{
    Proba_engine = Prob;
    
    return 0;
}


int Event_generator::find_sigma(double E)
{
	if(smearing_matrix!=SOL_YES)
	{
		double alpha = resolution[0];
		double beta = resolution[1];
		double gamma = resolution[2];
	
		
		Sigma = alpha*E+beta*sqrt(E)+ gamma ;

	}
	
    
	return 0;
	
}


int Event_generator::generate_events()
{
    
    Events.clear();
    

	if(res_stat==SOL_YES)
	{

		if(FAST_GEN!=SOL_YES)
		{
	
			for(int i=0;i<n_bins;i++)
			{
				if((log10((bin_center[i]+1.0)*1e-3)<=cross_max) || (bin_center[i]<=flux_max))
				{
			
		
					double Ei = bin_i[i]+1.0;
					double Ei_del = bin_f[i]+1.0;
		
					find_sigma(bin_center[i]+1.0);
		
					integration_reconstruc _inte_recons(Proba_engine,cross_interpolator,e_min,e_max,Sigma);

                    _inte_recons._channel = _channel;
                    _inte_recons._flavour = _flavour;
                    _inte_recons._particle = _particle;
		
					integration _inte;
		
					double ev = _inte.simpson_3_8<integration_reconstruc,double>(_inte_recons,Ei,Ei_del);
		
					Events.push_back(norm*ev*efficiency);
				}
				else
				{
					Events.push_back(0);
				}
		
		
			}
		}
		if(FAST_GEN == SOL_YES)
		{
			for(int i =0;i<n_bins;i++)
			{
				double sum = 0;
            
				for(int j=0;j<n_samplings;j++)
				{
					double term = 0;
//					if(smear_mat[i][j]>1e-5)
//					{
                        if((log10(samplings[j])*1e-3)<=cross_max || (samplings[j]<=flux_max))
                        {
    			    		term = Proba_engine.Calculate_decayed_flux((samplings[j]+1.3),_flavour,_particle,_channel)*cross_interpolator.interpolate(log10(samplings[j])*1e-3)*(samplings[j])*smear_mat[i][j]*sampling_space;


    			         }
                        else
                        {
                            term = 0;
                        }
//					}
					sum  = sum + term;

				}
            
				Events.push_back(norm*sum*efficiency);
			}
        
        
		}
	}
	if(res_stat!=SOL_YES)
	{
		double ev;
		
		for(int i =0;i<n_bins;i++)
		{
		
			if((log10(bin_center[i]*1e-3)<=cross_max) || (bin_center[i]<=flux_max))
			{
		
			

				double ev1 = Proba_engine.Calculate_decayed_flux(bin_i[i]+1.3,_flavour,_particle,_channel)*cross_interpolator.interpolate(log10((bin_i[i]+1.3)*1e-3))*(bin_i[i]+1.3);
				
				double ev2 = Proba_engine.Calculate_decayed_flux(bin_f[i]+1.3,_flavour,_particle,_channel)*cross_interpolator.interpolate(log10((bin_f[i]+1.3)*1e-3))*(bin_f[i]+1.3);

                

				ev = 0.5*(bin_f[i]-bin_i[i])*(ev1+ev2);

			}
			else
			{
				ev = 0;
			}
	
//            std::cout<<norm*ev*efficiency<<"\n"; 
		
			Events.push_back(norm*ev*efficiency);
		}
		
		
	}
	
	return 0;
	
}



double Event_generator::Calc_prob(double E)
{
	/*
	 *  calculation of probability
	 * 
	 * 
	 * 
	 * 
	 * */
	
	
	
	
	
	return 1.0;
}




int Event_generator::create_bins()
{
	if(Man_bins!=SOL_YES)
	{
		bin_w = (e_max-e_min)/n_bins;
		


		for(int i=0;i<n_bins;i++)
		{
			bin_i.push_back(e_min+i*bin_w);
			bin_f.push_back(e_min+(i+1)*bin_w);
			
			bin_center.push_back(e_min+(0.5+i)*bin_w);
			
			
		}
		
	}
	if(Man_bins == SOL_YES)
	{
		for(int i=0;i<manual_bins.size()-1;i++)
		{
			
			bin_i.push_back(manual_bins[i]);
			bin_f.push_back(manual_f[i]);
			bin_center.push_back((manual_bins[i]+manual_f[i])/2.0);
			

			
		}
		
		n_bins = manual_bins.size()-1;
		
		e_min = manual_bins[0];
		e_max = manual_f[manual_bins.size()-2];
		
	}
	
	return 0;
}


int Event_generator::init_interpolate_flux()
{
	flux.read_flux("flux/b8spectrum.txt");
	
	vec flux_x,flux_y;
	
	for(int i =0;i<flux.Flux.data[0].size();i++)
	{
		flux_x.push_back(flux.Flux.data[0][i]);
		flux_y.push_back(flux.Flux.data[1][i]);
		
	}
	
	flux_max = flux_x[flux_x.size()-2];

	
	flux_interpolator.set_cubic_spline(flux_x,flux_y,false);
	
    return 0;	
	
}



int Event_generator::init_interpolate_cross()
{
	cross.read_cross("cross_section/XCCreactor.dat");
	
	
	
	
	vec cross_x,cross_y;
	
	for(int i =0;i<cross.Cross.data[0].size();i++)
	{
		cross_x.push_back(cross.Cross.data[0][i]);
		cross_y.push_back(cross.Cross.data[4][i]);
		
		
	}
	
	cross_max = cross_x[cross_x.size()-2];

	
	cross_interpolator.set_cubic_spline(cross_x,cross_y,false);

    return 0;
	
}




/*********************************************************************************************************************
 * 
 * 
 * 
 * ******************************************************************************************************************/

int conv_event::Init_evgen()
{
	smearing_matrix=SOL_NO;
	eff_vector=SOL_NO;
//	Man_bins=SOL_NO;
	sys_stat=SOL_YES;
	res_stat=SOL_YES;

	norm = Normalization*Exposure;
	
    create_bins();
    //~ find_sigma();
    //~ create_smearing_matrix();
    
    init_interpolate_flux();
    init_interpolate_cross();
    
	
	return 0;
}

int conv_event::init_interpolate_flux()
{
	flux.read_flux("flux/b8spectrum.txt");
	
	vec flux_x,flux_y;
	
	for(int i =0;i<flux.Flux.data[0].size();i++)
	{
		flux_x.push_back(flux.Flux.data[0][i]);
		flux_y.push_back(flux.Flux.data[1][i]);
		
	}
	
	flux_interpolator.set_cubic_spline(flux_x,flux_y,false);
	
	flux_max = flux_x[flux_x.size()-2];

	
	
	
	return 0;
}

int conv_event::init_interpolate_cross()
{
	
	cross.read_cross("cross_section/XCCreactor.dat");
	
	
	
	
	vec cross_x,cross_y;
	
	for(int i =0;i<cross.Cross.data[0].size();i++)
	{
		cross_x.push_back(cross.Cross.data[0][i]);
		cross_y.push_back(cross.Cross.data[4][i]);
		
		
	}
	
	cross_max = cross_x[cross_x.size()-2];


	
	cross_interpolator.set_cubic_spline(cross_x,cross_y,false);
	
	
	return 0;
}


int conv_event::create_bins()
{
	if(Man_bins!=SOL_YES)
	{
		bin_w = (e_max-e_min)/n_bins;
		
		for(int i=0;i<n_bins;i++)
		{
			bin_i.push_back(e_min+i*bin_w);
			bin_f.push_back(e_min+(i+1)*bin_w);
			
			bin_center.push_back(e_min+(0.5+i)*bin_w);
			
			
		}
		
	}
	if(Man_bins == SOL_YES)
	{
        
		for(int i=0;i<manual_bins.size()-1;i++)
		{
			
			bin_i.push_back(manual_bins[i]);
			bin_f.push_back(manual_f[i]);
			bin_center.push_back((manual_bins[i]+manual_f[i])/2.0);
			

			
		}
		
		n_bins = manual_bins.size();
		
		e_min = manual_bins[0];
		e_max = manual_f[manual_bins.size()-2];

	}
	
	
	return 0;
}

int conv_event::Set_probability_engine(converse Prob)
{
    Proba_engine = Prob;

	return 0;
}

int conv_event::find_sigma(double E)
{
	if(smearing_matrix!=SOL_YES)
	{
		double alpha = resolution[0];
		double beta = resolution[1];
		double gamma = resolution[2];
	
		
		Sigma = alpha*E+beta*sqrt(E)+ gamma ;

	}
	
    
	
	
	return 0;
}

int conv_event::generate_events()
{
	
	
	 for(int i=0;i<n_bins-1;i++)
     {
		
        double Ei = bin_i[i]+1.0;
        double Ei_del = bin_f[i]+1.0;
		
		
		if(res_stat==SOL_YES)
		{
			
			if((log10((bin_center[i]+1.0)*1e-3)<=cross_max) || ((bin_center[i]+1.0)<=flux_max))
			{
				find_sigma(bin_center[i]+1.0);
		
				integration_conv _inte_conv(Proba_engine,flux_interpolator,cross_interpolator,Ei,Ei_del,Sigma);
		
				integration _inte;
		
				double ev = _inte.composite_simpson_3_8<integration_conv,double,int>(_inte_conv,e_min,e_max,50);
			
				Events.push_back(norm*ev*efficiency);
			}
			else
			{
				Events.push_back(0);
			}

        }
        if(res_stat==SOL_NO)
        {
			if((log10((bin_center[i]+0.7)*1e-3)<=cross_max) && (bin_center[i]<=flux_max))
			{
			
			
				//double ev1 = 8.5e-4*Proba_engine.Calculate_probability(bin_i[i]+1.3)*flux_interpolator.interpolate(bin_i[i]+1.3)*cross_interpolator.interpolate(log10((bin_i[i]+1.3)*1e-3))*(bin_i[i]+1.3);
				
				double ev1 = Proba_engine.Calculate_probability(bin_i[i]+0.7)*flux_interpolator.interpolate(bin_i[i]+0.7)*cross_interpolator.interpolate(log10((bin_i[i]+1.3)*1e-3))*(bin_i[i]+0.7);
				
                double ev2 = Proba_engine.Calculate_probability(bin_f[i]+0.7)*flux_interpolator.interpolate(bin_f[i]+0.7)*cross_interpolator.interpolate(log10((bin_f[i]+0.7)*1e-3))*(bin_f[i]+0.7);

				double ev = 0.5*(bin_f[i]-bin_i[i])*(ev1+ev2);
		

				Events.push_back(norm*ev*efficiency);
			}
			else
			{
				Events.push_back(0);
			}
				
		}

		
		
		
		
     }
		
	
	
	return 0;
}



/*********************************************************************************************************************
 * 
 * 
 * 
 * ******************************************************************************************************************/

integration_true::integration_true(dec_prob _Prob,Cubic_interpolator _cross,double e_p,double _sigma,int flav,int parti,int channel)
{
    
    _flavour = flav;
    _particle = parti;
    _channel = channel;

    
	Prob = _Prob;
	Cross = _cross;
	
	x0 = e_p;
	
	
	
	Sigma = _sigma;
	
}

double integration_true::integrand(double x)
{
    double res = Prob.Calculate_decayed_flux(x,_flavour,_particle,_channel)*Cross.interpolate(log10(x*1e-3))*x*gauss(x,x0,Sigma);
    
    
    //~ std::cout<<Prob.Calculate_decayed_flux(x)<<"\t"<<gauss(x,x0,Sigma)<<"\n";
	
	
	
	return res;
}

/*********************************************************************************************************************
 * 
 * 
 * 
 * ******************************************************************************************************************/

integration_reconstruc::integration_reconstruc(dec_prob _Prob, Cubic_interpolator _cross, double _emin,double _emax,double _sigma)
{
	
	Prob = _Prob;
	Cross = _cross;
	
	Sigma = _sigma;
	
	emin = _emin;
	emax = _emax;
	
	
}

double integration_reconstruc::integrand(double x)
{
	
	
	integration_true _inte_true(Prob,Cross,x,Sigma,_flavour,_particle,_channel);
	
//	emin = x-3.0*Sigma;
//	emax = x+3.0*Sigma;
	
    int numbers = int((emax-emin)/Sigma);
    
    if(numbers>50)
    {
        numbers=50;
    }
    
	integration _inte;
	double res = _inte.composite_simpson_3_8<integration_true,double,int>(_inte_true,emin,emax,100);
	return res;
}

/*********************************************************************************************************************
 * 
 * 
 * 
 * ******************************************************************************************************************/

integration_conv::integration_conv(converse _Prob,Cubic_interpolator _flux, Cubic_interpolator _cross, double _emin,double _emax,double _sigma)
{
	
	Prob = _Prob;
	Cross = _cross;
	Flux = _flux;
	
	Sigma = _sigma;
	
	emin = _emin;
	emax = _emax;
	
	
}

double integration_conv::integrand(double x)
{
	
	
	integration_smear _inte_smear(x,Sigma);
	
	integration _inte;
	double K_e = _inte.composite_simpson_3_8<integration_smear,double,int>(_inte_smear,emin,emax,51);
	//~ std::cout<<x<<"\n";
	double P = Prob.Calculate_probability(x);


	double res = P*Flux.interpolate(x)*Cross.interpolate(log10(x*1e-3))*x*K_e;
	return res;
}




/*********************************************************************************************************************
 *
 *
 *
 * ******************************************************************************************************************/

integration_smear::integration_smear(double _E, double _Sigma)
{
    Energy = _E;
    Sigma = _Sigma;
}

double integration_smear::integrand(double E_p)
{
    
    
    double res = gauss(Energy,E_p,Sigma);
    return res;
    
}

double gauss(double x,double x0,double sigma)
{
        double f = 1/(sqrt(2*M_PI)*sigma)*exp(-0.5*square((x-x0)/sigma));
        //~ std::cout<<x<<"\t"<<x0<<"\n";
		return f;
}


template <class T, class Real>
Real integration::simpson_3_8(T obj,Real a,Real b)
{
	
	Real h = (b-a)/3.0;
	
	
	Real f = 3.0*h/8.0*(obj.integrand(a)+3.0*obj.integrand((2.0*a+b)/3.0)+3.0*obj.integrand((a+2.0*b)/3.0)+obj.integrand(b));
	
	return f;
	
	
}


template <class T, class Real>
Real integration::trapezoidal(T obj,Real a,Real b)
{
	
	Real h = (b-a)/2.0;
	
	
	Real f = h*(obj.integrand(a)+obj.integrand(b));
	
	return f;
	
	
}

template <class T, class Real, class integer>
Real integration::composite_simpson_3_8(T obj,Real a,Real b, integer n)
{
	
	Real h = (b-a)/n;
	
	integer N = n-1;
	integer M = n/3-1;
	
	
	
	Real f1 = obj.integrand(a);
	Real f4 = obj.integrand(b);
	
	Real sum1 = 0;
	Real sum2 = 0;

	for(int i=1;i<=N;i++)
	{
		if(i%3!=0)
		{
			Real term = obj.integrand(a+i*h);
			
			sum1 = sum1 + term;
		}
		
	}
	
	for(int i=1;i<=M;i++)
	{
		Real term = obj.integrand(a+3.0*i*h);
		
		sum2 = sum2 + term;
	}
	
	double res = 3.0*h/8.0*(f1 + 3.0*sum1 + 2.0*sum2 + f4);

	return res;
}

