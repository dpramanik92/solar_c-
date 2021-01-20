//
//  event.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 21/11/20.
//

#include "probability.hpp"
#include "KamLAND_anti.hpp"
#include "read_files.hpp"
#include "numerical.hpp"
#include <fstream>
#include <cmath>



int Event_generator::Init_evgen()
{
		
	
	smearing_matrix=SOL_NO;
	eff_vector=SOL_NO;
	Man_bins=SOL_NO;
	sys_stat=SOL_YES;
	res_stat=SOL_YES;

	
	create_bins();
	//~ find_sigma();
	//~ create_smearing_matrix();
	
	init_interpolate_flux();
	init_interpolate_cross();
	
	
	
	std::cout<<"The Event generator is initialized...\n";
	
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


	
	for(int i=0;i<n_bins;i++)
	{
		
		double Ei = bin_i[i];
		double Ei_del = bin_f[i];
		
		find_sigma(bin_center[i]);
		
		integration_reconstruc _inte_recons(flux_interpolator,cross_interpolator,e_min,e_max,Sigma);
		
		integration _inte;
		
		double ev = _inte.simpson_3_8<integration_reconstruc,double>(_inte_recons,Ei,Ei_del);
		
		Events.push_back(ev);
		
		
		
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
			bin_f.push_back(manual_bins[i+1]);
			bin_center.push_back((manual_bins[i]+manual_bins[i+1])/2.0);
			
			
			
		}
		
		n_bins = manual_bins.size();
		
		
	}
	
	return 0;
}


int Event_generator::init_interpolate_flux()
{
	flux.read_flux("flux/b8spec-2006.dat");
	
	vec flux_x,flux_y;
	
	for(int i =0;i<flux.Flux.data[0].size();i++)
	{
		flux_x.push_back(flux.Flux.data[0][i]);
		flux_y.push_back(flux.Flux.data[1][i]);
		
	}
	
	flux_interpolator.set_cubic_spline(flux_x,flux_y,false);
	
	
	
}



int Event_generator::init_interpolate_cross()
{
	cross.read_cross("cross_section/reactor.dat");
	
	
	
	
	vec cross_x,cross_y;
	
	for(int i =0;i<cross.Cross.data[0].size();i++)
	{
		cross_x.push_back(cross.Cross.data[0][i]);
		cross_y.push_back(cross.Cross.data[1][i]);
		
		
	}
	
	
	
	cross_interpolator.set_cubic_spline(cross_x,cross_y,false);
	
}




/*********************************************************************************************************************
 * 
 * 
 * 
 * ******************************************************************************************************************/

integration_true::integration_true(Cubic_interpolator _flux,Cubic_interpolator _cross,double e_p,double _sigma)
{
	Flux = _flux;
	Cross = _cross;
	
	x0 = e_p;
	
	
	
	Sigma = _sigma;
	
}

double integration_true::integrand(double x)
{
	double res = Flux.interpolate(x)*Cross.interpolate(x)*gauss(x,x0,Sigma);
	
	
	
	return res;
}

/*********************************************************************************************************************
 * 
 * 
 * 
 * ******************************************************************************************************************/

integration_reconstruc::integration_reconstruc(Cubic_interpolator _flux, Cubic_interpolator _cross, double _emin,double _emax,double _sigma)
{
	
	Flux = _flux;
	Cross = _cross;
	
	Sigma = _sigma;
	
	emin = _emin;
	emax = _emax;
	
	
}

double integration_reconstruc::integrand(double x)
{
	
	
	integration_true _inte_true(Flux,Cross,x,Sigma);
	
	emin = x-3.0*Sigma;
	emax = x+3.0*Sigma;
	
	integration _inte;
	double res = _inte.composite_simpson_3_8<integration_true,double,int>(_inte_true,emin,emax,9);
	return res;
}




double square(double x)
{
	double s = x*x;
	return s;
}


double gauss(double x,double x0,double sigma)
{
        double f = 1/(sqrt(2*M_PI)*sigma)*exp(-0.5*square((x-x0)/sigma));
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

