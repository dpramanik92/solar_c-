//
//  KamLAND_anti.hpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 21/11/20.
//

#ifndef event_hpp
#define event_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include "read_files.hpp"
#include "probability.hpp"
#include "numerical.hpp"
#include "visible_anti.hpp"
#include "convers_prob.hpp"


#define SOL_YES 1
#define SOL_NO 0

typedef std::vector<double> vec;

double gauss(double,double,double);


class integration_smear
{
private:
    double Energy,Sigma;
    
public:
    integration_smear(double,double);
    double integrand(double);

};

class integration_conv
{
	private:
	double emin,emax,Sigma,x0;
	Cubic_interpolator Cross,Flux;
	converse Prob;

	public:
	integration_conv(converse,Cubic_interpolator,Cubic_interpolator,double,double,double);
	double integrand(double);
};

class integration_true
{
	private:
	double emin,emax,Sigma,x0;
	Cubic_interpolator Cross;
    visible_anti_prob Prob;
	public:
	integration_true(visible_anti_prob,Cubic_interpolator,double,double);
	double integrand(double);
	
	
	
	
};

class integration_reconstruc
{
	private:
	double emin,emax,Sigma;
	Cubic_interpolator Cross;
    visible_anti_prob Prob;
	
	public:
	integration_reconstruc(visible_anti_prob,Cubic_interpolator,double,double,double);
	double integrand(double);
	
	
	
};

class integration
{

	public:
	
	template <class T, class Real>
	Real simpson_3_8(T,Real,Real);
	template <class T, class Real>
	Real trapezoidal(T,Real,Real);
	template <class T, class Real, class integer>
	Real composite_simpson_3_8(T,Real,Real,integer);

	
};

class Event_generator
{
	private:
		double smearing(double,double);
		//~ int create_smearing_matrix();
		//~ int interpolate_flux();
		//~ int interpolate_cross();
		int init_interpolate_flux();
		int init_interpolate_cross();
		int create_bins();
		int find_sigma(double);
		double Sigma;
		double Prob_E;
		neutrino_data flux;
		neutrino_data cross;
		//~ Probability prob;
		Cubic_interpolator flux_interpolator,cross_interpolator;
		double create_bin_events(int);
        int create_lookup_matrix();
        double sampling_space;
        visible_anti_prob Proba_engine;
        vec samplings;
        double **lookup;
        int FAST_GEN;
        double cross_max,flux_max;
        double norm;
	
	public:
		std::string flux_file,cross_file;
		std::string file_path;
		vec true_osc_params,test_osc_params;
		int smearing_matrix,eff_vector,Man_bins,sys_stat,res_stat;
        int Set_probability_engine(visible_anti_prob);
        double efficiency, resolution[3],e_min,e_max,bin_w,samp_min,samp_max;
        int n_bins,n_samplings;
        double **smear_mat;
        double Normalization,Exposure;
		vec eff;
		vec bin_center,bin_i,bin_f,manual_bins,manual_f;
		vec Events;
		vec systematics;
		int Init_evgen();
        int Init_fast_generator();
        int Set_fast_event_generator(int what,double,double,int);
		int generate_events();
		double Calc_prob(double);
			
	
};


class conv_event
{
	private:
		double smearing(double,double);
		int init_interpolate_flux();
		int init_interpolate_cross();
		int create_bins();
		int find_sigma(double);
		double Sigma;
		double Prob_E;
		neutrino_data flux,cross;
		Cubic_interpolator flux_interpolator,cross_interpolator;
		double create_bin_events(int);
		converse Proba_engine;
		vec samplings;
		double cross_max,flux_max;
		
	public:
		double Exposure,norm,Normalization;
		std::string flux_file,cross_file;
		std::string file_path;
		vec true_osc_params;
		int smearing_matrix,eff_vector,Man_bins,sys_stat,res_stat;
		int Set_probability_engine(converse);
		double efficiency, resolution[3],e_min,e_max,bin_w,samp_min,samp_max;
		int n_bins,n_samplings;
		double **smear_mat;
		vec eff;
		vec bin_center,bin_i,bin_f,manual_bins,manual_f;
		vec Events;
		vec systematics;
		int Init_evgen();
		int generate_events();
		
	
	
};

#endif /* event_hpp */


