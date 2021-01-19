//
//  read_files.hpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 20/11/20.
//

#ifndef read_files_hpp
#define read_files_hpp

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

typedef std::vector<double> vec;

class file_reader
{
private:
    int count_column(std::string);
    
    
public:
    int n_row;
    int n_col;
    vec *data;
    int col_count;
    file_reader();
    ~file_reader();
    int clean_data();
    int read_file(std::string);
    
    
};

class neutrino_data
{
public:
    std::string objname;
    file_reader Flux,Cross;
    std::string exp;
    int n_row,n_col;
    neutrino_data();
    ~neutrino_data();
    int read_flux(std::string);
    int read_cross(std::string);
    int print_flux();
    int print_cross();
    
    
    
};

class solar_density
{
private:
    
public:
    std::string objname;
    file_reader density;
    int n_row,n_col;
    int read_density(std::string);
    int clean_memory();
    int print_density();
    
};

class production
{
private:
    
public:
    ~production();
    file_reader prod;
    int n_row,n_col;
    int clean_memory();
    int read_prod_dist(std::string);
    int print_prod_dist();
    
};


class regeneration
{
public:
    ~regeneration();
    file_reader regen;
    int n_row,n_col;
    int clean_memory();
    int read_regeneration(std::string);
    int print_regeneration();
    
};

#endif /* read_files_hpp */
