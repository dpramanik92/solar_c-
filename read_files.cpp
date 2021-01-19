//
//  read_files.cpp
//  solar_neutrinos
//
//  Created by Dipyaman Pramanik on 20/11/20.
//

#include "read_files.hpp"

/*xxxxxxx**********************************************************************************
*  MEMBER FUNCTIONS OF THE CLASS FILE READER
*
******************************************************************************************* */

/*
 Constructor of the class file reader
 */
/* ***********************************************************************************************************************/


file_reader::file_reader()
{
    
}

file_reader::~file_reader()
{

}

int file_reader::clean_data()
{
    delete[] data;
    return 0;
}


/*
  This function counts the number of columns in a data file. And stores
  the number to the private variable col_count.
 */

/* ***********************************************************************************************************************/


int file_reader::count_column(std::string file_name)
{

    std::ifstream ifile;
    ifile.open(file_name);
    std::string line;
    
    std::getline(ifile,line);
    std::stringstream ss;
    ss<<line;
    
    int count = 0;
    double value;
    while(ss>>value) count++;
    
    file_reader::col_count=count;
    if(count==0)
    {
        std::cout<<"ERROR!! Number of columns is zero.\n";
        exit(-1);
        return -1;
        
    }
    
    return 0;
}

/* This function reads the data file and stores the information in the vector col. */

/* ***********************************************************************************************************************/


int file_reader::read_file(std::string file_name)
{
    file_reader::count_column(file_name);
    n_col = col_count;
    if(file_reader::col_count==0)
    {
        std::cout<<"ERROR!!Invalid file. zero columns encountered\n";
        exit(-1);
        return -1;
    }
  
    data = new vec[col_count];
  
    double *col;
    col=new double[col_count];
    
    std::ifstream ifile;
    ifile.open(file_name);

    while(ifile)
    {

        
        for(int i=0;i<col_count;i++)
        {
            
            ifile>>col[i];
            data[i].push_back(col[i]);

        
        }
        

    }

    n_row = int(data[0].size());
    
    col=NULL;
    delete[] col;
    
    ifile.close();
    
    for(int i =0;i<col_count;i++)
    {
        
        data[i][data[0].size()-1]=0;
    }
    
    if(data==NULL)
    {
        std::cout<<"Error!! File reading failed\n";
        exit(-1);
    }
    

    return 0;
}

/*xxxxxxx**********************************************************************************
*  MEMBER FUNCTIONS OF THE CLASS NEUTRINO DATA
*
******************************************************************************************* */



/*
  Constructor of Class neutrino data
 
 */

neutrino_data::neutrino_data()
{
    
    
}

neutrino_data::~neutrino_data()
{
    
    
}



/*
 This function reads the neutrino flux from a file
 and stores into a vector array
 
 */

/* ***********************************************************************************************************************/

int neutrino_data::read_flux(std::string  file_name)
{
    std::cout<<"Reading flux\n"; // of "<<neutrino_data::exp<<std::endl;
    Flux.read_file(file_name);
    
    neutrino_data::n_row = int(Flux.data[0].size());
    neutrino_data::n_col = Flux.col_count;
    
    return 0;
}

/*
 This function prints the flux stored in the vector array.
 */

/* ***********************************************************************************************************************/

int neutrino_data::print_flux()
{
    
    for(int i=0;i<n_row;i++)
    {
        for(int j=0;j<n_col;j++)
        {
            std::cout<<Flux.data[j][i];
            
        }
        std::cout<<"\n";
        
    }
  
    
    return 0;
}

/*
 This function reads the cross section fromt the file and stores
 in a vector array.
 
 */
/* ***********************************************************************************************************************/


int neutrino_data::read_cross(std::string file_name)
{
    std::cout<<"Reading cross section of\n";
    Cross.read_file(file_name);
    
    neutrino_data::n_row = int(Cross.data[0].size());
    neutrino_data::n_col = Cross.col_count;
    
    return 0;
}

/*
 This function prints the cross section stored in the vector array.
 
 */

/* ***********************************************************************************************************************/

int neutrino_data::print_cross()
{
    
    for(int i=0;i<n_row;i++)
    {
        for(int j=0;j<n_col;j++)
        {
            std::cout<<Cross.data[j][i]<<"\t";
            
        }
        std::cout<<"\n";
        
    }
  
    
    
    return 0;
}

/***********************************************************************************************\
    MEMBER FUNCTIONS OF THE CLASS SOLAR_DENSITY
 
 *********************************************************************************************/

int solar_density::clean_memory()
{
    density.clean_data();
    
    return 0;
}


/*
 
 This function reads the density of the sun and stores in a
 vector array
 
 */
/* ***********************************************************************************************************************/

int solar_density::read_density(std::string file_name)
{
 
    density.read_file(file_name);

    n_row = int(density.data[0].size());
    n_col = density.col_count;
    
    return 0;
    
}

/*
 This function prints the data stored in the array
 
 */
/* ***********************************************************************************************************************/

int solar_density::print_density()
{
    for(int i=0;i<n_row;i++)
    {
        for(int j=0;j<n_col;j++)
        {
            std::cout<<density.data[j][i]<<"\t";
            
        }
        std::cout<<"\n";
        
    }
  
    
    
    return 0;
    
}

/***********************************************************************************************\
    MEMBER FUNCTIONS OF THE CLASS PRODUCTION
 
 *********************************************************************************************/

int production::clean_memory()
{
    prod.clean_data();
    
    return 0;
}

production::~production()
{
    
}

/*
 
 This function reads the Production distribution in the sun and
 stores in a vector array
 
 */

/* ***********************************************************************************************************************/

int production::read_prod_dist(std::string file_name)
{
 
    prod.read_file(file_name);

    n_row = int(prod.data[0].size());
    n_col = prod.col_count;
    
    return 0;
    
}

/*
 This function prints the data stored in the array
 
 */
/* ***********************************************************************************************************************/

int production::print_prod_dist()
{
    for(int i=0;i<n_row;i++)
    {
        for(int j=0;j<n_col;j++)
        {
            std::cout<<prod.data[j][i]<<"\t";
            
        }
        std::cout<<"\n";
        
    }
  
    
    
    return 0;
    
}


/***********************************************************************************************\
    MEMBER FUNCTIONS OF THE CLASS PRODUCTION
 
 *********************************************************************************************/

/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/

int regeneration::clean_memory()
{
    regen.clean_data();
    
    return 0;
}

regeneration::~regeneration()
{
    
    
}

/* ***********************************************************************************************************************/


/* ***********************************************************************************************************************/


int regeneration::read_regeneration(std::string file_name)
{
 
    regen.read_file(file_name);

    n_row = int(regen.data[0].size());
    n_col = regen.col_count;
    
    return 0;
    
}

/* ***********************************************************************************************************************/

int regeneration::print_regeneration()
{
   for(int i=0;i<n_row;i++)
   {
       for(int j=0;j<n_col;j++)
       {
           std::cout<<regen.data[j][i]<<"\t";
           
       }
       std::cout<<"\n";
       
   }
 
   
   
   return 0;
   
}
