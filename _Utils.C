#include <TGraphErrors.h>
#include <TLatex.h>
#include <TString.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <dirent.h>
#include <string.h>



//____________________________________________
void PrintTime(time_t t0, time_t t1)
{
    double tDay = int((t1-t0)/(24*60*60.));
    double tH = int((t1-t0)/(60*60.) - tDay*24);
    double tMin = int((t1-t0)/(60.) - tDay*24*60 - tH*60);
    double tSec = (t1-t0) - tDay*24*60*60 - tH*60*60 - tMin*60;
    std::cout << "Duration of the simulation = " << tDay << " days " << tH << "h " << tMin << " min " << tSec << "s." << std::endl;
    
    std::cout << std::endl;
}


int GetNumberOfFiles(TString path, TString name) {
    Int_t num = 0;
    struct dirent **namelist;
    Int_t n = scandir(path, &namelist, 0, alphasort);
    if (n < 1) {std::cout << "empty folder" << std::endl;}
    else { while (n--) { if (strstr(namelist[n]->d_name, name) != NULL) num++;} }
	return num;
}

int GetNumberOfFiles(TString path) {
	TString name = "";
	return GetNumberOfFiles(path, name);
}

//____________________________________________
void PrintVector(
                 TString name,
                 Double_t* values,
                 Int_t size,
                 TString format = "%f" )
{
    //std::cout << "const Double_t " << name << "["<<size<< "] = {";
    std::cout << "const std::vector<Double_t> " << name << " = {";
    for( Int_t i=0; i<size; ++i )
    {
        std::cout << Form( format.Data(), values[i] );
        if( i != size-1 ) std::cout << ", ";
        else std::cout << "};";
    }
    std::cout << std::endl;
}

//____________________________________________
void PrintVector(
    TString name,
    std::vector<Double_t> values,
    TString format = "%f" )
{ PrintVector( name, &values[0], values.size(), format ); }


//____________________________________________
void PrintList(
                 TString name,
                 Double_t* list,
                 Int_t num,
                 TString format = "%f" )
{
    //std::cout << "const Double_t " << name << "["<<size<< "] = {";
    std::cout << "Double_t " << name << "[" << num << "]" << " = {";
    for( Int_t i=0; i<num; ++i )
    {
        std::cout << Form( format.Data(), list[i] );
        if( i != num-1 ) std::cout << ", ";
        else std::cout << "};";
    }
    std::cout << std::endl;
}
