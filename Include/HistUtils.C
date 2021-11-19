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
#include <fstream>
#include <cmath>
#include <dirent.h>
#include <iostream>     // std::cout
#include <sstream>      // std::istringstream
#include <string>       // std::string

using namespace std;

//____________________________________________
Int_t GetMaximumBin( TH1*h, Int_t firstBin = 1, Int_t lastBin = 0 )
{
    if( lastBin <= 0 ) lastBin = h->GetNbinsX();
    
    Int_t iBinMax = firstBin;
    Double_t maxValue = h->GetBinContent( firstBin );
    
    for( Int_t iBin = firstBin+1; iBin < lastBin; ++iBin )
    {
        const Double_t value = h->GetBinContent( iBin );
        if( value > maxValue )
        {
            maxValue = value;
            iBinMax = iBin;
        }
    }
    
    return iBinMax;
}

Int_t GetMaximumPos( TH1*h, Double_t firstPos = 1, Double_t lastPos = 0 )
{
    if( lastPos <= 0 ) lastPos = h->GetXaxis()->GetBinCenter( h->GetNbinsX() );
    
    Int_t firstBin = h->GetXaxis()->FindBin(firstPos);
    Int_t lastBin = h->GetXaxis()->FindBin(lastPos);
    
    return GetMaximumBin(h, firstBin, lastBin);
}

Int_t GetMinimumBin( TH1*h, Int_t firstBin = 1, Int_t lastBin = 0 )
{
    if( lastBin <= 0 ) lastBin = h->GetNbinsX();
    
    Int_t iBinMin = firstBin;
    Double_t minValue = h->GetBinContent( firstBin );
    
    for( Int_t iBin = firstBin+1; iBin < lastBin; ++iBin )
    {
        const Double_t value = h->GetBinContent( iBin );
        if( value < minValue )
        {
            minValue = value;
            iBinMin = iBin;
        }
    }
    
    return iBinMin;
}

//____________________________________________
Double_t GetMaximum( TH1*h, Int_t firstBin = 1, Int_t lastBin = 0 )
{
    if( lastBin <= 0 ) lastBin = h->GetNbinsX();
    
    Int_t iBinMax = firstBin;
    Double_t maxValue = h->GetBinContent( firstBin );
    
    for( Int_t iBin = firstBin+1; iBin < lastBin; ++iBin )
    {
        const Double_t value = h->GetBinContent( iBin );
        if( value > maxValue )
        {
            maxValue = value;
        }
    }
    return maxValue;
}


//____________________________________________
TH1* ReadData( const TString file)
{
    // histogram
    TH1* h = new TH1F( "mca", "", 2048, 0, 2048. );
    //TH1* h = new TH1F( "mca", "mca", 1000, 0, 50 );
    //TH1* h = new TH1F( "mca", "mca", 100, 0, 5 );
    
    std::string line;
    ifstream in( file.Data() );
    if (in) {
        Int_t iBin = 0;
        while( getline( in, line ) )
        {
            istringstream stream( line );
            Int_t value;
            stream >> value;
            if( stream )
            {
                ++iBin;
                //h->SetBinContent( iBin/coarseGain, value );
                h->Fill( iBin, value );
            }
            
        }
        /*
        h2->SetXTitle("MCA channels / Coarse Gain");
        for (int i = 1; i < h->GetNbinsX()+1; i++) {
            h2->SetBinContent(i, h->GetBinContent(i));
        }
         */
        
        //Double_t histMax = h->GetMaximum();
        while (h->GetMaximum()<100) h->Rebin(2);
        //h->Scale(1/histMax);
        h->SetOption("hist");
        h->SetXTitle("MCA channels");
        h->SetYTitle("normalised # counts");
    }
    else {cout << "\n\nImpossible d'ouvrir le fichier MCA: " << file << endl << endl << endl;}
    
    //return h;
    return h;
}

//____________________________________________
TH1* ReadData( const TString file, Int_t coarseGain)
{
    // histogram
    TH1* h = ReadData(file);
    //h->GetXaxis()->SetLimits(0, 1024./coarseGain );
    h->GetXaxis()->SetLimits(0, 2048./coarseGain );
    h->SetXTitle("MCA channels / Coarse Gain");
    //h->GetXaxis()->SetRangeUser(0,50);
    return h;
}


//____________________________________________
Double_t Square( Double_t x ) { return x*x; }

//____________________________________________
Double_t FitFunctionExp( Double_t* x, Double_t* par )
{ return TMath::Exp( par[0] + par[1]*x[0] ); }

//____________________________________________
Double_t FitFunctionLin( Double_t* x, Double_t* par )
{ return  par[0] + par[1]*x[0] ; }

//____________________________________________
Double_t FitFunctionConst( Double_t* par )
{ return ( par[0]); }


Double_t FitGauss( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
    return  par[2]*TMath::Gaus( x[0], par[0], par[1]); }

Double_t Fit3Gauss( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
    return  par[2]*TMath::Gaus( x[0], par[0], par[1]) + par[5]*TMath::Gaus( x[0], par[3], par[4]) + par[8]*TMath::Gaus( x[0], par[6], par[7]) ;}

Double_t Fit3GaussAndExp( Double_t* x, Double_t* par ) { //(Double_t x, Double_t mean = 0, Double_t sigma = 1, Bool_t norm = kFALSE)
    return  par[2]*TMath::Gaus( x[0], par[0], par[1]) + par[5]*TMath::Gaus( x[0], par[3], par[4]) + par[8]*TMath::Gaus( x[0], par[6], par[7]) + TMath::Exp( par[9] + par[10]*x[0] ); }

Double_t FitPoisson( Double_t* x, Double_t* par ) {
    return par[1]*par[2]*TMath::Poisson(x[0],par[0]);
}

//____________________________________________
Double_t CrystalBall2( Double_t x, Double_t mean, Double_t sigma, Double_t alpha1, Double_t n1, Double_t alpha2, Double_t n2 )
{
    Double_t t = (x-mean)/sigma;
    
    if( t < -alpha1 )
    {
        Double_t a = TMath::Power( n1/alpha1, n1 )*TMath::Exp( -Square( alpha1 )/2 );
        Double_t b = n1/alpha1 - alpha1;
        return a/TMath::Power( b - t, n1 );
        
    } else if( t > alpha2 ) {
        
        Double_t a = TMath::Power( n2/alpha2, n2 )*TMath::Exp( -Square( alpha2 )/2 );
        Double_t b = n2/alpha2 - alpha2;
        return a/TMath::Power( b + t, n2 );
        
    } else return TMath::Exp( -Square( t )/2 );
}

Double_t FitFunctionCrystalBall2( Double_t* x, Double_t* par )
//{ return par[2]*(CrystalBall2( x[0], par[0], par[1], par[3], par[4], par[5], par[6] )+ (0.15/0.85)*CrystalBall2( x[0], par[0]*0.5, par[1]*0.7, par[3], par[4], par[5], par[6] ) ); }
{ return par[2]*(CrystalBall2( x[0], par[0], par[1], par[3], par[4], par[5], par[6] )+ (0.12)*CrystalBall2( x[0], par[0]*2.9/5.9, par[1]*0.57, par[3], par[4], par[5], par[6] ) ); }

Double_t FitFunctionCrystalBall2Calib( Double_t* x, Double_t* par )
{ return par[2]*(CrystalBall2( x[0], par[0], par[1], par[3], par[4], par[5], par[6] )+ (0.15)*CrystalBall2( x[0], par[0]*0.5, par[1]*0.7, par[3], par[4], par[5], par[6] ) ) + par[9]*TMath::Gaus( x[0], par[7], par[8]) ; }


Double_t GetResolution( TF1* f )
{
    // calculate relative width
    Double_t relWidth = f->GetParameter(1)/f->GetParameter(0);
    Double_t relFWHM = relWidth * 2*sqrt( 2*TMath::Log(2) ); // 2.sqrt(2.ln(2)) = 2.35
    //std::cout << "GetResolution - Energy resolution: " << relWidth << " (sigma), " << relFWHM << " (FWHM)" << std::endl;
    return relWidth;
}


Double_t GetResolutionError( TF1* f )
{ return GetResolution(f)*TMath::Sqrt((Square(f->GetParError(0)/f->GetParameter(0)) + Square(f->GetParError(1)/f->GetParameter(1)) )); }

Double_t GetFWHM( TF1* f )
{
    // calculate relative width
    Double_t relWidth = f->GetParameter(1)/f->GetParameter(0);
    Double_t relFWHM = relWidth * 2*sqrt( 2*TMath::Log(2) ); // 2.sqrt(2.ln(2)) = 2.35
    //std::cout << "GetResolution - Energy resolution: " << relWidth << " (sigma), " << relFWHM << " (FWHM)" << std::endl;
    return relFWHM;
}


Double_t GetFWHMError( TF1* f )
{ return GetFWHM(f)*TMath::Sqrt((Square(f->GetParError(0)/f->GetParameter(0)) + Square(f->GetParError(1)/f->GetParameter(1)) )); }


TH1* SmoothHist(TH1* hist, Int_t halfwidth) {
    // width defines on how many bins you want to average the hist
    Int_t n = hist->GetNbinsX();
    double xmax = hist->GetXaxis()->GetBinCenter(n);
    double xmin = hist->GetXaxis()->GetBinCenter(0);
    TString histName = Form("%s-smooth", hist->GetName());
    TH1* hSmooth = new TH1F(histName, histName, n, xmin, xmax);
    
    //std::cout << "xmax = " << xmax << " and xmin = " << xmin << std::endl;
    //std::cout << "number of bins = " << n << std::endl;
    
    // add a cut on spikes
    
    for (int i = 0; i < halfwidth; i++) {
        double meanValue = 0.;
        for (int k = 0; k < i+halfwidth; k++) meanValue += hist->GetBinContent(k);
        meanValue /= (halfwidth+i);
        hSmooth->SetBinContent(i, meanValue);
    }
    for (int i = halfwidth; i < n-halfwidth; i++) {
        double meanValue = 0;
        for (int k = 0; k < 2*halfwidth; k++) meanValue += hist->GetBinContent(i-halfwidth+k);
        meanValue /= (2*halfwidth);
        hSmooth->SetBinContent(i, meanValue);
    }
    for (int i = n-halfwidth; i < n; i++) {
        double meanValue = 0;
        for (int k = 0; k < halfwidth + (n-i); k++) meanValue += hist->GetBinContent(i-halfwidth+k);
        meanValue /= (halfwidth+(n-i));
        hSmooth->SetBinContent(i, meanValue);
    }
    
    return hSmooth;
    
}
