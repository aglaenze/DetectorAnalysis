#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "HistUtils.C"

using namespace std;

//____________________________________________
//void DrawSpectrum( TString path, string detectorName, vector <int> hvList, Int_t coarseGain, Double_t& gain, Double_t& gainError, Double_t& fwhm, Double_t& fwhmError, Double_t& resolution, Double_t& resolutionError, bool drawFit=true)
void DrawSpectrum(string detectorName, string date, TString fileName, Int_t coarseGain, Double_t& gain, Double_t& gainError, Double_t& fwhm, Double_t& fwhmError, Double_t& resolution, Double_t& resolutionError, bool drawFit=true)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    /*
    TString fileName = path + "spectrum";
    for (int k = 0; k<(int) hvList.size(); k++) {fileName += Form("-%d", hvList[k]);}
    fileName += Form("-%d.mca", coarseGain);
     */
    
    // histogram
    TH1* h = ReadData( fileName, coarseGain );
    h->Rebin(4);
    
    //Set up fit
    Int_t iBinMax = GetMaximumBin( h, 20);
    if (coarseGain > 120) iBinMax = GetMaximumBin( h, 100 );
    iBinMax = GetMaximumPos(h, 1);
    Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
    Double_t maximum = GetMaximum(h, iBinMax );
    
    /*
    while (maximum < 20) {
        h->Rebin(2);
        maximum = GetMaximum(h, iBinMax );
    }
     */
    
    //const Double_t maximum = h->GetMaximum();
    std::cout << "xMax " << xMax << std::endl;
    std::cout << "iMax " << iBinMax << std::endl;
    std::cout << "Maximum " << maximum << std::endl;
    
    h->Scale(1/maximum);
    h->SetMaximum(1.1);
    Double_t scaledRMS = h->GetRMS()/coarseGain;
    cout << "scaledRMS = " << scaledRMS << endl;
    
    
    //if (xMax < 0.9) {h->Draw(); return;}
    
    Double_t fitRangeMin;
    Double_t fitRangeMax;
    fitRangeMin = xMax - 0.7*scaledRMS ;
    fitRangeMax = xMax + 0.8*scaledRMS ;
    if (coarseGain > 60) {
        fitRangeMin = xMax - 5*scaledRMS ;
        fitRangeMax = xMax + 7*scaledRMS ;
    }
    TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
    f->SetParNames("Mean", "Sigma", "Amplitude");
    f->SetParameters(xMax, scaledRMS , 1);
    
    h->Draw();
    if (drawFit) {
        h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
        
        h->GetXaxis()->SetRangeUser(0, 3*xMax);
        f->Draw("same");
    }
    
    
    // Compute gain
    const Double_t nPrimary = 228.403; // Ar-iC4H10 95/5
    Double_t alpha = GetCalibrationAlpha(detectorName, date)*nPrimary;
    gain = f->GetParameter(0)/alpha;
    gainError = f->GetParError(0)/alpha;
    
    fwhm = 100.0 * GetFWHM(f) ;
    fwhmError = 100.0 * GetFWHMError(f) ;
    
    resolution = 100.0 * GetResolution(f) ;
    resolutionError = 100.0 * GetResolutionError(f) ;
    
}
