#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "_HistUtils.C"

using namespace std;

//____________________________________________
void DrawSpectrum( TString path, string detectorName, vector <int> hvList, Int_t coarseGain, Double_t& gain, Double_t& gainError, Double_t& fwhm, Double_t& fwhmError, Double_t& resolution, Double_t& resolutionError)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	TString fileName = path + "spectrum";
	for (int k = 0; k<(int) hvList.size(); k++) {fileName += Form("-%d", hvList[k]);}
	fileName += Form("-%d.mca", coarseGain);
	
	// histogram
	TH1* h = ReadData( fileName, coarseGain );
	
	//Set up fit
	Int_t iBinMax = GetMaximumBin( h, 20 );
	if (coarseGain > 120) iBinMax = GetMaximumBin( h, 50 );
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	const Double_t maximum = GetMaximum(h, iBinMax );
	//const Double_t maximum = h->GetMaximum();
	std::cout << "xMax " << xMax << std::endl;
	std::cout << "iMax " << iBinMax << std::endl;
	std::cout << "Maximum " << maximum << std::endl;
	
	h->Scale(1/maximum);
	h->SetMaximum(1.1);
	Double_t scaledRMS = h->GetRMS()/coarseGain;
	
	if (xMax < 0.9) {h->Draw(); return;}
	
	Double_t fitRangeMin = xMax - 0.8*scaledRMS ;
	Double_t fitRangeMax = xMax + scaledRMS ;
	TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("Mean", "Sigma", "Amplitude");
	f->SetParameters(xMax, scaledRMS , maximum);
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	
	h->GetXaxis()->SetRangeUser(0, 3*xMax);
	h->Draw();
	f->Draw("same");
	
	// Compute gain
	const Double_t nPrimary = 228.403; // Ar-iC4H10 95/5
	Double_t alpha = GetCalibrationAlpha(detectorName)*nPrimary;
	gain = f->GetParameter(0)/alpha;
	gainError = f->GetParError(0)/alpha;

	fwhm = 100.0 * GetFWHM(f) ;
	fwhmError = 100.0 * GetFWHMError(f) ;
	
	resolution = 100.0 * GetResolution(f) ;
	resolutionError = 100.0 * GetResolutionError(f) ;
}
