#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "HistUtils.C"

using namespace std;

//____________________________________________
void DrawSpectrum( TString path, string detectorName, vector <int> hvList, Int_t coarseGain, Double_t& gain, Double_t& gainError, Double_t& fwhm, Double_t& fwhmError, Double_t& resolution, Double_t& resolutionError)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	bool simple = false;
	int arbitraryHv = hvList[0];
	for (int k = 1; k<(int)hvList.size(); k++) {
		if (hvList[k] == arbitraryHv) {simple = true; arbitraryHv = hvList[k];}
	}
	
	TString fileName = path + "spectrum";
	for (int k = 0; k<(int) hvList.size(); k++) {fileName += Form("-%d", hvList[k]);}
	fileName += Form("-%d.mca", coarseGain);
	
	// histogram
	TH1* h = ReadData( fileName, coarseGain );
	
	//Set up fit
	Int_t iBinMax = GetMaximumBin( h, 100 );
	if (coarseGain > 120) iBinMax = GetMaximumBin( h, 300 );
	Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
	
	const Double_t maximum = GetMaximum(h, iBinMax );
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
	TF1* f;
	simple = true;
	if (simple) {
		fitRangeMin = xMax - 0.5*scaledRMS ;
		fitRangeMax = xMax + 1.*scaledRMS ;
		f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
		f->SetParNames("Mean", "Sigma", "Amplitude");
		f->SetParameters(xMax, scaledRMS , maximum);
	}
	else {
        /*
		fitRangeMin = xMax - scaledRMS ;
		fitRangeMax = xMax + scaledRMS ;
         */
        fitRangeMin = 1.7;
        fitRangeMax = 16;

        double d1 = GetzElectrode(detectorName, "GEM bottom") - GetzElectrode(detectorName, "mesh");
        double d2 = GetzElectrode(detectorName, "mesh top") - GetzElectrode(detectorName, "GEM bottom");
        double d3 = GetzElectrode(detectorName, "drift") - GetzElectrode(detectorName, "mesh top");
        //double zDrift = GetzElectrode(detectorName, "drift");
        //cout << "\n\n\n\nzDrift = " << zDrift << endl;
		f = new TF1( "FitFunction", Fit3GaussAndExp, fitRangeMin, fitRangeMax, 11);
		f->SetParNames("MeanTotal", "SigmaTotal", "AmplitudeTotal", "MeanGem", "SigmaGem", "AmplitudeGem", "MeanMm", "SigmaMm", "AmplitudeMm", "BkgAmp", "BkgSlope");
		//f->SetParameters(xMax, scaledRMS , maximum, xMax*1.5, scaledRMS*1.5, maximum*d2/d3, xMax/2, scaledRMS/2, maximum*d1/d3, 10, -0.2);
        f->SetParameters(6, 1, 1, 11, 5, 0.2, 2, 0.5, 0.5, 10, -0.2);
        /*
        f = new TF1( "FitFunction", Fit3GaussAndExp, fitRangeMin, fitRangeMax, 9);
        f->SetParNames("MeanTotal", "SigmaTotal", "AmplitudeTotal", "MeanGem", "SigmaGem", "AmplitudeGem", "MeanMm", "SigmaMm", "AmplitudeMm");
        f->SetParameters(xMax, scaledRMS , maximum, xMax*2, scaledRMS*3, maximum*d2/d3, xMax/2, scaledRMS/3, maximum*d1/d3);
        */
        /*
        f->SetParLimits(0, xMax*0.8, xMax*1.3);
        f->SetParLimits(1, 0., 3*scaledRMS);
        //f->SetParLimits(3, 2, 10*xMax);
        f->SetParLimits(4, 0., 2*scaledRMS);
        //f->SetParLimits(6, 0., 10);
        f->SetParLimits(7, 0., 10*scaledRMS);
         */
        f->SetParLimits(0, 5, 8);
        f->SetParLimits(1, 0.1, 3);
        f->SetParLimits(2, 0, 3);
        f->SetParLimits(3, 7, 15);
        f->SetParLimits(4, 2, 5);
        f->SetParLimits(5, 0., 0.8);
        f->SetParLimits(6, 1.5, 2.5);
        f->SetParLimits(7, 0.2, 1.5);
        f->SetParLimits(8, 0.3, 0.7);

	}
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
