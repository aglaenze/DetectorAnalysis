#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>

#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooHist.h"
#include "RooFitResult.h"

#include "HistUtils.C"
#include "DetectorInfo.C"

using namespace RooFit;
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
    h->Rebin(10);
    
    // Update RooFit here
    RooRealVar gainVar("gain","MCA/Coarse Gain",0, 100);
    
    //RooArgSet variables(gainVar);
    RooDataHist* rooHistGain = new RooDataHist("rooHistGain","Gain spectrum", RooArgList(gainVar), h);
    
    // MM peak
    RooRealVar mmMean("mmMean","mmMean",2,1,3);
    RooRealVar mmSigma("mmSigma","mmSigma",0.2,0,2);
    RooGaussian *mmPeak = new RooGaussian("mmPeak","MM gain", gainVar, mmMean, mmSigma);
    
    // GEM peak
    RooRealVar gemMean("gemMean","gemMean",10,6,16);
    RooRealVar gemSigma("gemSigma","gemSigma",2,0,7);
    RooGaussian *gemPeak = new RooGaussian("gemPeak","GEM gain", gainVar, gemMean, gemSigma);
    
    // mesh top
    RooRealVar meshTopMean("meshTopMean","meshTopMean",6,5,9);
    RooRealVar meshTopSigma("meshTopSigma","meshTopSigma",0.5,0,3);
    RooGaussian *meshTopPeak = new RooGaussian("meshTopPeak","Mesh top gain", gainVar, meshTopMean, meshTopSigma);
    
    // Additional background
    RooRealVar b("b","b", -1, -20, 1.);
    RooExponential* bkg = new RooExponential("bkg","bkg", gainVar, b);
    
    RooRealVar mmYield("mmYield","mmYield", 3, 0, 200);
    RooRealVar gemYield("gemYield","gemYield", 1, 0, 100);
    RooRealVar meshTopYield("meshTopYield","meshTopYield", 1, 0, 500);
    RooRealVar bkgYield("bkgYield","bkgYield", 3, 0, 200);
    
    RooArgList* pdfList = new RooArgList(*mmPeak, *gemPeak, *meshTopPeak, *bkg);
    RooArgList yieldList = RooArgList(mmYield, gemYield, meshTopYield, bkgYield);
    // Create fit model
    RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
    RooFitResult* r = fitModel->fitTo(*rooHistGain, Range(1.7, 20), Extended(), Minos(true), Strategy(1), Save());
    
    RooPlot* frame = gainVar.frame(Title("Gain spectrum frame"));
    //rooHistGain->plotOn(frame, DataError(RooAbsData::SumW2), Binning(200));
    rooHistGain->plotOn(frame, DataError(RooAbsData::None), Binning(200));

    fitModel->plotOn(frame, Name("sum"), LineColor(kRed));
    fitModel->plotOn(frame,Name("mmPeak"),Components(*mmPeak),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    fitModel->plotOn(frame,Name("gemPeak"),Components(*gemPeak),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    fitModel->plotOn(frame,Name("meshTopPeak"),Components(*meshTopPeak),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    fitModel->plotOn(frame,Name("bkg"),Components(*bkg),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    
    frame->Draw();
    
    
    // End of RooFit
	
    /*
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
     */
	
	
	// Compute gain
	const Double_t nPrimary = 228.403; // Ar-iC4H10 95/5
	Double_t alpha = GetCalibrationAlpha(detectorName)*nPrimary;
    /*
	gain = f->GetParameter(0)/alpha;
	gainError = f->GetParError(0)/alpha;
	
	fwhm = 100.0 * GetFWHM(f) ;
	fwhmError = 100.0 * GetFWHMError(f) ;
	
	resolution = 100.0 * GetResolution(f) ;
	resolutionError = 100.0 * GetResolutionError(f) ;
     */
    gain = meshTopMean.getVal()/alpha;
    gainError = meshTopMean.getError()/alpha;
    
    resolution = 100.0 * meshTopSigma.getVal()/meshTopMean.getVal();
    resolutionError = resolution *TMath::Sqrt((Square(meshTopMean.getError()/meshTopMean.getVal()) + Square(meshTopSigma.getError()/meshTopSigma.getVal()) ));
    
    fwhm = resolution * 2*sqrt( 2*TMath::Log(2) );
    fwhmError = resolutionError * 2*sqrt( 2*TMath::Log(2) );
}


void RooGain() {
    
    TString path = "/Users/aglaenzer/Documents/Detector-Analysis/MCA/MGEM1/2021-02-17/";
    string detectorName = "MGEM1";
    vector <int> hvList = {360, 700, 950, 1200, 1400};
    Int_t coarseGain = 50;
    Double_t gain, gainError, fwhm, fwhmError, resolution, resolutionError;
    
    TCanvas* cv = new TCanvas();
    DrawSpectrum(path, detectorName, hvList, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError);
    
    cv->SaveAs("RooGain-test.pdf");
    
    cout << endl << endl << endl;
    cout << "Gain = " << gain << " pm " << gainError << endl;
    cout << "Resolution = " << resolution << " pm " << resolutionError << endl;
    cout << "FWHM = " << fwhm << " pm " << fwhmError << endl;
}
