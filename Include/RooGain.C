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

/*
 To do: ajuster les yields en fonction de la proportion de photons à chaque étage
 */

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
    //h->Rebin(2);
    
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
    
    double fitRangeMin = 0;
    for (int l = 0; l < h->GetNbinsX(); l++) {
        if (h->GetBinContent(l) > 0.1) {
            fitRangeMin = h->GetXaxis()->GetBinCenter(l-1);
            break;
        }
    }
    cout << "fitRangeMin = " << fitRangeMin << endl;
    double d1 = GetzElectrode(detectorName, "GEM bottom") - GetzElectrode(detectorName, "mesh");
    double d2 = GetzElectrode(detectorName, "mesh top") - GetzElectrode(detectorName, "GEM bottom");
    double d3 = GetzElectrode(detectorName, "drift") - GetzElectrode(detectorName, "mesh top");
    
    // Update RooFit here
    RooRealVar gainVar("gain","MCA/Coarse Gain",0, 100);
    
    RooDataHist* rooHistGain = new RooDataHist("rooHistGain","Gain spectrum", RooArgList(gainVar), h);
    
    // MM peak
    RooRealVar mmMean("mmMean","mmMean", 2, 1., 2.5);
    RooRealVar mmSigma("mmSigma","mmSigma", 0.2, 0.01, 2.5);
    RooGaussian *mmPeak = new RooGaussian("mmPeak","MM gain", gainVar, mmMean, mmSigma);
    
    // GEM peak
    RooRealVar gemMean("gemMean","gemMean", 10, 4, 16);
    RooRealVar gemSigma("gemSigma","gemSigma",2, 0.5, 7);
    RooGaussian *gemPeak = new RooGaussian("gemPeak","GEM gain", gainVar, gemMean, gemSigma);
    
    // mesh top
    RooRealVar meshTopMean("meshTopMean","meshTopMean", xMax, 5, 9);
    RooRealVar meshTopSigma("meshTopSigma","meshTopSigma", scaledRMS, 0.1, 3);
    RooGaussian *meshTopPeak = new RooGaussian("meshTopPeak","Mesh top gain", gainVar, meshTopMean, meshTopSigma);
    
    // Additional background
    RooRealVar b("b","b", 0, -2, 1.);
    RooExponential* bkg = new RooExponential("bkg","bkg", gainVar, b);
    
    RooRealVar mmYield("mmYield","mmYield", 3, 0, 200);
    RooRealVar gemYield("gemYield","gemYield", 1, 0, 100);
    RooRealVar meshTopYield("meshTopYield","meshTopYield", xMax*scaledRMS, 0, 500);
    RooRealVar bkgYield("bkgYield","bkgYield", 3, 0, 1000);
    
    RooArgList* pdfList = new RooArgList(*mmPeak, *gemPeak, *meshTopPeak, *bkg);
    RooArgList yieldList = RooArgList(mmYield, gemYield, meshTopYield, bkgYield);
    // Create fit model
    RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, yieldList, kFALSE);
    RooFitResult* r = fitModel->fitTo(*rooHistGain, Range(fitRangeMin, 20), Extended(), Minos(true), Strategy(1), Save());
    
    RooPlot* frame = gainVar.frame(Title("Gain spectrum frame"));
    //rooHistGain->plotOn(frame, DataError(RooAbsData::SumW2), Binning(200));
    rooHistGain->plotOn(frame, DataError(RooAbsData::None));

    fitModel->plotOn(frame, Name("sum"), LineColor(kRed));
    fitModel->plotOn(frame,Name("mmPeak"),Components(*mmPeak),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    fitModel->plotOn(frame,Name("gemPeak"),Components(*gemPeak),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    fitModel->plotOn(frame,Name("meshTopPeak"),Components(*meshTopPeak),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    fitModel->plotOn(frame,Name("bkg"),Components(*bkg),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    
    frame->Draw();
    
    
	// Compute gain
	const Double_t nPrimary = 228.403; // Ar-iC4H10 95/5
	Double_t alpha = GetCalibrationAlpha(detectorName)*nPrimary;

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
