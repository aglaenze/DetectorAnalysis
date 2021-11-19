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

/*
 #include "HistUtils.C"
 #include "DetectorInfo.C"
 */

using namespace RooFit;
using namespace std;

/*
 To do: ajuster les yields en fonction de la proportion de photons à chaque étage
 */


//____________________________________________
//void DrawRooSpectrum( TString path, string detectorName, vector <int> hvList, Int_t coarseGain, Double_t& gain, Double_t& gainError, Double_t& fwhm, Double_t& fwhmError, Double_t& resolution, Double_t& resolutionError, int numGauss = 3)
RooPlot* DrawRooSpectrum(string detectorName, string date, TString fileName, Int_t coarseGain, Double_t& gain, Double_t& gainError, Double_t& gainMm, Double_t& gainMmError, Double_t& fwhm, Double_t& fwhmError, Double_t& resolution, Double_t& resolutionError, int numGauss = 3, bool pileUp = true, bool escape = false, bool alpha = false)
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
    h->Rebin(2);

    //Set up fit
    //Int_t iBinMax = GetMaximumBin( h, 10, 2OO);
    Int_t iBinMax = GetMaximumPos( h, 2);   // use position on x axis instead of bin number
    //Int_t iBinMax = GetMaximumBin(h);
    //Int_t iBinMax = GetMaximumPos( h, 2);
    //if (coarseGain > 120) iBinMax = GetMaximumBin( h, 100 );
    Double_t xMax = h->GetXaxis()->GetBinCenter( iBinMax );
    //xMax = 2.4;
    
    const Double_t maximum = GetMaximum(h, iBinMax );
    //const Double_t maximum = h->GetMaximum();
    cout << endl << endl << endl << "xMax " << xMax << endl;
    cout << "iMax " << iBinMax << endl;
    cout << "Maximum " << maximum << endl;
    h->Scale(1./maximum);
    //h->SetMaximum(1.1);
    Double_t scaledRMS = h->GetRMS()/coarseGain;
    cout << "scaledRMS = " << scaledRMS << endl;
    
    bool noisy = false;
    if (h->GetMaximum()>1.5*maximum) noisy = true;
    
    //return nullptr;
    
    double fitRangeMin = 0;
    for (int l = 0; l < h->GetNbinsX(); l++) {
        if (h->GetBinContent(l) > 0.1) {
            fitRangeMin = h->GetXaxis()->GetBinCenter(l);
            break;
        }
    }
    if (noisy) fitRangeMin = 0.2;
    cout << "fitRangeMin = " << fitRangeMin << endl << endl << endl << endl;
    //return nullptr;
    double d1 = GetzElectrode(detectorName, "pileUp bottom") - GetzElectrode(detectorName, "mesh");
    double d2 = GetzElectrode(detectorName, "mesh top") - GetzElectrode(detectorName, "pileUp bottom");
    double d3 = GetzElectrode(detectorName, "drift") - GetzElectrode(detectorName, "mesh top");
    
    // Update RooFit here
    //RooRealVar gainVar("gain","MCA/Coarse Gain", 1, 100);
    RooRealVar gainVar("gain","MCA/Coarse Gain", 1, 30);
    
    RooDataHist* rooHistGain = new RooDataHist("rooHistGain","Gain spectrum", RooArgList(gainVar), h);
    
    // MM peak
    RooRealVar mmMean("mmMean","mmMean", 5., 1, 9.);
    //RooRealVar mmMean("mmMean","mmMean", 1.4, 1.2, 1.6);
    //RooRealVar mmMean("mmMean","mmMean", 1.2);
    RooRealVar mmSigma("mmSigma","mmSigma", 0.05, 0.001, 5);
    RooGaussian *mmPeak = new RooGaussian("mmPeak","MM gain", gainVar, mmMean, mmSigma);
    
    if (numGauss < 3) {
        mmMean.setRange(xMax*0.9, xMax*1.1);
        mmMean.setVal(xMax);
    }
    
    // gem peak
    RooRealVar gemMean("gemMean","gemMean", 4, 1, 15);
    RooRealVar gemSigma("gemSigma","gemSigma", 2, 0.1, 10);
    RooGaussian *gemPeak = new RooGaussian("gemPeak","gem gain", gainVar, gemMean, gemSigma);
    
    // mesh top
    RooRealVar meshTopMean("meshTopMean","meshTopMean", xMax, xMax*0.9, xMax*1.1);
    RooRealVar meshTopSigma("meshTopSigma","meshTopSigma", 2, 0.001, 5);
    RooGaussian *meshTopPeak = new RooGaussian("meshTopPeak","Mesh top gain", gainVar, meshTopMean, meshTopSigma);
    
    // pile-up peak
    //RooRealVar pileUpMean("pileUpMean","pileUpMean", 10, 3, 16);
    RooRealVar pileUpScaling("pileUpScaling", "pileUpScaling", 1.8, 1.5, 2.2);
    RooFormulaVar pileUpMean("pileUpMean","meshTopMean*pileUpScaling", RooArgSet(meshTopMean, pileUpScaling));
    //RooFormulaVar pileUpMean("pileUpMean","mmMean*pileUpScaling",RooArgSet(mmMean, pileUpScaling));
    RooRealVar pileUpSigma("pileUpSigma","pileUpSigma", 1, 0.05, 10);
    RooGaussian *pileUpPeak = new RooGaussian("pileUpPeak","pileUp gain", gainVar, pileUpMean, pileUpSigma);
    
    /*
    // escape peak for the GEM
    //RooRealVar escapeScaling("escapeScaling", "escapeScaling", 0.5, 0.42, 0.62);
    RooRealVar escapeScaling("escapeScaling", "escapeScaling", 1.0847457);
    RooFormulaVar escapeMean("escapeMean","gemMean*escapeScaling", RooArgSet(gemMean, escapeScaling));
    //RooFormulaVar pileUpMean("pileUpMean","mmMean*pileUpScaling",RooArgSet(mmMean, pileUpScaling));
    RooRealVar escapeSigma("escapeSigma","escapeSigma", 0.05, 0.02, 10);
    RooGaussian *escapePeak = new RooGaussian("escapePeak","escape peak gain", gainVar, escapeMean, escapeSigma);
     */
    
    // escape peak for the mesh
    //RooRealVar escapeScaling("escapeScaling", "escapeScaling", 0.5, 0.42, 0.62);
    RooRealVar escapeScaling("escapeScaling", "escapeScaling", 0.4915);
    RooFormulaVar escapeMean("escapeMean","mmMean*escapeScaling", RooArgSet(mmMean, escapeScaling));
    //RooFormulaVar pileUpMean("pileUpMean","mmMean*pileUpScaling",RooArgSet(mmMean, pileUpScaling));
    RooRealVar escapeSigma("escapeSigma","escapeSigma", 0.3, 0.001, 3);
    RooGaussian *escapePeak = new RooGaussian("escapePeak","escape peak gain", gainVar, escapeMean, escapeSigma);
    
    // alpha peak at 6.4 keV
    //RooRealVar escapeScaling("escapeScaling", "escapeScaling", 0.5, 0.42, 0.62);
    RooRealVar alphaScaling("alphaScaling", "alphaScaling", 1.102);
    RooFormulaVar alphaMean("alphaMean","mmMean*alphaScaling", RooArgSet(mmMean, alphaScaling));
    //RooFormulaVar pileUpMean("pileUpMean","mmMean*pileUpScaling",RooArgSet(mmMean, pileUpScaling));
    RooRealVar alphaSigma("alphaSigma","alphaSigma", 0.1, 0.03, 5);
    RooGaussian *alphaPeak = new RooGaussian("alphaPeak","alpha peak gain", gainVar, alphaMean, alphaSigma);
    
    // Additional background
    RooRealVar b("b","b", -1, -10, -0.00001);
    RooExponential* bkg = new RooExponential("bkg","bkg", gainVar, b);
    RooRealVar b2("b2","b2", -50, -100, -0.1);
    RooExponential* bkg2 = new RooExponential("bkg2","bkg2", gainVar, b2);
    
    RooRealVar mmYield("mmYield","mmYield", 200, 10, 1000);
    //RooRealVar mmYield("mmYield","mmYield", 0);
    RooRealVar gemYield("gemYield","gemYield", 100, 0, 500);
    //RooRealVar gemYield("gemYield","gemYield", 0);
    RooRealVar meshTopYield("meshTopYield","meshTopYield", 300, 0, 1000);
    RooRealVar pileUpYield("pileUpYield","pileUpYield", 100, 0, 5000);
    RooRealVar escapeYield("escapeYield","escapeYield", 20, 1, 500);
    RooRealVar alphaYield("alphaYield","alphaYield", 10, 1, 20);
    /*
    RooRealVar mmYield("mmYield","mmYield", 0);
    RooRealVar gemYield("gemYield","gemYield", 0);
    RooRealVar meshTopYield("meshTopYield","meshTopYield", 0);
    RooRealVar pileUpYield("pileUpYield","pileUpYield", 0);
    RooRealVar escapeYield("escapeYield","escapeYield", 0);
    RooRealVar alphaYield("escapeYield","escapeYield", 0);
    */
    RooRealVar bkgYield("bkgYield","bkgYield", 20, 1, 10000);
    RooRealVar bkgYield2("bkgYield2","bkgYield2", 1000, 0, 100000000);
    
    //if (numGauss == 3 && coarseGain < 90) {mmYield.setVal(0); mmYield.setConstant();}
    //mmYield.setVal(0); mmYield.setConstant();
    
    RooArgList* pdfList = new RooArgList(*mmPeak, *bkg);
    RooArgList* yieldList = new RooArgList(mmYield, bkgYield);
    
    if (numGauss > 1) {
        pdfList->add(*gemPeak);
        yieldList->add(gemYield);
    }
    if (numGauss > 2) {
        pdfList->add(*meshTopPeak);
        yieldList->add(meshTopYield);
    }
    if (pileUp) {
        pdfList->add(*pileUpPeak);
        yieldList->add(pileUpYield);
    }
    if (escape) {
        pdfList->add(*escapePeak);
        yieldList->add(escapeYield);
    }
    if (alpha) {
        pdfList->add(*alphaPeak);
        yieldList->add(alphaYield);
    }
    if (noisy && gainVar.getMin() < 0.7) {
        pdfList->add(*bkg2);
        yieldList->add(bkgYield2);
    }
    
    
    RooPlot* frame = gainVar.frame(Title("Gain spectrum frame"));
    //frame->SetMaximum(1.0);
    //rooHistGain->plotOn(frame, DataError(RooAbsData::SumW2), Binning(200));
    rooHistGain->plotOn(frame, DataError(RooAbsData::None), MarkerStyle(1));
    //frame->GetXaxis()->SetRangeUser(0, 3*xMax);
    
    // Create fit model
    Double_t fitRangeMax = gainVar.getMax()*0.98;
    fitRangeMax = 20;
    //if (noisy && gainVar.getMin() > 0.2)
    //fitRangeMin = gainVar.getMin()*1.2;
    fitRangeMin=1.;
    RooAbsPdf* fitModel = new RooAddPdf("model", "model", *pdfList, *yieldList, kFALSE);
    RooFitResult* r = fitModel->fitTo(*rooHistGain, Range(fitRangeMin, fitRangeMax), Extended(), Minos(true), Strategy(1), Save());
    
    fitModel->plotOn(frame, Name("sum"), LineColor(kRed));
    fitModel->plotOn(frame,Name("mmPeak"),Components(*mmPeak),LineStyle(kDashed), LineColor(kGreen+2), LineWidth(1));
    
    if (numGauss > 1) {
        fitModel->plotOn(frame, Name("gemPeak"), Components(*gemPeak), LineStyle(kDashed), LineColor(kPink-7), LineWidth(1));
    }
    if (numGauss > 2) {
        fitModel->plotOn(frame, Name("meshTopPeak"), Components(*meshTopPeak), LineStyle(kDashed), LineColor(kBlue-4), LineWidth(1));
    }
    if (escape) {
        fitModel->plotOn(frame, Name("escapePeak"), Components(*escapePeak), LineStyle(kDashed), LineColor(kBlue), LineWidth(1));
    }
    if (pileUp) {
        fitModel->plotOn(frame,Name("pileUpPeak"), Components(*pileUpPeak), LineStyle(kDashed), LineColor(kRed-4), LineWidth(1));
    }
    if (alpha) {
        fitModel->plotOn(frame,Name("alphaPeak"), Components(*alphaPeak), LineStyle(kDashed), LineColor(kRed-2), LineWidth(1));
    }
    fitModel->plotOn(frame,Name("bkg"), Components(*bkg), LineStyle(kDashed), LineColor(kBlack), LineWidth(1));

    /*
     frame->Draw();
     lgd->Draw();
     */
    
    // Compute gain
    const Double_t nPrimary = 228.403; // Ar-iC4H10 95/5
    Double_t alphaParam = GetCalibrationAlpha(detectorName, date)*nPrimary;
    
    if (numGauss == 1) {    // one gaussian when there is only MM
        gain = mmMean.getVal()/alphaParam;
        gainError = mmMean.getError()/alphaParam;
        
        resolution = 100.0 * mmSigma.getVal()/mmMean.getVal();
        resolutionError = resolution *TMath::Sqrt((Square(mmMean.getError()/mmMean.getVal()) + Square(mmSigma.getError()/mmSigma.getVal()) ));
        
        fwhm = resolution * 2*sqrt( 2*TMath::Log(2) );
        fwhmError = resolutionError * 2*sqrt( 2*TMath::Log(2) );
    }
    else if (numGauss == 2) {   // 2 gaussian for MM + GEM
        gainMm = mmMean.getVal()/alphaParam;
        gainMmError = mmMean.getError()/alphaParam;
        
        gain = gemMean.getVal()/alphaParam;
        gainError = gemMean.getError()/alphaParam;
        
        resolution = 100.0 * gemSigma.getVal()/gemMean.getVal();
        resolutionError = resolution *TMath::Sqrt((Square(gemMean.getError()/gemMean.getVal()) + Square(gemSigma.getError()/gemSigma.getVal()) ));
        
        fwhm = resolution * 2*sqrt( 2*TMath::Log(2) );
        fwhmError = resolutionError * 2*sqrt( 2*TMath::Log(2) );
    }
    else {                  // Finally, 3 or 4 (because of pile-up) gaussians for all amplification stages
        gain = meshTopMean.getVal()/alphaParam;
        gainError = meshTopMean.getError()/alphaParam;
        
        resolution = 100.0 * meshTopSigma.getVal()/meshTopMean.getVal();
        resolutionError = resolution *TMath::Sqrt((Square(meshTopMean.getError()/meshTopMean.getVal()) + Square(meshTopSigma.getError()/meshTopSigma.getVal()) ));
        
        fwhm = resolution * 2*sqrt( 2*TMath::Log(2) );
        fwhmError = resolutionError * 2*sqrt( 2*TMath::Log(2) );
    }
    return frame;
}

TLegend* GetLegend(RooPlot* frame, int numGauss = 3, bool pileUp = true, bool escape = false, bool alpha = false) {
    TLegend* lgd = new TLegend(0.6, 0.6, 0.9, 0.9);
    lgd->AddEntry(frame->findObject("sum"), "Sum","L");
    lgd->AddEntry(frame->findObject("bkg"), "Background","L");
    lgd->AddEntry(frame->findObject("mmPeak"), "Micromegas","L");
    if (numGauss > 1) {
        lgd->AddEntry(frame->findObject("gemPeak"), "GEM","L");
    }
    if (numGauss > 2) {
        lgd->AddEntry(frame->findObject("meshTopPeak"), "Total","L");
    }
    if (escape) {
        lgd->AddEntry(frame->findObject("escapePeak"), "Escape peak","L");
    }
    if (pileUp) {
        lgd->AddEntry(frame->findObject("pileUpPeak"), "Pile-up","L");
    }
    if (alpha) {
        lgd->AddEntry(frame->findObject("alphaPeak"), "Alpha peak","L");
    }
    return lgd;
}



void RooGain() {
    
    TString path = "/Users/aglaenzer/Documents/Detector-Analysis/MCA/MpileUp1/2021-02-17/";
    string detectorName = "MpileUp1";
    vector <int> hvList = {360, 700, 950, 1200, 1400};
    Int_t coarseGain = 50;
    Double_t gain, gainError, fwhm, fwhmError, resolution, resolutionError;
    Double_t gainMm, gainMmError;
    
    string date = "";
    TObjArray* matches = TPRegexp("(.*)(\\d+)-(\\d+)-(\\d+)(.*)").MatchS( path );
    for (int k = 1; k < 4; k++) {
        date += (static_cast<TObjString*>(matches->At(k)))->GetString();
        if (k<3) date += "-";
    }
    
    TCanvas* cv = new TCanvas();
    //DrawRooSpectrum(path, detectorName, hvList, coarseGain, gain, gainError, gainMm, gainMmError, fwhm, fwhmError, resolution, resolutionError);
    TString fileName = path + "spectrum";
    for (int k = 0; k<(int) hvList.size(); k++) {fileName += Form("-%d", hvList[k]);}
    fileName += Form("-%d.mca", coarseGain);
    DrawRooSpectrum(detectorName, date, fileName, coarseGain, gain, gainError, gainMm, gainMmError, fwhm, fwhmError, resolution, resolutionError);
    
    cv->SaveAs("RooGain-test.pdf");
    
    cout << endl << endl << endl;
    cout << "Gain = " << gain << " pm " << gainError << endl;
    cout << "Resolution = " << resolution << " pm " << resolutionError << endl;
    cout << "FWHM = " << fwhm << " pm " << fwhmError << endl;
}
