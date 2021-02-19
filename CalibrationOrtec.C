#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "Include/Utils.C"

using namespace std;

//____________________________________________
Double_t FitFunctionLin( Double_t* x, Double_t* par )
{ return (par[0] + par[1] * x[0]); }

//____________________________________________
void CalibrationOrtec( string detectorName = "MGEM1" )
{
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	const Double_t Gain = 2000;
	//const Double_t nPrimary = 157;  //in Ne-CF4
	//const Double_t nPrimary = 222.68; // Ar-CO2
	const Double_t nPrimary = 228.403; // Ar-iC4H10 (95/5)
	
    /*
	const Int_t n = 9;
	//const Double_t inputListmV[n] = {47, 33, 26, 25, 17, 13, 11.5};
    //const Double_t inputListmV[n] = {44.5, 35, 31.2, 25.4, 22.9, 18.3, 16.1, 12.9, 11.4};
    
	// For MM
	//const Double_t mcaList[n] = {737, 512, 404, 359, 520, 410, 364};
    const Double_t mcaList[n] = {800, 629, 555, 435, 387, 303, 267, 206, 184};
	//const Double_t coarseGainList[n] = {100, 100, 100, 100, 200, 200, 200};
    const Double_t coarseGainList[n] = {100, 100, 100, 100, 100, 100, 100, 100, 100};
    */
    
    /*
    const Int_t n = 9;
    const Double_t inputListmV[n] = {147, 115, 101, 81.2, 73.7, 58.3, 50, 40, 35};
    const Double_t mcaList[n] = {793, 632, 556, 442, 819, 654, 577, 460, 404};
    const Double_t coarseGainList[n] = {100, 100, 100, 100, 200, 200, 200, 200, 200};
     */
    
    /*
    // 17/02/2021
    const Int_t n = 4;
    const Double_t inputListmV[n] = {33.3, 23.3, 16.5, 11.5};
    const Double_t mcaList[n] = {912, 631, 464, 331};
    const Double_t coarseGainList[n] = {200, 200, 200, 200};
     */
    
    // 18/02/2021
    const Int_t n = 4;
    const Double_t inputListmV[n] = {33.3, 23.3, 16.5, 11.5};
    const Double_t mcaList[n] = {495, 347, 509, 370};
    const Double_t coarseGainList[n] = {100, 100, 200, 200};
	
	Double_t mcaListNorm[n] = {};
	Double_t mcaListNorm2[n] = {};
	Double_t eNumber[n] = {};
	
	//const Int_t coarseGain = 200;
	const Double_t qe = 1.602e-19;
	const Double_t capa = 2.e-12;		// 2pF = capa de la mesh du détecteur
	const Double_t conv = 1.e-3;        // conversion factor from mV to V
	
	for (Int_t i=0; i<n; i++) {
		//mcaListNorm[i] = mcaList[i] / coarseGain;
		mcaListNorm[i] = mcaList[i] / coarseGainList[i];
		eNumber[i] = inputListmV[i]*(capa * conv) /qe;
	}
	
	gStyle->SetPadGridX( kTRUE );
	gStyle->SetPadGridY( kTRUE );
	
	// MCA channel vs input mV
	TCanvas* cv1 = new TCanvas( "cv1", "cv1", 900, 600 );
	TGraph* gr = new TGraph(n, eNumber, mcaListNorm);
	//TGraph* gr = new TGraph(n, inputListmV, mcaListNorm);
	gr->Draw("AP*");
	gr->SetTitle( "Ortec calibration" );
	gr->GetXaxis()->SetTitle( "Number of electrons" );
	//gr->GetXaxis()->SetTitle( "mV" );
	gr->GetYaxis()->SetTitle( "MCA channels/ Coarse Gain" );
	
	Double_t xMin = 100;
	Double_t xMax = 0.;
	for (Int_t i =0; i<n; i++){
		if (xMax < eNumber[i]) xMax= eNumber[i];
		if (xMin > eNumber[i]) xMin= eNumber[i];
	}
	xMax *= 1.1;
	xMin *= 0.9;
	Int_t paramNum = 2;
	
	//xMin = inputListmV[1];
	//xMax = inputListmV[n-3];
	
	TF1* f = new TF1( "FitFunctionLin", FitFunctionLin, xMin, xMax, paramNum);
	f->SetParName( 0, "const" );
	f->SetParName( 1, "slope" );
	
	f->SetParameter(0,1e-6);
	f->SetParameter(1,1e-6);
	
	gr->Fit( f , "0", "", xMin, xMax);
	
	f->SetLineColor(2);
	f->Draw("same");
	
	double yMax = gr->GetHistogram()->GetMaximum();
	double yMin = gr->GetHistogram()->GetMinimum();
	double xLine = Gain*nPrimary/(capa * conv) *qe;
	//double xLine = Gain*nPrimary;
	TLine *l = new TLine(xLine, yMin, xLine, yMax);
	l->SetLineColor(kBlue);
	l->Draw("same");
	
	std::cout<< "\n\nFor a coarse gain of 200, you have a peak at MCA : " << 200*f->Eval(xLine) << std::endl << std::endl << std::endl;
	
	// PutText( 0.2, 0.8, Form( "y =  %.3g + %.3g x", f->GetParameter(0), f->GetParameter(1) ) );
	
	TLegend* legend = new TLegend(0.55,0.25,0.9,0.5);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legend->AddEntry(gr, "Data with coarse gain = 200", "P*");
	legend->AddEntry(l, "Gain = 2000", "l");
	legend->AddEntry(f, Form( "y =  %.3g + %.3g x", f->GetParameter(0), f->GetParameter(1) ), "l");
	legend->Draw();
	
	cv1->SaveAs(Form("Figures/%s/Calibration_Ortec-new2.pdf", detectorName.c_str()));
	
	
}
