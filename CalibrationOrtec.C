#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "Include/Utils.C"

using namespace std;

//____________________________________________
Double_t FitFunctionLin( Double_t* x, Double_t* par )
{ return (par[0] + par[1] * x[0]); }

//____________________________________________
void CalibrationOrtec( string detectorName = "MGEM3" )
{
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	const Double_t Gain = 2000;
	//const Double_t nPrimary = 157;  //in Ne-CF4
	//const Double_t nPrimary = 222.68; // Ar-CO2
	//const Double_t nPrimary = 228.403; // Ar-iC4H10 (95/5)
    const Double_t nPrimary = 224; // Ar-iC4H10 (95/5) d'apres Garfield++
    
    bool fEl = true;   // if true, draws as a function of number of electrons
                        // if false, draws as a function of mV
	
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
    
    /*
    // 18/02/2021
    const Int_t n = 4;
    const Double_t inputListmV[n] = {33.3, 23.3, 16.5, 11.5};
    const Double_t mcaList[n] = {495, 347, 509, 370};
    const Double_t coarseGainList[n] = {100, 100, 200, 200};
     */
    
    /*
    // 18/02/2021
    const Int_t n = 7;
    const Double_t inputListmV[n] = {33.3, 26, 23.3, 16.5, 13.3, 11.5, 9.3};
    const Double_t mcaList[n] = {684, 545, 485, 346, 277, 246, 197};
    const Double_t coarseGainList[n] = {200, 200, 200, 200, 200, 200, 200};
     */
    
    /*
    // 25/02/2021 MGEM3
    const Int_t n = 8;
    const Double_t inputListmV[n] = {131, 93.3, 66.2, 46.4, 33.3, 23.2, 16.3, 11.8};
    const Double_t mcaList[n] = {653, 466, 663, 476, 338, 514, 375, 280};
    const Double_t coarseGainList[n] = {50, 50, 100, 100, 100, 200, 200, 200};
     */
    
    /*
    // 08/04/2021 MGEM1
    const Int_t n = 4;
    const Double_t inputListmV[n] = {65.4, 46.8, 32.7, 23.1};
    const Double_t mcaList[n] = {650, 460, 326, 231};
    const Double_t coarseGainList[n] = {100, 100, 100, 100};
     */
    
    /*
    // 15/04/2021 MGEM3 mesh down
    const Int_t n = 4;
    const Double_t inputListmV[n] = {65.4, 46.8, 32.7, 23.1};
    const Double_t mcaList[n] = {722, 505, 354, 249};
    const Double_t coarseGainList[n] = {100, 100, 100, 100};
    */
    
    /*
    // 15/04/2021 MGEM3 mesh top
    const Int_t n = 4;
    const Double_t inputListmV[n] = {65.4, 46.8, 32.7, 23.1};
    const Double_t mcaList[n] = {675, 479, 341, 242};
    const Double_t coarseGainList[n] = {100, 100, 100, 100};
     */
    
    /*
    // 6/10/2021 MGEM3 mesh down
    const Int_t n = 4;
    const Double_t inputListmV[n] = {75.0, 53.5, 38.0, 27.2};
    const Double_t mcaList[n] = {734, 521, 373, 268};
    const Double_t coarseGainList[n] = {200, 200, 200, 200};
     */
    
    // 7/10/2021 MGEM1 mesh down
    const Int_t n = 4;
    const Double_t inputListmV[n] = {75.0, 53.5, 38.0, 27.2};
    const Double_t mcaList[n] = {349, 249, 180, 129};
    const Double_t coarseGainList[n] = {100, 100, 100, 100};
	
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
    
    Double_t xMin = 100;
    Double_t xMax = 0.;
	
	// MCA channel vs input mV
	TCanvas* cv1 = new TCanvas( "cv1", "cv1", 900, 600 );
    TGraph* gr = nullptr;
    double xLine;
    if (fEl) {
        gr = new TGraph(n, eNumber, mcaListNorm);
        gr->GetXaxis()->SetTitle( "Number of electrons" );
        for (Int_t i =0; i<n; i++){
            if (xMax < eNumber[i]) xMax= eNumber[i];
            if (xMin > eNumber[i]) xMin= eNumber[i];
        }
        xLine = Gain*nPrimary;
    }
    else {
        gr = new TGraph(n, inputListmV, mcaListNorm);
        gr->GetXaxis()->SetTitle( "mV" );
        for (Int_t i =0; i<n; i++){
            if (xMax < inputListmV[i]) xMax= inputListmV[i];
            if (xMin > inputListmV[i]) xMin= inputListmV[i];
        }
        xLine = Gain*nPrimary/(capa * conv) *qe;
    }
    xMax *= 1.1;
    xMin *= 0.9;
	gr->Draw("AP*");
	gr->SetTitle( "Ortec calibration" );
	gr->GetYaxis()->SetTitle( "MCA channels/ Coarse Gain" );

	Int_t paramNum = 2;
	

	TF1* f = new TF1( "FitFunctionLin", FitFunctionLin, xMin, xMax, paramNum);
	f->SetParName( 0, "const" );
	f->SetParName( 1, "slope" );
	
	f->SetParameter(0,1e-6);
	f->SetParameter(1,1e-6);
    //f->FixParameter( 0, 0);
	
	gr->Fit( f , "0", "", xMin, xMax);
	
	f->SetLineColor(2);
	f->Draw("same");
	
	double yMax = gr->GetHistogram()->GetMaximum();
	double yMin = gr->GetHistogram()->GetMinimum();
	TLine *l = new TLine(xLine, yMin, xLine, yMax);
	l->SetLineColor(kBlue);
	l->Draw("same");

	std::cout<< "\n\nFor a coarse gain of 100, you have a gain of 2000 at MCA : " << 100*f->Eval(xLine) << std::endl << std::endl << std::endl;
	
	// PutText( 0.2, 0.8, Form( "y =  %.3g + %.3g x", f->GetParameter(0), f->GetParameter(1) ) );
	
	TLegend* legend = new TLegend(0.55,0.25,0.9,0.5);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legend->AddEntry(gr, "Data with coarse gain = 200", "P*");
	legend->AddEntry(l, "Gain = 2000", "l");
	legend->AddEntry(f, Form( "y =  %.3g + %.3g x", f->GetParameter(0), f->GetParameter(1) ), "l");
	legend->Draw();
	
	cv1->SaveAs(Form("Figures/%s/Calibration_Ortec-mesh-down.pdf", detectorName.c_str()));
	
	
}
