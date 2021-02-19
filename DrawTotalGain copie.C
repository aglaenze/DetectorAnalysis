#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "_Utils.C"
#include "_HistUtils.C"

using namespace std;


//____________________________________________
int DrawTotalGain( string detectorName = "MGEM1" )
{
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetLabelSize(0.05, "X");
	gStyle->SetLabelSize(0.05, "Y");
	
	
	// Gain
	const int num = 4;
	int dV1[num] = {250, 280, 310, 340}; 	// GEM and mesh top
	int dV2[num] = {200, 250, 280, 310, 340}; 	// in GEM
	
	Double_t gainList[11] = {750.006, 1581.879, 1881.088, 545.404, 1182.419, 578.932, 1245.877, 1493.786, 607.657, 1693.861, 2004.780};
	Double_t gainErrorList[11] = {89.617, 193.627, 171.350, 111.128, 169.138, 110.095, 144.929, 193.373, 91.624, 187.809, 234.869};
	Double_t hvMeshList[11] = {390, 390, 390, 390, 390, 390, 390, 390, 390, 390, 390};
	Double_t hvGemDownList[11] = {600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600};
	Double_t hvGemTopList[11] = {850, 850, 850, 880, 880, 910, 910, 910, 940, 940, 940};
	Double_t hvMeshTopList[11] = {1180, 1200, 1200, 1190, 1220, 1220, 1250, 1250, 1250, 1280, 1280};
	Double_t dvMeshTopList[11] = {330, 350, 350, 310, 340, 310, 340, 340, 310, 340, 340};
	Double_t coarseGainList[11] = {100, 100, 50, 100, 100, 100, 100, 50, 100, 100, 50};
	
	
	
	TCanvas* cv = new TCanvas();
	cv->Divide(2,2);
	
	for (int i = 0; i< num; i++) {
		
		cv->cd(i+1);
		gPad->SetLeftMargin(0.15);
		gPad->SetGrid();
		gPad->SetLogy();
		
		vector<float> gainVec[num] = {};
		vector<float> gainErrorVec[num] = {};
		vector<int> dvMeshTopVec[num] = {};
		
		for (int l = 0; l<num; l++) {
			for (int k = 0; k < 30; k++) {
				if ( (hvGemDownList[k] - hvMeshList[k] == dV1[l]) && (hvGemTopList[k]-hvGemDownList[k] == dV2[i]) ) {
					gainVec[l].push_back(gainList[k]);
					gainErrorVec[l].push_back(gainErrorList[k]);
					dvMeshTopVec[l].push_back(dvMeshTopList[k]);
					
				}
			}
		}
		
		
		TLegend* legend1 = new TLegend(0.55,0.1,0.9,0.4);
		legend1->SetTextSize(0.045);
		
		TLegend* legend2 = new TLegend(0.1,0.1,0.9,0.9);
		legend2->SetTextSize(0.05);
		
		for (int l = 0; l<num; l++) {
			// create list from vector
			const int nn = (int)dvMeshTopVec[l].size();
			//TGraphErrors* gr1 = CreateTGraphError(hvMeshTopVec[l], gainVec[l], gainErrorVec[l]);
			Double_t li1[nn], li2[nn], li3[nn], li4[nn], li5[nn];
			for (int k = 0; k< nn; k++) {
				li1[k] = dvMeshTopVec[l][k];
				li3[k] = gainVec[l][k];
				li4[k] = gainErrorVec[l][k];
			}
			TGraphErrors* gr1 = new TGraphErrors(nn, li1, li3, 0, li4);
			gr1->SetTitle(Form("Gain = f(HV) for dV(GEM) = %dV", dV1[i]));
			gr1->GetXaxis()->SetTitle( "HV (V)" );
			gr1->GetXaxis()->SetLimits(140, 270);
			gr1->GetYaxis()->SetTitle( "Gain" );
			gr1->SetMarkerStyle(20);
			gr1->GetHistogram()->SetMinimum(1);   // along Y axis
			gr1->GetHistogram()->SetMaximum(15000);   // along Y axis
			
			legend1->AddEntry(gr1, Form("dV(tr) = %d V", dV2[l]), "lp");
			
			/*
			// Exp fit function
			TF1* f = new TF1( "FitFunctionExp", FitFunctionExp, li1[0]-5, li1[nn-1]+5, 2 );
			f->SetParameter(0, -7 );
			f->SetParameter(1, 0.04 );
			f->SetParName( 0, "const" );
			f->SetParName( 1, "slope" );
			gr1->Fit( f, "0" );
			f->SetLineColor(l+1);
			f->Draw("same");
			
			legend2->AddEntry(f, Form("dV = %dV: y = exp(%.4f + %.4f x)", dV1[l], f->GetParameter(0), f->GetParameter(1) ), "l");
			 */
			
			gr1->SetLineColor(l+1);
			gr1->SetMarkerColor(l+1);
			if (l == 0) gr1->Draw("ALP");
			else gr1->Draw("LP");
			
		}
		
		legend1->Draw();
	}
	
	/*
	 cv->cd(3);
	 legend2->Draw();
	 */
	
	
	cv->SaveAs(Form("Figures/%s/GainTotal-test.pdf", detectorName.c_str()));
	
	return 0;
	
}

