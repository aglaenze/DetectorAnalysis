#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "_Utils.C"
#include "_HistUtils.C"

using namespace std;


//____________________________________________
int DrawTotal1Curve( string detectorName = "MGEM1" )
{
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetLabelSize(0.05, "X");
	gStyle->SetLabelSize(0.05, "Y");
	
	
	// Gain
	Double_t gainList[21] = {752.718, 828.799, 710.361, 762.010, 670.914, 716.349, 648.440, 660.789, 618.227, 674.999, 624.699, 659.180, 343.753, 619.266, 670.848, 757.449, 737.441, 757.230, 838.239, 745.389, 785.140};
	Double_t gainErrorList[21] = {24.472, 22.348, 24.694, 19.093, 23.897, 18.988, 27.449, 17.543, 24.464, 20.872, 27.067, 19.069, 25.647, 24.468, 20.029, 24.207, 18.692, 22.990, 19.623, 23.444, 19.445};
	Double_t hvMeshList[21] = {320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320, 320};
	Double_t hvGemDownList[21] = {600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600};
	Double_t hvGemTopList[21] = {900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900};
	Double_t hvMeshTopList[21] = {1000, 1000, 1020, 1020, 1050, 1050, 1070, 1070, 1100, 1100, 1150, 1150, 900, 920, 920, 940, 940, 960, 960, 980, 980};
	Double_t dvMeshTopList[21] = {100, 100, 120, 120, 150, 150, 170, 170, 200, 200, 250, 250, 0, 20, 20, 40, 40, 60, 60, 80, 80};
	Double_t coarseGainList[21] = {100, 200, 100, 200, 100, 200, 100, 200, 100, 200, 100, 200, 200, 100, 200, 100, 200, 100, 200, 100, 200};
	

	vector<double> gain1Vec = {}, gain2Vec = {}, gainError1Vec = {}, gainError2Vec = {};
	vector<int> hv1Vec = {}, hv2Vec = {};
	
	for (int l = 0; l<21; l++) {
		if (coarseGainList[l] == 100) {
			gain1Vec.push_back(gainList[l]);
			gainError1Vec.push_back(gainErrorList[l]);
			hv1Vec.push_back(dvMeshTopList[l]);
		}
		else if (coarseGainList[l] == 200) {
			gain2Vec.push_back(gainList[l]);
			gainError2Vec.push_back(gainErrorList[l]);
			hv2Vec.push_back(dvMeshTopList[l]);
		}
		else return -1;
	}
	
	
	
	TCanvas* cv = new TCanvas("","", 800, 500);
	cv->Divide(2);
	
	
	cv->cd(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetGrid();
	gPad->SetLogy();

	
	// create list from vector
	const int nn = hv1Vec.size();
	//TGraphErrors* gr1 = CreateTGraphError(hvMeshTopVec[l], gainVec[l], gainErrorVec[l]);
	Double_t li1[nn], li2[nn], li3[nn], li4[nn], li5[nn];
	Double_t li12[nn], li32[nn], li42[nn];
	for (int k = 0; k< nn; k++) {
		li1[k] = hv1Vec[k];
		li3[k] = gain1Vec[k];
		li4[k] = gainError1Vec[k];
	}
	const int nn2 =  hv2Vec.size();
	for (int k = 0; k< nn2; k++) {
		li12[k] = hv2Vec[k];
		li32[k] = gain2Vec[k];
		li42[k] = gainError2Vec[k];
	}
	TGraphErrors* gr1 = new TGraphErrors(nn, li1, li3, 0, li4);
	gr1->SetTitle(Form("Gain = f(#DeltaV_{mesh top})"));
	gr1->GetXaxis()->SetTitle( "#DeltaV (V)" );
	gr1->GetXaxis()->SetLimits(-10, 270);
	gr1->GetYaxis()->SetTitle( "Gain" );
	gr1->SetMarkerStyle(20);
	gr1->GetHistogram()->SetMinimum(50);   // along Y axis
	gr1->GetHistogram()->SetMaximum(3000);   // along Y axis
	
	int l = 0;
	gr1->SetLineColor(l+1);
	gr1->SetMarkerColor(l+1);
	gr1->Draw("AP");
	
	TGraphErrors* gr12 = new TGraphErrors(nn2, li12, li32, 0, li42);
	gr12->SetMarkerStyle(20);
	gr12->SetMarkerColor(kRed);
	gr12->Draw("P same");
	
	TLegend* lgd = new TLegend(0.3,0.2, 0.9, 0.4);
	lgd->SetTextSize(0.05);
	lgd->AddEntry(gr1, "Coarse Gain = 100", "P");
	lgd->AddEntry(gr12, "Coarse Gain = 200", "P");
	lgd->Draw();
	
	cv->cd(2);
	gPad->SetLeftMargin(0.15);
	gPad->SetGrid();
	//gPad->SetLogy();
	
	TGraphErrors* gr2 = new TGraphErrors(nn, li1, li4, 0, 0);
	gr2->SetTitle(Form("Resolution = f(#DeltaV_{mesh top})"));
	gr2->GetXaxis()->SetTitle( "#DeltaV (V)" );
	gr2->GetXaxis()->SetLimits(-10, 270);
	gr2->GetYaxis()->SetTitle( "Resolution (%)" );
	gr2->SetMarkerStyle(20);
	gr2->GetHistogram()->SetMinimum(1);   // along Y axis
	gr2->GetHistogram()->SetMaximum(50);   // along Y axis
	
	gr2->SetLineColor(l+1);
	gr2->SetMarkerColor(l+1);
	gr2->Draw("AP");
	
	TGraphErrors* gr22 = new TGraphErrors(nn2, li12, li42, 0, 0);
	gr22->SetMarkerStyle(20);
	gr22->SetMarkerColor(kRed);
	gr22->Draw("P same");
	
	cv->SaveAs(Form("Figures/%s/GainTotal-mesh-top.pdf", detectorName.c_str()));
	
	return 0;
	
}

