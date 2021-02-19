#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "_Utils.C"
#include "_HistUtils.C"

using namespace std;


//____________________________________________
int DrawTotalGain( string detectorName = "MGEM3" )
{
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetLabelSize(0.05, "X");
	gStyle->SetLabelSize(0.05, "Y");
	
	
	// Gain
	const int num1 = 2;
	const int num2 = 4;
	int dV1[num2] = {250, 280, 310, 340}; 	// in GEM
	int dV2[num1] = {310, 340}; 	// GEM and mesh top
	
	Double_t gainList[11] = {750.006, 1581.879, 1881.088, 545.404, 1182.419, 578.932, 1245.877, 1493.786, 607.657, 1693.861, 2004.780};
	Double_t gainErrorList[11] = {89.617, 193.627, 171.350, 111.128, 169.138, 110.095, 144.929, 193.373, 91.624, 187.809, 234.869};
	Double_t hvMeshList[11] = {390, 390, 390, 390, 390, 390, 390, 390, 390, 390, 390};
	Double_t hvGemDownList[11] = {600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600};
	Double_t hvGemTopList[11] = {850, 850, 850, 880, 880, 910, 910, 910, 940, 940, 940};
	Double_t hvMeshTopList[11] = {1180, 1200, 1200, 1190, 1220, 1220, 1250, 1250, 1250, 1280, 1280};
	Double_t dvMeshTopList[11] = {330, 350, 350, 310, 340, 310, 340, 340, 310, 340, 340};
	Double_t coarseGainList[11] = {100, 100, 50, 100, 100, 100, 100, 50, 100, 100, 50};
	
	TCanvas* cv = new TCanvas();
	gPad->SetLeftMargin(0.15);
	gPad->SetGrid();
	gPad->SetLogy();
	
	TLegend* legend1 = new TLegend(0.55,0.1,0.9,0.3);
	legend1->SetTextSize(0.045);
	
	vector <TGraphErrors*> graphVec = {};
	
	for (int i = 0; i< num1; i++) {
		//int i = 0;
		
		vector<float> gainVec = {};
		vector<float> gainErrorVec= {};
		vector<int> dvMeshTopVec = {};
		vector<int> dvGem = {};
		
		for (int k = 0; k < 11; k++) {
			//if ( (hvGemTopList[k] - hvGemDownList[k] == dV1[i])) {
			if ( (hvMeshTopList[k] - hvGemTopList[k] == dV2[i])) {
				gainVec.push_back(gainList[k]);
				gainErrorVec.push_back(gainErrorList[k]);
				dvMeshTopVec.push_back(dvMeshTopList[k]);
				dvGem.push_back(hvGemTopList[k] - hvGemDownList[k]);
			}
		}
		
		// create list from vector
		const int nn = (int)dvMeshTopVec.size();
		//TGraphErrors* gr1 = CreateTGraphError(hvMeshTopVec[l], gainVec[l], gainErrorVec[l]);
		Double_t li1[nn], li2[nn], li3[nn], li4[nn], li5[nn];
		for (int k = 0; k< nn; k++) {
			//li1[k] = dvMeshTopVec[k];
			li1[k] = hvMeshTopList[k] - hvGemTopList[k];
			li3[k] = gainVec[k];
			li4[k] = gainErrorVec[k];
			//cout << "ok" << endl;
		}
		for (int k = 0; k< nn; k++) { cout << li1[k] << " " << li3[k] << " " << li4[k] << endl;}
		TGraphErrors* gr1 = new TGraphErrors(nn, li1, li3, 0, li4);
		gr1->SetTitle("Gain = f(HV_{GEM})");
		gr1->GetXaxis()->SetTitle( "HV (V)" );
		gr1->GetXaxis()->SetLimits(300, 360);
		gr1->GetYaxis()->SetTitle( "Gain" );
		gr1->SetMarkerStyle(20);
		gr1->GetHistogram()->SetMinimum(90);   // along Y axis
		gr1->GetHistogram()->SetMaximum(11000);   // along Y axis
		
		legend1->AddEntry(gr1, Form("dV_{mesh top} = %d V", dV2[i]), "lp");
		
		gr1->SetLineColor(i+1);
		gr1->SetMarkerColor(i+1);

		if (i == 0) gr1->Draw("ALP");
		else gr1->Draw("LP");

	}
	
	legend1->Draw();


	
	cv->SaveAs(Form("Figures/%s/GainTotal-test2.pdf", detectorName.c_str()));
	
	return 0;
	
}

