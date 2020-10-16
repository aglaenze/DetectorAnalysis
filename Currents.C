#include "TNtuple.h"
#include "TString.h"
#include "TDatime.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <dirent.h>

#include <TF1.h>
#include <TH1.h>


using namespace std;

TH1F* GetCurrentHist(TString fileName, string electrodeName, string detectorName) {
	
	TH1F* hImon = new TH1F("hImon", Form("Currents on %s", electrodeName.c_str()), 4000, -100, 100);
	
	map <string, int, NoSorting> electrodeMap;
	int ok = LoadElectrodeMap(detectorName, electrodeMap);
	if (ok != 0) {cout << "\n\n\nMap not loaded\n\n\n" << endl; return hImon;}
	
	int elIdx = -1;
	for(auto it = electrodeMap.cbegin(); it != electrodeMap.cend(); ++it) {
		if (it->first == electrodeName) {
			//std::cout << "\n\n\n" << it->first << " " << it->second << "\n";
			elIdx = it->second;
		}
	}
	if (elIdx == -1) {cout << "\n\n\nElectrode " << electrodeName << " does not exist\n\n\n" << endl; return hImon;}
	
	Double_t imon[13];
	time_t ftime, time = 0, time0 = 0;
	ifstream myfile(fileName.Data());
	if (myfile) {
	string line;
	int nlines = 0;
	while (myfile.is_open() && myfile.good() ) {
		getline (myfile,line); nlines++;
		if(!(nlines%2)) continue; //avoid errors in the data stream
								  //	8/29/2019 6:30:38 PM
		int day = 0, month = 0, year = 0; int hour = 0, min = 0, sec = 0;
		char ampm[20];
		int ncols = sscanf(line.c_str(),"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d/%d/%d\t%d:%d:%d %s"
						   ,&imon[1], &imon[2], &imon[3], &imon[4], &imon[5], &imon[6], &imon[7], &imon[8], &imon[9], &imon[10], &imon[11], &imon[12]
						   ,&month, &day, &year, &hour, &min, &sec, ampm);
		if (ncols != 19) continue;
		
		if(strcmp(ampm,"PM")==0) hour+=12;
		TDatime date(year, month, day, hour, min, sec);
		if (time0 == 0) time0 = date.Convert();
		ftime = date.Convert() - time0;

		hImon->Fill(imon[elIdx]);
	}
	hImon->SetXTitle(Form("i_{%s} [nA]", electrodeName.c_str()));
	while (hImon->GetMaximum()<30) hImon->Rebin(2);
	hImon->Scale(1./hImon->GetMaximum());
	hImon->SetMaximum(1.1);
	}
	else {cout << "\n\nImpossible d'ouvrir le fichier log\n\n" << endl;}
	return hImon;
}

TF1* GetFitCurve(TH1F* h) {
	//Set up fit
	Int_t iBinMax = h->GetMaximumBin();
	Double_t xMax = h->GetXaxis()->GetBinCenter(iBinMax);

	Double_t fitRangeMin = xMax - 2 * h->GetRMS() ;
	Double_t fitRangeMax = xMax + 2.5 * h->GetRMS() ;
	/*
	if (h->GetRMS() > 0.5) {fitRangeMin = xMax - 20; fitRangeMax = xMax + 0.5;}
	 */
	TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
	f->SetParNames("Mean", "Sigma", "Amplitude");
	f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
	
	h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
	return f;
}


void DrawCurrents(TString path, vector<int> hvList, string detectorName, string electrodeName) {
	
	/* With draw on one plot currents with and without a Fe source */
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);


	TString fileNameWith = path+ "log";
	TString fileNameWithout = path+ "log";
	for (int k = 0; k<(int) hvList.size(); k++) {
		fileNameWith += Form("-%d", hvList[k]);
		fileNameWithout += Form("-%d", hvList[k]);
	}
	fileNameWith += "-with.lvm";
	fileNameWithout += "-without.lvm";
	
	TH1F* hImonWith = GetCurrentHist(fileNameWith, electrodeName, detectorName);
	TH1F* hImonWithout = GetCurrentHist(fileNameWithout, electrodeName, detectorName);
	
	Double_t xMin = min(hImonWith->GetMean()-6*hImonWith->GetRMS(), hImonWithout->GetMean()-6*hImonWithout->GetRMS());
	Double_t xMax = max(hImonWith->GetMean()+8*hImonWith->GetRMS(), hImonWithout->GetMean()+8*hImonWithout->GetRMS());

	hImonWith->GetXaxis()->SetRangeUser(xMin, xMax);
	//hImonWith->GetXaxis()->SetRangeUser(-1, 1);
	hImonWith->SetLineColor(kBlue);
	hImonWithout->SetLineColor(kBlack);
	
	hImonWith->Draw("hist");
	hImonWithout->Draw("same hist");
	
	// Draw fit curves
	TF1* fIwith = GetFitCurve(hImonWith);
	fIwith->SetLineColor(kBlue);
	TF1* fIwithout = GetFitCurve(hImonWithout);
	fIwithout->SetLineColor(kBlack);
	fIwith->Draw("same");
	fIwithout->Draw("same");
	
	// Draw legend
	TLegend* lgd = new TLegend(0.5, 0.7, 0.9, 0.9);
	lgd->SetTextSize(0.05);
	lgd->AddEntry(hImonWith, "with source", "l");
	lgd->AddEntry(hImonWithout, "without source", "l");
	lgd->Draw();
}

void ComputeIbf(TString path, vector<int> hvList, string detectorName, double& ibf, double& ibfError) {
	
	TString fileNameWith = path+ "log";
	TString fileNameWithout = path+ "log";
	for (int k = 0; k<(int) hvList.size(); k++) {
		fileNameWith += Form("-%d", hvList[k]);
		fileNameWithout += Form("-%d", hvList[k]);
	}
	fileNameWith += "-with.lvm";
	fileNameWithout += "-without.lvm";
	
	TH1F* hImonMeshWith = GetCurrentHist(fileNameWith, "mesh", detectorName);
	TH1F* hImonMeshWithout = GetCurrentHist(fileNameWithout, "mesh", detectorName);
	
	TH1F* hImonDriftWith = GetCurrentHist(fileNameWith, "drift", detectorName);
	TH1F* hImonDriftWithout = GetCurrentHist(fileNameWithout, "drift", detectorName);

	// Fit curves
	TF1* fIMeshWith = GetFitCurve(hImonMeshWith);
	TF1* fIMeshWithout = GetFitCurve(hImonMeshWithout);
	
	TF1* fIDriftWith = GetFitCurve(hImonDriftWith);
	TF1* fIDriftWithout = GetFitCurve(hImonDriftWithout);
	
	double iMeshWith = fIMeshWith->GetParameter(0);
	double iMeshWithout = fIMeshWithout->GetParameter(0);
	double iDriftWith = fIDriftWith->GetParameter(0);
	double iDriftWithout = fIDriftWithout->GetParameter(0);
	
	double iErrorMeshWith = fIMeshWith->GetParError(0);
	double iErrorMeshWithout = fIMeshWithout->GetParError(0);
	double iErrorDriftWith = fIDriftWith->GetParError(0);
	double iErrorDriftWithout = fIDriftWithout->GetParError(0);
	
	
	ibf = abs( (iDriftWith-iDriftWithout)/(iMeshWith-iMeshWithout) );
	ibfError = ibf * sqrt( Square(iErrorMeshWith/iMeshWith) + Square(iErrorMeshWithout/iMeshWithout) + Square(iErrorDriftWith/iDriftWith) + Square(iErrorDriftWithout/iDriftWithout));
	
	//cout << endl << endl << endl << iMeshWith-iMeshWithout << endl;
	//cout << endl << endl << endl << iDriftWith-iDriftWithout << endl;

	return;
}
