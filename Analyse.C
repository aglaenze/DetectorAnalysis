#include <iostream>
#include <fstream>
#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include <TCanvas.h>
#include <TROOT.h>
#include <TF1.h>
#include <TF2.h>
#include <TObjString.h>
#include <TPRegexp.h>
#include <TFile.h>
#include <TLegend.h>
#include <TMath.h>

#include "_DetectorInfo.C"
#include "_Utils.C"
#include "Gain.C"
#include "Currents.C"

using namespace std;

vector<int> IdentifyHV(string detectorName, TString fileName) {
	int nElectrodes = GetElectrodeNum(detectorName);
	
	TString regExpr = "";
	for (int k = 0; k < nElectrodes; k++) {regExpr += "-(\\d+)";}
	regExpr += ".*-(\\d+)(?:\\.mca|_)";
	TObjArray* matches = TPRegexp( regExpr ).MatchS( fileName );
	
	vector<int> hvList = {};
	for (int k = 0; k < nElectrodes; k++) {hvList.push_back( (static_cast<TObjString*>(matches->At(k+1)))->GetString().Atoi());}
	return hvList;
}

Int_t GetCoarseGain(string detectorName, TString fileName) {
	int nElectrodes = GetElectrodeNum(detectorName);
	//cout << "There are " << nElectrodes << " electrodes" << endl;
	
	TString regExpr = "";
	for (int k = 0; k < nElectrodes; k++) {regExpr += "-(\\d+)";}
	regExpr += ".*-(\\d+)(?:\\.mca|_)";
	TObjArray* matches = TPRegexp( regExpr ).MatchS( fileName );
	
	Int_t coarseGain = (static_cast<TObjString*>(matches->At(matches->GetLast())))->GetString().Atoi();
	return coarseGain;
}

int Analyse(string detectorName, string folder) {
	
	TString mcaPath = Form("MCA/%s/%s/", detectorName.c_str(), folder.c_str());
	TString logPath = Form("logFiles/%s/%s/", detectorName.c_str(), folder.c_str());
	
	vector<TString> mcaFiles;
	struct dirent **namelist;
	Int_t n = scandir(mcaPath, &namelist, 0, alphasort);
	//cout << mcaPath << endl;
	if (n < 1) {cout << "empty folder" << endl; return 0;}
	else {
		while (n--) {
			if (strstr(namelist[n]->d_name, "mca") != NULL) {
				mcaFiles.insert(mcaFiles.begin(), mcaPath + namelist[n]->d_name);
			}
		}
	}
	mcaFiles.push_back(TString());
	Int_t nFiles = GetNumberOfFiles(mcaPath, "mca");
	cout << "Number of files: " << nFiles << endl;
	
	int nLogFound = 0;
	for( Int_t iFile = 0; iFile < nFiles; ++iFile ) {
		TString file  = mcaFiles[iFile];
		cout << "Processing file " << file << endl;
		
		Int_t coarseGain = GetCoarseGain(detectorName, file);
		vector<int> hvList = IdentifyHV(detectorName, file);
		
		// Chercher si le log correspondant existe
		bool logExist = false;
		struct dirent **namelist2;
		Int_t n2 = scandir(logPath, &namelist2, 0, alphasort);
		
		TString logName = "log";
		for (int k = 0; k< (int)hvList.size(); k++) {logName += Form("-%d", hvList[k]);}
		if (n2 > 0) {
			while (n2--) {
				if (strstr(namelist2[n2]->d_name, logName) != NULL) {
					logExist = true;
					nLogFound++;
					break;
				}
			}
		}
		if (logExist) cout << "YES!!" << endl;
		
		// Commencer les plots ici
		
		TString figureName = Form("Figures/%s/%s/gain", detectorName.c_str(), folder.c_str());
		if (logExist) figureName += "-ibf";
		for (int k = 0; k< (int)hvList.size(); k++) {figureName += Form("-%d", hvList[k]);}
		figureName += ".pdf";
		
		TCanvas* cv = new TCanvas();
		if (logExist) cv->Divide(3, 2);
		else {cv->Divide(3); cv->SetCanvasSize(900, 300);}
		
		double gain = 0, gainError = 0, fwhm = 0, fwhmError = 0, resolution = 0, resolutionError = 0;
		cv->cd(1);
		//gPad->SetLeftMargin(0.15);
		DrawDetector(detectorName, hvList);
		
		cv->cd(2);
		DrawSpectrum(mcaPath, detectorName, hvList, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError);
		
		// Écrire le gain + la résolution
		gStyle->SetTextSize(0.07);
		cv->cd(3);
		TLatex* txtGain = new TLatex(0.1, 0.7, Form("Gain = %.1f #pm %.1f", gain, gainError));
		TLatex* txtResolution = new TLatex(0.1, 0.6, Form("Resolution = %.1f #pm %.1f %s", resolution, resolutionError, "%"));
		TLatex* txtFwhm = new TLatex(0.1, 0.5, Form("FWHM = %.1f #pm %.1f %s", fwhm, fwhmError, "%"));
		txtGain->Draw();
		txtResolution->Draw();
		txtFwhm->Draw();
		
		
		// Draw current here
		if (logExist) {
			cv->cd(4);
			cout << "hello" << endl;
			DrawCurrents(logPath, hvList, detectorName, "mesh");
			cv->cd(5);
			DrawCurrents(logPath, hvList, detectorName, "drift");
			
			// Écrire l'IBF
			cv->cd(6);
			TLatex* ibfDef = new TLatex(0.1, 0.7, "#bf{IBF = #||{#frac{ i_{drift}^{with} - i_{drift}^{without}}{ i_{mesh}^{with} - i_{mesh}^{without}}}}");
			ibfDef->Draw();
			
			double ibf = 0, ibfError = 0;
			ibf *= 100;
			ibfError *= 100;
			ComputeIbf(logPath, hvList, detectorName, ibf, ibfError);
			TLatex* ibfResult = new TLatex(0.1, 0.5, Form("#bf{IBF = %.3f #pm %.3f %s}", ibf, ibfError, "%"));
			ibfResult->Draw();
			
		}
		
		cv->SaveAs(figureName);
	}
	
	cout << "\n\nNumber of log files found corresponding to MCA files = " << nLogFound << endl;
	return 0;
}


/*
 const TString path = "MCA/";
 const TString folder = "February_05_2020";
 
 // conversion factor (from fit to calibrated data)
 const Double_t nPrimary = 228.403; // Ar-iC4H10 95/5
 Double_t alpha = 5.45217e-06*nPrimary;  // Calibration de l'ORTEC
 
 const TString completePath = path+folder+"/";
 
 std::vector<TString> files;
 struct dirent **namelist;
 Int_t n = scandir(completePath, &namelist, 0, alphasort);
 if (n < 1) {std::cout << "empty folder" << std::endl; return 0;}
 else {
 while (n--) {
 if (strstr(namelist[n]->d_name, "mca") != NULL) files.insert(files.begin(), completePath + namelist[n]->d_name);
 //std::cout << path + namelist[n]->d_name << std::endl;
 }
 }
 //std::cout << "done" << std::endl;
 files.push_back(TString());
 
 // count files
 Int_t nFiles = 0;
 for( ; !files[nFiles].IsNull(); ++nFiles )
 {}
 
 std::cout << "DrawFeSignal - nFiles: " << nFiles << std::endl;
 */
