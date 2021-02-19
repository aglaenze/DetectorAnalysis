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

#include "Include/DetectorInfo.C"
#include "Include/Utils.C"
#include "Include/Gain.C"
#include "Include/Currents.C"

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
	
	gStyle->SetTitleFontSize(.05);
	gStyle->SetTitleXSize(.05);
	gStyle->SetTitleYSize(.05);
	gStyle->SetLabelSize(0.05, "X");
	gStyle->SetLabelSize(0.05, "Y");
	
	TString mcaPath = Form("MCA/%s/%s/", detectorName.c_str(), folder.c_str());
	TString logPath = Form("logFiles/%s/%s/", detectorName.c_str(), folder.c_str());
    TString picoPath = Form("/Users/aglaenzer/Softwares/pico-env/picopy/data/%s/%s/", detectorName.c_str(), folder.c_str());
	
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
	
	vector<double> gainList = {};
	vector<double> gainErrorList = {};
	vector<int> hvMeshList = {};
	vector<int> hvGemDownList = {};
	vector<int> hvGemTopList = {};
	vector<int> hvMeshTopList = {};
	vector<int> hvDvList = {};
	vector<int> coarseGainList = {};
	
	int nCurrentsFound = 0;
	for( Int_t iFile = 0; iFile < nFiles; ++iFile ) {
		TString file  = mcaFiles[iFile];
		cout << "Processing file " << file << endl;
		
		Int_t coarseGain = GetCoarseGain(detectorName, file);
        //if (coarseGain != 50) continue;
		vector<int> hvList = IdentifyHV(detectorName, file);
		
		int hvMesh = hvList[0];
		int hvGemDown = hvList[1];
		int hvGemTop = hvList[2];
		int hvMeshTop = hvList[3];
		int dVMeshTop = hvList[3]-hvList[2];
		
		// Chercher si le log correspondant existe
        bool logExist = FileExist(logPath, "log", hvList, nCurrentsFound);
        bool picoExist = FileExist(picoPath, "pico", hvList, nCurrentsFound);
        
		if (logExist || picoExist) cout << "YES!!" << endl;
		
		// Commencer les plots ici
		
		TString figureName = Form("Figures/%s/%s/gain", detectorName.c_str(), folder.c_str());
		if (logExist || picoExist) figureName += "-ibf";
		for (int k = 0; k< (int)hvList.size(); k++) {figureName += Form("-%d", hvList[k]);}
		figureName += Form("-%d.pdf", coarseGain);
		
		TCanvas* cv = new TCanvas();
		//if (logExist)
		cv->Divide(3, 2);
		//else {cv->Divide(3); cv->SetCanvasSize(900, 300);}
		
		double gain = 0, gainError = 0, fwhm = 0, fwhmError = 0, resolution = 0, resolutionError = 0;
		cv->cd(1);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		DrawDetector(detectorName, hvList);
		
		cv->cd(2);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		DrawSpectrum(mcaPath, detectorName, hvList, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError);
		
		// Draw gain but in log scale
		cv->cd(3);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetLogy();
		DrawSpectrum(mcaPath, detectorName, hvList, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError);
		
		// Écrire le gain + la résolution
		gStyle->SetTextSize(0.07);
		cv->cd(6);
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		TLatex* txtGain = new TLatex(0.1, 0.8, Form("#bf{Gain = %.1f #pm %.1f}", gain, gainError));
		TLatex* txtResolution = new TLatex(0.1, 0.7, Form("#bf{Resolution = %.1f #pm %.1f %s}", resolution, resolutionError, "%"));
		TLatex* txtFwhm = new TLatex(0.1, 0.6, Form("#bf{FWHM = %.1f #pm %.1f %s}", fwhm, fwhmError, "%"));
		txtGain->Draw();
		txtResolution->Draw();
		txtFwhm->Draw();
		
		
		// Draw current here
		if (logExist || picoExist) {
            TString tag, currentsPath;
            if (logExist) {tag = "log"; currentsPath = logPath;}
            else {tag = "pico"; currentsPath = picoPath;}
            string meshName = "mesh";
            string driftName = "drift";
            //string driftName = "GEM bottom";
			cv->cd(4);
			gPad->SetLeftMargin(0.15);
			gPad->SetBottomMargin(0.15);
			cout << "hello" << endl;
			DrawCurrents(currentsPath, tag, hvList, detectorName, meshName);
			cv->cd(5);
			gPad->SetLeftMargin(0.15);
			gPad->SetBottomMargin(0.15);
            DrawCurrents(currentsPath, tag, hvList, detectorName, driftName);
			
			// Écrire l'IBF
			cv->cd(6);
			gPad->SetLeftMargin(0.15);
			gPad->SetBottomMargin(0.15);
			TLatex* ibfDef = new TLatex(0.1, 0.3, "#bf{IBF = #||{#frac{ i_{drift}^{with} - i_{drift}^{without}}{ i_{mesh}^{with} - i_{mesh}^{without}}} - #frac{1}{gain}}");
			ibfDef->Draw();
			
			double ibf = 0, ibfError = 0;
            ComputeIbf(currentsPath, tag, hvList, detectorName, meshName, driftName, ibf, ibfError);
			ibf -= 1/gain;
			ibf *= 100;
			ibfError *= 100;
			TLatex* ibfResult = new TLatex(0.1, 0.1, Form("#bf{IBF = %.3f #pm %.3f %s}", ibf, ibfError, "%"));
			ibfResult->Draw();
			
		}
		
		cv->SaveAs(figureName);
		
		gainList.push_back(gain);
		gainErrorList.push_back(gainError);
		hvMeshList.push_back(hvMesh);
		hvGemDownList.push_back(hvGemDown);
		hvGemTopList.push_back(hvGemTop);
		hvMeshTopList.push_back(hvMeshTop);
		hvDvList.push_back(dVMeshTop);
		coarseGainList.push_back(coarseGain);
	}
	
	PrintList( "gainList", gainList, "%.3f" );
	PrintList( "gainErrorList", gainErrorList, "%.3f" );
	PrintList( "hvMeshList", hvMeshList, "%d" );
	PrintList( "hvGemDownList", hvGemDownList, "%d" );
	PrintList( "hvGemTopList", hvGemTopList, "%d" );
	PrintList( "hvMeshTopList", hvMeshTopList, "%d" );
	PrintList( "dvMeshTopList", hvDvList, "%d" );
	PrintList( "coarseGainList", coarseGainList, "%d" );
	
	cout << "\n\nNumber of current files found corresponding to MCA files = " << nCurrentsFound << endl;
	return 0;
}


