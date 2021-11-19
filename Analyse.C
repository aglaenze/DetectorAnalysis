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
#include "Include/RooGain.C"

using namespace std;

TObjArray* GetMatches(string detectorName, TString fileName) {
    int nElectrodes = GetElectrodeNum(detectorName);
    TString regExpr = "";
    for (int k = 0; k < nElectrodes; k++) {regExpr += "-(\\d+)";}
    //regExpr += ".*-(\\d+)(?:\\.mca|_)";
    //if (detectorName == "RD3SP1") regExpr += "-(\\d+)(_|-)(\\.mca)";
    regExpr += "-(\\d+)(_|-)?(.*)(\\.mca)";       // OK
    if (detectorName == "RD3SP1") regExpr="-(\\d+)-(\\d+)-(\\d+)-(\\d+)-(\\d+)-(\\d+)-(\\d+)(.mca)";
    TObjArray* matches = TPRegexp( regExpr ).MatchS( fileName );
    //cout << endl << endl << fileName << endl << regExpr << endl << endl;
    return matches;
}

vector<int> IdentifyHV(string detectorName, TString fileName) {
    int nElectrodes = GetElectrodeNum(detectorName);
    TObjArray* matches = GetMatches(detectorName, fileName);
    
    vector<int> hvList = {};
    for (int k = 0; k < nElectrodes; k++) {hvList.push_back( (static_cast<TObjString*>(matches->At(k+1)))->GetString().Atoi());}
    return hvList;
}

Int_t GetCoarseGain(string detectorName, TString fileName) {
    int nElectrodes = GetElectrodeNum(detectorName);
    TObjArray* matches = GetMatches(detectorName, fileName);
    
    //Int_t coarseGain = (static_cast<TObjString*>(matches->At(matches->GetLast())))->GetString().Atoi();
    Int_t coarseGain = (static_cast<TObjString*>(matches->At(nElectrodes+1)))->GetString().Atoi();
    return coarseGain;
}

Bool_t IsCollimated(string detectorName, TString fileName) {
    if (detectorName == "RD3SP1") return false;
    int nElectrodes = GetElectrodeNum(detectorName);
    TObjArray* matches = GetMatches(detectorName, fileName);
    //cout << "HERE " << (static_cast<TObjString*>(matches->At(nElectrodes+3)))->GetString() << endl;
    if ((static_cast<TObjString*>(matches->At(nElectrodes+3)))->GetString() == "collimated") return true;
    else return false;
}

string FindDate(string folder) {
    string date = "";
    TObjArray* matches = TPRegexp("(\\d+)-(\\d+)-(\\d+)(.*)").MatchS(folder);
    for (int k = 1; k < 4; k++) {
        date += (static_cast<TObjString*>(matches->At(k)))->GetString();
        if (k<3) date += "-";
    }
    return date;
}

int Analyse(string detectorName, string folder) {
    
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetLabelSize(0.05, "X");
    gStyle->SetLabelSize(0.05, "Y");
    
    string date = FindDate(folder);
        
    // variables here
    bool drawFit = true;
    bool roo = true;   // for complicated energy spectra
    int numGauss = 1;
    bool pileUp = false;
    bool escape = true;
    bool alpha = true;
    string meshName = "mesh";
    //string meshName = "mesh top";
    //string driftName = "GEM bottom";     // to read current: drift, if MM only: GEM bottom
    //string driftName = "mesh top";     // to read current: drift, if MM +GEM only: mesh top
    string driftName = "drift";     // to read current: drift
    
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
    vector<int> dvGemList = {};
    vector<int> hvMeshTopList = {};
    vector<int> hvDriftList = {};
    vector<int> hvDvList = {};
    vector<int> coarseGainList = {};
    
    vector<double> resolutionList = {};
    vector<double> resolutionErrorList = {};
    
    vector<double> gainMmList = {};
    vector<double> gainMmErrorList = {};
    
    vector<double> ibfList = {};
    vector<double> ibfErrorList = {};
    
    int nCurrentsFound = 0;
    for( Int_t iFile = 0; iFile < nFiles; ++iFile ) {
        TString file  = mcaFiles[iFile];
        cout << "Processing file " << file << endl;
        
        string suffix = "";
        bool collimated = IsCollimated(detectorName, file);
        if (collimated) {suffix= "-collimated";}
        vector<int> hvList = IdentifyHV(detectorName, file);
        Int_t coarseGain = GetCoarseGain(detectorName, file);
        
        
        int hvMesh = hvList[0];
        int hvGemDown = hvList[1];
        int hvGemTop = hvList[2];
        int hvMeshTop = hvList[3];
        int dvMeshTop = hvList[3]-hvList[2];
        int hvDrift = hvList[4];
        
        /*
         cout << "hvMesh = " << hvMesh << endl;
         cout << "hvGemDown = " << hvGemDown << endl;
         cout << "hvGemTop = " << hvGemTop << endl;
         cout << "hvMeshTop = " << hvMeshTop << endl;
         cout << "hvDrift = " << hvDrift << endl;
         cout << "coarseGain = " << coarseGain << endl;
         continue;
         */
        
        /*
        if (coarseGain != 50) {
            cout << "coarse gain = " << coarseGain << endl;
            continue;
        }
         */

        //if (hvMesh !=  320) continue;
        //if (hvGemTop !=  930) continue;
        if (hvDrift !=  1540) continue;
        //if (dvMeshTop != 130) continue;
        
        /*
         if (hvMeshTop-hvGemTop < 250) pileUp = true;
         else pileUp = false;
         */
        /*
         if (coarseGain>50) roo = false;
         else roo = true;
         */
        // Chercher si le log correspondant existe
        bool logExist = FileExist(logPath, "log", hvList, nCurrentsFound);
        bool picoExist = FileExist(picoPath, "pico", hvList, nCurrentsFound);
        
        if (logExist || picoExist) cout << "YES!!" << endl;
        //continue;
        
        // Commencer les plots ici
        
        TString figureName = Form("Figures/%s/%s/gain", detectorName.c_str(), folder.c_str());
        if (logExist || picoExist) figureName += "-ibf";
        for (int k = 0; k< (int)hvList.size(); k++) {figureName += Form("-%d", hvList[k]);}
        figureName += Form("-%d%s.pdf", coarseGain, suffix.c_str());
        cout << figureName << endl;
        
        TCanvas* cv = new TCanvas();
        //if (logExist)
        cv->Divide(3, 2);
        //else {cv->Divide(3); cv->SetCanvasSize(900, 300);}
        
        double gain = 0, gainError = 0, fwhm = 0, fwhmError = 0, resolution = 0, resolutionError = 0;
        double gainMm = 0, gainMmError = 0;
        cv->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        DrawDetector(detectorName, hvList);
        string gasName = "Ar-isobutane (95/5)";
        TText* txtGas = new TText(.3,.95,Form("Gas: %s", gasName.c_str()));
        txtGas->Draw("same");
        
        cv->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        
        // Draw gain but in log scale
        cv->cd(3);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetLogy();
        
        if (drawFit && roo) {
            cv->cd(2);
            //RooPlot* frame = DrawRooSpectrum(mcaPath, detectorName, hvList, coarseGain, gain, gainError, gainMm, gainMmError, fwhm, fwhmError, resolution, resolutionError, numGauss, pileUp, escape);
            RooPlot* frame = DrawRooSpectrum(detectorName, date, file, coarseGain, gain, gainError, gainMm, gainMmError, fwhm, fwhmError, resolution, resolutionError, numGauss, pileUp, escape, alpha);
            TLegend* lgd = GetLegend(frame, numGauss, pileUp, escape, alpha);
            frame->SetMaximum(1.2);
            //frame->GetXaxis()->SetRangeUser(2, 25);
            frame->GetXaxis()->SetRangeUser(0, 15);
            frame->Draw();
            lgd->Draw();
            cv->cd(3);
            frame->SetMinimum(1.e-3);
            frame->Draw();
            lgd->Draw();
        }
        else {
            cv->cd(2);
            //DrawSpectrum(mcaPath, detectorName, hvList, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError, drawFit);
            DrawSpectrum(detectorName, date, file, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError, drawFit);
            cv->cd(3);
            //DrawSpectrum(mcaPath, detectorName, hvList, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError, drawFit);
            DrawSpectrum(detectorName, date, file, coarseGain, gain, gainError, fwhm, fwhmError, resolution, resolutionError, drawFit);
        }
        
        // Écrire le gain + la résolution
        gStyle->SetTextSize(0.07);
        cv->cd(6);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        TLatex* txtGain = new TLatex(0.1, 0.8, Form("#bf{Gain = %.1f #pm %.1f}", gain, gainError));
        TLatex* txtResolution = new TLatex(0.1, 0.7, Form("#bf{Resolution = (%.1f #pm %.1f) %s}", resolution, resolutionError, "%"));
        TLatex* txtFwhm = new TLatex(0.1, 0.6, Form("#bf{FWHM = (%.1f #pm %.1f) %s}", fwhm, fwhmError, "%"));
        txtGain->Draw();
        txtResolution->Draw();
        txtFwhm->Draw();
        
        
        // Draw current here
        if (logExist || picoExist) {
            TString tag, currentsPath;
            if (logExist) {tag = "log"; currentsPath = logPath;}
            else {tag = "pico"; currentsPath = picoPath;}
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
            //ibfError = sqrt(Square(ibfError) + Square(gainError/gain) );
            ibfError *= 100;
            TLatex* ibfResult = new TLatex(0.1, 0.1, Form("#bf{IBF = (%.3f #pm %.3f) %s}", ibf, ibfError, "%"));
            ibfResult->Draw();
            
            ibfList.push_back(ibf);
            ibfErrorList.push_back(ibfError);
            
        }
        
        cv->SaveAs(figureName);
        
        gainList.push_back(gain);
        gainErrorList.push_back(gainError);
        hvMeshList.push_back(hvMesh);
        hvGemDownList.push_back(hvGemDown);
        hvGemTopList.push_back(hvGemTop);
        dvGemList.push_back(hvGemTop-hvGemDown);
        hvMeshTopList.push_back(hvMeshTop);
        hvDriftList.push_back(hvDrift);
        hvDvList.push_back(dvMeshTop);
        coarseGainList.push_back(coarseGain);
        
        resolutionList.push_back(resolution);
        resolutionErrorList.push_back(resolutionError);
        
        gainMmList.push_back(gainMm);
        gainMmErrorList.push_back(gainMmError);
    }
    
    PrintList( "gainList", gainList, "%.3f" );
    PrintList( "gainErrorList", gainErrorList, "%.3f" );
    PrintList( "hvMeshList", hvMeshList, "%d" );
    PrintList( "hvGemDownList", hvGemDownList, "%d" );
    PrintList( "hvGemTopList", hvGemTopList, "%d" );
    PrintList( "dvGemList", dvGemList, "%d" );
    PrintList( "hvMeshTopList", hvMeshTopList, "%d" );
    PrintList( "hvDriftList", hvDriftList, "%d" );
    PrintList( "dvMeshTopList", hvDvList, "%d" );
    PrintList( "coarseGainList", coarseGainList, "%d" );
    cout << endl;
    PrintList( "gainMmList", gainMmList, "%.3f" );
    PrintList( "gainMmErrorList", gainMmErrorList, "%.3f" );
    cout << endl;
    PrintList( "resolutionList", resolutionList, "%.2f" );
    PrintList( "resolutionErrorList", resolutionErrorList, "%.2f" );
    
    cout << endl;
    PrintList( "ibfList", ibfList, "%.3f" );
    PrintList( "ibfErrorList", ibfErrorList, "%.3f" );
    
    //for (int i = 0; i < ibfList.size(); i++) cout << ibfList[i] << endl;
    
    cout << "\n\nNumber of current files found corresponding to MCA files = " << nCurrentsFound << endl;
    return 0;
}


