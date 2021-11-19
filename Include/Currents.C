#include "TNtuple.h"
#include "TString.h"
#include "TDatime.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <dirent.h>

#include <TF1.h>
#include <TH1.h>
#include <TStyle.h>
#include <TROOT.h>


using namespace std;

TH1F* GetCurrentHist(TString fileName, string electrodeName, string detectorName) {
    
    //TH1F* hImon = new TH1F("hImon", Form("Currents on %s", electrodeName.c_str()), 4000, -100, 100);
    TH1F* hImon = new TH1F("hImon", Form("Currents on %s", electrodeName.c_str()), 20000, -500, 500);
    
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
            //    8/29/2019 6:30:38 PM
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
        //while (hImon->GetMaximum()<30) hImon->Rebin(2);
        //while (hImon->GetMaximum()<10) hImon->Rebin(2);
        hImon->Scale(1./hImon->GetMaximum());
        hImon->SetMaximum(1.25);
    }
    else {cout << "\n\nImpossible d'ouvrir le fichier log\n\n" << endl;}
    return hImon;
}

TH1F* GetCurrentHistCsv(TString fileName, string electrodeName, string detectorName) {
    
    //TH1F* hImon = new TH1F("hImon", Form("Currents on %s", electrodeName.c_str()), 4000, -100, 100);
    TH1F* hImon = new TH1F("hImon", Form("Currents on %s", electrodeName.c_str()), 20000, -500, 500);
    
    /* Load electrode map */
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
    
    /* start reading currents file */
    Double_t imon[13];
    time_t ftime, time = 0, time0 = 0;
    ifstream myfile(fileName.Data());
    if (myfile) {
        string line;
        int nlines = 0;
        getline (myfile,line);  // first line is just number
        while (myfile.is_open() && myfile.good() ) {
            getline (myfile,line);
            nlines++;
            int day = 0, month = 0, year = 0; int hour = 0, min = 0;
            double sec = 0;
            int idk1 = 3, idk2 = 2;
            
            //      2021-02-17 13:32:40.373927+00:00
            int ncols = sscanf(line.c_str(),"%d-%d-%d %d:%d:%lf+%d:%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                               &year, &month, &day, &hour, &min, &sec, &idk1, &idk2, &imon[1], &imon[2], &imon[3], &imon[4], &imon[5], &imon[6], &imon[7], &imon[8], &imon[9], &imon[10], &imon[11], &imon[12]);
            
            if (ncols != 20) continue;
            
            TDatime date(year, month, day, hour, min, sec);
            if (time0 == 0) time0 = date.Convert();
            ftime = date.Convert() - time0;
            //if (imon[elIdx] < 1000 && imon[elIdx] > - 1000)
            //hImon->Fill(imon[elIdx]+0.1);
            hImon->Fill(imon[elIdx]);
        }
        hImon->SetXTitle(Form("i_{%s} [nA]", electrodeName.c_str()));
        //while (hImon->GetMaximum()<30*hImon->GetRMS()) hImon->Rebin(2);
        //while (hImon->GetMaximum()<120) hImon->Rebin(2);
        hImon->Scale(1./hImon->GetMaximum());
        hImon->SetMaximum(1.4);
    }
    else {cout << "\n\nImpossible d'ouvrir le fichier csv\n\n" << endl;}
    return hImon;
}


TF1* GetFitCurve(TH1F* h, string electrodeName) {
    //Set up fit
    Int_t iBinMax = h->GetMaximumBin();
    //Int_t iBinMax = GetMaximumBin(h, 9000);
    Double_t xMax = h->GetXaxis()->GetBinCenter(iBinMax);
    
    cout << "xMax = " << xMax << endl;
    
    
    double prop = 0.0;
    if (h->GetRMS() > 5) prop = 0.5;
    else if (h->GetRMS() > 2) prop = 0.6;
    else if (h->GetRMS() > 1) prop = 1;
    else if (h->GetRMS() > 0.5) prop = 1.5;
    else prop = 0.3;
    if (abs(xMax) < 0.3 && h->GetRMS() < 1) prop = 0.3;
    //prop = 0.7;
    //if (h->GetRMS() < 1) prop = 2;
    prop = 1;
    //if (electrodeName == "drift") prop = 0.1;
    Double_t fitRangeMin = xMax - prop*h->GetRMS() ;
    Double_t fitRangeMax = xMax + prop*h->GetRMS() ;
    if (electrodeName == "drift") {
        fitRangeMin = xMax - 0.3 ;
        fitRangeMax = xMax + 0.3 ;
    }
    
    /*
     Double_t fitRangeMin = xMax - 0.1*h->GetRMS() ;
     Double_t fitRangeMax = xMax + 0.2*h->GetRMS() ;
     */
    
    cout << endl << "h->GetRMS() = " << h->GetRMS() << endl << endl;
    
    //if ( (h->GetRMS() > 10 || h->GetRMS() < 1) && abs(xMax) < 1.2) {fitRangeMin = xMax - 0.5; fitRangeMax = xMax + 0.5;}
    
    TF1* f = new TF1( "FitFunction", FitGauss, fitRangeMin, fitRangeMax, 3);
    f->SetParNames("Mean", "Sigma", "Amplitude");
    //f->SetParameters(xMax, h->GetRMS(), h->GetMaximum());
    f->SetParameters(xMax, h->GetRMS(), 1.);
    /*
    f->SetParLimits(3, 0.3, 10);
    if (electrodeName == "drift") {
        f->SetParLimits(1, -10, 1);
        f->SetParLimits(2, 0.0001, 3);
    }
    else {
        f->SetParLimits(1, -200, 20);
        f->SetParLimits(2, 0.002, 5);
    }
     */
    
    //if (xMax < -110) f->FixParameter(0, xMax);
    
    h->Fit(f, "0", "0", fitRangeMin, fitRangeMax);
    return f;
}


void DrawCurrents(TString path, TString tag, vector<int> hvList, string detectorName, string electrodeName) {
    
    /* With draw on one plot currents with and without a Fe source */
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    TH1F* hImonWith = nullptr, *hImonWithout = nullptr;
    
    TString fileNameWith = path+ tag;
    TString fileNameWithout = path+ tag;
    for (int k = 0; k<(int) hvList.size(); k++) {
        fileNameWith += Form("-%d", hvList[k]);
        fileNameWithout += Form("-%d", hvList[k]);
    }
    if (tag == "log") {
        fileNameWith += "-with.lvm";
        fileNameWithout += "-without.lvm";
        hImonWith = GetCurrentHist(fileNameWith, electrodeName, detectorName);
        hImonWithout = GetCurrentHist(fileNameWithout, electrodeName, detectorName);
    }
    else if (tag == "pico") {
        fileNameWith += "-with.csv";
        fileNameWithout += "-without.csv";
        hImonWith = GetCurrentHistCsv(fileNameWith, electrodeName, detectorName);
        hImonWithout = GetCurrentHistCsv(fileNameWithout, electrodeName, detectorName);
    }
    else {cout << "Tag " << tag << " does not exist" << endl; return;}
    
    
    // Get fit curves
    TF1* fIwith = GetFitCurve(hImonWith, electrodeName);
    TF1* fIwithout = GetFitCurve(hImonWithout, electrodeName);
    
    // Draw
    /*
    Double_t xMin = hImonWith->GetMean() - 1*hImonWith->GetRMS();
    Double_t xMax = hImonWith->GetMean() + 1*hImonWith->GetRMS();
    xMin = -135.;
    xMax = -125.;
     */
    
    Double_t prop = 10;
    if (electrodeName == "drift") prop = 3;
     Double_t xMin = min(fIwith->GetParameter(0) - prop*abs(fIwith->GetParameter(1)), fIwithout->GetParameter(0) - prop*abs(fIwithout->GetParameter(1)));
     Double_t xMax = max(fIwith->GetParameter(0) + prop*abs(fIwith->GetParameter(1)), fIwithout->GetParameter(0) + prop*abs(fIwithout->GetParameter(1)));

    hImonWith->GetXaxis()->SetRangeUser(xMin, xMax);
    //hImonWith->GetXaxis()->SetRangeUser(-3, 1);
    hImonWith->SetLineColor(kBlue);
    hImonWithout->SetLineColor(kBlack);
    
    hImonWith->Draw("hist");
    hImonWithout->Draw("same hist");
    
    fIwith->SetLineColor(kBlue);
    fIwithout->SetLineColor(kBlack);
    fIwith->Draw("same");
    fIwithout->Draw("same");
    
    // Draw legend
    TLegend* lgd = new TLegend(0.5, 0.75, 0.9, 0.9);
    lgd->SetTextSize(0.05);
    lgd->AddEntry(hImonWith, "with source", "l");
    lgd->AddEntry(hImonWithout, "without source", "l");
    lgd->Draw();
}

void ComputeIbf(TString path, TString tag, vector<int> hvList, string detectorName, string electrodeName, string driftName, double& ibf, double& ibfError) {
    
    TString fileNameWith = path+ tag;
    TString fileNameWithout = path+ tag;
    for (int k = 0; k<(int) hvList.size(); k++) {
        fileNameWith += Form("-%d", hvList[k]);
        fileNameWithout += Form("-%d", hvList[k]);
    }
    
    TH1F* hImonMeshWith = nullptr, *hImonMeshWithout = nullptr, *hImonDriftWith = nullptr, *hImonDriftWithout = nullptr;
    
    if (tag == "log") {
        fileNameWith += "-with.lvm";
        fileNameWithout += "-without.lvm";
        hImonMeshWith = GetCurrentHist(fileNameWith, electrodeName, detectorName);
        hImonMeshWithout = GetCurrentHist(fileNameWithout, electrodeName, detectorName);
        hImonDriftWith = GetCurrentHist(fileNameWith, driftName, detectorName);
        hImonDriftWithout = GetCurrentHist(fileNameWithout, driftName, detectorName);
    }
    else if (tag == "pico") {
        fileNameWith += "-with.csv";
        fileNameWithout += "-without.csv";
        hImonMeshWith = GetCurrentHistCsv(fileNameWith, electrodeName, detectorName);
        hImonMeshWithout = GetCurrentHistCsv(fileNameWithout, electrodeName, detectorName);
        hImonDriftWith = GetCurrentHistCsv(fileNameWith, driftName, detectorName);
        hImonDriftWithout = GetCurrentHistCsv(fileNameWithout, driftName, detectorName);
    }
    else {cout << "Tag " << tag << " does not exist" << endl; return;}
    
    // Fit curves
    TF1* fIMeshWith = GetFitCurve(hImonMeshWith, electrodeName);
    TF1* fIMeshWithout = GetFitCurve(hImonMeshWithout, electrodeName);
    
    TF1* fIDriftWith = GetFitCurve(hImonDriftWith, electrodeName);
    TF1* fIDriftWithout = GetFitCurve(hImonDriftWithout, electrodeName);
    
    double iMeshWith = fIMeshWith->GetParameter(0);
    double iMeshWithout = fIMeshWithout->GetParameter(0);
    
    double iErrorMeshWith = fIMeshWith->GetParError(0);
    double iErrorMeshWithout = fIMeshWithout->GetParError(0);
    /*
     double iErrorMeshWith = fIMeshWith->GetParameter(1);
     double iErrorMeshWithout = fIMeshWithout->GetParameter(1);
     */
    
    double iDriftWith, iDriftWithout, iErrorDriftWith, iErrorDriftWithout;
    if (true) {
        iDriftWith = fIDriftWith->GetParameter(0);
        iDriftWithout = fIDriftWithout->GetParameter(0);
        
        iErrorDriftWith = fIDriftWith->GetParError(0);
        iErrorDriftWithout = fIDriftWithout->GetParError(0);
        /*
         iErrorDriftWith = fIDriftWith->GetParameter(1);
         iErrorDriftWithout = fIDriftWithout->GetParameter(1);
         */
    }
    else {
        iDriftWith = hImonDriftWith->GetMean();
        iDriftWithout = hImonDriftWithout->GetMean();
        iErrorDriftWith = hImonDriftWith->GetMeanError();
        iErrorDriftWithout = hImonDriftWithout->GetMeanError();
    }
    
    
    
    ibf = abs( (iDriftWith-iDriftWithout)/(iMeshWith-iMeshWithout) );
    
    cout << endl << endl;
    cout << "iDriftWith = " << iDriftWith << endl;
    cout << "iDriftWithout = " << iDriftWithout << endl;
    cout << "iMeshWith = " << iMeshWith << endl;
    cout << "iMeshWithout = " << iMeshWithout << endl;
    cout << endl << endl;
    
    ibfError = ibf * sqrt( Square((iErrorMeshWith+iErrorMeshWithout)/(iMeshWith-iMeshWithout)) + Square((iErrorDriftWith+iErrorDriftWithout)/(iDriftWith-iDriftWithout)) );
    /*
    double firstTermSquared = (Square(iErrorDriftWith)+Square(iErrorDriftWithout))/Square(iDriftWith-iDriftWithout);
    double secondTermSquared = (Square(iErrorMeshWith)+Square(iErrorMeshWithout))/Square(iMeshWith-iMeshWithout);
    ibfError = sqrt( firstTermSquared + secondTermSquared );
     */
    
    //cout << endl << endl << endl << iMeshWith-iMeshWithout << endl;
    //cout << endl << endl << endl << iDriftWith-iDriftWithout << endl;
    
    return;
}
