#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "Include/Utils.C"
#include "Include/HistUtils.C"

using namespace std;


//____________________________________________
int DrawGainMM( string detectorName = "MGEM1" )
{
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleFontSize(.05);
    gStyle->SetTitleXSize(.05);
    gStyle->SetTitleYSize(.05);
    gStyle->SetLabelSize(0.05, "X");
    gStyle->SetLabelSize(0.05, "Y");
    
    
    // Gain
    /*
     // Pour MGEM3
     const int num = 2;
     int dV[num] = {250, 350};     // 250 300 350 400 450 500
     Double_t gainList[8] = {286.611, 255.277, 587.699, 523.402, 1749.067, 1557.851, 5642.467, 4999.005};
     Double_t gainErrorList[8] = {16.138, 20.507, 13.328, 20.580, 65.146, 39.770, 132.502, 167.821};
     Double_t hvMeshList[8] = {350, 350, 370, 370, 400, 400, 430, 430};
     Double_t hvGemDownList[8] = {600, 700, 620, 720, 650, 750, 680, 780};
     Double_t hvGemTopList[8] = {850, 950, 870, 970, 900, 1000, 930, 1030};
     Double_t hvMeshTopList[8] = {800, 900, 800, 900, 800, 900, 800, 900};
     */
    const int num = 1;
    int dV[num] = {290};     // 250 300 350 400 450 500
    
    Double_t gainList[6] = {1034.285, 2052.183, 4046.364, 8189.668, 17205.632, 25288.181};
    Double_t gainErrorList[6] = {12.972, 38.570, 72.416, 158.497, 478.747, 473.427};
    Double_t hvMeshList[6] = {340, 360, 380, 400, 420, 430};
    Double_t hvGemDownList[6] = {630, 650, 670, 690, 710, 720};
    Double_t hvGemTopList[6] = {600, 600, 600, 650, 650, 650};
    Double_t hvMeshTopList[6] = {600, 600, 600, 650, 650, 650};
    Double_t dvMeshTopList[6] = {0, 0, 0, 0, 0, 0};
    Double_t coarseGainList[6] = {500, 200, 100, 50, 20, 20};
    
    
    
    vector<float> gainVec[num] = {};
    vector<float> gainErrorVec[num] = {};
    vector<int> hvMeshVec[num] = {};
    vector<int> hvGemDownVec[num] = {};
    vector<float> hvRatio[num] = {};
    
    for (int l = 0; l<num; l++) {
        for (int k = 0; k < 30; k++) {
            if (hvGemDownList[k] - hvMeshList[k] == dV[l]) {
                gainVec[l].push_back(gainList[k]);
                gainErrorVec[l].push_back(gainErrorList[k]);
                hvMeshVec[l].push_back(hvMeshList[k]);
                hvGemDownVec[l].push_back(hvGemDownList[k]);
                hvRatio[l].push_back((double)hvMeshList[k]/(hvGemDownList[k]-hvMeshList[k])* 5./0.128 );
                
            }
        }
    }
    
    TCanvas* cv = new TCanvas();
    cv->Divide(2,2);
    
    cv->cd(1);
    gPad->SetGrid();
    gPad->SetLogy();
    
    cv->cd(2);
    gPad->SetGrid();
    gPad->SetLogy();
    
    TLegend* legend1 = new TLegend(0.6,0.1,0.9,0.4);
    legend1->SetTextSize(0.045);
    
    TLegend* legend2 = new TLegend(0.1,0.1,0.9,0.9);
    legend2->SetTextSize(0.05);
    
    for (int l = 0; l<num; l++) {
        // create list from vector
        const int nn = (int)hvMeshVec[l].size();
        //TGraphErrors* gr1 = CreateTGraphError(hvMeshVec[l], gainVec[l], gainErrorVec[l]);
        Double_t li1[nn], li2[nn], li3[nn], li4[nn], li5[nn];
        for (int k = 0; k< nn; k++) {
            li1[k] = hvMeshVec[l][k];
            li2[k] = hvGemDownVec[l][k];
            li3[k] = gainVec[l][k];
            li4[k] = gainErrorVec[l][k];
            li5[k] = hvRatio[l][k];
        }
        
        Double_t yMin = 5000, yMax = 100;
        int numberofelements = sizeof(gainList)/sizeof(gainList[0]);
        for (int k = 0; k < numberofelements; k++)
        {
            if (yMin > gainList[k]) yMin = gainList[k];
            if (yMax < gainList[k]) yMax = gainList[k];
        }
        yMin *= 0.5;
        yMax *= 2;
        TGraphErrors* gr1 = new TGraphErrors(nn, li1, li3, 0, li4);
        gr1->SetTitle("Gain = f(HV)");
        gr1->GetXaxis()->SetTitle( "HV (V)" );
        gr1->GetYaxis()->SetTitle( "Gain" );
        gr1->SetMarkerStyle(20);
        gr1->GetHistogram()->SetMinimum(yMin);   // along Y axis
        gr1->GetHistogram()->SetMaximum(yMax);   // along Y axis
        
        legend1->AddEntry(gr1, Form("dV = %d V", dV[l]), "lp");
        
        // Exp fit function
        TF1* f = new TF1( "FitFunctionExp", FitFunctionExp, li1[0]-5, li1[nn-1]+5, 2 );
        f->SetParameter(0, -7 );
        f->SetParameter(1, 0.04 );
        f->SetParName( 0, "const" );
        f->SetParName( 1, "slope" );
        gr1->Fit( f, "0" );
        f->SetLineColor(l+1);
        f->Draw("same");
        
        legend2->AddEntry(f, Form("dV = %dV: y = exp(%.4f + %.4f x)", dV[l], f->GetParameter(0), f->GetParameter(1) ), "l");
        
        cv->cd(1);
        gr1->SetLineColor(l+1);
        gr1->SetMarkerColor(l+1);
        if (l == 0) gr1->Draw("ALP");
        else gr1->Draw("LP");
        
        TGraphErrors* gr2 = new TGraphErrors(nn, li5, li3, 0, li4);
        gr2->SetTitle("Gain = f(HV)");
        gr2->GetXaxis()->SetTitle( "HV Ratio" );
        gr2->GetYaxis()->SetTitle( "Gain" );
        gr2->GetHistogram()->GetXaxis()->SetLimits(30, 120);
        gr2->GetHistogram()->SetMinimum(yMin);   // along Y axis
        gr2->GetHistogram()->SetMaximum(yMax);   // along Y axis
        gr2->SetMarkerStyle(20);
        cv->cd(2);
        gr2->SetLineColor(l+1);
        gr2->SetMarkerColor(l+1);
        if (l == 0) gr2->Draw("ALP");
        else gr2->Draw("LP");
    }
    
    cv->cd(1);
    legend1->Draw();
    cv->cd(2);
    legend1->Draw();
    
    cv->cd(3);
    legend2->Draw();
    
    
    cv->SaveAs(Form("Figures/%s/GainMM2.pdf", detectorName.c_str()));
    
    return 0;
    
}

