#include <cmath>
#include <dirent.h>
#include <stdio.h>

#include "_Utils.C"
#include "_HistUtils.C"


//____________________________________________
int DrawIbf( void )
{

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

     
    // Gain
    
    Double_t gainList[11] = {430.621, 601.575, 835.974, 1167.004, 1626.864, 2250.949, 3148.231, 4456.768, 6222.110, 8731.924, 12189.061};
    Double_t gainErrorList[11] = {13.817, 12.924, 11.392, 12.381, 13.891, 23.617, 26.841, 47.777, 57.497, 118.237, 140.382};
    //Double_t hvMmList[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
    Double_t resolutionList[11] = {35.433, 30.824, 28.267, 27.155, 26.954, 26.493, 26.466, 26.976, 27.315, 28.186, 29.380};
    Double_t resolutionErrorList[11] = {12.108, 7.574, 5.653, 4.522, 3.770, 4.498, 3.839, 4.628, 3.993, 5.553, 4.886};
    
    // log_2020_02_05
    Double_t hvMmList[11] = {340.000, 350.000, 360.000, 370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000};
    Double_t currentMeshListWith[11] = {-181797.344, -217663.190, -322755.399, -468123.690, -642788.442, -877367.507, -1230585.564, -1727479.197, -2413370.707, -3333273.896, -4661297.009};
    Double_t currentMeshErrorListWith[11] = {2821.869, 2767.678, 4126.040, 2493.752, 3665.984, 3814.833, 3645.747, 3777.000, 1894.379, 2914.395, 2887.029};
    Double_t currentDriftListWith[11] = {-2322.678, -4891.550, -4627.981, -6981.046, -8995.081, -10112.968, -12024.583, -15926.902, -20716.124, -25601.323, -34972.231};
    Double_t currentDriftErrorListWith[11] = {382.101, 316.579, 369.986, 248.967, 313.666, 417.592, 434.191, 337.505, 326.667, 493.164, 371.412};
    Double_t currentMeshListWithout[11] = {-6105.812, -12242.192, -7341.650, -11894.114, -17399.778, -11616.907, -13540.589, -1121.364, -15469.705, 1603.582, 388.781};
    Double_t currentMeshErrorListWithout[11] = {4960.335, 2700.575, 3244.432, 3165.233, 3889.531, 2898.654, 3236.111, 4524.517, 2970.527, 3222.727, 4456.962};
    Double_t currentDriftListWithout[11] = {-180.962, -1157.481, -263.850, -2519.345, -2495.841, -1416.220, 386.998, 461.031, -359.302, -1478.734, -765.492};
    Double_t currentDriftErrorListWithout[11] = {457.699, 350.644, 300.628, 319.786, 348.323, 411.511, 503.659, 362.432, 368.152, 343.907, 393.273};

    const int n = 11;
    
    
    Double_t ibf[n], ibfError[n];
    for (int i=0; i<n; i++) {
        double correctedCurrentDrift = currentDriftListWith[i] - currentDriftListWithout[i];
        double correctedCurrentMesh = currentMeshListWith[i] - currentMeshListWithout[i];
        /*
        double correctedCurrentDrift = currentDriftListWith[i];
        double correctedCurrentMesh = currentMeshListWith[i];
        */
        double correctedCurrentDriftError = currentDriftErrorListWith[i] + currentDriftErrorListWithout[i];
        double correctedCurrentMeshError = currentMeshErrorListWith[i] + currentMeshErrorListWithout[i];
        
        ibf[i] = correctedCurrentDrift/correctedCurrentMesh * 100.;
        //ibf14[i] = currentDriftList2[i]/currentMeshList2[i] * 100.;
        ibfError[i] = ibf[i] * TMath::Sqrt(Square(correctedCurrentDriftError/correctedCurrentDrift) + Square(correctedCurrentMeshError/correctedCurrentMesh) );
        if ( ibf[i] < 0 ) {ibf[i] = 15; ibfError[i] = 10;}
    }
    
    PrintList( "ibf", ibf, n, "%.3f" );
    PrintList( "ibfError", ibfError, n, "%.3f" );
    
    TCanvas* cv = new TCanvas("cv", "cv", 200, 400);
    cv->Divide(1,3);


    cv->cd(1);
    gPad->SetGrid();
    gPad->SetLogy();
    TGraphErrors* gr1 = new TGraphErrors(n, hvMmList, ibf, 0, ibfError);
    gr1->SetTitle("IBF = f(HV)");
    gr1->GetXaxis()->SetTitle( "HV (V)" );
    gr1->GetYaxis()->SetTitle( "Ibf (%)" );
    //gr->GetHistogram()->GetXaxis()->SetLimits(310, 450);
    gr1->GetHistogram()->SetMinimum(0.01);   // along Y axis
    gr1->GetHistogram()->SetMaximum(30);   // along Y axis
    //gr1->SetMarkerStyle(20);
    gr1->SetLineColor(kRed);
    gr1->Draw("ALP");

    
    TLegend* legend1 = new TLegend(0.2,0.7,0.6,0.9);
    //legend1->AddEntry(gr1,"Data (14/11)","lp");
    //legend1->Draw();
    
    cv->cd(2);
    gPad->SetGrid();
    gPad->SetLogx();
    gPad->SetLogy();
    TGraphErrors* gr2 = new TGraphErrors(n, gainList, ibf, gainErrorList, ibfError);
    gr2->SetTitle("IBF = f(Gain)");
    gr2->GetXaxis()->SetTitle( "Gain" );
    gr2->GetYaxis()->SetTitle( "Ibf (%)" );
    //gr2->GetHistogram()->GetXaxis()->SetLimits(8.e2, 3.e4);
    gr2->GetHistogram()->GetXaxis()->SetLimits(gainList[0]*0.8, gainList[n-1]*1.2);
    gr2->GetHistogram()->SetMinimum(0.01);   // along Y axis
    gr2->GetHistogram()->SetMaximum(30);   // along Y axis
    //gr2->SetMarkerStyle(20);
    gr2->SetLineColor(kRed);
    gr2->Draw("ALP");

    
    TLegend* legend2 = new TLegend(0.2,0.7,0.6,0.9);
    //legend2->AddEntry(gr2,"Data (14/11)","lp");
    //legend2->Draw();
    

    cv->cd(3);
    gPad->SetGrid();
    gPad->SetLogx();
    TGraphErrors* gr3 = new TGraphErrors(n, ibf, resolutionList, ibfError, resolutionErrorList);
    gr3->SetTitle("Resolution = f(IBF)");
    gr3->GetXaxis()->SetTitle( "IBF(%)" );
    gr3->GetYaxis()->SetTitle( "Resolution (%)" );
    //gr3->GetHistogram()->GetXaxis()->SetLimits(310, 450);
    gr3->GetHistogram()->GetXaxis()->SetLimits(0.01, 30);
    gr3->GetHistogram()->SetMinimum(0);   // along Y axis
    gr3->GetHistogram()->SetMaximum(100);   // along Y axis
    //gr3->SetMarkerStyle(20);
    gr3->SetLineColor(kRed);
    gr3->Draw("ALP");

    
    TLegend* legend3 = new TLegend(0.2,0.7,0.6,0.9);
    //legend3->AddEntry(gr3,"Data (14/11)","lp");
    //legend3->Draw();

    
    cv->SaveAs("Figures/Ibf.pdf");
    
    return 0;
    
}
