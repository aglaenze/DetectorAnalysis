#include <iostream>
#include <fstream>
#include <cmath>

#include <string.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TString.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

using namespace std;
// Set up detector geometry

struct NoSorting {
	bool operator()(string const& lhsIn, string const& rhsIn) const {
		if (lhsIn > rhsIn) return true;	// if condition to silence the warning
		else return true;
	}
};

int LoadElectrodeMap(string detectorName, map <string, int, NoSorting>& electrodeMap) {
	electrodeMap = {};
	if (detectorName=="ZZBOT") {
		electrodeMap["pad"] = 0;
		electrodeMap["mesh"] = 1;
		electrodeMap["drift"] = 2;
		electrodeMap["screen"] = 3;
		return 0;
	}
	else if (detectorName=="RD3SP1") {
		electrodeMap["pad"] = 0;
		electrodeMap["mesh"] = 1;
		electrodeMap["GEM1 bottom"] = 2;
		electrodeMap["GEM1 top"] = 3;
		electrodeMap["GEM2 bottom"] = 4;
		electrodeMap["GEM2 top"] = 5;
		electrodeMap["drift"] = 6;
		return 0;
	}
	else if (detectorName=="RD3SP3") {
		electrodeMap["pad"] = 0;
		electrodeMap["mesh"] = 1;
		electrodeMap["DM bottom"] = 2;
		electrodeMap["DM top"] = 3;
		electrodeMap["drift"] = 4;
		return 0;
	}
	else if (detectorName=="MGEM1" || detectorName=="MGEM3") {
		electrodeMap["pad"] = 0;
		electrodeMap["mesh"] = 1;
		electrodeMap["GEM bottom"] = 2;
		electrodeMap["GEM top"] = 3;
		electrodeMap["mesh top"] = 4;
		electrodeMap["drift"] = 5;
		return 0;
	}
	else {cout << "no info for this model" << endl; return -1;}
}


void LoadParameters(string detectorName, vector<double>& zElectrodes) {
	if (detectorName == "ZZBOT") {
		zElectrodes = {0.5, 0.0128, 0};
	}
	else if (detectorName == "RD3SP1") {
		double damp = 0.2+0.4+0.0125;
		zElectrodes = {0.8+damp, 0.0060+damp, damp, 0.0060+0.4+0.0125, 0.4+0.0125, 0.0125, 0};
	}
	else if (detectorName == "RD3SP3") {
		zElectrodes = {4.3+0.0128+0.19+0.0128, 0.0128+0.19+0.0128, 0.19+0.0128, 0.0128, 0};
	}
	else if (detectorName == "RD3SP4") {
		zElectrodes = {0.87+0.0128+0.4, 0.0128+0.4, 0.4, 0};
	}
	else if (detectorName == "MGEM1") {  // MM + GEM + MM
		double damp = 0.0128+0.722+0.0060+0.0660;
		zElectrodes = {damp+1.05, damp, 0.0128+0.722+0.0060, 0.0128+0.722, 0.0128, 0};
	}
	else if (detectorName == "MGEM3") {  // MM + GEM + MM
		double damp = 0.0128+0.52;
		zElectrodes = {damp+1.38, damp, 0.0128+0.52-0.0128, 0.0128+0.52-0.0128-0.0060, 0.0128, 0};
	}
	else if (detectorName == "LittleChinese") {   // MM + GEM
		zElectrodes = {0.3+0.0320+0.0220, 0.0320+0.0220, 0.0220, 0};
	}
	
	else {std::cout << "What detector?" << std::endl; return;}
}

int GetElectrodeNum(string detectorName) {
	vector<double> zElectrodes;
	LoadParameters(detectorName, zElectrodes);
	return (int)zElectrodes.size()-1;
}

void DrawDetector(string detectorName, vector<int> hvList) {
	// Draw geometry with voltages
	
	gStyle->SetTextSize(0.05);
	
	map <string, int, NoSorting> electrodeMap;
	LoadElectrodeMap(detectorName, electrodeMap);
	int electrodeNum = GetElectrodeNum(detectorName)+1;
	
	std::vector<double> zElectrodes = {};
	LoadParameters(detectorName, zElectrodes);
	
	// insert pad voltage in HV list
	hvList.insert(hvList.begin(),0);
	
	double yMax = 0.8; // coord of drift electrode
	double space = yMax/(electrodeNum-1);
	
	TBox* detectorBox = new TBox(0,0,1,1);
	detectorBox->Draw();
	
	int i = 0;	// i is the stage index
	
	// draw the electrodes from top to bottom (as declared in the function LoadElectrodeMap)
	
	map<string, int>::iterator it = electrodeMap.begin();
	// Iterate over the map using Iterator until end.
	while (it != electrodeMap.end()) {
		double yCoord = yMax-i*space;
		TLine* electrodeLine = new TLine(0, yCoord, 1, yCoord);
		electrodeLine->SetLineStyle(9);
		// 0 = drift
		if (i == 0 || i == electrodeNum-1) electrodeLine->SetLineStyle(1);
		TLatex* electrodeText = new TLatex(0.05, yCoord+0.03, Form("V_{%s} = -%d V", (it->first).c_str(), hvList[electrodeNum-1-i]));
		TText* zElectrode = new TText(0.6, yCoord+0.03, Form("z = %.3f mm", zElectrodes[i]*10));
		electrodeLine->Draw("same");
		electrodeText->Draw("same");
		zElectrode->Draw("same");
		// Computation of electric field
		if (i != electrodeNum-1) {
			double electricField = (hvList[electrodeNum-1-i]-hvList[electrodeNum-1-(i+1)])/(zElectrodes[i]-zElectrodes[i+1]);
			TString fieldtxt = Form("E = %.1f V/cm", electricField);
			if (electricField>1000.) fieldtxt = Form("E = %.2f kV/cm", electricField/1000.);
			TText* electricFieldtxt = new TText(0.05, yCoord-0.5*space, fieldtxt);
			electricFieldtxt->SetTextColor(kBlue);
			electricFieldtxt->Draw("same");
		}
		// Increment the Iterator to point to next entry
		it++;
		i++;
	}
}

double GetzElectrode(string detectorName, string electrodeName) {
    vector<double> zElectrodes = {};
    map <string, int, NoSorting> electrodeMap;
    LoadParameters(detectorName, zElectrodes);
    LoadElectrodeMap(detectorName, electrodeMap);
    int i = 0;
    map<string, int>::iterator it = electrodeMap.begin();
    while (it != electrodeMap.end()) {
        if (it->first == electrodeName) {
            return zElectrodes[i];
        }
        it++;
        i++;
    }
    return -1;
}

double GetCalibrationAlpha(string detectorName, string date) {
	if (detectorName == "ZZBOT") {
		return 5.45217e-06;
	}
	else if (detectorName == "RD3SP1") {
		return 5.07508e-06;
	}
	else if (detectorName == "RD3SP3") {
		return 5.07508e-06;
	}
	else if (detectorName == "RD3SP4") {
		return 5.07508e-06;
	}
	else if (detectorName == "MGEM1") {  // MM + GEM + MM
										 //return 1.08932e-05;
		//return 1.22298e-05;
        //return 1.51055e-05;
        //return 4.23323e-06;
        //return 1.06289e-05;     // data taken in TL, 17/02/2021
        //return 1.13753e-05/2;     // data taken in TL, 18/02/2021
        //return 1.0095e-05;      // data taken in TL, 24/02/2021
        //return 8.16756e-06;      // data taken in TL, 24/02/2021
        //return 4.23323e-06/1.5;
        //return 4.23323e-06*1.5;
        if (date == "2021-04-08" || date == "2021-04-12" || date == "2021-04-13" || date == "2021-10-07" || date == "2021-10-11" || date == "2021-10-14") return 7.9131e-06;
        return 3.67496e-06;
	}
	else if (detectorName == "MGEM3") {  // MM + GEM + MM
		//return 1.25197e-05;
        //return 7.81109e-06;
        //return 8.93988e-06;     // 15/04/2021 mesh down
        if (date == "2021-04-15" || date == "2021-04-21"  || date == "2021-08-30"  || date == "2021-08-31" ) return 8.17276e-06;     // 15/04/2021 mesh top
        return 3.90004e-06;       // 06/10/2021 mesh down
	}
	else if (detectorName == "LittleChinese") {   // MM + GEM
		return 5.65359e-06;
	}
	else {std::cout << "What detector?" << std::endl; return -1.;}
}
