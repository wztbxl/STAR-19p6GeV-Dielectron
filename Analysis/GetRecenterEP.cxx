#include <algorithm>
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TF2.h"
#include "TFile.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TFormula.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TProfile.h"

#include "/Users/wztbxl/Documents/macros/function.C"

void GetRecenterEP()
{
    TCanvas* c1 = new TCanvas("c1","c1",1600,900);
    
    for(int i = 0; i <= 30; i++)
    {
        TFile* inFile = new TFile(Form("%d.histo.root",i));
        auto his = (TH1D*)inFile->Get("hReCenterEventPlane");
        his->Draw();
        c1->SaveAs(Form("./plots/%d.png",i));
        // inFile->Delete();
        // his->Delete();
    }

    TFile* inFile = new TFile("newEP.root");
    TH3D* hQXvsQYvsRunIndex = (TH3D*)inFile->Get("hQXvsQYvsRunIndex");
    TH3D* hQXvsQYvsRunIndex_raw = (TH3D*)inFile->Get("hQXvsQYvsRunIndex_raw");
    TH3D* hQXvsQYvsRunIndex_recenter_west = (TH3D*)inFile->Get("hQXvsQYvsRunIndex_recenter_west");
    TH3D* hQXvsQYvsRunIndex_recenter_east = (TH3D*)inFile->Get("hQXvsQYvsRunIndex_recenter_east");

    hQXvsQYvsRunIndex->GetZaxis()->SetRange(2,10);
    hQXvsQYvsRunIndex_raw->GetZaxis()->SetRange(2,10);
    TH2D* hQXvsQY = (TH2D*)hQXvsQYvsRunIndex->Project3D("xy");
    TH2D* hQXvsQY_raw = (TH2D*)hQXvsQYvsRunIndex_raw->Project3D("xy");

    c1->SetLogz();
    hQXvsQY->Draw("colz");
    c1->SaveAs("./plots/recenterEP.png");

    hQXvsQY_raw->Draw("colz");
    c1->SaveAs("./plots/recenterEP_Raw.png");
}