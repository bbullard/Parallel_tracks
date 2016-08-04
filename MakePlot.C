#include "MakePlot.h"
#include "FastMC.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TString.h"
#include "TProfile.h"
#include "TCut.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "THStack.h"
#include "TLegend.h"

#include <iostream>

using namespace std;

MakePlot::MakePlot(const Double_t angle, const Bool_t r)
{
    fm = new FastMC();
    SetStyle();
    n_events = fm->fChain->GetEntriesFast();
    cutAngle = angle;
    rot = r;
}

MakePlot::~MakePlot()
{
    delete fm;
}

void MakePlot::PlotEv_vs_Angle(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw("Ev:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCQE(100,0,50,200,0,20)","(cc&&qel)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvsAngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Neutrino Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("Ev:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCRes(100,0,50,200,0,20)","(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Neutrino Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("Ev:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCDis(100,0,50,200,0,20)","(cc&&dis)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Neutrino Energy (GeV)");

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEl_vs_Angle(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCQE(100,0,50,200,0,20)","(cc&&qel)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvsAngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCRes(100,0,50,200,0,20)","(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("El:asin(abs(pxl/pl))/3.14159*180>>EvsAngleCCDis(100,0,50,200,0,20)","(cc&&dis)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvsAngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEl_vs_Length(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCQE(100,0,20,200,0,20)","(cc&&qel)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("El_vs_LengthCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(2);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCRes(100,0,20,200,0,20)","(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("El_vs_LengthCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    c1->cd(3);
    fm->fChain->Draw("El:0.4/sin(abs(asin(pxl/pl)))>>El_vs_LengthCCDis(100,0,20,200,0,20)","(cc&&dis)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("El_vs_LengthCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    h->GetYaxis()->SetTitle("Lepton Energy (GeV)");

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEv(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",800, 600);

    fm->fChain->Draw("Ev>>Ev(400,0.1,40)");
    TH1F *h = (TH1F*)gDirectory->Get("Ev");
    h->SetTitle("");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    gPad->SetLogy();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotOscWeight_vs_Ev(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",800, 600);

    fm->fChain->Draw("OSCWeight:Ev>>Osc(400,0.1,40)");
    TH1F *h = (TH1F*)gDirectory->Get("Osc");
    h->SetTitle(";Neutrino Energy (GeV);Oscillation Weight");
    gPad->SetLogx();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEl_vs_Ev(const char* saveAs)
{
    int nTotal = fm->fChain->Draw("","cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("El:Ev>>EvElCCQE(200,0,20,200,0,20)","(cc&&qel)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvElCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetYaxis()->SetTitle("Lepton energy (GeV)");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    c1->cd(4);
    TProfile *px = h->ProfileX("pxCCQE");
    TH1D *pxp = px->ProjectionX("pxpCCQE");
    for (int i=1; i<=200; i++) {
        pxp->SetBinContent(i, pxp->GetBinContent(i)/i*10) ;
        pxp->SetBinError(i, pxp->GetBinError(i)/i*10);
    }
    pxp->GetYaxis()->SetTitle("Lepton Energy / Neutrino Energy");
    pxp->Draw();

    c1->cd(2);
    fm->fChain->Draw("El:Ev>>EvElCCRes(200,0,20,200,0,20)","(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvElCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetYaxis()->SetTitle("Lepton energy (GeV)");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    c1->cd(5);
    px = h->ProfileX("pxCCRes");
    pxp = px->ProjectionX("pxpCCRes");
    for (int i=1; i<=200; i++) {
        pxp->SetBinContent(i, pxp->GetBinContent(i)/i*10) ;
        pxp->SetBinError(i, pxp->GetBinError(i)/i*10);
    }
    pxp->GetYaxis()->SetTitle("Lepton Energy / Neutrino Energy");
    pxp->Draw();

    c1->cd(3);
    fm->fChain->Draw("El:Ev>>EvElCCDis(200,0,20,200,0,20)","(cc&&dis)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvElCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetYaxis()->SetTitle("Lepton energy (GeV)");
    h->GetXaxis()->SetTitle("Neutrino Energy (GeV)");
    c1->cd(6);
    px = h->ProfileX("pxCCDis");
    pxp = px->ProjectionX("pxpCCDis");
    for (int i=1; i<=200; i++) {
        pxp->SetBinContent(i, pxp->GetBinContent(i)/i*10) ;
        pxp->SetBinError(i, pxp->GetBinError(i)/i*10);
    }
    pxp->GetYaxis()->SetTitle("Lepton Energy / Neutrino Energy");
    pxp->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotAngle_cut_E(const char* saveAs)
{
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", (ecut && "cc")*"OSCWeight");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCQE(100,0,50)", ecut && "cc&&qel");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCQE(180,0,90)", (ecut && "cc&&qel")*"OSCWeight");
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(2);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCRes(100,0,50)", ecut && "cc&&res");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCRes(180,0,90)", (ecut && "cc&&res")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(3);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCDis(100,0,50)", ecut && "cc&&dis");
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180>>AngleCCDis(180,0,90)", (ecut && "cc&&dis")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotAngle_cut_E_Perp(const char* saveAs)
{
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCQE(100,0,50)", ecut && "cc&&qel");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCQE(180,0,90)", (ecut && "cc&&qel")*"OSCWeight");
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(2);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCRes(100,0,50)", ecut && "cc&&res");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCRes(180,0,90)", (ecut && "cc&&res")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    c1->cd(3);
    // fm->fChain->Draw("acos(cthl)/3.14159*180>>AngleCCDis(100,0,50)", ecut && "cc&&dis");
    fm->fChain->Draw("abs(asin(pzl/pl))/3.14159*180>>AngleCCDis(180,0,90)", (ecut && "cc&&dis")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotLength_cut_E(const char* saveAs)
{
    TCut ecut("El<6");
    TH1F *h = 0;
    TH1 *hc = 0;

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCQE(100,0,20)", (ecut && "cc&&qel")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("LengthCCQE");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(4);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(2);
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCRes(100,0,20)", (ecut && "cc&&res")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("LengthCCRes");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(3);
    fm->fChain->Draw("0.4/sin(abs(asin(pxl/pl)))>>LengthCCDis(100,0,20)", (ecut && "cc&&dis")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("LengthCCDis");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotLength_cut_E_Perp(const char* saveAs)
{
    TCut ecut("El<6");
    TH1F *h = 0;
    TH1 *hc = 0;

    int nTotal = fm->fChain->Draw("", ecut && "cc");

    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCQE(100,0,20)", (ecut && "cc&&qel")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("LengthCCQE");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(4);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(2);
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCRes(100,0,20)", (ecut && "cc&&res")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("LengthCCRes");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    c1->cd(3);
    fm->fChain->Draw("0.4/sin(abs(asin(pzl/pl)))>>LengthCCDis(100,0,20)", (ecut && "cc&&dis")*"OSCWeight");
    h = (TH1F*)gDirectory->Get("LengthCCDis");
    cout << "average length: " << h->GetMean() << endl;
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Lepton length in beginning 4-mm slice [cm]");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,101));
    for (int i=1; i<=100; i++) {
        hc->SetBinContent(i, 1-hc->GetBinContent(i));
    }
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    // gPad->SetLogy();
    cout << hc->GetBinContent(3*5) << ", " << hc->GetBinContent(6*5) << ", " << hc->GetBinContent(10*5) << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEv_vs_Ef(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("Ev:Ef>>EvEfCCQE(50,0,5,200,0,20)","(cc&&qel&&pf>0)*OSCWeight","colz");
    TH2F *h = (TH2F*)gDirectory->Get("EvEfCCQE");
    h->SetTitle("CC-QE");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    c1->cd(4);
    h->ProjectionX()->Draw();


    c1->cd(2);
    fm->fChain->Draw("Ev:Ef>>EvEfCCRes(50,0,5,200,0,20)","(cc&&res&&pf>0)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvEfCCRes");
    h->SetTitle("CC-Res");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    c1->cd(5);
    h->ProjectionX()->Draw();

    c1->cd(3);
    fm->fChain->Draw("Ev:Ef>>EvEfCCDis(50,0,5,200,0,20)","(cc&&dis&&pf>0)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("EvEfCCDis");
    h->SetTitle("CC-Dis");
    h->GetYaxis()->SetTitle("Neutrino energy (GeV)");
    h->GetXaxis()->SetTitle("Hadron Energy (GeV)");
    gPad->SetLogz();
    c1->cd(6);
    h->ProjectionX()->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotNf_vs_Ratio(const char* particle, const char* title, const char* saveAs)
{
    TH2F *h = 0;
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    c1->Divide(3,1);
    c1->cd(1);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCQE(51,0,1.02,51,0,51)", particle), "(cc&&qel)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCQE");
    h->SetTitle("CC-QE");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    c1->cd(2);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCRes(51,0,1.02,51,0,51)", particle),"(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCRes");
    h->SetTitle("CC-Res");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    c1->cd(3);
    fm->fChain->Draw(TString::Format("nf:(%s)/nf>>Nf_RatioCCDis(51,0,1.02,51,0,51)", particle),"(cc&&res)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Nf_RatioCCDis");
    h->SetTitle("CC-Dis");
    h->GetYaxis()->SetTitle("# Final state particles (hadron side)");
    h->GetXaxis()->SetTitle(title);
    gPad->SetLogz();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotAnglel_vs_Angleh(const char* saveAs)
{
    TH2F *h = 0;
    TCanvas *c1 = new TCanvas("c1","",1200, 800);
    c1->Divide(3,2);
    c1->cd(1);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCQE(180,0,90,180,0,90)","(cc&&qel&&El<6&&pf>0&&pdgf!=2112)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCQE");
    h->SetTitle("CC-QE");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(4);
    h->ProjectionX()->Draw();

    c1->cd(2);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCRes(180,0,90,180,0,90)","(cc&&res&&El<6&&pf>0&&pdgf!=2112)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCRes");
    h->SetTitle("CC-Res");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(5);
    h->ProjectionX()->Draw();

    c1->cd(3);
    fm->fChain->Draw("abs(asin(pxl/pl))/3.14159*180:abs(asin(pxf/pf))/3.14159*180>>Angleh_vs_AnglelCCDis(180,0,90,180,0,90)","(cc&&dis&&El<6&&pf>0&&pdgf!=2112)*OSCWeight","colz");
    h = (TH2F*)gDirectory->Get("Angleh_vs_AnglelCCDis");
    h->SetTitle("CC-Dis");
    h->GetXaxis()->SetTitle("hadron angle to wire plane (deg)");
    h->GetYaxis()->SetTitle("lepton angle to wire plane (deg)");
    gPad->SetLogz();
    c1->cd(6);
    h->ProjectionX()->Draw();

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotMinAngle_lepton(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 800);

    TFile *f = new TFile("numu-nue-cc-hadron_angle.root");
    TTree *t = (TTree*)f->Get("T");
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", (ecut && "cc")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    cout<<"-----------"<<endl;
    cout<<"Total simulated events: "<<fm->fChain->GetEntries()<<endl;
    cout<<"Number of CC events: "<<nTotal<<endl;
    cout<<"Number of long hadron events: "<<t->Draw("",Form("(El<6&&cc&&h_min_angle<91)*OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr))<<endl;
    c1->Divide(3,2);
    c1->cd(1);
    t->Draw("l_angle>>LAngleCCQE(180,0,90)", (ecut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("LAngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum lepton angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    Double_t I = h->Integral();
    cout << "CCQE:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << "-----------" << endl;

    c1->cd(2);
    t->Draw("l_angle>>LAngleCCRes(180,0,90)", (ecut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("LAngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum lepton angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCRes:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << "-----------" << endl;

    c1->cd(3);
    t->Draw("l_angle>>LAngleCCDis(180,0,90)", (ecut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("LAngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum lepton angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCDis:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << "-----------" << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotMinAngle_hadron(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 800);

    TFile *f = new TFile("numu-nue-cc-hadron_angle.root");
    TTree *t = (TTree*)f->Get("T");
    TCut ecut("El<6");

    int nTotal = fm->fChain->Draw("", (ecut && "cc")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    cout<<"-----------"<<endl;
    cout<<"Total simulated events: "<<fm->fChain->GetEntries()<<endl;
    cout<<"Number of CC events: "<<nTotal<<endl;
    cout<<"Number of long hadron events: "<<t->Draw("",Form("(El<6&&cc&&h_min_angle<91)*OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr))<<endl;
    c1->Divide(3,2);
    c1->cd(1);
    t->Draw("h_min_angle>>AngleCCQE(180,0,90)", (ecut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("AngleCCQE");
    h->SetTitle(TString::Format("CC-QE (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(4);
    TH1* hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    Double_t I = h->Integral();
    cout << "CCQE:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*0.315 << ", " << hc->GetBinContent(7.5*2)*0.453 << ", " << hc->GetBinContent(10*2)*0.563 << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << hc->GetBinContent(5*2)*0.315*I << ", " << hc->GetBinContent(7.5*2)*0.453*I << ", " << hc->GetBinContent(10*2)*0.563*I << endl;
    cout << "-----------" << endl;

    c1->cd(2);
    t->Draw("h_min_angle>>AngleCCRes(180,0,90)", (ecut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("AngleCCRes");
    h->SetTitle(TString::Format("CC-Res (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(5);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCRes:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*0.280 << ", " << hc->GetBinContent(7.5*2)*0.401 << ", " << hc->GetBinContent(10*2)*0.507 << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << hc->GetBinContent(5*2)*0.280*I << ", " << hc->GetBinContent(7.5*2)*0.401*I << ", " << hc->GetBinContent(10*2)*0.507*I << endl;
    cout << "-----------" << endl;

    c1->cd(3);
    t->Draw("h_min_angle>>AngleCCDis(180,0,90)", (ecut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("AngleCCDis");
    h->SetTitle(TString::Format("CC-Dis (%.1f%%)", h->GetEntries()/nTotal*100));
    h->GetXaxis()->SetTitle("Minimum hadron angle to wire plane (deg)");
    c1->cd(6);
    hc = h->GetCumulative();
    hc->Scale(1./h->Integral(0,181));
    hc->GetYaxis()->SetTitle("Cumulative ratio");
    hc->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    I = h->Integral();
    cout << "CCDis:" << endl;
    cout << I <<endl;
    cout << hc->GetBinContent(5*2) << ", " << hc->GetBinContent(7.5*2) << ", " << hc->GetBinContent(10*2) << endl;
    cout << hc->GetBinContent(5*2)*0.176 << ", " << hc->GetBinContent(7.5*2)*0.270 << ", " << hc->GetBinContent(10*2)*0.364 << endl;
    cout << hc->GetBinContent(5*2)*I << ", " << hc->GetBinContent(7.5*2)*I << ", " << hc->GetBinContent(10*2)*I << endl;
    cout << hc->GetBinContent(5*2)*0.176*I << ", " << hc->GetBinContent(7.5*2)*0.270*I << ", " << hc->GetBinContent(10*2)*0.364*I << endl;
    cout << "-----------" << endl;

    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::CreateTree(const char* filename)
{
    TFile f(filename,"recreate");
    TTree *T = new TTree("T", "New Tree from Fast MC");
    float h_min_angle = 91;
    int h_min_angle_pdg = 0;
    float h_min_angle_Ef = 0;
    float h_min_angle_eside = 91;
    int h_min_angle_eside_pdg = 0;
    float h_min_angle_eside_Ef = 0;
    // float h_sec_min_angle = 91;
    float l_angle = 91;
    float a;

    T->Branch("h_min_angle", &h_min_angle);
    T->Branch("h_min_angle_pdg", &h_min_angle_pdg);
    T->Branch("h_min_angle_Ef", &h_min_angle_Ef);
    T->Branch("h_min_angle_eside", &h_min_angle_eside);
    T->Branch("h_min_angle_eside_pdg", &h_min_angle_eside_pdg);
    T->Branch("h_min_angle_eside_Ef", &h_min_angle_eside_Ef);
    // T->Branch("h_sec_min_angle", &h_sec_min_angle);
    T->Branch("l_angle", &l_angle);
    T->Branch("cc", &(fm->cc) );
    T->Branch("qel", &(fm->qel) );
    T->Branch("res", &(fm->res) );
    T->Branch("dis", &(fm->dis) );
    T->Branch("El", &(fm->El));
    T->Branch("nf", &(fm->nf));
    T->Branch("POTWeight",&(fm->POTWeight));
    T->Branch("OSCWeight",&(fm->OSCWeight));
    T->Branch("n_events",&n_events);
    T->Branch("DetMass",&DetMass);
    T->Branch("exposure",&exposure);
    T->Branch("POTperYr",&POTperYr);

    int nEntries = fm->fChain->GetEntries();
    for (int i=0; i<nEntries; i++) {
        fm->fChain->GetEntry(i);

        if(!rot) l_angle = asin(fm->pxl/fm->pl)/3.14159*180;
	if(rot) l_angle = asin(fm->pzl/fm->pl)/3.14159*180;

        for (int j=0; j<fm->nf; j++) {
            if (fm->pf[j]<=0 || fm->pdgf[j] == 2112 || fm->pdgf[j] == 111 || fm->pdgf[j] == 22) continue;
	    //Cut hadrons below threshold
	    if (fm->pdgf[j] == 2212 && fm->Ef[j] < 0.039) continue; //cut protons below 39 MeV
	    if ((fm->pdgf[j] == 321 || fm->pdgf[j] == -321) && fm->Ef[j] < 0.03) continue; //cut kaons below 30 MeV
	    if ((fm->pdgf[j] == 211 || fm->pdgf[j] == -211) && fm->Ef[j] < 0.0175) continue; //cut pions below 17.5 MeV
            if(!rot) a = asin(fm->pxf[j]/fm->pf[j])/3.14159*180;
	    if(rot) a = asin(fm->pzf[j]/fm->pf[j])/3.14159*180;
            if (abs(a)<h_min_angle_eside && l_angle*a>0) {
                h_min_angle_eside = abs(a);
                h_min_angle_eside_pdg = fm->pdgf[j];
                h_min_angle_eside_Ef = fm->Ef[j];
            }
	    if (abs(a)<h_min_angle) {
                h_min_angle = abs(a);
                h_min_angle_pdg = fm->pdgf[j];
                h_min_angle_Ef = fm->Ef[j];
            }
            // if (a < h_sec_min_angle_wp && a>h_min_angle_wp) h_sec_min_angle_wp = a;
        }
	l_angle = abs(l_angle);
        T->Fill();

        h_min_angle = 91;
        h_min_angle_pdg = 0;
        h_min_angle_Ef = 0;
	h_min_angle_eside = 91;
        h_min_angle_eside_pdg = 0;
        h_min_angle_eside_Ef = 0;
        l_angle = 91;
    }
    T->Write();
    f.Close();
}

void MakePlot::CreateCCEMTree(const char* filename)
{
	TFile f(filename,"recreate");
    	TTree *T = new TTree("T", "New Tree from Fast MC");
	
	float l_angle = 91;
	float ph_min_angle = 91;
	float l_ph_convert = 0;
	float Eph = 0;
	float a;
	TRandom3 r;
	
	bool isE = 0;

    	T->Branch("l_angle", &l_angle);
	T->Branch("EM_min_angle", &ph_min_angle);
	T->Branch("l_convert", &l_ph_convert);
    	T->Branch("cc", &(fm->cc) );
	T->Branch("nc", &(fm->nc) );
    	T->Branch("qel", &(fm->qel) );
    	T->Branch("res", &(fm->res) );
    	T->Branch("dis", &(fm->dis) );
    	T->Branch("El", &(fm->El));
    	T->Branch("EEM", &Eph);
	T->Branch("POTWeight",&(fm->POTWeight));
	T->Branch("OSCWeight",&(fm->OSCWeight));
	T->Branch("n_events",&n_events);
	T->Branch("DetMass",&DetMass);
	T->Branch("exposure",&exposure);
	T->Branch("POTperYr",&POTperYr);
	
	int nEntries = fm->fChain->GetEntries();

    	for (int i=0; i<nEntries; i++) {
        	fm->fChain->GetEntry(i);
		if(!rot) l_angle = abs(asin(fm->pxl/fm->pl))/3.14159*180;
		if(rot) l_angle = abs(asin(fm->pzl/fm->pl))/3.14159*180;
        	for (int j=0; j<fm->nEMf; j++) {
			if(fm->EMEf[j]>6) continue;
			if(!rot) a = abs(asin(fm->EMpxf[j]/fm->EMpf[j]))/3.14159*180;
			if(rot) a = abs(asin(fm->EMpzf[j]/fm->EMpf[j]))/3.14159*180;
			if(a < ph_min_angle) {
				if(fm->EMpdgf[j] == 22) isE = 0;
				else isE = 1;
				float l = r.Exp(18);
				if(isE == 0 && l>3) continue;
				l_ph_convert = l;
				ph_min_angle = a;
				Eph = fm->EMEf[j];
			}
        	}
		if(isE) l_ph_convert = 0;
		T->Fill();
		ph_min_angle = 91;
		l_angle = 91;
    	}	
	T->Write();
	f.Close();
}

void MakePlot::CreateNCEMTree(const char* filename)
{
	TFile f(filename,"recreate");
    	TTree *T = new TTree("T", "New Tree from Fast MC");
	
	float l_angle = 91;
	float ph_min_angle = 91;
	float l_ph_convert = 0;
	float Eph = 0;
	float a, ap, ah;
	TRandom3 r;
	bool isE = 0, hasPi0 = 0;
	
    	T->Branch("l_angle", &l_angle);
	T->Branch("EM_min_angle", &ph_min_angle);
	T->Branch("l_convert", &l_ph_convert);
    	T->Branch("cc", &(fm->cc) );
	T->Branch("nc", &(fm->nc) );
    	T->Branch("qel", &(fm->qel) );
    	T->Branch("res", &(fm->res) );
    	T->Branch("dis", &(fm->dis) );
    	T->Branch("El", &(fm->El));
    	T->Branch("EEM", &Eph);
	T->Branch("isE", &isE);
	T->Branch("hasPi0", &hasPi0);
	T->Branch("POTWeight",&(fm->POTWeight));
	T->Branch("OSCWeight",&(fm->OSCWeight));
	T->Branch("n_events",&n_events);
	T->Branch("DetMass",&DetMass);
	T->Branch("exposure",&exposure);
	T->Branch("POTperYr",&POTperYr);
	
	int nEntries = fm->fChain->GetEntries();

    	for (int i=0; i<nEntries; i++) {
        	fm->fChain->GetEntry(i);
		hasPi0 = 0;
		for (int k=0; k<fm->nf; k++) {
           		if (fm->pf[k]>0 && fm->pdgf[k] == 111) hasPi0 = 1;
        	}
		l_angle = asin(fm->pxl/fm->pl)/3.14159*180;
        	for (int j=0; j<fm->nEMf; j++) {
			if(fm->EMEf[j]>6) continue;
			if(!rot) a = asin(fm->EMpxf[j]/fm->EMpf[j])/3.14159*180;
			if(rot) a = asin(fm->EMpzf[j]/fm->EMpf[j])/3.14159*180;
			if(abs(a) < ph_min_angle) {
				if(fm->EMpdgf[j] != 22) 
				{
					for (int k=0; k<fm->nf; k++) {
						if (fm->pf[k]<=0 || fm->pdgf[k] == 2112 || fm->pdgf[k] == 111 || fm->pdgf[k] == 22) continue;
            					if(!rot) ah = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
						if(rot) ah = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
            					if (abs(ah)>cutAngle) continue;
						if (ah*a < 0) continue;
						ph_min_angle = abs(a);
						isE = 1;
						Eph = fm->EMEf[j];
        				}
				}
				else 
				{
					for (int k=0; k<fm->nf; k++) {
           					if (fm->pf[k]<=0 || (fm->pdgf[k] != 211 && fm->pdgf[k] != -211)) continue;
            					if(!rot) ap = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
						if(rot) ap = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
						if (abs(ap)>cutAngle) continue;
						if (ap*a<0) continue;
						float l = r.Exp(18);
						if (l <1 || l>5 ) continue;
						l_ph_convert = l;	
						ph_min_angle = abs(a);
						isE = 0;
						Eph = fm->EMEf[j];
					}
				}
			}
        	}
		if(isE) l_ph_convert = 0;
		T->Fill();
		ph_min_angle = 91;
		l_angle = 91;
    	}	
	T->Write();
	f.Close();
}

void MakePlot::PlotEvents_vs_El_l(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    gStyle->SetLegendBorderSize(0);

    TFile *fccs = new TFile(Form("numu-nue-cc-hadron_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccs = (TTree*)fccs->Get("T");
    TFile *fccb = new TFile(Form("numu-nue-cc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccb = (TTree*)fccb->Get("T");
    TFile *fncb = new TFile(Form("numu-numu-nc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tncb = (TTree*)fncb->Get("T");
    TCut phcut2("isE==0&&hasPi0");
    TCut ecut("isE==1&&hasPi0");
    TCut EMacut(Form("EM_min_angle<%f",cutAngle));
    TCut lgood(Form("l_angle>%f",cutAngle));
    TCut lbad(Form("l_angle<%f",cutAngle));
    
    cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    cout<<"Lepton cut"<<endl;
    cout<<"Is rotated: "<<rot<<endl;
    cout<<"Angle cut: "<<cutAngle<<endl;
    cout<<"QE"<<endl;
    cout<<"----------------------------"<<endl;
    c1->Divide(3,1);
    THStack *hs = new THStack("hs","QE");
    TLegend *leg = new TLegend(.52,.68,.89,.89);
    c1->cd(1);
    tncb->Draw("EEM>>GNCbQE(40,0,6)", (phcut2 && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("GNCbQE");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow#gamma+X NC [#theta_{#gamma}<%.1f#circ]",cutAngle),"f");
    tncb->Draw("EEM>>ENCbQE(40,0,6)", (ecut && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbQE");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e+X NC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    tccb->Draw("EEM>>ECCbQE(40,0,6)", (EMacut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbQE");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e/#gamma CC [#theta_{EM}<%.1f#circ]",cutAngle),"f");
    tccs->Draw("El>>ECCsQE(40,0,6)", (lbad && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQE");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsQEg(40,0,6)", (lgood && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQEg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e}>%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    
    cout<<"Res"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Res");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(2);
    tncb->Draw("EEM>>GNCbRes(40,0,6)", (phcut2 && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbRes");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbRes(40,0,6)", (ecut && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbRes");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbRes(40,0,6)", (EMacut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbRes");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsRes(40,0,6)", (lbad && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsRes");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsResg(40,0,6)", (lgood && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsResg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Res;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    
    cout<<"Dis"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Dis");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(3);
    tncb->Draw("EEM>>GNCbDis(40,0,6)", (phcut2 && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbDis");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbDis(40,0,6)", (ecut && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbDis");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbDis(40,0,6)", (EMacut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbDis");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsDis(40,0,6)", (lbad && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDis");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsDisg(40,0,6)", (lgood && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDisg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Dis;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    cout<<"----------------------------------------------------------------------------"<<endl<<endl;
    
    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEvents_vs_El_h(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    gStyle->SetLegendBorderSize(0);

    TFile *fccs = new TFile(Form("numu-nue-cc-hadron_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccs = (TTree*)fccs->Get("T");
    TFile *fccb = new TFile(Form("numu-nue-cc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccb = (TTree*)fccb->Get("T");
    TFile *fncb = new TFile(Form("numu-numu-nc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tncb = (TTree*)fncb->Get("T");
    TCut phcut2("isE==0&&hasPi0");
    TCut ecut("isE==1&&hasPi0");
    TCut EMacut(Form("EM_min_angle<%f",cutAngle));
    TCut hgood(Form("h_min_angle>%f",cutAngle));
    TCut hbad(Form("h_min_angle<%f",cutAngle));
    
    cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    cout<<"Hadron cut"<<endl;
    cout<<"Is rotated: "<<rot<<endl;
    cout<<"Angle cut: "<<cutAngle<<endl;
    cout<<"QE"<<endl;
    cout<<"----------------------------"<<endl;
    c1->Divide(3,1);
    THStack *hs = new THStack("hs","QE");
    TLegend *leg = new TLegend(.52,.68,.89,.89);
    c1->cd(1);
    tncb->Draw("EEM>>GNCbQE(40,0,6)", (phcut2 && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("GNCbQE");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow#gamma+X NC [#theta_{#gamma}<%.1f#circ]",cutAngle),"f");
    tncb->Draw("EEM>>ENCbQE(40,0,6)", (ecut && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbQE");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e+X NC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    tccb->Draw("EEM>>ECCbQE(40,0,6)", (EMacut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbQE");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e/#gamma CC [#theta_{EM}<%.1f#circ]",cutAngle),"f");
    tccs->Draw("El>>ECCsQE(40,0,6)", (hbad && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQE");
    h->SetFillColor(kRed);
    cout<<"CC Problem (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{h}<%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsQEg(40,0,6)", (hgood && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQEg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{h}>%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    
    cout<<"Res"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Res");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(2);
    tncb->Draw("EEM>>GNCbRes(40,0,6)", (phcut2 && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbRes");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbRes(40,0,6)", (ecut && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbRes");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbRes(40,0,6)", (EMacut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbRes");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsRes(40,0,6)", (hbad && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsRes");
    h->SetFillColor(kRed);
    cout<<"CC Problem (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsResg(40,0,6)", (hgood && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsResg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Res;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    
    cout<<"Dis"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Dis");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(3);
    tncb->Draw("EEM>>GNCbDis(40,0,6)", (phcut2 && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbDis");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbDis(40,0,6)", (ecut && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbDis");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbDis(40,0,6)", (EMacut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbDis");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsDis(40,0,6)", (hbad && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDis");
    h->SetFillColor(kRed);
    cout<<"CC Problem (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsDisg(40,0,6)", (hgood && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDisg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Dis;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    cout<<"----------------------------------------------------------------------------"<<endl<<endl;
    
    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::PlotEvents_vs_El_lh(const char* saveAs)
{
    TCanvas *c1 = new TCanvas("c1","",1200, 400);
    gStyle->SetLegendBorderSize(0);

    TFile *fccs = new TFile(Form("numu-nue-cc-hadron_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccs = (TTree*)fccs->Get("T");
    TFile *fccb = new TFile(Form("numu-nue-cc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tccb = (TTree*)fccb->Get("T");
    TFile *fncb = new TFile(Form("numu-numu-nc-EM_angle_%.1f_%i.root",cutAngle,rot));
    TTree *tncb = (TTree*)fncb->Get("T");
    TCut phcut2("isE==0&&hasPi0");
    TCut ecut("isE==1&&hasPi0");
    TCut EMacut(Form("EM_min_angle<%f",cutAngle));
    TCut lhgood(Form("l_angle>%f||h_min_angle>%f",cutAngle,cutAngle));
    TCut lhbad(Form("l_angle<%f&&h_min_angle<%f",cutAngle,cutAngle));
    
    cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    cout<<"Lepton and Hadron cut"<<endl;
    cout<<"Is rotated: "<<rot<<endl;
    cout<<"Angle cut: "<<cutAngle<<endl;
    cout<<"QE"<<endl;
    cout<<"----------------------------"<<endl;
    c1->Divide(3,1);
    THStack *hs = new THStack("hs","QE");
    TLegend *leg = new TLegend(.52,.68,.89,.89);
    c1->cd(1);
    tncb->Draw("EEM>>GNCbQE(40,0,6)", (phcut2 && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    TH1F *h = (TH1F*)gDirectory->Get("GNCbQE");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow#gamma+X NC [#theta_{#gamma}<%.1f#circ]",cutAngle),"f");
    tncb->Draw("EEM>>ENCbQE(40,0,6)", (ecut && EMacut && "nc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbQE");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e+X NC [#theta_{e}<%.1f#circ]",cutAngle),"f");
    tccb->Draw("EEM>>ECCbQE(40,0,6)", (EMacut && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbQE");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#pi^{0}#rightarrow e/#gamma CC [#theta_{EM}<%.1f#circ]",cutAngle),"f");
    tccs->Draw("El>>ECCsQE(40,0,6)", (lhbad && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQE");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e},#theta_{h}<%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsQEg(40,0,6)", (lhgood && "cc&&qel")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsQEg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton and hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    leg->AddEntry(h,Form("#nu_{e} CC [#theta_{e} or #theta_{h}>%.1f#circ]",cutAngle),"f");
    hs->SetTitle("QE;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    
    cout<<"Res"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Res");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(2);
    tncb->Draw("EEM>>GNCbRes(40,0,6)", (phcut2 && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbRes");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbRes(40,0,6)", (ecut && EMacut && "nc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbRes");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbRes(40,0,6)", (EMacut && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbRes");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsRes(40,0,6)", (lhbad && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsRes");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsResg(40,0,6)", (lhgood && "cc&&res")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsResg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton and hadron angle): "<<h->Integral()<<endl<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Res;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    
    cout<<"Dis"<<endl;
    cout<<"----------------------------"<<endl;
    hs = new THStack("hs","Dis");
    //leg = new TLegend(.52,.7,.89,.89);
    c1->cd(3);
    tncb->Draw("EEM>>GNCbDis(40,0,6)", (phcut2 && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("GNCbDis");
    h->SetFillColor(kViolet+3);
    cout<<"NC gamma Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{#mu} NC [#theta_{EM}<7.5#circ]","f");
    tncb->Draw("EEM>>ENCbDis(40,0,6)", (ecut && EMacut && "nc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ENCbDis");
    h->SetFillColor(kViolet);
    cout<<"NC e Background: "<<h->Integral()<<endl;
    hs->Add(h);
    tccb->Draw("EEM>>ECCbDis(40,0,6)", (EMacut && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCbDis");
    h->SetFillColor(kOrange);
    cout<<"CC Background: "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"EM #nu_{e} CC [#theta_{EM}<7.5#circ]","f");
    tccs->Draw("El>>ECCsDis(40,0,6)", (lhbad && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDis");
    h->SetFillColor(kRed);
    cout<<"CC Problem (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}<7.5#circ]","f");
    hs->SetTitle("QE;Angle to wire plane (deg);Events/150 MeV/34 kT/3#times10^{20} POT");
    tccs->Draw("El>>ECCsDisg(40,0,6)", (lhgood && "cc&&dis")*Form("OSCWeight*POTWeight*%f*%f/%f*%f",DetMass,exposure,n_events,POTperYr));
    h = (TH1F*)gDirectory->Get("ECCsDisg");
    h->SetFillColor(kBlue);
    cout<<"CC Good (lepton and hadron angle): "<<h->Integral()<<endl;
    hs->Add(h);
    //leg->AddEntry(h,"#nu_{e} CC [#theta_{e}>7.5#circ]","f");
    hs->SetTitle("Dis;E_{e/#gamma} [GeV];Events/150 MeV/34 kT/3.3#times10^{21} POT");
    hs->Draw("HIST");
    leg->Draw();
    cout<<"----------------------------------------------------------------------------"<<endl<<endl;
    
    if (saveAs) {
        c1->SaveAs(saveAs);
    }
}

void MakePlot::TabulateNCEvents()
{
	TFile f("NCcuts.root","recreate");
    	TTree *T = new TTree("T", "New Tree from Fast MC");
	
	float l_ph_convert = 0;
	float a, ap, ah;
	TRandom3 r;
	
	bool hasPi0 = 0, hasShallowE = 0, EhasH = 0, hasShallowG = 0, GdecR = 0, GhasPi = 0;
	int t = 1;
	
	T->Branch("t",&t);
    	T->Branch("cc", &(fm->cc) );
	T->Branch("nc", &(fm->nc) );
    	T->Branch("qel", &(fm->qel) );
    	T->Branch("res", &(fm->res) );
    	T->Branch("dis", &(fm->dis) );
	T->Branch("hasPi0", &hasPi0);
	T->Branch("hasShallowE", &hasShallowE);
	T->Branch("EhasH", &EhasH);
	T->Branch("hasShallowG", &hasShallowG);
	T->Branch("GdecR", &GdecR);
	T->Branch("GhasPi", &GhasPi);
	T->Branch("POTWeight",&(fm->POTWeight));
	T->Branch("OSCWeight",&(fm->OSCWeight));
	T->Branch("n_events",&n_events);
	T->Branch("DetMass",&DetMass);
	T->Branch("exposure",&exposure);
	T->Branch("POTperYr",&POTperYr);
	
	int nEntries = fm->fChain->GetEntries();

    	for (int i=0; i<nEntries; i++) {
        	fm->fChain->GetEntry(i);
		for (int k=0; k<fm->nf; k++) {
           		if (fm->pf[k]>0 && fm->pdgf[k] == 111) hasPi0 = 1;
        	}
        	for (int j=0; j<fm->nEMf; j++) {
			if(fm->EMEf[j]>6) continue;
			if(!rot) a = asin(fm->EMpxf[j]/fm->EMpf[j])/3.14159*180;
			if(rot) a = asin(fm->EMpzf[j]/fm->EMpf[j])/3.14159*180;
			if(fm->EMpdgf[j] != 22 && abs(a)<cutAngle) 
			{
				hasShallowE = 1;
				for (int k=0; k<fm->nf; k++) {
					if (fm->pf[k]<=0 || fm->pdgf[k] == 2112 || fm->pdgf[k] == 111 || fm->pdgf[k] == 22) continue;
            				if(!rot) ah = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
					if(rot) ah = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
            				if (abs(ah)>cutAngle) continue;
					if (ah*a < 0) continue;
					EhasH = 1;
        				}
			}
			if(fm->EMpdgf[j] == 22 && abs(a)<cutAngle){
				hasShallowG = 1;
				float l = r.Exp(18);
				if (l >=1 && l<=5 ) GdecR = 1;	
				for (int k=0; k<fm->nf; k++) {
            				if (fm->pf[k]<=0 || (fm->pdgf[k] != 211 && fm->pdgf[k] != -211)) continue;
            				if(!rot) ap = asin(fm->pxf[k]/fm->pf[k])/3.14159*180;
					if(rot) ap = asin(fm->pzf[k]/fm->pf[k])/3.14159*180;
					if (abs(ap)>cutAngle) continue;
					if (ap*a<0) continue;
					if (GdecR == 0 ) continue;
					GhasPi = 1;
				}
			}	
        	}
		T->Fill();
		hasPi0 = 0; hasShallowE = 0; EhasH = 0; hasShallowG = 0; GdecR = 0; GhasPi = 0; l_ph_convert = 0;
    	}	
	T->Write();
	
	TCanvas can("can","can",500,500);
	cout<<endl<<"----------------------------------------------------------------------------"<<endl;
    	cout<<"Is rotated: "<<rot<<endl;
    	cout<<"Angle cut: "<<cutAngle<<endl;
	T->Draw("t>>h1(40,0,6)","(dis||res||qel)*nc*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	TH1F *h = (TH1F*)gDirectory->Get("h1");
	cout<<"Total NC: "<<h->Integral()<<endl;
	T->Draw("t>>h2(40,0,6)","(dis||res||qel)*nc*hasPi0*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h2");
	cout<<"NC with pi0: "<<h->Integral()<<endl;
	T->Draw("t>>h3(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowE*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h3");
	cout<<"NC with pi0 and shallow e: "<<h->Integral()<<endl;
	T->Draw("t>>h3p1(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowE*EhasH*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h3p1");
	cout<<"NC with pi0 and shallow e and hadron: "<<h->Integral()<<endl;
	T->Draw("t>>h4(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4");
	cout<<"NC with pi0 and shallow gamma: "<<h->Integral()<<endl;
	T->Draw("t>>h4p1(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*GdecR*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4p1");
	cout<<"NC with pi0 and shallow gamma in decay range: "<<h->Integral()<<endl;
	T->Draw("t>>h4p2(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*GdecR*GhasPi*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4p2");
	cout<<"NC with pi0 and shallow gamma in decay range and pion: "<<h->Integral()<<endl;
	T->Draw("t>>h4p2a3p1(40,0,6)","(dis||res||qel)*nc*hasPi0*hasShallowG*GdecR*GhasPi*hasShallowE*EhasH*OSCWeight*POTWeight*DetMass*exposure/n_events*POTperYr");
	h = (TH1F*)gDirectory->Get("h4p2a3p1");
	cout<<"NC with pi0 and shallow gamma in decay range and pion, and shallow e and hadron: "<<h->Integral()<<endl;
	f.Close();
}

void MakePlot::SetStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleFont(42, "x");
    gStyle->SetTitleFont(42, "y");
    gStyle->SetTitleFont(42, "z");
    gStyle->SetLabelFont(42, "x");
    gStyle->SetLabelFont(42, "y");
    gStyle->SetLabelFont(42, "z");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.35, "y");
    gStyle->SetHistLineWidth(2);
    gStyle->SetHistLineColor(kBlack);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetPalette(1);
    gROOT->ForceStyle();
    gROOT->UseCurrentStyle();
    fm->fChain->UseCurrentStyle();

}
