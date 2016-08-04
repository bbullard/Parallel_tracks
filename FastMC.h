//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 20 17:23:24 2016 by ROOT version 5.34/32
// from TTree gst/GENIE Summary Event Tree
// found on file: fastmcNtp_20140711_lbne_g4lbnev3r2p4_nuflux_numuflux_nue_LAr_1_g280_Ar40_5000_Default.root
//////////////////////////////////////////////////////////

#ifndef FastMC_h
#define FastMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class FastMC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Char_t          FluxVersion;
   Int_t           iev;
   Int_t           neu;
   Int_t           fspl;
   Int_t           tgt;
   Int_t           Z;
   Int_t           A;
   Int_t           hitnuc;
   Int_t           hitqrk;
   Int_t           resid;
   Bool_t          sea;
   Bool_t          qel;
   Bool_t          mec;
   Bool_t          res;
   Bool_t          dis;
   Bool_t          coh;
   Bool_t          cohe;
   Bool_t          dfr;
   Bool_t          imd;
   Bool_t          imda;
   Bool_t          nuel;
   Bool_t          em;
   Bool_t          cc;
   Bool_t          nc;
   Bool_t          charm;
   Int_t           neut_code;
   Int_t           nuance_code;
   Double_t        wght;
   Double_t        POTWeight;
   Double_t        POTperYear;
   Double_t        nomPOTperYear;
   Double_t        OSCWeight;
   Double_t        xs;
   Double_t        ys;
   Double_t        ts;
   Double_t        Q2s;
   Double_t        Ws;
   Double_t        x;
   Double_t        y;
   Double_t        t;
   Double_t        Q2;
   Double_t        P;
   Double_t        ET;
   Double_t        W;
   Double_t        W2;
   Double_t        Ev;
   Double_t        Ehad;
   Double_t        Evtx;
   Double_t        p0v;
   Double_t        pxv;
   Double_t        pyv;
   Double_t        pzv;
   Double_t        x_reco;
   Double_t        y_reco;
   Double_t        t_reco;
   Double_t        Q2_reco;
   Double_t        QE_Q2_reco;
   Double_t        QE_Ev_reco;
   Double_t        W_reco;
   Double_t        W2_reco;
   Double_t        Ev_reco;
   Double_t        EvB_reco;
   Double_t        EvUB_reco;
   Double_t        Ehad_reco;
   Double_t        EhadB_reco;
   Double_t        EhadUB_reco;
   Double_t        Ehad_unbias_default;
   Double_t        Ehad_unbias_numu;
   Double_t        Ehad_unbias_numubar;
   Double_t        Ehad_unbias_nue;
   Double_t        Ehad_unbias_nuebar;
   Double_t        Ehad_unbias_nc;
   Double_t        Ehad_unbias_ncbar;
   Double_t        p0v_reco;
   Double_t        pxv_reco;
   Double_t        pyv_reco;
   Double_t        pzv_reco;
   Int_t           EvClass_reco;
   Int_t           isQE_reco;
   Bool_t          IsContv_reco;
   Int_t           hasElikeEM;
   Double_t        En;
   Double_t        p0n;
   Double_t        pxn;
   Double_t        pyn;
   Double_t        pzn;
   Int_t           Pdgl;
   Double_t        El;
   Double_t        p0l;
   Double_t        pxl;
   Double_t        pyl;
   Double_t        pzl;
   Double_t        pl;
   Double_t        cthl;
   Int_t           PdgL;
   Double_t        EL;
   Double_t        p0L;
   Double_t        pxL;
   Double_t        pyL;
   Double_t        pzL;
   Double_t        pL;
   Double_t        cthL;
   Int_t           ParentPdgL;
   Double_t        ParentEL;
   Double_t        p0LParent;
   Double_t        pxLParent;
   Double_t        pyLParent;
   Double_t        pzLParent;
   Double_t        pLParent;
   Double_t        cthLParent;
   Int_t           PdgL_reco;
   Double_t        EL_reco;
   Double_t        p0L_reco;
   Double_t        pxL_reco;
   Double_t        pyL_reco;
   Double_t        pzL_reco;
   Double_t        pL_reco;
   Double_t        cthL_reco;
   Double_t        TrkL_reco;
   Bool_t          IsContL_reco;
   Int_t           Pdgl_reco;
   Double_t        El_reco;
   Double_t        p0l_reco;
   Double_t        pxl_reco;
   Double_t        pyl_reco;
   Double_t        pzl_reco;
   Double_t        pl_reco;
   Double_t        Trkl_reco;
   Double_t        cthl_reco;
   Int_t           nfp;
   Int_t           nfn;
   Int_t           nfpip;
   Int_t           nfpim;
   Int_t           nfpi0;
   Int_t           nfkp;
   Int_t           nfkm;
   Int_t           nfk0;
   Int_t           nfem;
   Int_t           nfother;
   Int_t           nip;
   Int_t           nin;
   Int_t           nipip;
   Int_t           nipim;
   Int_t           nipi0;
   Int_t           nikp;
   Int_t           nikm;
   Int_t           nik0;
   Int_t           niem;
   Int_t           niother;
   Int_t           ni;
   Int_t           pdgi[20];   //[ni]
   Int_t           resc[20];   //[ni]
   Double_t        Ei[20];   //[ni]
   Double_t        p0i[20];   //[ni]
   Double_t        pxi[20];   //[ni]
   Double_t        pyi[20];   //[ni]
   Double_t        pzi[20];   //[ni]
   Int_t           nf;
   Int_t           pdgf[50];   //[nf]
   Double_t        Ef[50];   //[nf]
   Double_t        p0f[50];   //[nf]
   Double_t        pxf[50];   //[nf]
   Double_t        pyf[50];   //[nf]
   Double_t        pzf[50];   //[nf]
   Double_t        pf[50];   //[nf]
   Double_t        cthf[50];   //[nf]
   Double_t        Ef_reco[50];   //[nf]
   Double_t        p0f_reco[50];   //[nf]
   Double_t        pxf_reco[50];   //[nf]
   Double_t        pyf_reco[50];   //[nf]
   Double_t        pzf_reco[50];   //[nf]
   Double_t        Trkf_reco[50];   //[nf]
   Double_t        cthf_reco[50];   //[nf]
   Double_t        pT_Imbalance;
   Int_t           nEMf;
   Int_t           EMpdgf[22];   //[nEMf]
   Double_t        EMp0f[22];   //[nEMf]
   Double_t        EMpxf[22];   //[nEMf]
   Double_t        EMpyf[22];   //[nEMf]
   Double_t        EMpzf[22];   //[nEMf]
   Double_t        EMpf[22];   //[nEMf]
   Double_t        EMcthf[22];   //[nEMf]
   Double_t        EMEf[22];   //[nEMf]
   Double_t        EMEf_reco[22];   //[nEMf]
   Double_t        EMp0f_reco[22];   //[nEMf]
   Double_t        EMpxf_reco[22];   //[nEMf]
   Double_t        EMpyf_reco[22];   //[nEMf]
   Double_t        EMpzf_reco[22];   //[nEMf]
   Double_t        EMTrkf_reco[22];   //[nEMf]
   Double_t        EMcthf_reco[22];   //[nEMf]
   Double_t        EMPizMass_reco[22];   //[nEMf]
   Double_t        EMElPizM_reco[22];   //[nEMf]
   Double_t        vtxx;
   Double_t        vtxy;
   Double_t        vtxz;
   Double_t        vtxt;
   Double_t        calresp0;
   Double_t        caltotal;
   Double_t        calneutron;
   Double_t        calproton;
   Double_t        calpion;
   Double_t        calmeson;
   Double_t        calbaryon;
   Double_t        calhadron;
   Double_t        calem;
   Double_t        calother;
   Double_t        calbinding;
   Double_t        Tau_Prob_nue;
   Double_t        Tau_Prob_numu;
   Int_t           rw7;
   Double_t        rw7_qelCCvalenciaRPA[7];   //[rw7]
   Double_t        rw7_qelNCvalenciaRPA[7];   //[rw7]
   Double_t        rw7_reslowq2[7];   //[rw7]
   Int_t           frw7;
   Int_t           fNsyst;
   Double_t        frw7_TargetYTilt0_5mm[7];   //[frw7]
   Double_t        frw7_TargetYOffset0_5mm[7];   //[frw7]
   Double_t        frw7_TargetXTilt0_5mm[7];   //[frw7]
   Double_t        frw7_TargetXOffset0_5mm[7];   //[frw7]
   Double_t        frw7_TargetBe[7];   //[frw7]
   Double_t        frw7_ProtonP90GeV[7];   //[frw7]
   Double_t        frw7_ProtonP80GeV[7];   //[frw7]
   Double_t        frw7_ProtonP70GeV[7];   //[frw7]
   Double_t        frw7_ProtonP60GeV[7];   //[frw7]
   Double_t        frw7_ProtonP50GeV[7];   //[frw7]
   Double_t        frw7_ProtonP40GeV[7];   //[frw7]
   Double_t        frw7_ProtonP30GeV[7];   //[frw7]
   Double_t        frw7_ProtonP20GeV[7];   //[frw7]
   Double_t        frw7_ProtonP150GeV[7];   //[frw7]
   Double_t        frw7_ProtonP140GeV[7];   //[frw7]
   Double_t        frw7_ProtonP130GeV[7];   //[frw7]
   Double_t        frw7_ProtonP110GeV[7];   //[frw7]
   Double_t        frw7_ProtonP100GeV[7];   //[frw7]
   Double_t        frw7_Horn2YOffset0_5mm[7];   //[frw7]
   Double_t        frw7_Horn2XOffset0_5mm[7];   //[frw7]
   Double_t        frw7_Horn1YOffset0_5mm[7];   //[frw7]
   Double_t        frw7_Horn1XOffset0_5mm[7];   //[frw7]
   Double_t        frw7_DecayPipeRadius3m[7];   //[frw7]
   Double_t        frw7_DecayPipeRadius1_9m[7];   //[frw7]
   Double_t        frw7_DecayPipeLength250m[7];   //[frw7]
   Double_t        frw7_BeamSigmaY1_4mm[7];   //[frw7]
   Double_t        frw7_BeamSigmaX1_4mm[7];   //[frw7]
   Double_t        frw7_700kW[7];   //[frw7]
   Double_t        frw7_230kA[7];   //[frw7]

   // List of branches
   TBranch        *b_FluxVersion;   //!
   TBranch        *b_iev;   //!
   TBranch        *b_neu;   //!
   TBranch        *b_fspl;   //!
   TBranch        *b_tgt;   //!
   TBranch        *b_Z;   //!
   TBranch        *b_A;   //!
   TBranch        *b_hitnuc;   //!
   TBranch        *b_hitqrk;   //!
   TBranch        *b_resid;   //!
   TBranch        *b_sea;   //!
   TBranch        *b_qel;   //!
   TBranch        *b_mec;   //!
   TBranch        *b_res;   //!
   TBranch        *b_dis;   //!
   TBranch        *b_coh;   //!
   TBranch        *b_cohe;   //!
   TBranch        *b_dfr;   //!
   TBranch        *b_imd;   //!
   TBranch        *b_imda;   //!
   TBranch        *b_nuel;   //!
   TBranch        *b_em;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_nc;   //!
   TBranch        *b_charm;   //!
   TBranch        *b_neut_code;   //!
   TBranch        *b_nuance_code;   //!
   TBranch        *b_wght;   //!
   TBranch        *b_POTWeight;   //!
   TBranch        *b_POTperYear;   //!
   TBranch        *b_nomPOTperYear;   //!
   TBranch        *b_OSCWeight;   //!
   TBranch        *b_xs;   //!
   TBranch        *b_ys;   //!
   TBranch        *b_ts;   //!
   TBranch        *b_Q2s;   //!
   TBranch        *b_Ws;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_t;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_P;   //!
   TBranch        *b_ET;   //!
   TBranch        *b_W;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_Ev;   //!
   TBranch        *b_Ehad;   //!
   TBranch        *b_Evtx;   //!
   TBranch        *b_p0v;   //!
   TBranch        *b_pxv;   //!
   TBranch        *b_pyv;   //!
   TBranch        *b_pzv;   //!
   TBranch        *b_x_reco;   //!
   TBranch        *b_y_reco;   //!
   TBranch        *b_t_reco;   //!
   TBranch        *b_Q2_reco;   //!
   TBranch        *b_QE_Q2_reco;   //!
   TBranch        *b_QE_Ev_reco;   //!
   TBranch        *b_W_reco;   //!
   TBranch        *b_W2_reco;   //!
   TBranch        *b_Ev_reco;   //!
   TBranch        *b_EvB_reco;   //!
   TBranch        *b_EvUB_reco;   //!
   TBranch        *b_Ehad_reco;   //!
   TBranch        *b_EhadB_reco;   //!
   TBranch        *b_EhadUB_reco;   //!
   TBranch        *b_Ehad_unbias_default;   //!
   TBranch        *b_Ehad_unbias_numu;   //!
   TBranch        *b_Ehad_unbias_numubar;   //!
   TBranch        *b_Ehad_unbias_nue;   //!
   TBranch        *b_Ehad_unbias_nuebar;   //!
   TBranch        *b_Ehad_unbias_nc;   //!
   TBranch        *b_Ehad_unbias_ncbar;   //!
   TBranch        *b_p0v_reco;   //!
   TBranch        *b_pxv_reco;   //!
   TBranch        *b_pyv_reco;   //!
   TBranch        *b_pzv_reco;   //!
   TBranch        *b_EvClass_reco;   //!
   TBranch        *b_isQE_reco;   //!
   TBranch        *b_IsContv_reco;   //!
   TBranch        *b_hasElikeEM;   //!
   TBranch        *b_En;   //!
   TBranch        *b_p0n;   //!
   TBranch        *b_pxn;   //!
   TBranch        *b_pyn;   //!
   TBranch        *b_pzn;   //!
   TBranch        *b_Pdgl;   //!
   TBranch        *b_El;   //!
   TBranch        *b_p0l;   //!
   TBranch        *b_pxl;   //!
   TBranch        *b_pyl;   //!
   TBranch        *b_pzl;   //!
   TBranch        *b_pl;   //!
   TBranch        *b_cthl;   //!
   TBranch        *b_PdgL;   //!
   TBranch        *b_EL;   //!
   TBranch        *b_p0L;   //!
   TBranch        *b_pxL;   //!
   TBranch        *b_pyL;   //!
   TBranch        *b_pzL;   //!
   TBranch        *b_pL;   //!
   TBranch        *b_cthL;   //!
   TBranch        *b_ParentPdgL;   //!
   TBranch        *b_ParentEL;   //!
   TBranch        *b_p0LParent;   //!
   TBranch        *b_pxLParent;   //!
   TBranch        *b_pyLParent;   //!
   TBranch        *b_pzLParent;   //!
   TBranch        *b_pLParent;   //!
   TBranch        *b_cthLParent;   //!
   TBranch        *b_PdgL_reco;   //!
   TBranch        *b_EL_reco;   //!
   TBranch        *b_p0L_reco;   //!
   TBranch        *b_pxL_reco;   //!
   TBranch        *b_pyL_reco;   //!
   TBranch        *b_pzL_reco;   //!
   TBranch        *b_pL_reco;   //!
   TBranch        *b_cthL_reco;   //!
   TBranch        *b_TrkL_reco;   //!
   TBranch        *b_IsContL_reco;   //!
   TBranch        *b_Pdgl_reco;   //!
   TBranch        *b_El_reco;   //!
   TBranch        *b_p0l_reco;   //!
   TBranch        *b_pxl_reco;   //!
   TBranch        *b_pyl_reco;   //!
   TBranch        *b_pzl_reco;   //!
   TBranch        *b_pl_reco;   //!
   TBranch        *b_Trkl_reco;   //!
   TBranch        *b_cthl_reco;   //!
   TBranch        *b_nfp;   //!
   TBranch        *b_nfn;   //!
   TBranch        *b_nfpip;   //!
   TBranch        *b_nfpim;   //!
   TBranch        *b_nfpi0;   //!
   TBranch        *b_nfkp;   //!
   TBranch        *b_nfkm;   //!
   TBranch        *b_nfk0;   //!
   TBranch        *b_nfem;   //!
   TBranch        *b_nfother;   //!
   TBranch        *b_nip;   //!
   TBranch        *b_nin;   //!
   TBranch        *b_nipip;   //!
   TBranch        *b_nipim;   //!
   TBranch        *b_nipi0;   //!
   TBranch        *b_nikp;   //!
   TBranch        *b_nikm;   //!
   TBranch        *b_nik0;   //!
   TBranch        *b_niem;   //!
   TBranch        *b_niother;   //!
   TBranch        *b_ni;   //!
   TBranch        *b_pdgi;   //!
   TBranch        *b_resc;   //!
   TBranch        *b_Ei;   //!
   TBranch        *b_p0i;   //!
   TBranch        *b_pxi;   //!
   TBranch        *b_pyi;   //!
   TBranch        *b_pzi;   //!
   TBranch        *b_nf;   //!
   TBranch        *b_pdgf;   //!
   TBranch        *b_Ef;   //!
   TBranch        *b_p0f;   //!
   TBranch        *b_pxf;   //!
   TBranch        *b_pyf;   //!
   TBranch        *b_pzf;   //!
   TBranch        *b_pf;   //!
   TBranch        *b_cthf;   //!
   TBranch        *b_Ef_reco;   //!
   TBranch        *b_p0f_reco;   //!
   TBranch        *b_pxf_reco;   //!
   TBranch        *b_pyf_reco;   //!
   TBranch        *b_pzf_reco;   //!
   TBranch        *b_Trkf_reco;   //!
   TBranch        *b_cthf_reco;   //!
   TBranch        *b_pT_Imbalance;   //!
   TBranch        *b_nEMf;   //!
   TBranch        *b_EMpdgf;   //!
   TBranch        *b_EMp0f;   //!
   TBranch        *b_EMpxf;   //!
   TBranch        *b_EMpyf;   //!
   TBranch        *b_EMpzf;   //!
   TBranch        *b_EMpf;   //!
   TBranch        *b_EMcthf;   //!
   TBranch        *b_EMEf;   //!
   TBranch        *b_EMEf_reco;   //!
   TBranch        *b_EMp0f_reco;   //!
   TBranch        *b_EMpxf_reco;   //!
   TBranch        *b_EMpyf_reco;   //!
   TBranch        *b_EMpzf_reco;   //!
   TBranch        *b_EMTrkf_reco;   //!
   TBranch        *b_EMcthf_reco;   //!
   TBranch        *b_EMPizMass_reco;   //!
   TBranch        *b_EMElPizM_reco;   //!
   TBranch        *b_vtxx;   //!
   TBranch        *b_vtxy;   //!
   TBranch        *b_vtxz;   //!
   TBranch        *b_vtxt;   //!
   TBranch        *b_calresp0;   //!
   TBranch        *b_caltotal;   //!
   TBranch        *b_calneutron;   //!
   TBranch        *b_calproton;   //!
   TBranch        *b_calPion;   //!
   TBranch        *b_calmeson;   //!
   TBranch        *b_calbaryon;   //!
   TBranch        *b_calhadron;   //!
   TBranch        *b_calem;   //!
   TBranch        *b_calother;   //!
   TBranch        *b_calbinding;   //!
   TBranch        *b_Tau_Prob_nue;   //!
   TBranch        *b_Tau_Prob_numu;   //!
   TBranch        *b_rw7;   //!
   TBranch        *b_rw7_qelCCvalenciaRPA;   //!
   TBranch        *b_rw7_qelNCvalenciaRPA;   //!
   TBranch        *b_rw7_reslowq2;   //!
   TBranch        *b_frw7;   //!
   TBranch        *b_fNsyst;   //!
   TBranch        *b_frw7_TargetYTilt0_5mm;   //!
   TBranch        *b_frw7_TargetYOffset0_5mm;   //!
   TBranch        *b_frw7_TargetXTilt0_5mm;   //!
   TBranch        *b_frw7_TargetXOffset0_5mm;   //!
   TBranch        *b_frw7_TargetBe;   //!
   TBranch        *b_frw7_ProtonP90GeV;   //!
   TBranch        *b_frw7_ProtonP80GeV;   //!
   TBranch        *b_frw7_ProtonP70GeV;   //!
   TBranch        *b_frw7_ProtonP60GeV;   //!
   TBranch        *b_frw7_ProtonP50GeV;   //!
   TBranch        *b_frw7_ProtonP40GeV;   //!
   TBranch        *b_frw7_ProtonP30GeV;   //!
   TBranch        *b_frw7_ProtonP20GeV;   //!
   TBranch        *b_frw7_ProtonP150GeV;   //!
   TBranch        *b_frw7_ProtonP140GeV;   //!
   TBranch        *b_frw7_ProtonP130GeV;   //!
   TBranch        *b_frw7_ProtonP110GeV;   //!
   TBranch        *b_frw7_ProtonP100GeV;   //!
   TBranch        *b_frw7_Horn2YOffset0_5mm;   //!
   TBranch        *b_frw7_Horn2XOffset0_5mm;   //!
   TBranch        *b_frw7_Horn1YOffset0_5mm;   //!
   TBranch        *b_frw7_Horn1XOffset0_5mm;   //!
   TBranch        *b_frw7_DecayPipeRadius3m;   //!
   TBranch        *b_frw7_DecayPipeRadius1_9m;   //!
   TBranch        *b_frw7_DecayPipeLength250m;   //!
   TBranch        *b_frw7_BeamSigmaY1_4mm;   //!
   TBranch        *b_frw7_BeamSigmaX1_4mm;   //!
   TBranch        *b_frw7_700kW;   //!
   TBranch        *b_frw7_230kA;   //!

   FastMC(TTree *tree=0);
   virtual ~FastMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

