#include "FastMC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

FastMC::FastMC(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("fastmcNtp_20140711_lbne_g4lbnev3r2p4_nuflux_numuflux_numu_LAr_1_g280_Ar40_5000_Default.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("data/fastmcNtp_20140711_lbne_g4lbnev3r2p4_nuflux_numuflux_numu_LAr_1_g280_Ar40_5000_Default.root");
      }
      f->GetObject("gst",tree);

   }
   Init(tree);
}

FastMC::~FastMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FastMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FastMC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FastMC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("FluxVersion", &FluxVersion, &b_FluxVersion);
   fChain->SetBranchAddress("iev", &iev, &b_iev);
   fChain->SetBranchAddress("neu", &neu, &b_neu);
   fChain->SetBranchAddress("fspl", &fspl, &b_fspl);
   fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
   fChain->SetBranchAddress("Z", &Z, &b_Z);
   fChain->SetBranchAddress("A", &A, &b_A);
   fChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
   fChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
   fChain->SetBranchAddress("resid", &resid, &b_resid);
   fChain->SetBranchAddress("sea", &sea, &b_sea);
   fChain->SetBranchAddress("qel", &qel, &b_qel);
   fChain->SetBranchAddress("mec", &mec, &b_mec);
   fChain->SetBranchAddress("res", &res, &b_res);
   fChain->SetBranchAddress("dis", &dis, &b_dis);
   fChain->SetBranchAddress("coh", &coh, &b_coh);
   fChain->SetBranchAddress("cohe", &cohe, &b_cohe);
   fChain->SetBranchAddress("dfr", &dfr, &b_dfr);
   fChain->SetBranchAddress("imd", &imd, &b_imd);
   fChain->SetBranchAddress("imda", &imda, &b_imda);
   fChain->SetBranchAddress("nuel", &nuel, &b_nuel);
   fChain->SetBranchAddress("em", &em, &b_em);
   fChain->SetBranchAddress("cc", &cc, &b_cc);
   fChain->SetBranchAddress("nc", &nc, &b_nc);
   fChain->SetBranchAddress("charm", &charm, &b_charm);
   fChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
   fChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
   fChain->SetBranchAddress("wght", &wght, &b_wght);
   fChain->SetBranchAddress("POTWeight", &POTWeight, &b_POTWeight);
   fChain->SetBranchAddress("POTperYear", &POTperYear, &b_POTperYear);
   fChain->SetBranchAddress("nomPOTperYear", &nomPOTperYear, &b_nomPOTperYear);
   fChain->SetBranchAddress("OSCWeight", &OSCWeight, &b_OSCWeight);
   fChain->SetBranchAddress("xs", &xs, &b_xs);
   fChain->SetBranchAddress("ys", &ys, &b_ys);
   fChain->SetBranchAddress("ts", &ts, &b_ts);
   fChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
   fChain->SetBranchAddress("Ws", &Ws, &b_Ws);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("P", &P, &b_P);
   fChain->SetBranchAddress("ET", &ET, &b_ET);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("Ev", &Ev, &b_Ev);
   fChain->SetBranchAddress("Ehad", &Ehad, &b_Ehad);
   fChain->SetBranchAddress("Evtx", &Evtx, &b_Evtx);
   fChain->SetBranchAddress("p0v", &p0v, &b_p0v);
   fChain->SetBranchAddress("pxv", &pxv, &b_pxv);
   fChain->SetBranchAddress("pyv", &pyv, &b_pyv);
   fChain->SetBranchAddress("pzv", &pzv, &b_pzv);
   fChain->SetBranchAddress("x_reco", &x_reco, &b_x_reco);
   fChain->SetBranchAddress("y_reco", &y_reco, &b_y_reco);
   fChain->SetBranchAddress("t_reco", &t_reco, &b_t_reco);
   fChain->SetBranchAddress("Q2_reco", &Q2_reco, &b_Q2_reco);
   fChain->SetBranchAddress("QE_Q2_reco", &QE_Q2_reco, &b_QE_Q2_reco);
   fChain->SetBranchAddress("QE_Ev_reco", &QE_Ev_reco, &b_QE_Ev_reco);
   fChain->SetBranchAddress("W_reco", &W_reco, &b_W_reco);
   fChain->SetBranchAddress("W2_reco", &W2_reco, &b_W2_reco);
   fChain->SetBranchAddress("Ev_reco", &Ev_reco, &b_Ev_reco);
   fChain->SetBranchAddress("EvB_reco", &EvB_reco, &b_EvB_reco);
   fChain->SetBranchAddress("EvUB_reco", &EvUB_reco, &b_EvUB_reco);
   fChain->SetBranchAddress("Ehad_reco", &Ehad_reco, &b_Ehad_reco);
   fChain->SetBranchAddress("EhadB_reco", &EhadB_reco, &b_EhadB_reco);
   fChain->SetBranchAddress("EhadUB_reco", &EhadUB_reco, &b_EhadUB_reco);
   fChain->SetBranchAddress("Ehad_unbias_default", &Ehad_unbias_default, &b_Ehad_unbias_default);
   fChain->SetBranchAddress("Ehad_unbias_numu", &Ehad_unbias_numu, &b_Ehad_unbias_numu);
   fChain->SetBranchAddress("Ehad_unbias_numubar", &Ehad_unbias_numubar, &b_Ehad_unbias_numubar);
   fChain->SetBranchAddress("Ehad_unbias_nue", &Ehad_unbias_nue, &b_Ehad_unbias_nue);
   fChain->SetBranchAddress("Ehad_unbias_nuebar", &Ehad_unbias_nuebar, &b_Ehad_unbias_nuebar);
   fChain->SetBranchAddress("Ehad_unbias_nc", &Ehad_unbias_nc, &b_Ehad_unbias_nc);
   fChain->SetBranchAddress("Ehad_unbias_ncbar", &Ehad_unbias_ncbar, &b_Ehad_unbias_ncbar);
   fChain->SetBranchAddress("p0v_reco", &p0v_reco, &b_p0v_reco);
   fChain->SetBranchAddress("pxv_reco", &pxv_reco, &b_pxv_reco);
   fChain->SetBranchAddress("pyv_reco", &pyv_reco, &b_pyv_reco);
   fChain->SetBranchAddress("pzv_reco", &pzv_reco, &b_pzv_reco);
   fChain->SetBranchAddress("EvClass_reco", &EvClass_reco, &b_EvClass_reco);
   fChain->SetBranchAddress("isQE_reco", &isQE_reco, &b_isQE_reco);
   fChain->SetBranchAddress("IsContv_reco", &IsContv_reco, &b_IsContv_reco);
   fChain->SetBranchAddress("hasElikeEM", &hasElikeEM, &b_hasElikeEM);
   fChain->SetBranchAddress("En", &En, &b_En);
   fChain->SetBranchAddress("p0n", &p0n, &b_p0n);
   fChain->SetBranchAddress("pxn", &pxn, &b_pxn);
   fChain->SetBranchAddress("pyn", &pyn, &b_pyn);
   fChain->SetBranchAddress("pzn", &pzn, &b_pzn);
   fChain->SetBranchAddress("Pdgl", &Pdgl, &b_Pdgl);
   fChain->SetBranchAddress("El", &El, &b_El);
   fChain->SetBranchAddress("p0l", &p0l, &b_p0l);
   fChain->SetBranchAddress("pxl", &pxl, &b_pxl);
   fChain->SetBranchAddress("pyl", &pyl, &b_pyl);
   fChain->SetBranchAddress("pzl", &pzl, &b_pzl);
   fChain->SetBranchAddress("pl", &pl, &b_pl);
   fChain->SetBranchAddress("cthl", &cthl, &b_cthl);
   fChain->SetBranchAddress("PdgL", &PdgL, &b_PdgL);
   fChain->SetBranchAddress("EL", &EL, &b_EL);
   fChain->SetBranchAddress("p0L", &p0L, &b_p0L);
   fChain->SetBranchAddress("pxL", &pxL, &b_pxL);
   fChain->SetBranchAddress("pyL", &pyL, &b_pyL);
   fChain->SetBranchAddress("pzL", &pzL, &b_pzL);
   fChain->SetBranchAddress("pL", &pL, &b_pL);
   fChain->SetBranchAddress("cthL", &cthL, &b_cthL);
   fChain->SetBranchAddress("ParentPdgL", &ParentPdgL, &b_ParentPdgL);
   fChain->SetBranchAddress("ParentEL", &ParentEL, &b_ParentEL);
   fChain->SetBranchAddress("p0LParent", &p0LParent, &b_p0LParent);
   fChain->SetBranchAddress("pxLParent", &pxLParent, &b_pxLParent);
   fChain->SetBranchAddress("pyLParent", &pyLParent, &b_pyLParent);
   fChain->SetBranchAddress("pzLParent", &pzLParent, &b_pzLParent);
   fChain->SetBranchAddress("pLParent", &pLParent, &b_pLParent);
   fChain->SetBranchAddress("cthLParent", &cthLParent, &b_cthLParent);
   fChain->SetBranchAddress("PdgL_reco", &PdgL_reco, &b_PdgL_reco);
   fChain->SetBranchAddress("EL_reco", &EL_reco, &b_EL_reco);
   fChain->SetBranchAddress("p0L_reco", &p0L_reco, &b_p0L_reco);
   fChain->SetBranchAddress("pxL_reco", &pxL_reco, &b_pxL_reco);
   fChain->SetBranchAddress("pyL_reco", &pyL_reco, &b_pyL_reco);
   fChain->SetBranchAddress("pzL_reco", &pzL_reco, &b_pzL_reco);
   fChain->SetBranchAddress("pL_reco", &pL_reco, &b_pL_reco);
   fChain->SetBranchAddress("cthL_reco", &cthL_reco, &b_cthL_reco);
   fChain->SetBranchAddress("TrkL_reco", &TrkL_reco, &b_TrkL_reco);
   fChain->SetBranchAddress("IsContL_reco", &IsContL_reco, &b_IsContL_reco);
   fChain->SetBranchAddress("Pdgl_reco", &Pdgl_reco, &b_Pdgl_reco);
   fChain->SetBranchAddress("El_reco", &El_reco, &b_El_reco);
   fChain->SetBranchAddress("p0l_reco", &p0l_reco, &b_p0l_reco);
   fChain->SetBranchAddress("pxl_reco", &pxl_reco, &b_pxl_reco);
   fChain->SetBranchAddress("pyl_reco", &pyl_reco, &b_pyl_reco);
   fChain->SetBranchAddress("pzl_reco", &pzl_reco, &b_pzl_reco);
   fChain->SetBranchAddress("pl_reco", &pl_reco, &b_pl_reco);
   fChain->SetBranchAddress("Trkl_reco", &Trkl_reco, &b_Trkl_reco);
   fChain->SetBranchAddress("cthl_reco", &cthl_reco, &b_cthl_reco);
   fChain->SetBranchAddress("nfp", &nfp, &b_nfp);
   fChain->SetBranchAddress("nfn", &nfn, &b_nfn);
   fChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
   fChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
   fChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
   fChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
   fChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
   fChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
   fChain->SetBranchAddress("nfem", &nfem, &b_nfem);
   fChain->SetBranchAddress("nfother", &nfother, &b_nfother);
   fChain->SetBranchAddress("nip", &nip, &b_nip);
   fChain->SetBranchAddress("nin", &nin, &b_nin);
   fChain->SetBranchAddress("nipip", &nipip, &b_nipip);
   fChain->SetBranchAddress("nipim", &nipim, &b_nipim);
   fChain->SetBranchAddress("nipi0", &nipi0, &b_nipi0);
   fChain->SetBranchAddress("nikp", &nikp, &b_nikp);
   fChain->SetBranchAddress("nikm", &nikm, &b_nikm);
   fChain->SetBranchAddress("nik0", &nik0, &b_nik0);
   fChain->SetBranchAddress("niem", &niem, &b_niem);
   fChain->SetBranchAddress("niother", &niother, &b_niother);
   fChain->SetBranchAddress("ni", &ni, &b_ni);
   fChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
   fChain->SetBranchAddress("resc", resc, &b_resc);
   fChain->SetBranchAddress("Ei", Ei, &b_Ei);
   fChain->SetBranchAddress("p0i", p0i, &b_p0i);
   fChain->SetBranchAddress("pxi", pxi, &b_pxi);
   fChain->SetBranchAddress("pyi", pyi, &b_pyi);
   fChain->SetBranchAddress("pzi", pzi, &b_pzi);
   fChain->SetBranchAddress("nf", &nf, &b_nf);
   fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
   fChain->SetBranchAddress("Ef", Ef, &b_Ef);
   fChain->SetBranchAddress("p0f", p0f, &b_p0f);
   fChain->SetBranchAddress("pxf", pxf, &b_pxf);
   fChain->SetBranchAddress("pyf", pyf, &b_pyf);
   fChain->SetBranchAddress("pzf", pzf, &b_pzf);
   fChain->SetBranchAddress("pf", pf, &b_pf);
   fChain->SetBranchAddress("cthf", cthf, &b_cthf);
   fChain->SetBranchAddress("Ef_reco", Ef_reco, &b_Ef_reco);
   fChain->SetBranchAddress("p0f_reco", p0f_reco, &b_p0f_reco);
   fChain->SetBranchAddress("pxf_reco", pxf_reco, &b_pxf_reco);
   fChain->SetBranchAddress("pyf_reco", pyf_reco, &b_pyf_reco);
   fChain->SetBranchAddress("pzf_reco", pzf_reco, &b_pzf_reco);
   fChain->SetBranchAddress("Trkf_reco", Trkf_reco, &b_Trkf_reco);
   fChain->SetBranchAddress("cthf_reco", cthf_reco, &b_cthf_reco);
   fChain->SetBranchAddress("pT_Imbalance", &pT_Imbalance, &b_pT_Imbalance);
   fChain->SetBranchAddress("nEMf", &nEMf, &b_nEMf);
   fChain->SetBranchAddress("EMpdgf", EMpdgf, &b_EMpdgf);
   fChain->SetBranchAddress("EMp0f", EMp0f, &b_EMp0f);
   fChain->SetBranchAddress("EMpxf", EMpxf, &b_EMpxf);
   fChain->SetBranchAddress("EMpyf", EMpyf, &b_EMpyf);
   fChain->SetBranchAddress("EMpzf", EMpzf, &b_EMpzf);
   fChain->SetBranchAddress("EMpf", EMpf, &b_EMpf);
   fChain->SetBranchAddress("EMcthf", EMcthf, &b_EMcthf);
   fChain->SetBranchAddress("EMEf", EMEf, &b_EMEf);
   fChain->SetBranchAddress("EMEf_reco", EMEf_reco, &b_EMEf_reco);
   fChain->SetBranchAddress("EMp0f_reco", EMp0f_reco, &b_EMp0f_reco);
   fChain->SetBranchAddress("EMpxf_reco", EMpxf_reco, &b_EMpxf_reco);
   fChain->SetBranchAddress("EMpyf_reco", EMpyf_reco, &b_EMpyf_reco);
   fChain->SetBranchAddress("EMpzf_reco", EMpzf_reco, &b_EMpzf_reco);
   fChain->SetBranchAddress("EMTrkf_reco", EMTrkf_reco, &b_EMTrkf_reco);
   fChain->SetBranchAddress("EMcthf_reco", EMcthf_reco, &b_EMcthf_reco);
   fChain->SetBranchAddress("EMPizMass_reco", EMPizMass_reco, &b_EMPizMass_reco);
   fChain->SetBranchAddress("EMElPizM_reco", EMElPizM_reco, &b_EMElPizM_reco);
   fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
   fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
   fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
   fChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
   fChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);
   fChain->SetBranchAddress("caltotal", &caltotal, &b_caltotal);
   fChain->SetBranchAddress("calneutron", &calneutron, &b_calneutron);
   fChain->SetBranchAddress("calproton", &calproton, &b_calproton);
   fChain->SetBranchAddress("calpion", &calpion, &b_calPion);
   fChain->SetBranchAddress("calmeson", &calmeson, &b_calmeson);
   fChain->SetBranchAddress("calbaryon", &calbaryon, &b_calbaryon);
   fChain->SetBranchAddress("calhadron", &calhadron, &b_calhadron);
   fChain->SetBranchAddress("calem", &calem, &b_calem);
   fChain->SetBranchAddress("calother", &calother, &b_calother);
   fChain->SetBranchAddress("calbinding", &calbinding, &b_calbinding);
   fChain->SetBranchAddress("Tau_Prob_nue", &Tau_Prob_nue, &b_Tau_Prob_nue);
   fChain->SetBranchAddress("Tau_Prob_numu", &Tau_Prob_numu, &b_Tau_Prob_numu);
   fChain->SetBranchAddress("rw7", &rw7, &b_rw7);
   fChain->SetBranchAddress("rw7_qelCCvalenciaRPA", rw7_qelCCvalenciaRPA, &b_rw7_qelCCvalenciaRPA);
   fChain->SetBranchAddress("rw7_qelNCvalenciaRPA", rw7_qelNCvalenciaRPA, &b_rw7_qelNCvalenciaRPA);
   fChain->SetBranchAddress("rw7_reslowq2", rw7_reslowq2, &b_rw7_reslowq2);
   fChain->SetBranchAddress("frw7", &frw7, &b_frw7);
   fChain->SetBranchAddress("fNsyst", &fNsyst, &b_fNsyst);
   fChain->SetBranchAddress("frw7_TargetYTilt0.5mm", frw7_TargetYTilt0_5mm, &b_frw7_TargetYTilt0_5mm);
   fChain->SetBranchAddress("frw7_TargetYOffset0.5mm", frw7_TargetYOffset0_5mm, &b_frw7_TargetYOffset0_5mm);
   fChain->SetBranchAddress("frw7_TargetXTilt0.5mm", frw7_TargetXTilt0_5mm, &b_frw7_TargetXTilt0_5mm);
   fChain->SetBranchAddress("frw7_TargetXOffset0.5mm", frw7_TargetXOffset0_5mm, &b_frw7_TargetXOffset0_5mm);
   fChain->SetBranchAddress("frw7_TargetBe", frw7_TargetBe, &b_frw7_TargetBe);
   fChain->SetBranchAddress("frw7_ProtonP90GeV", frw7_ProtonP90GeV, &b_frw7_ProtonP90GeV);
   fChain->SetBranchAddress("frw7_ProtonP80GeV", frw7_ProtonP80GeV, &b_frw7_ProtonP80GeV);
   fChain->SetBranchAddress("frw7_ProtonP70GeV", frw7_ProtonP70GeV, &b_frw7_ProtonP70GeV);
   fChain->SetBranchAddress("frw7_ProtonP60GeV", frw7_ProtonP60GeV, &b_frw7_ProtonP60GeV);
   fChain->SetBranchAddress("frw7_ProtonP50GeV", frw7_ProtonP50GeV, &b_frw7_ProtonP50GeV);
   fChain->SetBranchAddress("frw7_ProtonP40GeV", frw7_ProtonP40GeV, &b_frw7_ProtonP40GeV);
   fChain->SetBranchAddress("frw7_ProtonP30GeV", frw7_ProtonP30GeV, &b_frw7_ProtonP30GeV);
   fChain->SetBranchAddress("frw7_ProtonP20GeV", frw7_ProtonP20GeV, &b_frw7_ProtonP20GeV);
   fChain->SetBranchAddress("frw7_ProtonP150GeV", frw7_ProtonP150GeV, &b_frw7_ProtonP150GeV);
   fChain->SetBranchAddress("frw7_ProtonP140GeV", frw7_ProtonP140GeV, &b_frw7_ProtonP140GeV);
   fChain->SetBranchAddress("frw7_ProtonP130GeV", frw7_ProtonP130GeV, &b_frw7_ProtonP130GeV);
   fChain->SetBranchAddress("frw7_ProtonP110GeV", frw7_ProtonP110GeV, &b_frw7_ProtonP110GeV);
   fChain->SetBranchAddress("frw7_ProtonP100GeV", frw7_ProtonP100GeV, &b_frw7_ProtonP100GeV);
   fChain->SetBranchAddress("frw7_Horn2YOffset0.5mm", frw7_Horn2YOffset0_5mm, &b_frw7_Horn2YOffset0_5mm);
   fChain->SetBranchAddress("frw7_Horn2XOffset0.5mm", frw7_Horn2XOffset0_5mm, &b_frw7_Horn2XOffset0_5mm);
   fChain->SetBranchAddress("frw7_Horn1YOffset0.5mm", frw7_Horn1YOffset0_5mm, &b_frw7_Horn1YOffset0_5mm);
   fChain->SetBranchAddress("frw7_Horn1XOffset0.5mm", frw7_Horn1XOffset0_5mm, &b_frw7_Horn1XOffset0_5mm);
   fChain->SetBranchAddress("frw7_DecayPipeRadius3m", frw7_DecayPipeRadius3m, &b_frw7_DecayPipeRadius3m);
   fChain->SetBranchAddress("frw7_DecayPipeRadius1.9m", frw7_DecayPipeRadius1_9m, &b_frw7_DecayPipeRadius1_9m);
   fChain->SetBranchAddress("frw7_DecayPipeLength250m", frw7_DecayPipeLength250m, &b_frw7_DecayPipeLength250m);
   fChain->SetBranchAddress("frw7_BeamSigmaY1.4mm", frw7_BeamSigmaY1_4mm, &b_frw7_BeamSigmaY1_4mm);
   fChain->SetBranchAddress("frw7_BeamSigmaX1.4mm", frw7_BeamSigmaX1_4mm, &b_frw7_BeamSigmaX1_4mm);
   fChain->SetBranchAddress("frw7_700kW", frw7_700kW, &b_frw7_700kW);
   fChain->SetBranchAddress("frw7_230kA", frw7_230kA, &b_frw7_230kA);
   Notify();
}

Bool_t FastMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FastMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FastMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void FastMC::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L FastMC.C
//      Root > FastMC t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
