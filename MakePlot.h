#ifndef MakePlot_h
#define MakePlot_h

class FastMC;
class TCanvas;

class MakePlot {
public :
   MakePlot(const Double_t angle=7.5, const Bool_t r=0);
   virtual ~MakePlot();

   FastMC *fm;
   TCanvas *c1;
   Double_t DetMass = 34.0;
   Double_t exposure = 3.0;
   Double_t POTperYr = 1.1e21;
   Double_t n_events;
   Double_t cutAngle;
   Bool_t rot;

   void PlotEv_vs_Angle(const char* saveAs=0);
   void PlotEl_vs_Angle(const char* saveAs=0);
   void PlotEl_vs_Length(const char* saveAs=0);


   void PlotEv(const char* saveAs=0);
   void PlotOscWeight_vs_Ev(const char* saveAs=0);
   void PlotEl_vs_Ev(const char* saveAs=0);

   void PlotAngle_cut_E(const char* saveAs=0);
   void PlotAngle_cut_E_Perp(const char* saveAs=0);
   void PlotLength_cut_E(const char* saveAs=0);
   void PlotLength_cut_E_Perp(const char* saveAs=0);

   // hadron side
   void PlotEv_vs_Ef(const char* saveAs=0);
   void PlotNf_vs_Ratio(const char* particle, const char* title, const char* saveAs=0);
   void PlotAnglel_vs_Angleh(const char* saveAs=0);
   void PlotMinAngle_lepton(const char* saveAs=0);
   void PlotMinAngle_hadron(const char* saveAs=0);
   void PlotMinAngle_photon(const char* saveAs=0, const float l_max = 3);
   void PlotMinAngle_photonNC(const char* saveAs=0, const float l_max = 3);
   
   void PlotEvents_vs_El_l(const char* saveAs=0);
   void PlotEvents_vs_El_h(const char* saveAs=0);
   void PlotEvents_vs_El_lh(const char* saveAs=0);
   
   // create new tree
   void CreateTree(const char* filename);

   void CreateCCEMTree(const char* saveAs=0);
   void CreateNCEMTree(const char* saveAs=0);

   void TabulateNCEvents();
   
   void SetStyle();

};

#endif
