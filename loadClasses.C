void loadClasses()
{
    TString include = ".include ";
    TString load = ".L ";
    TString pwd = gSystem->pwd();
    pwd = pwd + '/';
//    gSystem->Load("libGLEW");
//    gSystem->Load("libFTGL");

    TString prefix;
    prefix = ".";
    gROOT->ProcessLine( include + pwd + prefix );
    gROOT->ProcessLine( load + prefix + "/FastMC.C+" );
    gROOT->ProcessLine( load + prefix + "/MakePlot.C+" );
}
