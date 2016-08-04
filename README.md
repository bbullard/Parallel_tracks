#Code to study parallel tracks in DUNE
Compile using 
```
root -b -q loadClasses.C
```
Run using 
```
root -b loadClasses.C run.C
```
Edit `run.C` and `MakePlot.C` to change analysis. FastMC files are linked in `FastMC.h`, and can be found [here] (http://www.phy.bnl.gov/bviren/data/fmcfiles/). This code is adapted from Chao's, which can be found [here] (https://github.com/czczc/DUNE-Parallel-Events).
