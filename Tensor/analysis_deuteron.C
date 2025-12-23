#include "Tensor_SIDIS_ND3.h"

using namespace std;

int main(int argc, char * argv[]){

  int nEvents = atoi(argv[1]);

  printf("nEvents = %i \n", nEvents);

  int RndSeed = atoi(argv[2]);

  printf("RndSeed = %i \n", RndSeed);
  
  gRandom->SetSeed(2);
  //gRandom->SetSeed(RndSeed);

  //set ouptupt file name with variable name based on random seed
  char out_file_pip[40];
  sprintf(out_file_pip, "output/coverage_11_pip_%i.root", RndSeed);

  
  //number of events:
  //default from original file: 100000000

  MakeKinematicCoveragePlots(nEvents, 11.0, "pi+", out_file_pip);
  //MakeKinematicCoveragePlots(nEvents, 11.0, "K+", "output/coverage_11_Kp.root");
  //MakeKinematicCoveragePlots(1000000000, 11.0, "pi-", "output/coverage_11m.root");
  //MakeKinematicCoveragePlots(100000, 8.8, "pi+", "output/coverage_8.8p.root");
  //MakeKinematicCoveragePlots(100000, 8.8, "pi-", "output/coverage_8.8m.root");
  


  return 0;
}

