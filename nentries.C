#define prt__beam
#include "prttools.C"

#include "../mz-unpacker-BarrelDirc/TPrtHit.h"
#include "../mz-unpacker-BarrelDirc/TPrtEvent.h"

void nentries(TString inFile = "M.root"){  
  // std::ofstream lStream( "/dev/null" );
  // std::cout.rdbuf( lStream.rdbuf() );
  PrtInit(inFile,0);
  // lStream.close();
}


