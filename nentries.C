#include "../mz-unpacker-BarrelDirc/TPrtHit.h"
#include "../mz-unpacker-BarrelDirc/TPrtEvent.h"

void nentries(TString inFile = "T.root"){  
  // std::ofstream lStream( "/dev/null" );
  // std::cout.rdbuf( lStream.rdbuf() );
  
  TChain *ch = new TChain("T");
  ch->Add(inFile);
  std::cout<<"Entries in chain:  "<< ch->GetEntries()<<std::endl;

  // lStream.close();
}


