// tcalibration - routine for the prtdirc data calibration 
// original author: Roman Dzhygadlo - GSI Darmstad


#define TTSelector_cxx
#include "prttools.C"

#include "TStyle.h"
#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGStatusBar.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGProgressBar.h>
#include <TThread.h>
#include <TProof.h>
#include <TGSplitter.h>
#include <TChainElement.h>
#include <TKey.h>

#include "tcalibration.h"

const Int_t maxfiles(200);
const Int_t maxChannel(3000);
const Int_t nmcp(15), npix(64);
TString fileList[maxfiles];


TH1F *hCh;

TString ginFile="";
Int_t gTrigger(0), gMode(0);;

const Int_t tdcnum(88);
const Int_t tdcmax(10000);
TString trbsid[tdcnum] = 
  {"0010","0011","0012","0013","0110","0111","0112","0113","0210","0211","0212","0213","0310","0311","0312","0313","0410","0411","0412","0413"
   ,"0510","0511","0512","0513","0610","0611","0612","0613","0710","0711","0712","0713","0810","0811","0812","0813","0910","0911","0912","0913"
   ,"1010","1011","1012","1013","1110","1111","1112","1113","1210","1211","1212","1213","1310","1311","1312","1313","1410","1411","1412","1413"
   ,"1510","1511","1512","1513","1610","1611","1612","1613","1710","1711","1712","1713","1810","1811","1812","1813","1910","1911","1912","1913"
   ,"2010","2011","2012","2013","2110","2111","2112","2113"};


Int_t tdcid[tdcnum];
Double_t trbRefTime[tdcnum];

Double_t timeTe0[tdcmax][50];
Int_t mult[tdcmax];

Int_t tdcmap[tdcmax];
Int_t mcpmap[tdcmax];
Int_t pixmap[tdcmax];

Int_t gComboId=0;
TGraph *gGr[nmcp][npix];

void CreateMap(){
  Int_t seqid =0;
  for(Int_t i=0; i<tdcmax; i++){
    tdcmap[i]=-1;
    mcpmap[i]=-1;
    for(Int_t j=0; j<tdcnum; j++){
      if(i==TString::BaseConvert(trbsid[j],16,10).Atoi()){
	tdcmap[i]=seqid++;
	mcpmap[i]=j/4;
	pixmap[i]=j-j/32;
	break;
      }
    }
  }
}

void TTSelector::Begin(TTree *){
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  gTrigger = ((TObjString*)strobj->At(0))->GetString().Atoi();
  gMode = ((TObjString*)strobj->At(1))->GetString().Atoi();
  CreateMap();
  TString filedir=ginFile;
  filedir.Remove(filedir.Last('.'));
  fFile = new TFile(filedir+"C.root","RECREATE");
  fTree = new TTree("M","Tree for GSI Prt Analysis");  
  fEvent = new TPrtEvent();
  fTree->Branch("TPrtEvent", "TPrtEvent", &fEvent, 64000, 2);
}

Bool_t TTSelector::Process(Long64_t entry){
  Int_t trbSeqId,ch;
  Double_t timeTot(0), grTime0=0, grTime1=0,timeLe=0, timeTe=0;
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  
  GetEntry(entry);
  
  fEvent = new TPrtEvent();
  fEvent->SetReferenceChannel(1920);
  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = 32*trbSeqId+Hits_nTdcChannel[i];
    if(++mult[ch]>50) continue;
    timeTe0[ch][mult[ch]]=Hits_fTime[i];
    if(Hits_nTdcChannel[i]==0 && Hits_bIsRefChannel[i]) {
      trbRefTime[trbSeqId] = Hits_fTime[i];
      if((ch-gTrigger)<64 && (ch-gTrigger)>=0) grTime0 = Hits_fTime[i];
    }
    if(ch==gTrigger+1) grTime1 = Hits_fTime[i];
  }

  if((grTime0>0 && grTime1>0) || gTrigger==0){
    for(Int_t i=0; i<Hits_; i++){
      // Double_t fHitTimeCoarse = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);
      trbSeqId = tdcmap[Hits_nTrbAddress[i]];
      ch = 32*trbSeqId+Hits_nTdcChannel[i];
      Int_t mcp = ch/128;
      Int_t pix = (ch - mcp*128)/2;
      Int_t col = 7-(pix/2 - 8*(pix/16));
      Int_t row = pix%2 + 2*(pix/16);

      if(ch%2==0) continue; // go away trailing edge
      if(ch<3000 && !Hits_bIsRefChannel[i]) {
	pix = col*8+row;
	// bad pixels
	// if(mcp==2  && pix==55) continue;
	// if(mcp==2  && pix==62) continue;
	// if(mcp==13 && pix==62) continue;
	// if(mcp==14 && pix==28) continue;
	// if(mcp==10 && pix==46) continue;
	//if(mcp<15)
	{

	  timeLe = Hits_fTime[i]-trbRefTime[trbSeqId];
	  timeTe = timeTe0[ch+1][0]-trbRefTime[trbSeqId];
	  
	  timeLe = timeLe - (grTime1-grTime0);
          timeTot = timeTe0[ch+1][0] - timeTe0[ch][0]; // timeTe - (grTime1-grTime0)
	  TPrtHit hit(Hits_nTrbAddress[i],Hits_nTdcChannel[i],ch,mcp,pix+1,timeLe,timeTot); // tdcId,tdcCh,ch,mcp,pix,timeLe,timeTot
	  fEvent->AddHit(hit);
	}
      }
    }
  }

  for(Int_t i=0; i<Hits_; i++){
    Int_t trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    Int_t ch = 32*trbSeqId+Hits_nTdcChannel[i];
    mult[ch]=-1;
    for(Int_t j=0; j<50; j++){
      timeTe0[ch][j]=0; 
    }
  }
  fTree->Fill();
  fEvent->Clear();
  delete fEvent;

  return kTRUE;
}

void TTSelector::Terminate(){
  fFile->Write();
  fFile->Close();
}


void tcalibration(TString inFile= "../../data/cc2.hld.root", Int_t trigger=1920, Int_t mode =0){
  ginFile = inFile;
  gTrigger = trigger;
  gMode=mode;


  TFile f("calib.root");
  TIter nextkey(f.GetListOfKeys());
  TKey *key;

  const Int_t nmcp = 15, npix = 64;
  TGraph *gGr[nmcp][npix];
  while (key = (TKey*)nextkey()) {
    TGraph *gr = (TGraph*)key->ReadObj();
    TString name = gr->GetName();
    TObjArray *sarr = name.Tokenize("_");
        
    Int_t mcp = ((TObjString*)sarr->At(0))->GetString().Atoi();
    Int_t pix = ((TObjString*)sarr->At(1))->GetString().Atoi();
    
    // gr->Draw("AL");
    // can->Modified();
    // can->Update();
    // can->WaitPrimitive();

    gGr[mcp][pix]=(TGraph*)key->ReadObj();
  }
  f.Close();

  TChain* ch = new TChain("T");
  ch->Add(ginFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
 
  TTSelector *selector = new TTSelector();
  TString option = Form("%d %d",gTrigger,gMode);
  std::cout<<"1111 " <<std::endl;
  
  ch->Process(selector,option,entries);
}
