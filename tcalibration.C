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
const Int_t maxch =3000;
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
Int_t chmap[nmcp][npix];

Int_t gComboId=0;
TGraph *gGr[maxch];
TH1F  *hL = new TH1F("hL", "hL" , 500,-100,-50);


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
  for(Int_t ch=0; ch<maxch; ch++){
    Int_t mcp = ch/128;
    Int_t pix = (ch - mcp*128)/2;
    Int_t col = 7-(pix/2 - 8*(pix/16));
    Int_t row = pix%2 + 2*(pix/16);
    pix = col*8+row;
    //std::cout<<"ch  "<<ch <<"  m "<<mcp <<"  p  "<<pix <<std::endl;
    
    chmap[mcp][pix]=ch;
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

  TFile f("../../data/calib.root");
  TIter nextkey(f.GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)nextkey())) {
    TGraph *gr = (TGraph*)key->ReadObj();
    TString name = gr->GetName();
    Int_t channel = name.Atoi();
    
    gGr[channel]= new TGraph(*gr);
    
    // TCanvas *can = new TCanvas("can","can",800,500);
    // std::cout<<" gGr["<<channel<<"]  "<<gGr[channel]->GetN()  <<std::endl;   
    // gGr[channel]->Draw("AL");
    // can->Modified();
    // can->Update();
    // can->WaitPrimitive();
  }
  f.Close();


}

//Double_t time[50000];
Bool_t TTSelector::Process(Long64_t entry){
  Int_t trbSeqId,ch,mcp,pix,col,row;;
  Double_t timeTot(0), grTime0(0), grTime1(0),timeLe(0), timeTe(0);
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  Double_t time[50000];
  GetEntry(entry);
  
  fEvent = new TPrtEvent();
  fEvent->SetReferenceChannel(gTrigger);
  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = 32*trbSeqId+Hits_nTdcChannel[i];
    //std::cout<<"ch  "<<ch <<std::endl;
    
    if(++mult[ch]>50) continue;
    Double_t coarseTime = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);
 
    Int_t chid = ch;
    //if(ch%2==0 && ch!=0) chid = ch-1;
    
    //time[i] = Hits_fTime[i] ;//
    time[i] = coarseTime-(Hits_nFineTime[i]-31)*0.0102;//0.0102;
    //time[i] = coarseTime-gGr[chid]->Eval(Hits_nFineTime[i]);

    if(++mult[ch]>50) continue;
    timeTe0[ch][mult[ch]]=time[i];
    if(Hits_nTdcChannel[i]==0) {  // is ref channel
      trbRefTime[trbSeqId] = time[i];
      if((gTrigger-ch)<=32 && (gTrigger-ch)>0) grTime0 = time[i];
    }
    if(ch==gTrigger) grTime1 = time[i];
  }

  if((grTime0>0 && grTime1>0) || gTrigger==0){
    for(Int_t i=0; i<Hits_; i++){
      if(Hits_nTrbAddress[i]==0) continue;
      trbSeqId = tdcmap[Hits_nTrbAddress[i]];
      ch = 32*trbSeqId+Hits_nTdcChannel[i];
      mcp = ch/128;
      pix = (ch%128)/2;	
      col = pix/2 - 8*(pix/16);
      row = pix%2 + 2*(pix/16);
      pix = (7-col)*8+row;

      if(ch%2==0) continue; // go away trailing edge
      if(ch<3000) {
	  timeLe = time[i]-trbRefTime[trbSeqId];
	  timeTe = timeTe0[ch+1][0]-trbRefTime[trbSeqId];
	  
	  timeLe = timeLe - (grTime1-grTime0);
	  if(ch == 363) hL->Fill(timeLe);
          timeTot = timeTe0[ch+1][0] - timeTe0[ch][0]; // timeTe - (grTime1-grTime0)
	  TPrtHit hit(Hits_nTrbAddress[i],Hits_nTdcChannel[i],ch,mcp,pix+1,timeLe,timeTot);
	  fEvent->AddHit(hit);
      }
    }
  }

  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = 32*trbSeqId+Hits_nTdcChannel[i];
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

  hL->Draw();
  hL->Fit("gaus","V","E1",-80,-70);
  gGr[363]->Draw("AP");
}


void tcalibration(TString inFile= "../../data/cc2.hld.root", Int_t trigger=1921, Int_t mode =0){
  ginFile = inFile;
  gTrigger = trigger;
  gMode=mode;


  TChain* ch = new TChain("T");
  ch->Add(ginFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
 
  TTSelector *selector = new TTSelector();
  TString option = Form("%d %d",gTrigger,gMode);
  
  ch->Process(selector,option,entries);
}
