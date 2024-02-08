//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <map>
#include "sys/types.h"
#include "dirent.h"

#include "math.h"
#include "string.h"
//Add the data structure
#include "miniDst.h"
//#include <stdio.h> 

#ifndef __CINT__  
#include "TROOT.h"
#include "TFile.h" 
#include "TChain.h"
#include "TMath.h"
#include "TH1.h" 
#include "TH2.h"   
#include "TH3.h" 
#include "TF1.h" 
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;
#endif 

vector<Int_t> runNumber; 

Int_t main(Int_t argc, char** argv) 
{
	//Float_t testvalue=0;
	if(argc!=1&&argc!=3) return -1;
	TString inFile="test.list";
	char outFile[1024];
	sprintf(outFile,"test");
	if(argc==3){
		inFile = argv[1];
		sprintf(outFile,"%s",argv[2]);
	}

	//+---------------------------------+
	//| open files and add to the chain |
	//+---------------------------------+
	TChain *chain = new TChain("miniDst");

	Int_t ifile=0;
	char filename[512];
	ifstream *inputStream = new ifstream;
	inputStream->open(inFile.Data());
	if (!(inputStream)) {
		printf("can not open list file\n");
		return 0;
	}
	for(;inputStream->good();){
		inputStream->getline(filename,512);
		if(inputStream->good()) {
			TFile *ftmp = new TFile(filename);
			if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
				cout<<"something wrong"<<endl;
			} else {
				cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
				chain->Add(filename);
				ifile++;
			}
			delete ftmp;
		}
	}
	delete inputStream;


	miniDst *event = new miniDst(chain);
	Int_t nEvents = chain->GetEntries();
	cout<<nEvents<<" events"<<endl;

	runNumber.clear();

	//+-------------+
	//| loop events |
	//+-------------+
	for(Int_t i=0;i<nEvents;i++){ 

		if(i%(nEvents/10)==0) cout << "begin " << i << "th entry...." << endl;

		event->GetEntry(i);
		
		bool eventflag=kFALSE;
		for(Int_t i=0;i<event->mNTrigs;i++){
			if(event->mTrigId[i] == 640001) eventflag = kTRUE; //19.6 minbias
			if(event->mTrigId[i] == 640011) eventflag = kTRUE; //19.6 minbias
			if(event->mTrigId[i] == 640021) eventflag = kTRUE; //19.6 minbias
			if(event->mTrigId[i] == 640031) eventflag = kTRUE; //19.6 minbias
			if(event->mTrigId[i] == 640041) eventflag = kTRUE; //19.6 minbias
			if(event->mTrigId[i] == 640051) eventflag = kTRUE; //19.6 minbias
		}
		if(!eventflag) continue;

		Int_t runId = event->mRunId; 
	  // cout << runId << endl;
    vector<Int_t>::iterator iter = find(runNumber.begin(), runNumber.end(), runId);
		if( iter == runNumber.end() ) runNumber.push_back(runId);

	} 

	cout<<"size of vector runNumber: " << runNumber.size()<<endl;
	sort(runNumber.begin(),runNumber.end());

	ofstream outdata;
	if(argc == 1) outdata.open("test_runNumber.dat");
	else if(argc == 3) outdata.open(Form("%s.dat",argv[2]));

	for(vector<Int_t>::size_type i=0; i<runNumber.size(); i++){
		cout<<runNumber[i]<<endl;
		outdata<<runNumber[i]<<endl;
	}

	cout<<"end of program"<<endl; 
	return 0;
}
