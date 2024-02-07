#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

void sortRunNumber() 
{

	system("mkdir -p runList");

    vector<Int_t> runNumber; 
	vector<Int_t> dayNumber;
	runNumber.clear();
	dayNumber.clear();

	ifstream indata("runList.list");

	Int_t runId;
	while(indata>>runId){ 
		vector<Int_t>::iterator iter = find(runNumber.begin(), runNumber.end(), runId);
		if( iter == runNumber.end() ) runNumber.push_back(runId);
	} 
	cout<<"Total size of runNumber list: " << runNumber.size()<<endl;

	Int_t dayId;
	ofstream outdata("runList/mTotalRunList.dat");
	for(UInt_t i=0; i<runNumber.size(); i++){
		dayId = (runNumber[i]/1000)%1000;
		cout<<"RunId:"<<runNumber[i]<<" \t"<<"dayId:"<<dayId<<endl;
		outdata<<runNumber[i]<<endl;

		vector<Int_t>::iterator iter = find(dayNumber.begin(), dayNumber.end(), dayId);
		if( iter == dayNumber.end() ) dayNumber.push_back(dayId);
	}
	outdata.close();
	cout<<"Total size of dayNumber list: " << dayNumber.size()<<endl;

	outdata.open("runList/mTotalDayList.dat");
	for(UInt_t i=0; i<dayNumber.size(); i++){
		cout<<"dayId:"<<dayNumber[i]<<endl;
		outdata<<dayNumber[i]<<endl;
	}

	indata.close();
	outdata.close();
	cout<<"end of program"<<endl; 
}
