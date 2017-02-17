/*

*/
#include <signal.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TSystem.h"
#include "TString.h"
#include "TApplication.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TCutG.h"

#include "and.h"

using namespace std;

int main(int argc, char* argv[]){
    //TApplication theApp("MKFILE", &argc, argv);
    const char* outputdir = "./";
    const char* rawfileprefix = "raw_";
    int run;
    TString ifname, ofname, ofprefix;

    if(argc<2){
        printf("Wrong arguments\n");
        printf("Arguments: ./makeraw inputfile.ridf output.root\n");
        printf("output.root filename is optional, by default raw_ preffix is added and .ridf suffix is changed to .root\n");
	return 0;
   }

   if(argc>1){   
	ofprefix = outputdir;
	ofprefix = ofprefix + rawfileprefix;
	ifname = argv[1];
        ofname = ifname;
        ofname.ReplaceAll(".ridf",".root");
        ofname.Replace(0,ofname.Last('/')+1,ofprefix.Data());
   }
   if(argc>2){
       ofname = argv[2];
   }
   cout<<"output file = \033[1;32m"<<ofname<<"\033[0m"<<endl;
    /*  Open File & Set EventStore, RawEventObject  */
    TArtEventStore *estore = new TArtEventStore();
    TArtRawEventObject *rawevent = new TArtRawEventObject();
    rawevent = estore->GetRawEventObject();/* return rawevent */
    estore->Open(ifname); /*confirm whether the ridf file is opened or not. */

    /*  Define File, Tree, Histo, Canvas  */
    TFile *fout = new TFile(ofname, "RECREATE");
    TTree *tree = new TTree("tree","tree");
    
    /*  Set Branch  */
    raw.set_branch_raw(tree);
    
    event_loop(estore,tree);    
    fout->Write();
    fout->Close();
    cout<<"\033[1;32mdone\033[0m"<<endl;
    return 0;
}

