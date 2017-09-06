#ifndef andrej_h
#define andrej_h
#include "./src/bigrips_zds.h"
#include "./src/ribf123_rawvar.h"
#include <math.h>
#include <string.h>


#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TDOMParser.h"
#include "TXMLNode.h"
#include <iostream>

#include "TArtEventStore.hh"
#include "TArtRawEventObject.hh"


void unpack(TArtRawEventObject *rawevent);
void unpack_dipoles(TArtRawEventObject *rawevent);

rawdata raw;

void event_loop(TArtEventStore *estore, TTree *tree){
	TArtRawEventObject *rawevent = new TArtRawEventObject();
        rawevent = estore->GetRawEventObject();
	int Nevent=0;
	int neve=0;

	while(estore->GetNextEvent()){
		Nevent = rawevent->GetEventNumber();
		if(Nevent<1)continue;
		// clear variables
		raw.clear_variables();
		if(Nevent==1)
		    unpack_dipoles(rawevent);
		unpack(rawevent);

		tree->Fill();
		neve++;
      		if(Nevent%1000==0){
      			std::cout << " decoded " << Nevent << "events" << std::flush << "\r";
  		}
      		estore->ClearData();
		
	} 
}


void unpack_dipoles(TArtRawEventObject *rawevent){    
    TString status_data;
    Int_t runnumber = rawevent->GetRunNumber();

    status_data = *(rawevent->GetStatusData());
    TDOMParser *domParser = new TDOMParser();
    domParser->SetValidate(false);

    memset(raw.Dipole,0,sizeof(raw.Dipole));

    Int_t parsecode = domParser->ParseBuffer(status_data.Data(), status_data.Length());
    if (parsecode < 0) {
	std::cerr << "Parser Error:" << std::endl;
	std::cerr << domParser->GetParseCodeMessage(parsecode) << std::endl;
	return;
      }
	
    TXMLNode * node = domParser->GetXMLDocument()->GetRootNode();
    node = node->GetChildren();
	
    for (; node; node = node->GetNextNode())if (node->GetNodeType() == TXMLNode::kXMLElementNode) break; //is runstatus
    node = node->GetChildren(); 
    for (; node; node = node->GetNextNode())if (node->GetNodeType() == TXMLNode::kXMLElementNode) break; //is dipole
	
	
    TString name;
    double nmr, brho;
	
    for (; node; node = node->GetNextNode()){
      if(node->GetNodeType() == TXMLNode::kXMLElementNode){
	  if(strcmp(node->GetNodeName(), "dipole") == 0){
	    TXMLNode * cnode = node->GetChildren();
	    for (; cnode; cnode = cnode->GetNextNode()){
	      if(cnode->GetNodeType() == TXMLNode::kXMLElementNode){
		  if (strcmp(cnode->GetNodeName(), "name") == 0)
		    name = cnode->GetText();
		  if (strcmp(cnode->GetNodeName(), "nmr") == 0)
		    nmr = (Double_t)atof(cnode->GetText());
		  if (strcmp(cnode->GetNodeName(), "brho") == 0)
		    brho = (Double_t)atof(cnode->GetText());
		  }
	      }
	    if(name.EqualTo("D1"))raw.Dipole[1]=brho;
	    if(name.EqualTo("D2"))raw.Dipole[2]=brho;
	    if(name.EqualTo("D3"))raw.Dipole[3]=brho;
	    if(name.EqualTo("D4"))raw.Dipole[4]=brho;
	    if(name.EqualTo("D5"))raw.Dipole[5]=brho;
	    if(name.EqualTo("D6"))raw.Dipole[6]=brho;
	    if(name.EqualTo("D7"))raw.Dipole[7]=brho;
	    if(name.EqualTo("D8"))raw.Dipole[8]=brho;
	    }// end of dipole node
	  }
      }
    /*
    for(int i=1;i<9;i++){
      std::cout<<"D"<<i<<" = "<<Dipole[i]<<std::endl;
    }
    */ 
    delete domParser;
}


void unpack(TArtRawEventObject *rawevent){
	int numseg = rawevent->GetNumSeg();/* number of segment*/
	for(int i=0; i<numseg; i++){
		TArtRawSegmentObject *seg = rawevent->GetSegment(i);
		int dev=seg->GetDevice(); // ex. BigRIPS
		int fpl=seg->GetFP(); // F5
		int det=seg->GetDetector();  // PPAC-T
		int mod=seg->GetModule(); // C16 16bit Fixed (TDC)
		int num=seg->GetNumData(); // 20ch	

	for(int j=0; j<num; j++){	    
	  TArtRawDataObject *d = seg->GetData(j);
	  int buf = d->GetVal();
	  int ch = d->GetCh(); 

	  if(mod == P7166){  /* Mask for Phillips */
	    buf = buf & 0x0fff;
	  }
	  switch(fpl){
	  case F2:
	    break;
	  case F8:
	    break;
	  case F12:
	    break;
	  default:
	    switch(det){                  
	    case DUPPLAT:
	      break;
	    case ISCALER:
	      break;
	    case ICT:
	      raw.GSI1290Raw[ch][raw.GSIMHit[ch]]=buf;
	      raw.GSIMHit[ch]++;
	      /*
	      if(ch<8){
		GSIICTRaw[0][ch] = buf;
		}  
	     if(ch>7 && ch<16){
                GSIICTRaw[1][ch-8] = buf;
                } 
		*/ 
	    break;
	    case B2SCALER:
	      switch(fpl){
	      case F11VME:
		raw.Scaler[0][ch] = buf;
		break;
	      case F10:
		raw.Scaler[6][ch] = buf;	  
		break;
	      case F7:
		raw.Scaler[7][ch] = buf;
		break;
	      case F11:
		raw.Scaler[8][ch] = buf;
		break;
	      }
	      break;
	    case SCALER:
	      switch(fpl){
	      case F1:
		raw.Scaler[1][ch] = buf;
		break;
	      case F5:
		raw.Scaler[2][ch] = buf;	  
		break;
	      case F7:
		raw.Scaler[3][ch] = buf;	  
		break;
	      case F8:
		raw.Scaler[4][ch] = buf;	  
		break;
	      case F9:
		raw.Scaler[5][ch] = buf;	  
		break;
	      }
	      break;
	    case COIN:
	      for(int k=0; k<16; k++){
		if(fpl==3){
		  raw.tCOIN[0][k] = (k+1)*((buf & (Int_t)pow(2,k)) >> k);
		}else if(fpl==8){
		  raw.tCOIN[1][k] = (k+1)*((buf & (Int_t)pow(2,k)) >> k);
		}
	      }
	      break;
	    case RF:
	      raw.tRFRaw = buf;
	      break;
	    case PPACQ:
	      if(ch<4){
		if(fpl==3)
		  raw.PPAC3_AQRaw[ch] = buf;
		if(fpl==5)
		  raw.PPAC5_AQRaw[ch] = buf;
		if(fpl==7)
		  raw.PPAC7_AQRaw[ch] = buf;
		if(fpl==9)
		  raw.PPAC9_AQRaw[ch] = buf;
		if(fpl==11)
		  raw.PPAC11_AQRaw[ch] = buf;
		}
	      else {
		if(fpl==3)
		  raw.PPAC3_QRaw[ch/4-1][ch%4] = buf;
		if(fpl==5)
		  raw.PPAC5_QRaw[ch/4-1][ch%4] = buf;
		if(fpl==7)
		  raw.PPAC7_QRaw[ch/4-1][ch%4] = buf;
		if(fpl==9)
		  raw.PPAC9_QRaw[ch/4-1][ch%4] = buf;
		if(fpl==11)
		  raw.PPAC11_QRaw[ch/4-1][ch%4] = buf;
	      }
	      break;	    	    	    
	    case PPACT:
	      if(ch<4){
		if(fpl==3)
		  raw.PPAC3_ATRaw[ch] = buf;
		if(fpl==5)
		  raw.PPAC5_ATRaw[ch] = buf;
		if(fpl==7)
		  raw.PPAC7_ATRaw[ch] = buf;
		if(fpl==9)
		  raw.PPAC9_ATRaw[ch] = buf;
		if(fpl==11)
		  raw.PPAC11_ATRaw[ch] = buf;
		}
	      else {
		if(fpl==3)
		  raw.PPAC3_TRaw[ch/4-1][ch%4] = buf;
		if(fpl==5)
		  raw.PPAC5_TRaw[ch/4-1][ch%4] = buf;
		if(fpl==7)
		  raw.PPAC7_TRaw[ch/4-1][ch%4] = buf;
		if(fpl==9)
		  raw.PPAC9_TRaw[ch/4-1][ch%4] = buf;
		if(fpl==11)
		  raw.PPAC11_TRaw[ch/4-1][ch%4] = buf;
	      }
	      break;	    	    	    
	    case PLAQ: //F11 plastics
	      if(fpl==F11){
		break;
	      }else if(fpl==B3F&&ch==0){
		raw.PL11_QRaw[2] = buf;
	      }else if(fpl==B3F&&ch==2){
		raw.PL11_QRaw[3] = buf;
	      }else if(fpl==B3F&&ch==4){
		raw.PL11_QRaw[0]=buf;
	      }else if(fpl==B3F&&ch==6){
		raw.PL11_QRaw[1]=buf;
	      }else if(fpl==B3F&&ch==8){
		raw.PL11long_QRaw[0] = buf;
	      }else if(fpl==B3F&&ch==10){
		raw.PL11long_QRaw[1] = buf;		  
	      }else{
		if(fpl==3)
		  raw.PL3_QRaw[ch] = buf;
		if(fpl==5)
		  raw.PL5_QRaw[ch] = buf;
		if(fpl==7)
		  raw.PL7_QRaw[ch] = buf;
	      }
	      break;
	    case PLAT:
	      if(mod==V1290 && fpl==B3F){
		switch (ch){
		case 31:
		  raw.tV1290 = buf;
		  break;
		case 16:
		  if(raw.PL3_MHit[0]<10)
		    raw.PL3_MTRaw[0][raw.PL3_MHit[0]] = buf;
		  raw.PL3_MHit[0]++;
		  break;
		case 17:
		  if(raw.PL3_MHit[1]<10)
		    raw.PL3_MTRaw[1][raw.PL3_MHit[1]] = buf;
		  raw.PL3_MHit[1]++;		
		  break;
		case 18:
		  if(raw.PL5_MHit[0]<10)
		    raw.PL5_MTRaw[0][raw.PL5_MHit[0]] = buf;
		  raw.PL5_MHit[0]++;		
		  break;
		case 19:
		  if(raw.PL5_MHit[1]<10)
		    raw.PL5_MTRaw[1][raw.PL5_MHit[1]] = buf;
		  raw.PL5_MHit[1]++;		
		  break;
		case 20:
		  if(raw.PL7_MHit[0]<10)
		    raw.PL7_MTRaw[0][raw.PL7_MHit[0]] = buf;
		  raw.PL7_MHit[0]++;		
		  break;
		case 21:
		  if(raw.PL7_MHit[1]<10)
		    raw.PL7_MTRaw[1][raw.PL7_MHit[1]] = buf;
		  raw.PL7_MHit[1]++;		
		  break;
		case 22:// F11PL ULT
		  if(raw.PL11_MHit[0]<10)
		    raw.PL11_MTRaw[0][raw.PL11_MHit[0]] = buf;
		  raw.PL11_MHit[0]++;	
		  break;
		case 23:// F11PL URT
		  if(raw.PL11_MHit[1]<10)
		    raw.PL11_MTRaw[1][raw.PL11_MHit[1]] = buf;
		  raw.PL11_MHit[1]++;		
		  break;
		case 24:// PL11long-LT
		  if(raw.PL11long_MHit[0]<10)
		    raw.PL11long_MTRaw[0][raw.PL11long_MHit[0]] = buf;
		  raw.PL11long_MHit[0]++;		
		  break;
		case 25:// PL11long-RT
		  if(raw.PL11long_MHit[1]<10)
		    raw.PL11long_MTRaw[1][raw.PL11long_MHit[1]] = buf;
		  raw.PL11long_MHit[1]++;		
		  break;
		case 26:// F11PL DLT
		  if(raw.PL11_MHit[2]<10)
		    raw.PL11_MTRaw[2][raw.PL11_MHit[2]] = buf;
		  raw.PL11_MHit[2]++;		
		  break;
		case 27:// F11PL DRT
		  if(raw.PL11_MHit[3]<10)  
		    raw.PL11_MTRaw[3][raw.PL11_MHit[3]] = buf;
		  raw.PL11_MHit[3]++;		
		  break;
		}
	      }else if(fpl==F11){
		switch(ch){
		case 0:
		  raw.PL11_TRaw[1] = buf;
		  break;
		case 1:
		  raw.PL11_TRaw[3] = buf;
		  break;
		case 2:
		  raw.PL11_TRaw[0] = buf;
		  break;
		case 3:
		  raw.PL11_TRaw[2] = buf;
		  break;
		}
	      }else{
		if(fpl==3)
		  raw.PL3_TRaw[ch] = buf;
		if(fpl==5)
		  raw.PL5_TRaw[ch] = buf;
		if(fpl==7)
		  raw.PL7_TRaw[ch] = buf;
	      }
	      break;
	    case STOPPLA:
	      //	    if(ch<2)PL_CFDRaw[fpl][ch] = buf;
	      if(fpl==11){
		switch(ch){
		case 0:
		  raw.PL11_CFDRaw[1] = buf;
		  break;
		case 1:
		  raw.PL11_CFDRaw[3] = buf;
		  break;
		case 2:
		  raw.PL11_CFDRaw[0] = buf;
		  break;
		case 3:
		  raw.PL11_CFDRaw[2] = buf;
		  break;
		}
	      }else{
		if(fpl==3)
		  raw.PL3_CFDRaw[ch] = buf;
		if(fpl==5)
		  raw.PL5_CFDRaw[ch] = buf;
		if(fpl==7)
		  raw.PL7_CFDRaw[ch] = buf;
	      }	    
	      break;
	    case ICE:
	      switch (fpl){
	      case F3://F3VME for F3IC and F5IC
		if(ch<6){
		  raw.IC3Raw[ch] = buf;
		}else if(15<ch && ch<21){
		  raw.IC5Raw[ch-16] = buf;
		}else if(ch==6){
		  raw.IC_GasRaw[3] = buf;
		} else{
		  raw.IC_GasRaw[5] = buf;
		}
		break;
	      case F7:
		raw.IC7Raw[ch] = buf;
		break;
	      case F11VME://GSI IC
		if(ch<8){
		  raw.GSIICERaw[0][ch] = buf;
		}else if(7<ch&&ch<16){
		  raw.GSIICERaw[1][ch-8] = buf;
		}	
		break;
	      }
	      break;
	    case ICGAS:
	      raw.IC_GasRaw[fpl] = buf;
	      break; 
	    case TS:
	      if(ch==0){
		switch (fpl){
		case B3F:
		  raw.TSRaw_C1 = buf << 16;
		  raw.TSRaw_C1 = raw.TSRaw_C1 << 16;
		  break;
		case F11:
		  raw.TSRaw_C8 = buf << 16;
		  raw.TSRaw_C8 = raw.TSRaw_C8 << 16;
		  break;
		case F3:
		  raw.TSRaw_F3 = buf << 16;
		  raw.TSRaw_C1 = raw.TSRaw_F3 << 16;
		  break;
		}
	      }else{	      
		switch (fpl){
		case B3F:
		  raw.TSRaw_C1 = raw.TSRaw_C1  + buf;
		  raw.TS_C1[0] = (double)raw.TSRaw_C1 * 10.; // ns
		  raw.TS_C1[1] = raw.TS_C1[0] *0.000001; // ms
		  raw.TS_C1[2] = raw.TS_C1[1] *0.001; // s
		  break;
		case F11:
		  raw.TSRaw_C8 = raw.TSRaw_C8  + buf;
		  raw.TS_C8[0] = (double)raw.TSRaw_C1 * 10.; // ns
		  raw.TS_C8[1] = raw.TS_C8[0] *0.000001; // ms
		  raw.TS_C8[2] = raw.TS_C8[1] *0.001; // s
		  break;
		case F3:
		  raw.TSRaw_F3 = raw.TSRaw_F3 + buf;
		  raw.TS_F3[0] = (double)raw.TSRaw_C1 * 10.; // ns
		  raw.TS_F3[1] = raw.TS_F3[0] *0.000001; // ms
		  raw.TS_F3[2] = raw.TS_F3[1] *0.001; // s
		  break;
		}
	      }
	      break;

	    }	  // Detector number switch
	    break;
	  }
	  //		    if(mod!=21&&fpl!=63){
	} // number of data loop	  	


	}//segment loop
	
}

#endif
