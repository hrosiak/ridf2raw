#ifndef ribf123_rawvar_h
#define ribf123_rawvar_h

//===== Definition ==============================================
double inan=TMath::QuietNaN();

struct rawdata{
    int Nevent; // Event number
    double Dipole[9]; //Brhos
    int Scaler[9][2]; // Scaler ch0:Trigger, ch1:Clock
    int tCOIN[2][16]; // Pileup flag in CCNET01
    
    // PPACs
    int PPAC3_QRaw[4][4];
    int PPAC3_TRaw[4][4];
    int PPAC3_AQRaw[4];
    int PPAC3_ATRaw[4];
    
    int PPAC5_QRaw[4][4];
    int PPAC5_TRaw[4][4];
    int PPAC5_AQRaw[4];
    int PPAC5_ATRaw[4];

    int PPAC7_QRaw[4][4];
    int PPAC7_TRaw[4][4];
    int PPAC7_AQRaw[4];
    int PPAC7_ATRaw[4];

    int PPAC9_QRaw[4][4];
    int PPAC9_TRaw[4][4];
    int PPAC9_AQRaw[4];
    int PPAC9_ATRaw[4];

    int PPAC11_QRaw[4][4];
    int PPAC11_TRaw[4][4];
    int PPAC11_AQRaw[4];
    int PPAC11_ATRaw[4];
    
    //===== Plastics ==========================
    int PL3_QRaw[2]; // Plastic Q raw
    int PL3_TRaw[2]; // Plastic T raw
    int PL3_CFDRaw[2]; // Plastic T(CFD) raw
    int PL3_MTRaw[2][10]; // Plastic Multi-Hit T raw
    int PL3_MHit[2]; // V1290 Multiplicity

    int PL5_QRaw[2]; // Plastic Q raw
    int PL5_TRaw[2]; // Plastic T raw
    int PL5_CFDRaw[2]; // Plastic T(CFD) raw
    int PL5_MTRaw[2][10]; // Plastic Multi-Hit T raw
    int PL5_MHit[2]; // V1290 Multiplicity

    int PL7_QRaw[2]; // Plastic Q raw
    int PL7_TRaw[2]; // Plastic T raw
    int PL7_CFDRaw[2]; // Plastic T(CFD) raw
    int PL7_MTRaw[2][10]; // Plastic Multi-Hit T raw
    int PL7_MHit[2]; // V1290 Multiplicity

    int PL11_QRaw[4]; // Plastic Q raw
    int PL11_TRaw[4]; // Plastic T raw
    int PL11_CFDRaw[4]; // Plastic T(CFD) raw
    int PL11_MTRaw[4][10]; // Plastic Multi-Hit T raw
    int PL11_MHit[4]; // V1290 Multiplicity
    
    int PL11long_QRaw[2]; // VETO Q raw
    int PL11long_MTRaw[2][10]; // VETO Multi-Hit T raw
    int PL11long_MHit[2]; // VETO V1290 Multiplicity

    int tV1290; // V1290 Trigger  
    int tRFRaw; // RF
    
    //===== MUSIC ==========================
    int IC3Raw[6];
    int IC5Raw[5];
    int IC7Raw[6];
    int IC_GasRaw[8];

    int GSIICERaw[2][8]; 
    int GSIICTRaw[2][8]; 
    int GSI1290Raw[32][10];
    int GSIMHit[32];
    
    //===== Time Stamp =====================    
    unsigned long long int TSRaw_C1;
    unsigned long long int TSRaw_C8;
    unsigned long long int TSRaw_F3;
    double TS_C1[3];
    double TS_C8[3];
    double TS_F3[3];    
    
    void clear_variables();
    void set_branch_raw(TTree *tree);
    void set_branch_addresses(TTree *tree2);
};

#ifndef ribf123_rawvar_cpp
void rawdata::clear_variables(){
      Nevent = -1;
      memset(Scaler,-100,sizeof(Scaler));
      memset(tCOIN,0,sizeof(tCOIN));

      memset(PPAC3_QRaw,-100,sizeof(PPAC3_QRaw));
      memset(PPAC3_TRaw,-100,sizeof(PPAC3_TRaw));
      memset(PPAC3_AQRaw,-100,sizeof(PPAC3_AQRaw));
      memset(PPAC3_ATRaw,-100,sizeof(PPAC3_ATRaw));
      
      memset(PPAC5_QRaw,-100,sizeof(PPAC5_QRaw));
      memset(PPAC5_TRaw,-100,sizeof(PPAC5_TRaw));
      memset(PPAC5_AQRaw,-100,sizeof(PPAC5_AQRaw));
      memset(PPAC5_ATRaw,-100,sizeof(PPAC5_ATRaw));
      
      memset(PPAC7_QRaw,-100,sizeof(PPAC7_QRaw));
      memset(PPAC7_TRaw,-100,sizeof(PPAC7_TRaw));
      memset(PPAC7_AQRaw,-100,sizeof(PPAC7_AQRaw));
      memset(PPAC7_ATRaw,-100,sizeof(PPAC7_ATRaw));
      
      memset(PPAC9_QRaw,-100,sizeof(PPAC9_QRaw));
      memset(PPAC9_TRaw,-100,sizeof(PPAC9_TRaw));
      memset(PPAC9_AQRaw,-100,sizeof(PPAC9_AQRaw));
      memset(PPAC9_ATRaw,-100,sizeof(PPAC9_ATRaw));
      
      memset(PPAC11_QRaw,-100,sizeof(PPAC11_QRaw));
      memset(PPAC11_TRaw,-100,sizeof(PPAC11_TRaw));
      memset(PPAC11_AQRaw,-100,sizeof(PPAC11_AQRaw));
      memset(PPAC11_ATRaw,-100,sizeof(PPAC11_ATRaw));
      
      memset(PL3_QRaw,-100,sizeof(PL3_QRaw));
      memset(PL3_TRaw,-100,sizeof(PL3_TRaw));
      memset(PL3_CFDRaw,-100,sizeof(PL3_CFDRaw));
      memset(PL3_MTRaw,-100,sizeof(PL3_MTRaw));
      memset(PL3_MHit,0,sizeof(PL3_MHit));
      
      memset(PL5_QRaw,-100,sizeof(PL5_QRaw));
      memset(PL5_TRaw,-100,sizeof(PL5_TRaw));
      memset(PL5_CFDRaw,-100,sizeof(PL5_CFDRaw));
      memset(PL5_MTRaw,-100,sizeof(PL5_MTRaw));
      memset(PL5_MHit,0,sizeof(PL5_MHit));
      
      memset(PL7_QRaw,-100,sizeof(PL7_QRaw));
      memset(PL7_TRaw,-100,sizeof(PL7_TRaw));
      memset(PL7_CFDRaw,-100,sizeof(PL7_CFDRaw));
      memset(PL7_MTRaw,-100,sizeof(PL7_MTRaw));
      memset(PL7_MHit,0,sizeof(PL7_MHit));
      
      memset(PL11_QRaw,-100,sizeof(PL11_QRaw));
      memset(PL11_TRaw,-100,sizeof(PL11_TRaw));
      memset(PL11_CFDRaw,-100,sizeof(PL11_CFDRaw));
      memset(PL11_MTRaw,-100,sizeof(PL11_MTRaw));
      memset(PL11_MHit,0,sizeof(PL11_MHit));
      
      memset(PL11long_MTRaw,-100,sizeof(PL11long_MTRaw));
      tV1290 = -100;
      
      
      memset(PL11long_MHit,0,sizeof(PL11long_MHit));
      tRFRaw = -100;
      
      memset(IC3Raw,-100,sizeof(IC3Raw));
      memset(IC5Raw,-100,sizeof(IC5Raw));
      memset(IC7Raw,-100,sizeof(IC7Raw));
      memset(IC_GasRaw,-100,sizeof(IC_GasRaw));
      
      memset(GSIICERaw,-100,sizeof(GSIICERaw));
      memset(GSIICTRaw,-100,sizeof(GSIICTRaw));
      memset(GSIMHit,0,sizeof(GSIMHit));
      memset(GSI1290Raw,-100,sizeof(GSI1290Raw));
      
      TSRaw_C1 = 0;
      TSRaw_C8 = 0;
      TSRaw_F3 = 0;

}


void rawdata::set_branch_raw(TTree *tree){
  tree->Branch("Nevent",&Nevent,"Nevent/I");
  tree->Branch("Scaler",Scaler,"Scaler[9][2]/I");
  tree->Branch("tCOIN",tCOIN,"tCOIN[2][16]/I");
  
  tree->Branch("PPAC3_QRaw",PPAC3_QRaw,"PPAC3_QRaw[4][4]/I");
  tree->Branch("PPAC3_TRaw",PPAC3_TRaw,"PPAC3_TRaw[4][4]/I");
  tree->Branch("PPAC3_AQRaw",PPAC3_AQRaw,"PPAC3_AQRaw[4]/I");
  tree->Branch("PPAC3_ATRaw",PPAC3_ATRaw,"PPAC3_ATRaw[4]/I");
  
  tree->Branch("PPAC5_QRaw",PPAC5_QRaw,"PPAC5_QRaw[4][4]/I");
  tree->Branch("PPAC5_TRaw",PPAC5_TRaw,"PPAC5_TRaw[4][4]/I");
  tree->Branch("PPAC5_AQRaw",PPAC5_AQRaw,"PPAC5_AQRaw[4]/I");
  tree->Branch("PPAC5_ATRaw",PPAC5_ATRaw,"PPAC5_ATRaw[4]/I");
  
  tree->Branch("PPAC7_QRaw",PPAC7_QRaw,"PPAC7_QRaw[4][4]/I");
  tree->Branch("PPAC7_TRaw",PPAC7_TRaw,"PPAC7_TRaw[4][4]/I");
  tree->Branch("PPAC7_AQRaw",PPAC7_AQRaw,"PPAC7_AQRaw[4]/I");
  tree->Branch("PPAC7_ATRaw",PPAC7_ATRaw,"PPAC7_ATRaw[4]/I");
  
  tree->Branch("PPAC9_QRaw",PPAC9_QRaw,"PPAC9_QRaw[4][4]/I");
  tree->Branch("PPAC9_TRaw",PPAC9_TRaw,"PPAC9_TRaw[4][4]/I");
  tree->Branch("PPAC9_AQRaw",PPAC9_AQRaw,"PPAC9_AQRaw[4]/I");
  tree->Branch("PPAC9_ATRaw",PPAC9_ATRaw,"PPAC9_ATRaw[4]/I");
  
  tree->Branch("PPAC11_QRaw",PPAC11_QRaw,"PPAC11_QRaw[4][4]/I");
  tree->Branch("PPAC11_TRaw",PPAC11_TRaw,"PPAC11_TRaw[4][4]/I");
  tree->Branch("PPAC11_AQRaw",PPAC11_AQRaw,"PPAC11_AQRaw[4]/I");
  tree->Branch("PPAC11_ATRaw",PPAC11_ATRaw,"PPAC11_ATRaw[4]/I");
  
  tree->Branch("PL3_QRaw",PL3_QRaw,"PL3_QRaw[2]/I");
  tree->Branch("PL3_TRaw",PL3_TRaw,"PL3_TRaw[2]/I");
  tree->Branch("PL3_CFDRaw",PL3_CFDRaw,"PL3_CFDRaw[2]/I");
  tree->Branch("PL3_MTRaw",PL3_MTRaw,"PL3_MTRaw[2][10]/I");
  tree->Branch("PL3_MHit",PL3_MHit,"PL3_MHit[2]/I");
  
  tree->Branch("PL5_QRaw",PL5_QRaw,"PL5_QRaw[2]/I");
  tree->Branch("PL5_TRaw",PL5_TRaw,"PL5_TRaw[2]/I");
  tree->Branch("PL5_CFDRaw",PL5_CFDRaw,"PL5_CFDRaw[2]/I");
  tree->Branch("PL5_MTRaw",PL5_MTRaw,"PL5_MTRaw[2][10]/I");
  tree->Branch("PL5_MHit",PL5_MHit,"PL5_MHit[2]/I");
  
  tree->Branch("PL7_QRaw",PL7_QRaw,"PL7_QRaw[2]/I");
  tree->Branch("PL7_TRaw",PL7_TRaw,"PL7_TRaw[2]/I");
  tree->Branch("PL7_CFDRaw",PL7_CFDRaw,"PL7_CFDRaw[2]/I");
  tree->Branch("PL7_MTRaw",PL7_MTRaw,"PL7_MTRaw[2][10]/I");
  tree->Branch("PL7_MHit",PL7_MHit,"PL7_MHit[2]/I");
  
  tree->Branch("PL11_QRaw",PL11_QRaw,"PL11_QRaw[4]/I");
  tree->Branch("PL11_TRaw",PL11_TRaw,"PL11_TRaw[4]/I");
  tree->Branch("PL11_CFDRaw",PL11_CFDRaw,"PL11_CFDRaw[4]/I");
  tree->Branch("PL11_MTRaw",PL11_MTRaw,"PL11_MTRaw[4][10]/I");
  tree->Branch("PL11_MHit",PL11_MHit,"PL11_MHit[4]/I");
  
  
  tree->Branch("tV1290",&tV1290,"tV1290/I");
  tree->Branch("tRFRaw",&tRFRaw,"tRFRaw/I");
  tree->Branch("PL11long_QRaw",PL11long_QRaw,"PL11long_QRaw[2]/I");
  tree->Branch("PL11long_MTRaw",PL11long_MTRaw,"PL11long_MTRaw[2][10]/I");
  tree->Branch("PL11long_MHit",PL11long_MHit,"PL11long_MHit[2][2]/I");
  
  // IC
  tree->Branch("IC3Raw",IC3Raw,"IC3Raw[6]/I");
  tree->Branch("IC5Raw",IC5Raw,"IC5Raw[5]/I");
  tree->Branch("IC7Raw",IC7Raw,"IC7Raw[6]/I");
  tree->Branch("IC_GasRaw",IC_GasRaw,"IC_GasRaw[8]/I");
  
  tree->Branch("GSIICERaw",GSIICERaw,"GSIICERaw[2][8]/I");
  tree->Branch("GSI1290Raw",GSI1290Raw,"GSI1290Raw[32][10]/I");
  
  tree->Branch("TSRaw_C1",&TSRaw_C1,"TSRaw_C1/l");
  tree->Branch("TSRaw_C8",&TSRaw_C8,"TSRaw_C8/l");
  tree->Branch("TSRaw_F3",&TSRaw_F3,"TSRaw_F3/l");
  tree->Branch("TS_C1",TS_C1,"TS_C1[3]/D");
  tree->Branch("TS_C8",TS_C8,"TS_C8[3]/D");
  tree->Branch("TS_F3",TS_F3,"TS_F3[3]/D");
  tree->Branch("Dipole",Dipole,"Dipole[9]/D");
  }
  
  
void rawdata::set_branch_addresses(TTree *tree2){
  tree2->SetBranchAddress("Nevent",&Nevent);
  tree2->SetBranchAddress("Scaler",Scaler);
  tree2->SetBranchAddress("tCOIN",tCOIN);

  tree2->SetBranchAddress("PPAC3_QRaw",PPAC3_QRaw);
  tree2->SetBranchAddress("PPAC3_TRaw",PPAC3_TRaw);
  tree2->SetBranchAddress("PPAC3_AQRaw",PPAC3_AQRaw);
  tree2->SetBranchAddress("PPAC3_ATRaw",PPAC3_ATRaw);
  
  tree2->SetBranchAddress("PPAC5_QRaw",PPAC5_QRaw);
  tree2->SetBranchAddress("PPAC5_TRaw",PPAC5_TRaw);
  tree2->SetBranchAddress("PPAC5_AQRaw",PPAC5_AQRaw);
  tree2->SetBranchAddress("PPAC5_ATRaw",PPAC5_ATRaw);
  
  tree2->SetBranchAddress("PPAC7_QRaw",PPAC7_QRaw);
  tree2->SetBranchAddress("PPAC7_TRaw",PPAC7_TRaw);
  tree2->SetBranchAddress("PPAC7_AQRaw",PPAC7_AQRaw);
  tree2->SetBranchAddress("PPAC7_ATRaw",PPAC7_ATRaw);
  
  tree2->SetBranchAddress("PPAC9_QRaw",PPAC9_QRaw);
  tree2->SetBranchAddress("PPAC9_TRaw",PPAC9_TRaw);
  tree2->SetBranchAddress("PPAC9_AQRaw",PPAC9_AQRaw);
  tree2->SetBranchAddress("PPAC9_ATRaw",PPAC9_ATRaw);
  
  tree2->SetBranchAddress("PPAC11_QRaw",PPAC11_QRaw);
  tree2->SetBranchAddress("PPAC11_TRaw",PPAC11_TRaw);
  tree2->SetBranchAddress("PPAC11_AQRaw",PPAC11_AQRaw);
  tree2->SetBranchAddress("PPAC11_ATRaw",PPAC11_ATRaw);
  
  tree2->SetBranchAddress("PL3_TRaw",PL3_TRaw);
  tree2->SetBranchAddress("PL3_QRaw",PL3_QRaw);
  tree2->SetBranchAddress("PL3_CFDRaw",PL3_CFDRaw);
  tree2->SetBranchAddress("PL3_MTRaw",PL3_MTRaw);
  tree2->SetBranchAddress("PL3_MHit",PL3_MHit);
  
  tree2->SetBranchAddress("PL5_TRaw",PL5_TRaw);
  tree2->SetBranchAddress("PL5_QRaw",PL5_QRaw);
  tree2->SetBranchAddress("PL5_CFDRaw",PL5_CFDRaw);
  tree2->SetBranchAddress("PL5_MTRaw",PL5_MTRaw);
  tree2->SetBranchAddress("PL5_MHit",PL5_MHit);
  
  tree2->SetBranchAddress("PL7_TRaw",PL7_TRaw);
  tree2->SetBranchAddress("PL7_QRaw",PL7_QRaw);
  tree2->SetBranchAddress("PL7_CFDRaw",PL7_CFDRaw);
  tree2->SetBranchAddress("PL7_MTRaw",PL7_MTRaw);
  tree2->SetBranchAddress("PL7_MHit",PL7_MHit);
  
  tree2->SetBranchAddress("PL11_TRaw",PL11_TRaw);
  tree2->SetBranchAddress("PL11_QRaw",PL11_QRaw);
  tree2->SetBranchAddress("PL11_CFDRaw",PL11_CFDRaw);
  tree2->SetBranchAddress("PL11_MTRaw",PL11_MTRaw);
  tree2->SetBranchAddress("PL11_MHit",PL11_MHit);
  
  tree2->SetBranchAddress("tV1290",&tV1290);
  tree2->SetBranchAddress("tRFRaw",&tRFRaw);
  tree2->SetBranchAddress("PL11long_QRaw",PL11long_QRaw);
  tree2->SetBranchAddress("PL11long_MTRaw",PL11long_MTRaw);
  tree2->SetBranchAddress("PL11long_MHit",PL11long_MHit);
  
  tree2->SetBranchAddress("IC3Raw",IC3Raw);
  tree2->SetBranchAddress("IC5Raw",IC5Raw);
  tree2->SetBranchAddress("IC7Raw",IC7Raw);
  
  tree2->SetBranchAddress("IC_GasRaw",IC_GasRaw);
  tree2->SetBranchAddress("GSIICERaw",GSIICERaw);
  tree2->SetBranchAddress("GSI1290Raw",GSI1290Raw);
  
  tree2->SetBranchAddress("TSRaw_C1",&TSRaw_C1);    
  tree2->SetBranchAddress("TSRaw_C8",&TSRaw_C8);    
  tree2->SetBranchAddress("TSRaw_F3",&TSRaw_F3);    
  tree2->SetBranchAddress("TS_C1",TS_C1);
  tree2->SetBranchAddress("TS_C8",TS_C8);
  tree2->SetBranchAddress("TS_F3",TS_F3);
  tree2->SetBranchAddress("Dipole",Dipole);
                  
}

#endif

#endif
