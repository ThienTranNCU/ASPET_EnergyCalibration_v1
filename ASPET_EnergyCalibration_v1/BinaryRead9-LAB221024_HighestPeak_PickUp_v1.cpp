//221216
//Pick up the highest energy peak ang give ToT value, resolution of the peak, and area of the peak for each channel
//output file:   gmslID / sticID / channel / peak_area / ToT_ADC / resolution_sigma / real_Energy_keV
//Out put file in .txt file of the same directory of input .bin file

#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TRint.h"
#include "TObject.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TSpectrum.h"
#include "TPaveStats.h"

#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TBrowser.h"

#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <complex>
#include <TApplication.h>
#include <TMath.h>

///TIME
#include <chrono>

//#include <vector>


using namespace std;


int compare(const void *arg1, const void *arg2)
{
    //    if(arg1>arg2){return 1;}
    //    if(arg1<arg2){return -1;}
    //    return (arg1>arg2);
    const int *lhs = static_cast<const  int*>(arg1);
    const int *rhs = static_cast<const  int*>(arg2);
    if ( lhs[0] < rhs[0] ) return -1;
    if ( lhs[0] > rhs[0] ) return +1;
    if ( lhs[0] == rhs[0] ) return 0;
}

int compare2(const void *arg1, const void *arg2)   //Frame | Time
{
    //    if(arg1>arg2){return 1;}
    //    if(arg1<arg2){return -1;}
    //    return (arg1>arg2);
    const int *lhs = static_cast<const  int*>(arg1);
    const int *rhs = static_cast<const  int*>(arg2);
    if ( lhs[3] > rhs[3] ) return +1;
    if ( lhs[3] < rhs[3] ) return -1;
    if(lhs[3]==rhs[3])
    {
        if ( lhs[0] < rhs[0] ) return -1;
        if ( lhs[0] > rhs[0] ) return +1;
        if ( lhs[0] == rhs[0] ) return 0;
    }
}

int compare3(const void *arg1, const void *arg2)   //Frame | STiC ID
{
    //    if(arg1>arg2){return 1;}
    //    if(arg1<arg2){return -1;}
    //    return (arg1>arg2);
    const int *lhs = static_cast<const  int*>(arg1);
    const int *rhs = static_cast<const  int*>(arg2);
    if ( lhs[3] > rhs[3] ) return +1;
    if ( lhs[3] < rhs[3] ) return -1;
    if(lhs[3]==rhs[3])
    {
        if ( lhs[1] < rhs[1] ) return -1;
        if ( lhs[1] > rhs[1] ) return +1;
        if ( lhs[1] == rhs[1] )
        {
            if ( lhs[6] < rhs[6] ) return -1;
            if ( lhs[6] > rhs[6] ) return +1;
        }
    }
}

//0tcc,1sid,2board,3frame,4ch,5energy, 6eventRecNo, 7time, 8gmslid
int compare4(const void *arg1, const void *arg2)   //Frame | GMSL ID |EventID
{
    //    if(arg1>arg2){return 1;}
    //    if(arg1<arg2){return -1;}
    //    return (arg1>arg2);
    const int *lhs = static_cast<const  int*>(arg1);
    const int *rhs = static_cast<const  int*>(arg2);
    if ( lhs[3] > rhs[3] ) return +1;
    if ( lhs[3] < rhs[3] ) return -1;
    if(lhs[3]==rhs[3])
    {
        if ( lhs[8] < rhs[8] ) return -1;
        if ( lhs[8] > rhs[8] ) return +1;
        if ( lhs[8] == rhs[8] )
        {
            if ( lhs[6] < rhs[6] ) return -1;
            if ( lhs[6] > rhs[6] ) return +1;
        }
    }
}

int compareX(const void *arg1, const void *arg2)
{
    //    if(arg1>arg2){return 1;}
    //    if(arg1<arg2){return -1;}
    //    return (arg1>arg2);
    const int *lhs = static_cast<const  int*>(arg1);
    const int *rhs = static_cast<const  int*>(arg2);
    if ( lhs[0] < rhs[0] ) return -1;
    if ( lhs[0] > rhs[0] ) return +1;
    if ( lhs[0] == rhs[0] ) return 0;
}
int compareY(const void *arg1, const void *arg2)
{
    //    if(arg1>arg2){return 1;}
    //    if(arg1<arg2){return -1;}
    //    return (arg1>arg2);
    const int *lhs = static_cast<const  int*>(arg1);
    const int *rhs = static_cast<const  int*>(arg2);
    if ( lhs[1] < rhs[1] ) return -1;
    if ( lhs[1] > rhs[1] ) return +1;
    if ( lhs[1] == rhs[1] ) return 0;
}


int mapx1(int channel)
{
    float x=0;
    if((channel==58)||(channel==57)||(channel==54)||(channel==52)||(channel==46)||(channel==44)||(channel==38)||(channel==38)){x=0;}
    if((channel==60)||(channel==59)||(channel==56)||(channel==50)||(channel==48)||(channel==42)||(channel==40)||(channel==36)){x=1;}
    if((channel==61)||(channel==55)||(channel==53)||(channel==47)||(channel==45)||(channel==39)||(channel==32)||(channel==35)){x=2;}
    if((channel==63)||(channel==62)||(channel==51)||(channel==49)||(channel==43)||(channel==41)||(channel==34)||(channel==33)){x=3;}
    if((channel==0)||(channel==1)||(channel==12)||(channel==14)||(channel==20)||(channel==22)||(channel==29)||(channel==30)){x=4;}
    if((channel==2)||(channel==8)||(channel==10)||(channel==16)||(channel==18)||(channel==24)||(channel==31)||(channel==28)){x=5;}
    if((channel==3)||(channel==4)||(channel==7)||(channel==13)||(channel==15)||(channel==21)||(channel==23)||(channel==27)){x=6;}
    if((channel==5)||(channel==6)||(channel==9)||(channel==11)||(channel==17)||(channel==19)||(channel==25)||(channel==26)){x=7;}
    //else{x=-1;}
    return x;
}
int mapy1(int channel)
{
    float y=0;
    if((channel==58)||(channel==60)||(channel==61)||(channel==63)||(channel==0)||(channel==2)||(channel==3)||(channel==5)){y=0;}
    if((channel==57)||(channel==59)||(channel==55)||(channel==62)||(channel==1)||(channel==8)||(channel==4)||(channel==6)){y=1;}
    if((channel==54)||(channel==56)||(channel==53)||(channel==51)||(channel==12)||(channel==10)||(channel==7)||(channel==9)){y=2;}
    if((channel==52)||(channel==50)||(channel==47)||(channel==49)||(channel==14)||(channel==16)||(channel==13)||(channel==11)){y=3;}
    if((channel==46)||(channel==48)||(channel==45)||(channel==43)||(channel==20)||(channel==18)||(channel==15)||(channel==17)){y=4;}
    if((channel==44)||(channel==42)||(channel==39)||(channel==41)||(channel==22)||(channel==24)||(channel==21)||(channel==19)){y=5;}
    if((channel==38)||(channel==40)||(channel==32)||(channel==34)||(channel==29)||(channel==31)||(channel==23)||(channel==25)){y=6;}
    if((channel==37)||(channel==36)||(channel==35)||(channel==33)||(channel==30)||(channel==28)||(channel==27)||(channel==37)){y=7;}
    //else{y=-1;}
    return y;
}

int mapxAll(int channel, int sid, int gid)
{
    float x=0;
    if((channel==58)||(channel==57)||(channel==54)||(channel==52)||(channel==46)||(channel==44)||(channel==38)||(channel==38)){x=0;}
    if((channel==60)||(channel==59)||(channel==56)||(channel==50)||(channel==48)||(channel==42)||(channel==40)||(channel==36)){x=1;}
    if((channel==61)||(channel==55)||(channel==53)||(channel==47)||(channel==45)||(channel==39)||(channel==32)||(channel==35)){x=2;}
    if((channel==63)||(channel==62)||(channel==51)||(channel==49)||(channel==43)||(channel==41)||(channel==34)||(channel==33)){x=3;}
    if((channel==0)||(channel==1)||(channel==12)||(channel==14)||(channel==20)||(channel==22)||(channel==29)||(channel==30)){x=4;}
    if((channel==2)||(channel==8)||(channel==10)||(channel==16)||(channel==18)||(channel==24)||(channel==31)||(channel==28)){x=5;}
    if((channel==3)||(channel==4)||(channel==7)||(channel==13)||(channel==15)||(channel==21)||(channel==23)||(channel==27)){x=6;}
    if((channel==5)||(channel==6)||(channel==9)||(channel==11)||(channel==17)||(channel==19)||(channel==25)||(channel==26)){x=7;}
    //else{x=-1;}
    if(sid/2==0){x=x;}
    if(sid/2==1){x= 7-x;}
    x= (gid%2)*16+(sid%2)*8 + x;
    return x;
}
int mapyAll(int channel, int sid, int gid)
{
    float y=0;
    if((channel==58)||(channel==60)||(channel==61)||(channel==63)||(channel==0)||(channel==2)||(channel==3)||(channel==5)){y=0;}
    if((channel==57)||(channel==59)||(channel==55)||(channel==62)||(channel==1)||(channel==8)||(channel==4)||(channel==6)){y=1;}
    if((channel==54)||(channel==56)||(channel==53)||(channel==51)||(channel==12)||(channel==10)||(channel==7)||(channel==9)){y=2;}
    if((channel==52)||(channel==50)||(channel==47)||(channel==49)||(channel==14)||(channel==16)||(channel==13)||(channel==11)){y=3;}
    if((channel==46)||(channel==48)||(channel==45)||(channel==43)||(channel==20)||(channel==18)||(channel==15)||(channel==17)){y=4;}
    if((channel==44)||(channel==42)||(channel==39)||(channel==41)||(channel==22)||(channel==24)||(channel==21)||(channel==19)){y=5;}
    if((channel==38)||(channel==40)||(channel==32)||(channel==34)||(channel==29)||(channel==31)||(channel==23)||(channel==25)){y=6;}
    if((channel==37)||(channel==36)||(channel==35)||(channel==33)||(channel==30)||(channel==28)||(channel==27)||(channel==37)){y=7;}
    //else{y=-1;}
    if(((sid/2)^(gid/2))==0){y=y;}
    if(((sid/2)^(gid/2))==1){y= 7-y;}
    y =  ((sid/2)^(gid/2))*8+y;
    return y;
}


std::string toBinary(int n)
{
    std::string r; int count=15;
    //while(n!=0) {
    while(count>0) {
        r=(n%2==0 ?"0":"1")+r; n/=2;
        count--;
    }
    return r;
}

std::string toBinaryBig(int n)
{
    std::string r; int count=30;
    //while(n!=0) {
    while(count>0) {
        r=(n%2==0 ?"0":"1")+r; n/=2;
        count--;
    }
    return r;
}

std::string toBinarySmall(int n, int count)
{
    std::string r;
    //while(n!=0) {
    while(count>0) {
        r=(n%2==0 ?"0":"1")+r; n/=2;
        count--;
    }
    return r;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////MAIN PROGRAM
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    cout<<"Start1"<<endl;
    cout<<"argc="<<argc<<endl;
    //cout<<"argv[0]="<<argv[0]<<"\t argv[1]="<<argv[1]<<"\t argv[2]"<<argv[2]<<"\n"<<endl;
    
    ////Variables definition
    char argName[1000], fileName[1000], file[10][1000], fileNameRoot[1000];//file name string
    const int N=20000;
    double xmin=0, xmax=32800,realE=511.0,peakLT=100,peakHT=100000;
    int tDAQ=1,ChA=0,ChB=0;
    double sigP=0.4,sigB=7;                //// peak find sigma and BG find sigma sigP=0.4 for 511keV, sigP=0.6 for 306KeV, sigB=7 for 511keV and 306 keV, sigB = 3 for >1000keV
    int ELChA[4][4][64], EHChA[4][4][64], ELChB[4][4][64], EHChB[4][4][64]; //A is primary, B is secondary
    //int ELChA=200, EHChA=700, ELChB=200, EHChB=700;
    int pos=0; //POSITION TAG
    int  Ecut=2; //number of sigmas + and - about the mean photopeak value (e.g Ecut=2=>+/- 2sigma)
    int px1,py1;
    if(argc<2){
        string usage1= "Usage: FilePath Elo Ehi t(s) ChA ChB ELChA EHChA ELChB EHChB\n";
        string usage2= "./BinaryRead9 test/PET/test/2021-... 200 700 10 63 63 511.0";
        string usage3= "./BinaryRead9 ~/../../Volumes/ASPETcalib/2022-.. 200 700 10 63 63 ";
        cout<<usage1<<"\n"<<usage2<<"\n"<<usage3<<endl; return -1;}
    if(argc>=3){xmin=strtod(argv[2],NULL); xmax=strtod(argv[3],NULL);} //ENERGY LIMITS for ToT ADC
    if(argc>=4){sigP=strtod(argv[4],NULL);} //Peak finding parameter
    if(argc>=5){sigB=strtod(argv[5],NULL);} //Back ground fitting parameter
    if(argc>=6){realE=strtod(argv[6],NULL);}
	if(argc>=8){peakLT=strtod(argv[7],NULL); peakHT=strtod(argv[8],NULL);}  ////Lower and Higher threshold for peak count
    //if(argc>=8){ELChA=strtod(argv[7],NULL),EHChA=strtod(argv[8],NULL) ,ELChB=strtod(argv[9],NULL) ,EHChB=strtod(argv[10],NULL);}
//    if(argc>=12){pos=strtod(argv[11],NULL);}
    
    //cout<<strtod(argv[4],NULL)<<"\t tDAQ="<<tDAQ<<endl;
    sprintf(argName,"%s",argv[1]);
    sprintf(fileName,"%s.bin",argName);
    cout<<"Opening file from the argument..\t"<<fileName<<endl;
    //// ********************************6 arguments 1: file name; 2,3: min and max energy accepted in ToT ADC; 4: Peak finding parameter; 5: Background finding parameter; 6: real energy of the highest energy peak; 
    
    //////////////////////////////
    
    const int nevents=1;//100000
    unsigned char buffer[8*nevents]; //char is 1 byte long, so char[8] is 8 bytes long - buffer size for one event
    unsigned int a[8];
    
    unsigned int ecoarse;
    unsigned int tcoarse;
    unsigned int energy;
    unsigned int efine;
    unsigned int tfine;
    unsigned int time;
    unsigned int ebad;
    unsigned int tbad;
    
    unsigned int channel;
    unsigned int frame;
    unsigned int gmslID;
    unsigned int sticID;
    unsigned int boardID;
    unsigned int Data0_Sync1;
    unsigned int RSTcountcc;
    
    
    //HISTOGRAMS
    TH1F *hEc=new TH1F("hEc","E coarse",32770,-1,32769);
    //TH1F *hTc=new TH1F("hTc","time coarse",32770,-1,32769);//32768
    TH1F *hTc=new TH1F("hTc","time coarse",40000,-1,39999);//32768
    TH1F *hTf=new TH1F("hTf","time fine",32,0,32);
    TH1F *hEf=new TH1F("hEf","energy fine",32,0,32);
    TH1F *hT=new TH1F("hT","Time",10000000,-5000000,5000000);
    TH1F *hEbad=new TH1F("hEbad","E-bad",2,0,2);
    TH1F *hTbad=new TH1F("hTbad","T-bad",2,0,2);
    
    TH1F *hFrame=new TH1F("hFrame","Frames",256,0,256);
    TH1F *hCh=new TH1F("hCh","Channels",64,0,64);
    TH1F *hSticID=new TH1F("hSticID","STiC ID",8,0,8);
    TH1F *hGMSLID=new TH1F("hGMSLID","FPGA ID",8,0,8);
    TH1F *hBoardID=new TH1F("hBoardID","Board ID",8,0,8);
    TH2F *h2EventEnergy = new TH2F("h2EventEnergy","h2EventEnergy",1000,0,1000000,33000,0,33000);
    TH2F *h2SticIDGMSLID=new TH2F("hSticIDGMSLID2","SticID-GMSLID",4,0,4,4,0,4); //210822// more info on distribution
    
    TH1F *hRST=new TH1F("hRST","Reset counter",10000000,0,10000000);
    TH1F *hRSTG0=new TH1F("hRSTG0","Reset counter G0",10000,0,100000);
    TH1F *hRSTG1=new TH1F("hRSTG1","Reset counter G1",1000000,0,100000);
    TH1F *hRSTG2=new TH1F("hRSTG2","Reset counter G2",1000000,0,100000);
    TH1F *hRSTG3=new TH1F("hRSTG3","Reset counter G3",1000000,0,100000);
    TH1F *hRSTdeltaG0=new TH1F("hRSTdeltaG0","Reset counter Delta",10000,0,10000);
    TH1F *hRSTdeltaG1=new TH1F("hRSTdeltaG1","Reset counter Delta",10000,0,10000);
    TH1F *hRSTdeltaG2=new TH1F("hRSTdeltaG2","Reset counter Delta",10000,0,10000);
    TH1F *hRSTdeltaG3=new TH1F("hRSTdeltaG3","Reset counter Delta",10000,0,10000);
    
    TH1F *hE[4][4][64];
    TH1F *Bas[4][4][64];
    TH1F *Sig[4][4][64];
    TH1F *h1PPmean[4][4], *h1PPamp[4][4];
    char hname[50];
    for(int g=0; g<4; g++)
    {
        for(int s=0; s<4; s++)
        {
            sprintf(hname,"h1PPmeanG%dS%d",g,s);
            h1PPmean[g][s]=new TH1F(hname,"E",1000,0,1000);
            sprintf(hname,"h1PPampG%dS%d",g,s);
            h1PPamp[g][s]=new TH1F(hname,"E",1000,0,1000);
            for(int i=0; i<64; i++)
            {
				px1=mapxAll(i,s,g);
				py1=mapyAll(i,s,g);
                sprintf(hname,"G%dS%dCh%d",g,s,i);
                hE[g][s][i]=new TH1F(hname,"E",8192,0,32768); // until 210915-binning one in 4
                sprintf(hname,"G%dS%dCh%d",g,s,i);
                Bas[g][s][i]=new TH1F(hname,"B",8192,0,32768);
                sprintf(hname,"G%dS%dCh%d",g,s,i);
                Sig[g][s][i]=new TH1F(hname,"S",8192,0,32768);
            }
        }
    }
    
    // Defining TFile and tree
    sprintf(fileNameRoot,"%s_outputTree.root",argName);
    //TFile* fileRoot_tree = new TFile("data.root", "RECREATE");
    TFile* fileRoot_tree = new TFile(fileNameRoot, "RECREATE");
    TTree *tree= new TTree("tree", "tree");
    //Setting up the branches..
    tree->Branch("ecoarse", &ecoarse, "ecoarse/I");
    tree->Branch("tcoarse", &tcoarse, "tcoarse/I");
    tree->Branch("efine", &efine, "efine/I");
    tree->Branch("tfine", &tfine, "tfine/I");
    tree->Branch("channel", &channel, "channel/I");
    tree->Branch("frame", &frame, "frame/I");
    tree->Branch("tbad", &tbad, "tbad/I");
    tree->Branch("ebad", &ebad, "ebad/I");
    tree->Branch("sticID", &sticID, "sticID/I");
    tree->Branch("gmslID", &gmslID, "gmslID/I");
    tree->Branch("boardID", &boardID, "boardID/I");
    tree->Branch("Data0_Sync1", &Data0_Sync1, "Data0_Sync1/I");
    tree->Branch("energy", &energy, "energy/I");
    tree->Branch("time", &time, "time/I");
    tree->Branch("RSTcountcc", &RSTcountcc, "RSTcountcc/I");
    
    //char path[200],tag[300], date[50];
    FILE *ptr;
    
    //OPENING FILE; run a loop, get events, close file; OPEN AGAIN
    cout<<"Opening file.."<<fileName<<endl;
    ptr=fopen(fileName,"rb");
    if(!ptr){cout<<"No such file or directory!"<<endl; return 0;}
    else{cout<<"opening succesful.."<<fileName<<endl;}
    int eventsTOTAL=0;
    while(std::fread(buffer, sizeof(buffer),1,ptr)!='\0')
    {
        eventsTOTAL++;
    }
    const int eventsTOTALmain=eventsTOTAL;
    const int loops = eventsTOTALmain/15320; //ARBITRARY ASSIGNMENT OF ARRAY SIZE LEADING TO MANY LOOPS
    cout<<"TOTAL EVENTS="<<eventsTOTALmain<<endl;
    //if(eventsTOTAL>200000){eventsTOTAL=200000;}
    if(eventsTOTAL>20480){eventsTOTAL=20480;} // Assigned to coin loops later
    fclose(ptr);
    cout<<"ANALYZED EVENTS="<<eventsTOTAL<<endl;
    ptr=fopen(fileName,"rb");                                                               // Open and read binary file
    
    //with an event rate of 1MHz, you will get 1 event in 1us, or 100 events in 1000 events in 1ms on average. So the chance of a true coincidence lying in the next ms is extremely low.
    //with 20K, we further increase this to a 20ms avg window.
    
    //////////////////////////////////////////////////////////////////////////
    //OUTFILE SINGLES
 //   ofstream outfileSingles;
    char fileNameSingles[1000];
    sprintf(fileNameSingles,"%s_Singles.txt",argName);
    //sprintf(fileNameD,"DBOX/DELAY-VALUES-WRITE/%s_Summary.txt",argName);
 //   outfileSingles.open(fileNameSingles, ios::out);
 //   if(outfileSingles.fail()){
  //      printf("Failed opening the Singles write file %s\n",fileNameSingles);
 //       return 0;
  //  }
  //  printf("Opening file.. \t %s\n",fileNameSingles);
    
    
    //////////////////////////////////////////////////////////////////////////
    //OUTFILE LOR
//    ofstream outfileLOR;
    char fileNameLOR[1000];
    sprintf(fileNameLOR,"%s_LOR.txt",argName);
    //sprintf(fileNameD,"DBOX/DELAY-VALUES-WRITE/%s_Summary.txt",argName);
 //   outfileLOR.open(fileNameLOR, ios::out);
 //   if(outfileLOR.fail()){
//        printf("Failed opening the LOR write file %s\n",fileNameLOR);
 //       return 0;
 //   }
  //  printf("Opening file.. \t %s\n",fileNameLOR);
    
    
    //////////////////////////////////////////////////////////////////////////
    //LFSR DECODING****
    //////////////////////////////////////////////////////////////////////////
    int16_t m_lut[ 1 << 15 ], encodedLfsr[1 << 15];
    m_lut[ 0x7FFF ] = -1; // invalid state
    uint16_t lfsr = 0x0000;
    for ( int16_t n = 0; n < ( 1 << 15 ) - 1; ++n )
    {
        m_lut[ lfsr ] = n;
        encodedLfsr[n]=lfsr;
        const uint8_t bits13_14 = lfsr >> 13;
        uint8_t new_bit;
        switch ( bits13_14 )
        { // new_bit = !(bit13 ^ bit14)
            case 0x00:
            case 0x03:
                new_bit = 0x01;
                break;
            case 0x01:
            case 0x02:
                new_bit = 0x00;
                break;
        }// switch
        lfsr = ( lfsr << 1 ) | new_bit; // add the new bit to the right
        lfsr &= 0x7FFF; // throw away the 16th bit from the shift
    }// for
    /////////////////////////////////////
    /////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////
    ///// EVENT MAIN LOOP: SHORT COIN BUFFERS, RST COUNTERS ETC
    ////////////////////////////////////////////////////////////////////////////////
    //0tcc,1sid,2board,3frame,4ch,5energy, 6eventRecNo, 7time, 8gmslid
    //READING FILE THROUGH THE BUFFER
    int ncount=0, ncountc=0, ncountall=0; // ncounts: stores the number of events, ncountc is a local counter for cm[][]
    int nloopc=0; //number of coincidence loops
    int ncountGood=0; //Stores the number of useful events
    //eventsTOTAL has the total events being analyzed for this task
    int ncoin=0;//Stores the total opposite coincidences registered
    
    int ncountBrdGmslStic[1][4][4];
    int ncountGoodBrdGmslStic[1][4][4];
    
    int RSTcountAbscc=0; // Absolute value of the RST counter as readout qqq
    int RSTcountf=0; // fine value of RST count given by the rolling Efine bins qqq
    int RSTcountccStored[4]; // Counter value as maintained by the number of RST events read by the code
    //long int RSTcountccStoredG=0; //? qqq
    int RSTcountT0[4]; //Start value of the RSTcountAbs or the first stored RST event value
    int oldRSTcountcc0=0,oldRSTcountcc1=0,oldRSTcountcc2=0, oldRSTcountcc3=0;// Previous value
    int RST[4][10000];//1M GOOD FOR 2000S ~33MINS
    int rcount[4][10000];
    int evcount[4][10000];
    float RSTLOR = 0;
    int RSTEventcc[4]; //Holds the latest RSTcc value for the corresponding gID
    int RSTEventff[4]; //RSTff for the current event
    bool RSTflag;
    int RSTnew[4], RSTold[4];//old and current until RSTcountT0 is found for all gmsl.
    //cout<<"SNo\t"<<"frame\t"<<"bID\t"<<"sID\t"<<"ch\t"<<"t\t"<<"e\t"<<"tcc\t"<<"ecc\t"<<"dtcc\t"<<endl;
    
    
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////COINCIDENCE CONDITIONS, DECLARING VARIABLES, BUFFERS, HISTOGRAMS, ETC
    ////////////////////////////////////////////////////////////////////////////////
    //const int cmsize=eventsTOTAL/8; //cmsize=eventsTOTAL; cmsize=eventsTOTAL/2
    const int cmsize=eventsTOTAL;
    const int cmwidth=11;
    int cm[cmsize][cmwidth]; cout<<"cm[cmsize][cmwidth] array opened"<<endl;
    
    //PLL DATA
    //const int windowTcc=3, windowEvents=2048; //windowEvents=512 // 32736
    //const int windowcluster=2; //No of tcc within which events are considered as clusters
    //NA22 DATA
    const int windowTcc=10, windowEvents=1024; //210526 ans 2106xx dats |windowTcc=10, windowEvents=1024
    //const int windowTcc=3, windowEvents=2048; //windowEvents=512 // 2109xx..
    const int windowcluster=2; //No of tcc within which events are considered as clusters
    
    
    
    //Channel Map | Printing Channel Map
    float x1=0,x2=0,xmean=0,y1=0,y2=0,ymean=0;
    int lch=0, rch=0; //continuous channel map with GID, SID weighted  - 512ch each
    cout<<"Channel"<<"\tmapX"<<"\tmapY"<<endl;
    for(int i=0;i<63;i++)
    {
        cout<<i<<"\t"<<mapx1(i)<<"\t"<<mapy1(i)<<endl;
    }
    cout<<"first few events"<<cm[0][0]<<"\t"<<cm[1][0]<<"\t"<<cm[2][0]<<endl;

    int CFtrans=0, CStrans=0;
    int Cframe=cm[0][3], Cframeold=cm[0][3];
    int CgmslID=cm[0][8], CgmslIDold=cm[0][8];
    int Cevent=0;
    
    int changeS=cm[0][8], changeF=cm[0][3]; //not needed
    
    int CeventStart1=0, CeventStart2=0, CeventEnd=0;
    int oldtcc1=0,oldtcc2=0; int cluster=0;
    ////////////////////////MAPS TIMEFRAME MAPS TIMEFRAME MAPS TIMEFRAME//////////////////
    //DISPLAY OF MAPS BROKEN IN TIMEFRAMES qqq
    float frameTime=1200, totalTime=1200; //1200S TOTAL
    int nframes=totalTime/frameTime;
    
    /////////////////Histogram definition annd fill on the fly
    TH2F *h3Coin =  new TH2F("h3Coin","h3Coin",512,0,512,512,0,512);//This is the raw coincidence map from the detector
    TH2F *h2Coin =  new TH2F("h2Coin","h2Coin",64,0,64,64,0,64); //This is the raw coin for 2x64 channel level
    TH1F *h1Coin = new TH1F("h1Coin","h2Coin",1000,-500,500); //Time difference spectrum L and R
    TH1F *hdEvents = new TH1F("hdEvents","hdEvents",512,-32768,32768);
    //TH2F *h2Map =  new TH2F("h2Map","h2Map",16,0,8,16,0,8); //Mapped distribution - single CHIP
    //TH2F *h2Map =  new TH2F("h2Map","h2Map",32,0,32,16,0,16); //Mapped distribution - ALL CHIPS - SIMPLE BINNING
    TH2F *h2Map =  new TH2F("h2Map","h2Map",64,0,32,32,0,16); //Mapped distribution - ALL CHIPS - OVER BINNING
    TH2F *h2tY =  new TH2F("h2tY","h2MapBP",120,0,120,32,0,16); //Y vs time
    TH2F *h2tX =  new TH2F("h2tX","h2MapBP",120,0,120,64,0,32); //X vs time
    char  hnamecoin[50]; sprintf(hnamecoin,"hCoinCh%d-Ch%d",ChA,ChB);
    TH1F *hCoinAB=new TH1F(hnamecoin,hnamecoin,1000,-500,500);// No Energy Sel
    sprintf(hnamecoin,"hCoinECh%d-Ch%d",ChA,ChB);
    TH1F *hCoinABE=new TH1F(hnamecoin,hnamecoin,1000,-500,500);//With Energy Sel
    
    //ANALYSIS - PHOTOPEAKS
    TH2F *hPPamp0 =  new TH2F("hPPamp0","hPPamp0",64,0,32,32,0,16);//PhotoPeak Amp map 2D for the entire detector B0:32X16CH
    TH2F *hPPamp1 =  new TH2F("hPPamp1","hPPamp1",64,0,32,32,0,16);//PhotoPeak Amp map 2D for the entire detector B1:32X16CH
    TH2F *hPPmean0 =  new TH2F("hPPmean0","hPPmean0",64,0,32,32,0,16);//PhotoPeak Amp map 2D for the entire detector B0:32X16CH
    TH2F *hPPmean1 =  new TH2F("hPPmean1","hPPmean1",64,0,32,32,0,16);//PhotoPeak Amp map 2D for the entire detector B1:32X16CH
    TH2F *hSingle0 =  new TH2F("hSingle0","Single board 0",32,0,32,16,0,16);
    TH2F *hSingle1 =  new TH2F("hSingle1","Single board 1",32,0,32,16,0,16);
    
    TH2F *hPeak0 =  new TH2F("hPeak0","Peak counts board 0",32,0,32,16,0,16);
    TH2F *hPeak1 =  new TH2F("hPeak1","Peak counts board 1",32,0,32,16,0,16);
    
    TH2F *hPeakC0 =  new TH2F("hPeakC0","Selected peak counts board 0",32,0,32,16,0,16);
    TH2F *hPeakC1 =  new TH2F("hPeakC1","Selected peak counts board 1",32,0,32,16,0,16);
    
    TH2F *hRes0 =  new TH2F("hRes0","Resolution board 0",32,0,32,16,0,16);
    TH2F *hRes1 =  new TH2F("hRes1","Resolution board 1",32,0,32,16,0,16);
    
    TH1F *hRes =  new TH1F("hRes","Resolution",100,0.0,0.5);
    
    /////////////////Set of Histograms 2D TIMEFRAMES qqq
    TH2F *h2Mapt[nframes];
    char h2MaptTitle[100];
    for(int i=0; i<nframes; i++)
    {
        sprintf(h2MaptTitle,"h2Mapt-%d",i);
        h2Mapt[i] = new TH2F(h2MaptTitle,h2MaptTitle,64,0,32,32,0,16);
        //cout<<"FRAMES - SETTING HIST:"<<i<<endl;
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //FIRST SCAN OF EVENTS - ENERGY HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
  //  cout<<"ncount"<<"\t"<<"gmslID"<<"\t"<<"sticID"<<"\t"<<"efine"<<"\t"<<"tcoarse"<<"\t"<<"ecoarse"<<"\t"<<"energy"<<"\t"<<"energy"<<endl;
    
    while(std::fread(buffer, sizeof(buffer),1,ptr)!='\0')
    {
        //if(ncount==0){timeOld=0; frameOld=0;}
        for(int j=0; j<nevents; j++)
        {
            for(int i=0; i<8; i++)
            {
                //a[i]=int(buffer[j*8+i]); //TESTER PACK
                //if(j<1000&&i==7){cout<<endl;cout<<"Row"<<ncount<<"\t";}
                //if(j<1000){cout<<toBinary(int(buffer[i]))<<"\t";};
            }
            
            
            if((buffer[j*8+0]/4+buffer[j*8+7]*256/4+(buffer[j*8+6]%2)*256*256/4)!=32767)
            {
                ecoarse=m_lut[ buffer[j*8+3]/32+buffer[j*8+2]*256/32+(buffer[j*8+1]%16)*256*256/32 ];
                tcoarse=m_lut[ (buffer[j*8+0]/4+buffer[j*8+7]*256/4+(buffer[j*8+6]%2)*256*256/4) ];
                //            efine= buffer[j*8+3]%32;
                //            tfine= buffer[j*8+1]/32+(buffer[j*8+0]%4)*256/32;
                channel=buffer[j*8+6]/4;
                //            frame=buffer[j*8+5]%64;
                //            tbad=(buffer[j*8+6]/2)%2;
                //            ebad=(buffer[j*8+1]/16)%2;
                sticID=buffer[j*8+5]/64;
                gmslID=buffer[j*8+4]%4;
                boardID=(buffer[j*8+4]/4)%32;
                //Data0_Sync1=buffer[j*8+4]/128;
                //RSTcountcc=oldRSTcountcc;
                ////////TIME-ENERGY CALCULATION
                //energy=ecoarse-tcoarse;
                if(ecoarse>tcoarse){energy=ecoarse-tcoarse;}
                if(ecoarse<tcoarse){energy=ecoarse+32766-tcoarse;}
                if(energy<=xmin||energy>=xmax){energy=-1;}
                //time=32*tcoarse+tfine;
              if((energy>=xmin)&&(energy<=xmax)){  
                hE[gmslID][sticID][channel]->Fill(energy);
			
					 if((gmslID==0)||(gmslID==1)){
							hSingle0->Fill(mapxAll(channel,sticID,gmslID),mapyAll(channel,sticID,gmslID));
					}
					if((gmslID==2)||(gmslID==3)){
							hSingle1->Fill(mapxAll(channel,sticID,gmslID),mapyAll(channel,sticID,gmslID));
					}
				}
			
            } //END of if statemnt
        }//END OF J dummy loop
    } // END of WHILE - event scanning
    fclose(ptr);
  //  cout<<"ENERGY READING COMPLE, FILE WILL BE CLOSED AND OPENED AGAIN"<<eventsTOTAL<<endl;
 //   ptr=fopen(fileName,"rb");
    Int_t Gevent=0, Bevent=0;

 char fileNamePeak[1000];
    sprintf(fileNamePeak,"%s_HighestPeak.txt",argName);

ofstream out;

out.open(Form(fileNamePeak));
    ////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////
    //PEAK FINDING - ENERGY HISTOGRAMS
    //////////////////////////////////////////////////////////// ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    TCanvas *cAll = new TCanvas("cAll","cAll",500,500);
    double PPmean[4][4][64], PPamp[4][4][64], PPsigma[4][4][64]; //511keV
    double PPmean2[4][4][64], PPamp2[4][4][64]; //1.2MeV gamma
    double factorgamma=1.7; //Ratio of 511keV to 1.27MeV peaks//observed are 6-8X; but 1/2 should be thresh
    double gsigma; //window of photopeak: PPmean*12%/2.355
    int noiseflag[4][4][64]; //if 511keV is not the doiminant peak, the flag is set. Helps with threhsolding
    for(int g1=0; g1<4; g1++)
    {
        for(int s1=0; s1<4; s1++)
        {
            cout<<"G"<<g1<<"S"<<s1<<"\n";
            for(int i=0; i<64; i++)
            {

                TSpectrum *sp = new TSpectrum();
                sp->Search( hE[g1][s1][i] ,1 ,"",sigP); //option: gon //0.05 toomany peaks   ///0.6
                sp->SetResolution();
                int a = sp->GetNPeaks(); //No of peaks
                Double_t *b = sp->GetPositionX(); //x vector
                Double_t *c = sp->GetPositionY(); //y vector
                int parray[a][2]; //2D array to  store the x and y values for each of the peaks
                //cout<<"Before sorting: X,Y"<<endl;
                for(int j=0; j<a; j++) //assign x,y  to an array
                {
                    parray[j][0]=b[j];//b[j] is X
                    parray[j][1]=c[j];//c[j] is Y
                    //cout<< parray[j][0]<<"\t"<< parray[j][1]<<endl;
                }
                //Identifying the dominant peaks
                Double_t *C = std::max_element(c,c+a); //Highest Y Peak
                int q = std::distance(c,C); //position of the Highest peak among peaks
                double c2 = c[q-1]; //second highest peak
                double d2 = c[q]; //fist highest peak
                //PPamp[g1][s1][i] = c[q];
                //PPmean[g1][s1][i] = b[q];
                
                //SORT the parray - peaks array: make two sort functions and simply call them.
                //Highest X peak is 1.2MeV gamma; So X sorting gives this
                std::qsort(parray, a, 2*sizeof(int), compareX); //SORT: X: peak position
                int xdist=parray[a-1][0]-parray[0][0];//position of 1.2M peak (highestX)
                PPmean2[g1][s1][i] = parray[a-1][0]; //1274keV peak X
                PPamp2[g1][s1][i] = parray[a-1][1]; //1274keV peak Y
                for(int j=a-1; j>=0; j--)
                {
                    if(parray[j][1]>factorgamma*parray[a-1][1])
                    {PPamp[g1][s1][i]=parray[j][1];
                        PPmean[g1][s1][i]=parray[j][0];
                        goto peakend; //exits the loop once this item is found
                    }
                }
				peakend:
                //PEAKS FOUND FOR 511KEV AND 1274KEV
                //511keV in PPmean[g1][s1][i], PPamp[g1][s1][i]
                //1274keV in PPmean2[g1][s1][i], PPamp2[g1][s1][i]
                
                if(PPamp[g1][s1][i]!=d2){noiseflag[g1][s1][i]=1;} //+1 indicates too low thresh
                if(PPamp[g1][s1][i]==d2){noiseflag[g1][s1][i]=0;}
                if(PPamp[g1][s1][i]==0){noiseflag[g1][s1][i]=-1;} //-1 indicates too high thresh
                if(a==0){noiseflag[g1][s1][i]=0;} //unused/masked channels
                
                double PPleft = float(PPmean2[g1][s1][i])*(1-sigP/8);
                double PPright = float(PPmean2[g1][s1][i])*(1+sigP/8);
               
                //FITTING 1274KEV GAMMA
               // PPleft = float(PPmean2[g1][s1][i])*(1-sigP);
               // PPright = float(PPmean2[g1][s1][i])*(1+sigP);
                TF1 *fG = new TF1("fG","[0]*exp(-0.5*((x-[1])/[2])^2)",PPleft,PPright);
                gsigma=sigB;
                fG->SetParameters(10.0,PPmean2[g1][s1][i],gsigma);
                hE[g1][s1][i]->Fit(fG,"RQ");
                
        //        cout<<"Fitting parameters: "<< g1<<"\t"<<s1<<"\t"<<i<<"\t"<<"  #########   "<<fG->GetParameter(0)<<"\t"<<fG->GetParameter(1)<<"\t"<<fG->GetParameter(2)<<endl;

			Double_t   peak=0,peak2=0;                  
			Double_t R1=fG->GetParameter(1)-1.5*abs(fG->GetParameter(2));
			Double_t R2=fG->GetParameter(1)+1.5*abs(fG->GetParameter(2));
			///////////////Find Background
			Double_t sou[8192]={0};
			
			TSpectrum *s = new TSpectrum();
			for (Int_t t = 0; t< 8192; t++){ sou[t]=hE[g1][s1][i]->GetBinContent(t+1);}
			s->Background(sou,8192,sigB,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing3,kFALSE);
			for (Int_t t = 0; t <8192; t++){ Bas[g1][s1][i]->SetBinContent(t +1,sou[t]);}
			//////////////////////Sinal after subtracting BG
			Sig[g1][s1][i]->Add(hE[g1][s1][i]);
			Sig[g1][s1][i]->Add(Bas[g1][s1][i],-1);
      
			fG = new TF1("fG","[0]*exp(-0.5*((x-[1])/[2])^2)",PPleft,PPright);
            gsigma=PPmean2[g1][s1][i]*0.12/2.355;
            fG->SetParameters(10.0,PPmean2[g1][s1][i],gsigma);
            Sig[g1][s1][i]->Fit(fG,"RQ");
      
			R1=fG->GetParameter(1)-abs(3*fG->GetParameter(2));
			R2=fG->GetParameter(1)+abs(3*fG->GetParameter(2));
			Double_t st=abs(6*fG->GetParameter(2))/1000;
			Int_t rb1=R1/4;
			Int_t rb2=R2/4;
		
			for(Int_t t = rb1;t<=rb2;t++){
				peak2+=Sig[g1][s1][i]->GetBinContent(t+1);
			}
		
			for(Int_t t=0;t<1000;t++){
				peak+=(fG->GetParameter(0))*exp(-0.5*(((R1+st*t)-fG->GetParameter(1))*((R1+st*t)-fG->GetParameter(1)))/((fG->GetParameter(2))*(fG->GetParameter(2))))*st;
			}
			peak=peak/4;
          
          
          px1=mapxAll(i,s1,g1);
          py1=mapyAll(i,s1,g1);
          if((g1==0)||(g1==1)){
			hPeak0->SetBinContent(px1+1,py1+1,peak);
			hRes0->SetBinContent(px1+1,py1+1,abs(fG->GetParameter(2)*2.355/fG->GetParameter(1)));
			}else{
			hPeak1->SetBinContent(px1+1,py1+1,peak);
			hRes1->SetBinContent(px1+1,py1+1,abs(fG->GetParameter(2)*2.355/fG->GetParameter(1)));	
			}
			hRes->Fill(abs(fG->GetParameter(2)*2.355/fG->GetParameter(1)));
			
			if(((a>0)&&(a<=7))&&((abs(fG->GetParameter(2))<sigB*2))&&((peak>=peakLT)&&(peak<=peakHT))){                                                                                    //// condition for good spectrum and good fit
                out<<g1<<"\t"<<s1<<"\t"<<i<<"\t"<<peak<<"\t"<<fG->GetParameter(1)<<"\t"<<abs(fG->GetParameter(2))<<"\t"<<realE<<endl;
				if((g1==0)||(g1==1)){
					hPeakC0->SetBinContent(px1+1,py1+1,peak);
				}else{
					hPeakC1->SetBinContent(px1+1,py1+1,peak);	
				}
			}else{
				out<<g1<<"\t"<<s1<<"\t"<<i<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<"\t"<<0.0<<endl;
				if((g1==0)||(g1==1)){
					hPeakC0->SetBinContent(px1+1,py1+1,0);
				}else{
					hPeakC1->SetBinContent(px1+1,py1+1,0);	
				}
			
			}
         }
            cout<<"\n"<<endl;
        }
    }
    
    out.close();

    ///////////////////////////////////////
        TCanvas *cE[4][4][4];
        char hname2[100];
    for(int g=0; g<4; g++)
    {
        for(int s=0;  s<4; s++)
        {
            //            h1PPamp[g][s]->Write();
            //            h1PPmean[g][s]->Write();
            for(int j=0; j<4; j++)
            {				
                sprintf(hname,"EBOX/Spectrum/G%dS%dE-ChSet%d.pdf",g,s,j);
                cE[g][s][j] = new TCanvas(hname,"Energy",1000,1000);
                int nx=4, ny=4;
                int number=0;
                cE[g][s][j]->Divide(nx,ny,0,0);
                for(int i=0; i<nx*ny; i++)
                {
                    number++;
                    cE[g][s][j]->cd(number);
                    px1=mapxAll(16*j+i,s,g);
					py1=mapyAll(16*j+i,s,g);
                    //h1->FillRandom("gaus",1000);
                    //hE[47+i]->GetXaxis()->SetXLimits();
                    //double xmin=0, xmax=800;
                    sprintf(hname2,"G%dS%d-Ch%d_Xindex-%d_Yindex-%d",g,s,16*j+i,px1,py1);
                    hE[g][s][16*j+i]->SetTitle(hname2);
                    hE[g][s][16*j+i]->SetAxisRange(xmin, xmax,"X");
                    hE[g][s][16*j+i]->GetXaxis()->SetLabelFont(53);
                    hE[g][s][16*j+i]->GetXaxis()->SetLabelSize(10);
                    hE[g][s][16*j+i]->GetXaxis()->SetNdivisions(5);
                    hE[g][s][16*j+i]->GetYaxis()->SetLabelFont(53);
                    hE[g][s][16*j+i]->GetYaxis()->SetLabelSize(10);
                    hE[g][s][16*j+i]->GetYaxis()->SetNdivisions(5);
                    hE[g][s][16*j+i]->DrawCopy();
                    Bas[g][s][16*j+i]->SetLineColor(1);
                    Bas[g][s][16*j+i]->DrawCopy("same");
                    Sig[g][s][16*j+i]->SetLineColor(kGreen+2);
                    Sig[g][s][16*j+i]->DrawCopy("same");
                    
                }
            //    sprintf(hname,"%s.pdf",hname);
                cE[g][s][j]->Write(); cE[g][s][j]->SaveAs(hname);
                sprintf(hname,"EBOX/Spectrum/G%dS%dE-ChSet%d.png",g,s,j);
                cE[g][s][j]->Write(); cE[g][s][j]->SaveAs(hname);
            }
        }
    }
    
    //sprintf(fhistname,"EBOX/h2Map.pdf",windowTcc, windowEvents,windowcluster);
    
    
    TCanvas *cSall= new TCanvas("cSall","cSall",1000,1000);
    gStyle->SetPalette(55);
    cSall->Divide(1,2);
    cSall->cd(1);
    hSingle0->Draw("colz");
    hSingle0->GetXaxis()->SetTitle("X index");
    hSingle0->GetYaxis()->SetTitle("Y index");
     cSall->cd(2);
    hSingle1->Draw("colz");
    hSingle1->GetXaxis()->SetTitle("X index");
    hSingle1->GetYaxis()->SetTitle("Y index");
    cSall->SaveAs("EBOX/SingleAll_v1.pdf");
    cSall->SaveAs("EBOX/SingleAll_v1.png");
    cSall->Write();
    
    TCanvas *cpeak= new TCanvas("cpeak","cpeak",1000,1000);
    gStyle->SetPalette(55);
    cpeak->Divide(1,2);
    cpeak->cd(1);
    hPeak0->Draw("colz");
    hPeak0->GetXaxis()->SetTitle("X index");
    hPeak0->GetYaxis()->SetTitle("Y index");
     cpeak->cd(2);
    hPeak1->Draw("colz");
    hPeak1->GetXaxis()->SetTitle("X index");
    hPeak1->GetYaxis()->SetTitle("Y index");
    cpeak->SaveAs("EBOX/Peak_2D_v1.pdf");
    cpeak->SaveAs("EBOX/Peak_2D_v1.png");
    cpeak->Write();
    
    TCanvas *cpeakC= new TCanvas("cpeakC","cpeakC",1000,1000);
    gStyle->SetPalette(55);
    cpeakC->Divide(1,2);
    cpeakC->cd(1);
    hPeakC0->Draw("colz");
    hPeakC0->GetXaxis()->SetTitle("X index");
    hPeakC0->GetYaxis()->SetTitle("Y index");
     cpeakC->cd(2);
    hPeakC1->Draw("colz");
    hPeakC1->GetXaxis()->SetTitle("X index");
    hPeakC1->GetYaxis()->SetTitle("Y index");
    cpeakC->SaveAs("EBOX/Peak_Selected_2D_v1.pdf");
    cpeakC->SaveAs("EBOX/Peak_Selected_2D_v1.png");
    cpeakC->Write();
    
    TCanvas *cRe= new TCanvas("cRe","cRe",1000,1000);
    gStyle->SetPalette(55);
    cRe->Divide(1,2);
    cRe->cd(1);
    hRes0->Draw("colz");
    hRes0->SetMaximum(0.1);
    hRes0->SetMinimum(0.0);
    hRes0->GetXaxis()->SetTitle("X index");
    hRes0->GetYaxis()->SetTitle("Y index");
     cRe->cd(2);
    hRes1->SetMaximum(.1);
    hRes1->SetMinimum(0.0);
    hRes1->Draw("colz");
    hRes1->GetXaxis()->SetTitle("X index");
    hRes1->GetYaxis()->SetTitle("Y index");
    cRe->SaveAs("EBOX/Resolution_2D_v1.pdf");
    cRe->SaveAs("EBOX/Resolution_2D_v1.png");
    cRe->Write();
    
    TCanvas *cR= new TCanvas("cR","cR",1000,1000);
    hRes->Draw();
    cR->SaveAs("EBOX/Resolution_Distribution_v1.pdf");
    cR->SaveAs("EBOX/Resolution_Distribution_v1.png");
    cR->Write();
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //DELETING HISTOGRAMS
    ////////////////////////////////////////////////////////////////////////////////////////////////////
//    delete hEc;
//    delete  hTc; delete hTf; delete hEf; delete hT; delete hEbad; delete hTbad; delete hFrame; delete hCh; delete hSticID; delete hGMSLID; delete hBoardID; delete h2EventEnergy;  delete h2SticIDGMSLID;
    delete hE,Bas,Sig;
    //delete  h1PPamp, h1PPmean;
//    delete hPPamp0; delete hPPamp1; delete hPPmean0; delete hPPmean1;
     
    return 0;
}
///////////THE END







