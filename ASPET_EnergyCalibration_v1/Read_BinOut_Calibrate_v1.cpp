//Tran Cong Thien: 16th Dec 2022
// This program read at least 3 _HighestPeak.txt for energy calibration for each channel of the ASPET
//This program required 2 folders in the same dirrectory as: EnergyCal and EnergyCal/fitting

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


using namespace std;

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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////MAIN PROGRAM
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	
	char argName[1000], fileName[1000];
	ifstream in;
	double load[argc][4][4][64][4]={0};
	int g1,s1,ch1;
	double  peakArea, ToT, res, realE;    /////////////// loading features
	
	for(int f=0;f<argc;f++){
		sprintf(argName,"%s",argv[f+1]);
		sprintf(fileName,"%s.txt",argName);
		in.open(Form(fileName));
		 while(1){
			if (!in.good()) break;  
			in>>g1>>s1>>ch1>>peakArea>>ToT>>res>>realE;             											////// Loading all data in to array load[file][g][s][ch][figures]
			load[f][g1][s1][ch1][0] = peakArea;
			load[f][g1][s1][ch1][1] = ToT;
			load[f][g1][s1][ch1][2] = res;
			load[f][g1][s1][ch1][3] = realE;
		}
		in.close();
	}
	

	
	
	
	ofstream out;
	out.open(Form("EnergyCal/Energy_Calibration_v1.txt"));
	
	
	
	TH1D *h1=new TH1D("h1","",100,0,2.0);
	TH2D *hChi0=new TH2D("hChi0","Chi2 values of board 0",32,0,32,16,0,16);
	TH2D *hChi1=new TH2D("hChi1","Chi2 values of board 1",32,0,32,16,0,16);
	
	double p0,p1,p2,chi2;
	double pr0,pr1,pr2;
	TGraphErrors *gr1[4][4][64];
	int fn=0,zero=1;
	int count=0;
	int px1,py1;
	
 for(int g2=0;g2<4;g2++){
	for(int s2=0;s2<4;s2++){
		for(int ch2=0;ch2<64;ch2++){
			p0=0;p1=0;p2=0;
			pr0=0;pr1=0;pr2=0;
			chi2=1000;
			fn=0;
			gr1[g2][s2][ch2] = new TGraphErrors();
			for(int c=0;c<argc;c++){
					zero=1;
					if(load[c][g2][s2][ch2][0]==0){zero=0;}	
					if(zero==1){
							gr1[g2][s2][ch2]->SetPoint(fn,load[c][g2][s2][ch2][3],load[c][g2][s2][ch2][1]);
							gr1[g2][s2][ch2]->SetPointError(fn,1.0,load[c][g2][s2][ch2][2]);
							fn++;
					}	
			}
					
			if(fn>=2){
				TF1 *fG = new TF1("fG","[0]*(exp(x/[1]))+[2]",190,1400);
				fG->SetParameters(-380.78,-432.601,565.605);
				gr1[g2][s2][ch2]->Fit(fG,"RQ");
				p0=fG->GetParameter(0);
				p1=fG->GetParameter(1);
				p2=fG->GetParameter(2);
				pr0=fG->GetParError(0);
				pr1=fG->GetParError(1);
				pr2=fG->GetParError(2);
				chi2=fG->GetChisquare();		
				h1->Fill(chi2);
			}
			out<<g2<<"\t"<<s2<<"\t"<<ch2<<"\t"<<p0<<"\t"<<pr0<<"\t"<<p1<<"\t"<<pr1<<"\t"<<p2<<"\t"<<pr2<<"\t"<<chi2<<endl;
			
			if(chi2<1){count++;}	
			px1=mapxAll(ch2,s2,g2);
			py1=mapyAll(ch2,s2,g2);
			if((g2==0)||(g2==1)){
				hChi0->SetBinContent(px1+1,py1+1,chi2);
			}else{
				hChi1->SetBinContent(px1+1,py1+1,chi2);
			}
					
		}
	}
}

double meanChi=h1->GetMean();

 
 
//############################################################################################
out.close();

   ///////////////////////////////////////
        TCanvas *cE[4][4][4];
        char hname[50],hname2[50];
    for(Int_t g=0; g<4; g++)
    {
        for(Int_t s=0;  s<4; s++)
        {
            for(Int_t j=0; j<4; j++)
            {				
                sprintf(hname,"EnergyCal/fitting/G%dS%dE-ChSet%d.png",g,s,j);
                cE[g][s][j] = new TCanvas(hname,"Energy fitting",1000,1000);
                Int_t nx=4, ny=4;
                Int_t number=0;
                cE[g][s][j]->Divide(nx,ny,0,0);
                for(Int_t i=0; i<16; i++)
                {
                    number++;
                    cE[g][s][j]->cd(number);
                    cE[g][s][j]->cd(number)->SetGrid(1);
                    px1=mapxAll(16*j+i,s,g);
					py1=mapyAll(16*j+i,s,g);
                    sprintf(hname2,"G%dS%d-Ch%d---Xindex-%d_Yindex-%d",g,s,16*j+i,px1,py1);
                    gr1[g][s][16*j+i]->SetTitle(hname2);
                //    gr1[g][s][16*j+i]->SetAxisRange(0,1500,"X");
                //    gr1[g][s][16*j+i]->SetAxisRange(0,1000,"Y");
                    gr1[g][s][16*j+i]->GetXaxis()->SetLabelFont(53);
                    gr1[g][s][16*j+i]->GetXaxis()->SetLabelSize(10);
                    gr1[g][s][16*j+i]->GetXaxis()->SetTitle("Corresponding energy (keV)");
                    gr1[g][s][16*j+i]->GetYaxis()->SetTitle("ToT (ADC)");
                    gr1[g][s][16*j+i]->GetXaxis()->SetNdivisions(5);
                    gr1[g][s][16*j+i]->GetYaxis()->SetLabelFont(53);
                    gr1[g][s][16*j+i]->GetYaxis()->SetLabelSize(10);
                    gr1[g][s][16*j+i]->GetYaxis()->SetNdivisions(5);
                    gr1[g][s][16*j+i]->GetXaxis()->SetRangeUser(0,1500);
                    gr1[g][s][16*j+i]->GetYaxis()->SetRangeUser(200,700);
                    gr1[g][s][16*j+i]->SetMarkerStyle(22);
					gr1[g][s][16*j+i]->SetMarkerSize(1.1);
                    gr1[g][s][16*j+i]->Draw("AP");
                    
                }
            //    sprintf(hname,"%s.pdf",hname);
                cE[g][s][j]->Write(); cE[g][s][j]->SaveAs(hname);
            }
        }
    }
    sprintf(hname,"Chi2<1.0_%d-channels",count);
    
 TCanvas* c1 = new TCanvas("c1","histogram",100,10,1000,1000);//700-rong, 900-dai window hien thi
 gStyle->SetPalette(55);
 c1->SetGrid(1);
 h1->SetTitle(hname);
 h1->GetXaxis()->SetTitle("Chi2");
 h1->GetYaxis()->SetTitle("Number of channels");
 h1->Draw();
 c1->SaveAs("EnergyCal/Chi2.png");   
  c1->SaveAs("EnergyCal/Chi2.pdf");   
  
 TCanvas* c11 = new TCanvas("c11","histogram",100,10,1000,1000);//700-rong, 900-dai window hien thi
gStyle->SetPalette(55);
c11->Divide(1,2);
c11->cd(1);
hChi0->Draw("colz");
hChi0->SetMaximum(meanChi*5);
hChi0->SetMinimum(0);
c11->cd(2);
hChi1->Draw("colz");
hChi1->SetMaximum(meanChi*5);
hChi1->SetMinimum(0);
 c11->SaveAs("EnergyCal/Chi2_2D.png");   
 c11->SaveAs("EnergyCal/Chi2_2D.pdf"); 
  
/*
 TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
 leg->AddEntry(h1,"201keV","L");
leg->AddEntry(h2,"306keV","L");
leg->AddEntry(h3,"511keV","L");
 leg->AddEntry(h4,"1274keV","L");
//leg->Draw();
 */

return 0;

}
