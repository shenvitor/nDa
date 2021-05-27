#include <g4root.hh>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH2D.h>

TCanvas* nImgFill1 = new TCanvas("nImg")
TH2D* 5nImg= new TH2D("5nMov", //name
                       "5nMovingImaging", //title
                       200, 0.0, 200.0,  //x-dimension
                       200, 0.0, 200.0) //y-dimension

//total bin= 200 * 200=40000
//HisFill times for each move= 40000 / 25 = 1600 (= 40*40)
// 5nImg->Fill(x,y,n)

//TFile *f1 = new TFile ("1.root")
//TTree *t1 = (TTree*) f1->Get("neutron")
//int t=1;
//for (int i=0; i<1600;i++){
//t+=5
// get entry of InRowNb
// get entry of LevelNb

//*******************
//(TFile *) f1 = new TFile ("5n-LD22-2Dimage.root")
//(TCanvas *) Canvas_1_n2
//(TH2f *) htemp
//int x=1, y=1;
//for (x; x<1601; x++){
    //(y; y<1601; y++){
//double n = htemp->GetBinContent(x,y)
//}
//}