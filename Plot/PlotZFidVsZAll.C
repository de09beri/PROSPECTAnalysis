#include "PROSPECT_Style.cc"
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TPave.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TColor.h"
#include "TExec.h"
#include <sstream>
#include "TLatex.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Header.C"

//Function for creating a relative to mean graph
TGraphErrors *makeRelGr(TGraphErrors *gr){
	gr->Fit("pol0","Q0");
	double mean = gr->GetFunction("pol0")->GetParameter(0);
	double meanErr = gr->GetFunction("pol0")->GetParError(0);

	TGraphErrors *grRel = (TGraphErrors*)gr->Clone();
	int numPt = gr->GetN();

	double rel, relErr; 
	double grx, gry, gryErr;

	for(int i=0;i<numPt;i++){
		gr->GetPoint(i,grx,gry);
		gryErr = gr->GetErrorY(i);

		rel = gry/mean;
		relErr = rel * sqrt(pow(gryErr/gry,2) + pow(meanErr/mean,2));	

		grRel->SetPoint(i,grx,rel);
		grRel->SetPointError(i,0,relErr);
	}

	return grRel;
}

TGraphErrors *makeResidGr(TGraphErrors *gr0, TGraphErrors *gr1){
	TGraphErrors *grResid = (TGraphErrors*)gr0->Clone();
	int numPt = gr0->GetN();
	
	double diff, diffErr;
	double grx0, gry0, gry0Err;
	double grx1, gry1, gry1Err;

	for(int i=0;i<numPt;i++){
		gr0->GetPoint(i,grx0,gry0);
		gry0Err = gr0->GetErrorY(i);	

		gr1->GetPoint(i,grx1,gry1);
		gry1Err = gr1->GetErrorY(i);	

		diff = gry0-gry1;
		diffErr = sqrt(pow(gry0Err,2) + pow(gry1Err,2));

		grResid->SetPoint(i,grx0,diff);
		grResid->SetPointError(i,0,diffErr);
	}

	return grResid;
}

TGraphErrors *makeETGr(TGraphErrors *gr){
        const int numET = 46;
        int ET[numET] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,27,28,41,42,55,56,69,70,83,84,97,98,111,112,125,126,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153};

        int numPt = gr->GetN();

        vector<double> vGrx, vGry, vGryErr;
        double grx, gry, gryErr;
        bool boolET;

        for(int i=0;i<numPt;i++){
                gr->GetPoint(i,grx,gry);
                gryErr = gr->GetErrorY(i);

                boolET = find(begin(ET), end(ET), grx) != end(ET);
                if(boolET){
                        vGrx.push_back(grx);
                        vGry.push_back(gry);
                        vGryErr.push_back(gryErr);
                }
        }

        int size = vGrx.size();
        TGraphErrors *grET = new TGraphErrors(size);

        for(int i=0;i<size;i++){
                grET->SetPoint(i,vGrx[i],vGry[i]);
                grET->SetPointError(i,0,vGryErr[i]);
        }

        return grET;
}


TGraphErrors *makeNoETGr(TGraphErrors *gr){
        const int numET = 46;
        int ET[numET] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,27,28,41,42,55,56,69,70,83,84,97,98,111,112,125,126,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153};

        int numPt = gr->GetN();

        vector<double> vGrx, vGry, vGryErr;
        double grx, gry, gryErr;
        bool boolET;

        for(int i=0;i<numPt;i++){
                gr->GetPoint(i,grx,gry);
                gryErr = gr->GetErrorY(i);

                boolET = find(begin(ET), end(ET), grx) != end(ET);
                if(!boolET){
                        vGrx.push_back(grx);
                        vGry.push_back(gry);
                        vGryErr.push_back(gryErr);
                }
        }

        int size = vGrx.size();
        TGraphErrors *grET = new TGraphErrors(size);

        for(int i=0;i<size;i++){
                grET->SetPoint(i,vGrx[i],vGry[i]);
                grET->SetPointError(i,0,vGryErr[i]);
        }

        return grET;
}


int PlotZFidVsZAll(){

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	TFile *fz = new TFile(Form("%s/Ac227_GraphsPerCell_AllZ.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"))); 
	if(!fz){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}
	
	TGraphErrors *grRate_z = (TGraphErrors*)fz->Get("grRate");
	TGraphErrors *grRelRate_z = makeRelGr(grRate_z);
	TGraphErrors *grETRelRate_z = makeETGr(grRelRate_z);
	TGraphErrors *grPoLifetime_z = (TGraphErrors*)fz->Get("grLifetime");
	TGraphErrors *grETPoLifetime_z = makeETGr(grPoLifetime_z);
	TGraphErrors *grTotEff_z = (TGraphErrors*)fz->Get("grTotEff");
	TGraphErrors *grETTotEff_z = makeETGr(grTotEff_z);	

	//-----------------------
	TFile *fzFid = new TFile(Form("%s/Ac227_GraphsPerCell_zFid.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"))); 
	if(!fzFid){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}
	
	TGraphErrors *grRate_zFid = (TGraphErrors*)fzFid->Get("grRate");
	TGraphErrors *grRelRate_zFid = makeRelGr(grRate_zFid);
	TGraphErrors *grETRelRate_zFid = makeETGr(grRelRate_zFid);
	TGraphErrors *grPoLifetime_zFid = (TGraphErrors*)fzFid->Get("grLifetime");
	TGraphErrors *grETPoLifetime_zFid = makeETGr(grPoLifetime_zFid);
	TGraphErrors *grTotEff_zFid = (TGraphErrors*)fzFid->Get("grTotEff");
	TGraphErrors *grETTotEff_zFid = makeETGr(grTotEff_zFid);	

		
	TGraphErrors *grResid = makeResidGr(grRelRate_z,grRelRate_zFid);
	TGraphErrors *grETResid = makeETGr(grResid);

	TGraphErrors *grNoETRelRate_z = makeNoETGr(grRelRate_z);
	TGraphErrors *grNoETRelRate_zFid = makeNoETGr(grRelRate_zFid);
	TGraphErrors *grNoETResid = makeNoETGr(grResid);

	//-----------------------
	grResid->SetMarkerStyle(24);
	grResid->SetMarkerColor(9);
	grResid->SetLineColor(9);
	grETResid->SetMarkerStyle(20);
	grETResid->SetMarkerColor(9);
	grETResid->SetLineColor(9);
	grNoETResid->SetMarkerStyle(20);
	grNoETResid->SetMarkerColor(9);
	grNoETResid->SetLineColor(9);

	grRelRate_z->SetMarkerStyle(24);
	grRelRate_zFid->SetMarkerStyle(24);
	grNoETRelRate_z->SetMarkerStyle(24);
	grNoETRelRate_zFid->SetMarkerStyle(24);
	grPoLifetime_z->SetMarkerStyle(24);
	grPoLifetime_zFid->SetMarkerStyle(24);
	grTotEff_z->SetMarkerStyle(24);
	grTotEff_zFid->SetMarkerStyle(24);
	
	grRelRate_zFid->SetMarkerColor(kRed);
	grRelRate_zFid->SetLineColor(kRed);
	grRelRate_z->SetTitle(";Segment;Relative Rate [mHz]");
	grResid->SetTitle(";Segment;z-z_{Fid}");
	grNoETRelRate_zFid->SetMarkerColor(kRed);
	grNoETRelRate_zFid->SetLineColor(kRed);
	grNoETRelRate_z->SetTitle(";Segment;Relative Rate [mHz]");
	grNoETResid->SetTitle(";Segment;z-z_{Fid}");

	grPoLifetime_zFid->SetMarkerColor(kRed);
	grPoLifetime_zFid->SetLineColor(kRed);
	grTotEff_zFid->SetMarkerColor(kRed);
	grTotEff_zFid->SetLineColor(kRed);

	grETRelRate_z->SetMarkerStyle(20);
	grETRelRate_zFid->SetMarkerColor(kRed);
	grETRelRate_zFid->SetLineColor(kRed);
	grETRelRate_zFid->SetMarkerStyle(20);
	grETPoLifetime_z->SetMarkerStyle(20);
	grETPoLifetime_zFid->SetMarkerColor(kRed);
	grETPoLifetime_zFid->SetLineColor(kRed);
	grETPoLifetime_zFid->SetMarkerStyle(20);
	grETTotEff_z->SetMarkerStyle(20);
	grETTotEff_zFid->SetMarkerColor(kRed);
	grETTotEff_zFid->SetLineColor(kRed);
	grETTotEff_zFid->SetMarkerStyle(20);

	//-----------------------
	TCanvas *cRate = new TCanvas("cRate","Rate Per Cell",1000,700);
	TPad *pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,0);
        TPad *pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,0);
        pad1->Draw();
        pad2->Draw();

	pad1->cd();
	grRelRate_z->Draw("AP");
	grRelRate_z->GetYaxis()->SetRangeUser(0.94,1.04);
	grETRelRate_z->Draw("P");
	grRelRate_zFid->Draw("P");	
	grETRelRate_zFid->Draw("P");
	TLegend *leg = new TLegend(0.8,0.85,0.99,0.99);
	leg->AddEntry(grRelRate_z,"-1000 < z < 1000 mm","l");
	leg->AddEntry(grRelRate_zFid,"-444 < z < 444 mm","l");
	leg->Draw();

	pad2->cd();
	grResid->Draw("AP");
	grResid->Fit("pol0");
	grETResid->Draw("P");
	grResid->GetFunction("pol0")->SetLineStyle(2);
	grResid->GetYaxis()->SetTitleSize(0.07);
        grResid->GetYaxis()->SetLabelSize(0.07);
        grResid->GetYaxis()->SetTitleOffset(0.4);
        grResid->GetXaxis()->SetTitleSize(0.07);
        grResid->GetXaxis()->SetLabelSize(0.07);
        grResid->GetXaxis()->SetTitleOffset(0.5);
        TPaveText *pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
        pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",grResid->GetFunction("pol0")->GetChisquare(),grResid->GetFunction("pol0")->GetNDF()));
        pt->AddText(Form("Prob   %.3f",grResid->GetFunction("pol0")->GetProb()));
        pt->AddText(Form("p0   %f #pm %f",grResid->GetFunction("pol0")->GetParameter(0),grResid->GetFunction("pol0")->GetParError(0)));
        pt->Draw();
	cRate->SaveAs(Form("%s/RateVsCell_ZvsZFid.pdf",gSystem->Getenv("AD_AC227_PLOTS")));
 	
	//-----------------------
	TCanvas *cNoETRate = new TCanvas("cNoETRate","Rate no ET",1000,700);
	pad1->Draw();
	pad2->Draw();
	pad1->cd();
	grNoETRelRate_z->Draw("AP");
	grNoETRelRate_z->GetYaxis()->SetRangeUser(0.96,1.04);
	grNoETRelRate_zFid->Draw("P");
	leg = new TLegend(0.8,0.85,0.99,0.99);
	leg->AddEntry(grNoETRelRate_z,"-1000 < z < 1000 mm","l");
	leg->AddEntry(grNoETRelRate_zFid,"-444 < z < 444 mm","l");
	leg->Draw();

	pad2->cd();
	grNoETResid->Draw("AP");
	grNoETResid->Fit("pol0");	
	grNoETResid->GetFunction("pol0")->SetLineStyle(2);
	grNoETResid->GetYaxis()->SetTitleSize(0.07);
        grNoETResid->GetYaxis()->SetLabelSize(0.07);
        grNoETResid->GetYaxis()->SetTitleOffset(0.4);
        grNoETResid->GetXaxis()->SetTitleSize(0.07);
        grNoETResid->GetXaxis()->SetLabelSize(0.07);
        grNoETResid->GetXaxis()->SetTitleOffset(0.5);
        pt = new TPaveText(0.8,0.8,0.99,0.99,"NDCNB");
        pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",grNoETResid->GetFunction("pol0")->GetChisquare(),grNoETResid->GetFunction("pol0")->GetNDF()));
        pt->AddText(Form("Prob   %.3f",grNoETResid->GetFunction("pol0")->GetProb()));
        pt->AddText(Form("p0   %f #pm %f",grNoETResid->GetFunction("pol0")->GetParameter(0),grNoETResid->GetFunction("pol0")->GetParError(0)));
        pt->Draw();
	cNoETRate->SaveAs(Form("%s/RateVsCell_NoET_ZvsZFid.pdf",gSystem->Getenv("AD_AC227_PLOTS")));


	//-----------------------
	double grLifetimexStart, grLifetimexEnd, grLifetimeyStart, grLifetimeyEnd;
        grPoLifetime_z->GetPoint(0,grLifetimexStart,grLifetimeyStart);    
        grPoLifetime_z->GetPoint(grPoLifetime_z->GetN()-1,grLifetimexEnd,grLifetimeyEnd);      
        grLifetimexStart = grLifetimexStart - 3;
        grLifetimexEnd = grLifetimexEnd + 3;  

	TLine *lLifetime = new TLine(grLifetimexStart,2.569,grLifetimexEnd,2.569);
        lLifetime->SetLineStyle(2);
        lLifetime->SetLineColor(30);

        TGraph *grlLifetime = new TGraph(4);
        grlLifetime->SetPoint(0,grLifetimexStart,2.562);
        grlLifetime->SetPoint(1,grLifetimexEnd,2.562);
        grlLifetime->SetPoint(2,grLifetimexEnd,2.576);
        grlLifetime->SetPoint(3,grLifetimexStart,2.576);
        grlLifetime->SetFillStyle(3001);
        grlLifetime->SetFillColor(30);

	grPoLifetime_z->SetTitle(";Segment;Po^{215} #tau");

	TCanvas *cLifetime = new TCanvas("cLifetime","Po Lifetime",1000,400);
	grPoLifetime_z->Draw("AP");
	lLifetime->Draw("same");
	grlLifetime->Draw("f");
	grPoLifetime_z->Draw("P");
	grPoLifetime_zFid->Draw("P");
	grETPoLifetime_z->Draw("P");
	grETPoLifetime_zFid->Draw("P");
	grPoLifetime_z->GetYaxis()->SetRangeUser(2.48,2.67);
	grPoLifetime_z->Fit("pol0");
	grPoLifetime_zFid->Fit("pol0");
	grPoLifetime_z->GetFunction("pol0")->SetLineColor(kBlue+2);
	grPoLifetime_z->GetFunction("pol0")->SetLineStyle(3);
	grPoLifetime_z->GetFunction("pol0")->SetLineWidth(1);
	grPoLifetime_zFid->GetFunction("pol0")->SetLineStyle(10);
	grPoLifetime_zFid->GetFunction("pol0")->SetLineWidth(1);
	leg = new TLegend(0.8,0.83,0.99,0.99);
	leg->AddEntry(grPoLifetime_z,"-1000 < z < 1000 mm","l");
	leg->AddEntry(grPoLifetime_zFid,"-444 < z < 444 mm","l");
	leg->Draw();
        pt = new TPaveText(0.1,0.76,0.29,0.92,"NDCNB");
        pt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",grPoLifetime_z->GetFunction("pol0")->GetChisquare(),grPoLifetime_z->GetFunction("pol0")->GetNDF()));
        pt->AddText(Form("Prob   %.2f",grPoLifetime_z->GetFunction("pol0")->GetProb()));
        pt->AddText(Form("#tau   %.4f #pm %.4f",grPoLifetime_z->GetFunction("pol0")->GetParameter(0),grPoLifetime_z->GetFunction("pol0")->GetParError(0)));
	pt->SetFillColorAlpha(kBlue+2,0.4);
        pt->Draw();
        TPaveText *ptt = new TPaveText(0.3,0.76,0.49,0.92,"NDCNB");
        ptt->AddText(Form("#Chi^{2}/NDF   %.2f/%d",grPoLifetime_zFid->GetFunction("pol0")->GetChisquare(),grPoLifetime_zFid->GetFunction("pol0")->GetNDF()));
        ptt->AddText(Form("Prob   %.2f",grPoLifetime_zFid->GetFunction("pol0")->GetProb()));
        ptt->AddText(Form("#tau   %.4f #pm %.4f",grPoLifetime_zFid->GetFunction("pol0")->GetParameter(0),grPoLifetime_zFid->GetFunction("pol0")->GetParError(0)));
	ptt->SetFillColorAlpha(kRed,0.4);
        ptt->Draw();
	cLifetime->SaveAs(Form("%s/LifetimeVsCell_ZvsZFid.pdf",gSystem->Getenv("AD_AC227_PLOTS")));

	//-----------------------
	grTotEff_z->SetTitle(";Segment;Efficiency");

	TCanvas *cTotEff = new TCanvas("cTotEff","Total Eff.",1000,400);
	grTotEff_z->Draw("AP");
	grETTotEff_z->Draw("P");
	grTotEff_zFid->Draw("P");
	grETTotEff_zFid->Draw("P");	
	grTotEff_z->GetYaxis()->SetTitleOffset(0.9);
	leg = new TLegend(0.8,0.85,0.99,0.99);
	leg->AddEntry(grTotEff_z,"-1000 < z < 1000 mm","l");
	leg->AddEntry(grTotEff_zFid,"-444 < z < 444 mm","l");
	leg->Draw();
	cTotEff->SaveAs(Form("%s/TotEffVsCell_ZvsZFid.pdf",gSystem->Getenv("AD_AC227_PLOTS")));


	return 0;
}
