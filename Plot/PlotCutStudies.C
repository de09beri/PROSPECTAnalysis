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

/*TGraphErrors *makeResidGr(TGraphErrors *gr0, TGraphErrors *gr1){
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
		//grResid->SetPointError(i,0,diffErr);
		grResid->SetPointError(i,0,0);
	}

	return grResid;
} */

TGraph *makeResidGr(TGraphErrors *gr0, TGraphErrors *gr1){
//	TGraphErrors *grResid = (TGraphErrors*)gr0->Clone();
	int numPt = gr0->GetN();

	TGraph *grResid = new TGraph(numPt);
	
	double diff;
	double grx0, gry0;
	double grx1, gry1;

	for(int i=0;i<numPt;i++){
		gr0->GetPoint(i,grx0,gry0);
		gr1->GetPoint(i,grx1,gry1);
		diff = gry0-gry1;

		grResid->SetPoint(i,grx0,diff);
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


int PlotCutStudies(){

	setup_PROSPECT_style();
  	gROOT->ForceStyle();

	string var = "PE";	
	string pdfVar = "PromptE";

	TFile *f0 = new TFile(Form("%s/Ac227_GraphsPerCell_40SigmaPE.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"))); 
	if(!f0){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}

	TGraphErrors *grRate_0 = (TGraphErrors*)f0->Get("grRate");
	TGraphErrors *grRelRate_0 = makeRelGr(grRate_0);
	//TGraphErrors *grResid_0 = makeResidGr(grRelRate_0,grRelRate_0);
	TGraph *grRelResid_0 = makeResidGr(grRelRate_0,grRelRate_0);
	TGraph *grResid_0 = makeResidGr(grRate_0,grRate_0);

	TGraphErrors *grEff_0 = (TGraphErrors*)f0->Get("grTotEff");


	//-----------------------
	TFile *f1 = new TFile(Form("%s/Ac227_GraphsPerCell_38Sigma%s.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"),var.c_str())); 
	if(!f1){
		printf("Graph file not found. Exiting. \n");
		return -1;
	}

	TGraphErrors *grRate_1 = (TGraphErrors*)f1->Get("grRate");
	TGraphErrors *grRelRate_1 = makeRelGr(grRate_1);
	//TGraphErrors *grResid_1 = makeResidGr(grRelRate_0,grRelRate_1);
	TGraph *grRelResid_1 = makeResidGr(grRelRate_0,grRelRate_1);
	TGraph *grResid_1 = makeResidGr(grRate_0,grRate_1);

	TGraphErrors *grEff_1 = (TGraphErrors*)f1->Get("grTotEff");

	//-----------------------
	TFile *f2 = new TFile(Form("%s/Ac227_GraphsPerCell_35Sigma%s.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"),var.c_str())); 
	if(!f2){
		printf("Graph file not found. Exiting. \n");
		return -2;
	}

	TGraphErrors *grRate_2 = (TGraphErrors*)f2->Get("grRate");
	TGraphErrors *grRelRate_2 = makeRelGr(grRate_2);
	//TGraphErrors *grResid_2 = makeResidGr(grRelRate_0,grRelRate_2);
	TGraph *grRelResid_2 = makeResidGr(grRelRate_0,grRelRate_2);
	TGraph *grResid_2 = makeResidGr(grRate_0,grRate_2);

	TGraphErrors *grEff_2 = (TGraphErrors*)f2->Get("grTotEff");

	//-----------------------
	TFile *f3 = new TFile(Form("%s/Ac227_GraphsPerCell_30Sigma%s.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"),var.c_str())); 
	if(!f3){
		printf("Graph file not found. Exiting. \n");
		return -3;
	}

	TGraphErrors *grRate_3 = (TGraphErrors*)f3->Get("grRate");
	TGraphErrors *grRelRate_3 = makeRelGr(grRate_3);
	//TGraphErrors *grResid_3 = makeResidGr(grRelRate_0,grRelRate_3);
	TGraph *grRelResid_3 = makeResidGr(grRelRate_0,grRelRate_3);
	TGraph *grResid_3 = makeResidGr(grRate_0,grRate_3);

	TGraphErrors *grEff_3 = (TGraphErrors*)f3->Get("grTotEff");

	//-----------------------
	TFile *f4 = new TFile(Form("%s/Ac227_GraphsPerCell_25Sigma%s.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS"),var.c_str())); 
	if(!f4){
		printf("Graph file not found. Exiting. \n");
		return -4;
	}

	TGraphErrors *grRate_4 = (TGraphErrors*)f4->Get("grRate");
	TGraphErrors *grRelRate_4 = makeRelGr(grRate_4);
	//TGraphErrors *grResid_4 = makeResidGr(grRelRate_0,grRelRate_4);
	TGraph *grRelResid_4 = makeResidGr(grRelRate_0,grRelRate_4);
	TGraph *grResid_4 = makeResidGr(grRate_0,grRate_4);

	TGraphErrors *grEff_4 = (TGraphErrors*)f4->Get("grTotEff");
	
	//-----------------------
	double grx,gry;
	
	TH1F *hRelRate_0 = new TH1F("hRelRate_0",";R_{i}/#LTR#GT;Counts",30,0.97,1.03);
	for(int i=0;i<grRelRate_0->GetN();i++){
		grRelRate_0->GetPoint(i,grx,gry);
		hRelRate_0->Fill(gry);
	}

	TH1F *hRelRate_1 = new TH1F("hRelRate_1",";R_{i}/#LTR#GT;Counts",30,0.97,1.03);
	TH1F *hRelResid_1 = new TH1F("hRelResid_1",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-5e-3,5e-3);
	TH1F *hResid_1 = new TH1F("hResid_1",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-6e-3,30e-3);
	for(int i=0;i<grRelRate_1->GetN();i++){
		grRelRate_1->GetPoint(i,grx,gry);
		hRelRate_1->Fill(gry);
		grRelResid_1->GetPoint(i,grx,gry);
		hRelResid_1->Fill(gry);
		grResid_1->GetPoint(i,grx,gry);
		hResid_1->Fill(gry);
	}

	TH1F *hRelRate_2 = new TH1F("hRelRate_2",";R_{i}/#LTR#GT;Counts",30,0.97,1.03);
	TH1F *hRelResid_2 = new TH1F("hRelResid_2",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-5e-3,5e-3);
	TH1F *hResid_2 = new TH1F("hResid_2",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-6e-3,30e-3);
	for(int i=0;i<grRelRate_2->GetN();i++){
		grRelRate_2->GetPoint(i,grx,gry);
		hRelRate_2->Fill(gry);
		grRelResid_2->GetPoint(i,grx,gry);
		hRelResid_2->Fill(gry);
		grResid_2->GetPoint(i,grx,gry);
		hResid_2->Fill(gry);
	}

	TH1F *hRelRate_3 = new TH1F("hRelRate_3",";R_{i}/#LTR#GT;Counts",30,0.97,1.03);
	TH1F *hRelResid_3 = new TH1F("hRelResid_3",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-5e-3,5e-3);
	TH1F *hResid_3 = new TH1F("hResid_3",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-6e-3,30e-3);
	for(int i=0;i<grRelRate_3->GetN();i++){
		grRelRate_3->GetPoint(i,grx,gry);
		hRelRate_3->Fill(gry);
		grRelResid_3->GetPoint(i,grx,gry);
		hRelResid_3->Fill(gry);
		grResid_3->GetPoint(i,grx,gry);
		hResid_3->Fill(gry);
	}

	TH1F *hRelRate_4 = new TH1F("hRelRate_4",";R_{i}/#LTR#GT;Counts",30,0.97,1.03);
	TH1F *hRelResid_4 = new TH1F("hRelResid_4",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-5e-3,5e-3);
	TH1F *hResid_4 = new TH1F("hResid_4",";(R_{i}/#LTR#GT)_{4#sigma} - (R_{i}/#LTR#GT)_{j#sigma};Counts",40,-6e-3,30e-3);
	for(int i=0;i<grRelRate_4->GetN();i++){
		grRelRate_4->GetPoint(i,grx,gry);
		hRelRate_4->Fill(gry);
		grRelResid_4->GetPoint(i,grx,gry);
		hRelResid_4->Fill(gry);
		grResid_4->GetPoint(i,grx,gry);
		hResid_4->Fill(gry);
	}

	//-----------------------
	//int col0 = 1, col1 = 2, col2 = 8, col3 = 4, col4 = 6;
	int col0 = 1, col1 = 46, col2 = 41, col3 = 28, col4 = 38;

	hRelRate_0->SetLineWidth(2);
	hRelRate_1->SetLineWidth(2);
	hRelResid_1->SetLineWidth(2);
	hResid_1->SetLineWidth(2);
	hResid_1->SetLineStyle(2);
	hRelRate_2->SetLineWidth(2);
	hRelResid_2->SetLineWidth(2);
	hResid_2->SetLineWidth(2);
	hResid_2->SetLineStyle(2);
	hRelRate_3->SetLineWidth(2);
	hRelResid_3->SetLineWidth(2);
	hResid_3->SetLineWidth(2);
	hResid_3->SetLineStyle(2);
	hRelRate_4->SetLineWidth(2);
	hRelResid_4->SetLineWidth(2);
	hResid_4->SetLineWidth(2);
	hResid_4->SetLineStyle(2);

	grRate_0->SetMarkerColor(col0);
	grRate_0->SetMarkerColor(col0);
	grRelRate_0->SetLineColor(col0);
	grRelRate_0->SetLineColor(col0);
	hRelRate_0->SetLineColor(col0);

	grRate_1->SetMarkerColor(col1);	
	grRate_1->SetLineColor(col1);
	grRelRate_1->SetMarkerColor(col1);	
	grRelRate_1->SetLineColor(col1);
	hRelRate_1->SetLineColor(col1);
	hResid_1->SetLineColor(col1);
	hRelResid_1->SetLineColor(col1);

	grRate_2->SetMarkerColor(col2);	
	grRate_2->SetLineColor(col2);
	grRelRate_2->SetMarkerColor(col2);	
	grRelRate_2->SetLineColor(col2);
	hRelRate_2->SetLineColor(col2);
	hResid_2->SetLineColor(col2);
	hRelResid_2->SetLineColor(col2);

	grRate_3->SetMarkerColor(col3);	
	grRate_3->SetLineColor(col3);
	grRelRate_3->SetMarkerColor(col3);	
	grRelRate_3->SetLineColor(col3);
	hRelRate_3->SetLineColor(col3);
	hResid_3->SetLineColor(col3);
	hRelResid_3->SetLineColor(col3);

	grRate_4->SetMarkerColor(col4);	
	grRate_4->SetLineColor(col4);
	grRelRate_4->SetMarkerColor(col4);	
	grRelRate_4->SetLineColor(col4);
	hRelRate_4->SetLineColor(col4);
	hResid_4->SetLineColor(col4);
	hRelResid_4->SetLineColor(col4);


	grResid_0->Fit("pol0");
	grRelResid_0->Fit("pol0");
	printf("Fit 38SigmaPE \n");
	grResid_1->Fit("pol0");
	grRelResid_1->Fit("pol0");
	printf("Fit 35SigmaPE \n");
	grResid_2->Fit("pol0");
	grRelResid_2->Fit("pol0");
	printf("Fit 30SigmaPE \n");
	grResid_3->Fit("pol0");
	grRelResid_3->Fit("pol0");
	printf("Fit 25SigmaPE \n");
	grResid_4->Fit("pol0");
	grRelResid_4->Fit("pol0");

	//-----------------------
	const int n = 5;
	double sigma[n] = {2.5,3.0,3.5,3.8,4.0};
	double chi2[n], chi2rel[n];
	double p0[n], p0rel[n];

	chi2[4] = grResid_0->GetFunction("pol0")->GetChisquare()/(double)grResid_0->GetFunction("pol0")->GetNDF();
	chi2rel[4] = grRelResid_0->GetFunction("pol0")->GetChisquare()/(double)grRelResid_0->GetFunction("pol0")->GetNDF();
	p0[4] = grResid_0->GetFunction("pol0")->GetParameter(0);
	p0rel[4] = grRelResid_0->GetFunction("pol0")->GetParameter(0);

	chi2[3] = grResid_1->GetFunction("pol0")->GetChisquare()/(double)grResid_1->GetFunction("pol0")->GetNDF();
	chi2rel[3] = grRelResid_1->GetFunction("pol0")->GetChisquare()/(double)grRelResid_1->GetFunction("pol0")->GetNDF();
	p0[3] = grResid_1->GetFunction("pol0")->GetParameter(0);
	p0rel[3] = grRelResid_1->GetFunction("pol0")->GetParameter(0);

	chi2[2] = grResid_2->GetFunction("pol0")->GetChisquare()/(double)grResid_2->GetFunction("pol0")->GetNDF();
	chi2rel[2] = grRelResid_2->GetFunction("pol0")->GetChisquare()/(double)grRelResid_2->GetFunction("pol0")->GetNDF();
	p0[2] = grResid_2->GetFunction("pol0")->GetParameter(0);
	p0rel[2] = grRelResid_2->GetFunction("pol0")->GetParameter(0);

	chi2[1] = grResid_3->GetFunction("pol0")->GetChisquare()/(double)grResid_3->GetFunction("pol0")->GetNDF();
	chi2rel[1] = grRelResid_3->GetFunction("pol0")->GetChisquare()/(double)grRelResid_3->GetFunction("pol0")->GetNDF();
	p0[1] = grResid_3->GetFunction("pol0")->GetParameter(0);
	p0rel[1] = grRelResid_3->GetFunction("pol0")->GetParameter(0);

	chi2[0] = grResid_4->GetFunction("pol0")->GetChisquare()/(double)grResid_4->GetFunction("pol0")->GetNDF();
	chi2rel[0] = grRelResid_4->GetFunction("pol0")->GetChisquare()/(double)grRelResid_4->GetFunction("pol0")->GetNDF();
	p0[0] = grResid_4->GetFunction("pol0")->GetParameter(0);
	p0rel[0] = grRelResid_4->GetFunction("pol0")->GetParameter(0);

	TGraph *grChi2 = new TGraph(n,sigma,chi2);
	TGraph *grRelChi2 = new TGraph(n,sigma,chi2rel);
	TGraph *grP0 = new TGraphErrors(n,sigma,p0);
	TGraph *grRelP0 = new TGraphErrors(n,sigma,p0rel);

	grChi2->SetTitle(";#sigma;Chi^{2}/NDF");
	grRelChi2->SetTitle(";#sigma;Chi^{2}/NDF");
	grP0->SetTitle(";#sigma;#LTR_{4#sigma}- R_{i#sigma}#GT");
	grRelP0->SetTitle(";#sigma;#LTr_{4#sigma}- r_{i#sigma}#GT");

	//-----------------------
	TCanvas *cRate = new TCanvas("cRate","Rate Per Cell",1000,400);
	grRate_0->Draw("AP");
	grRate_1->Draw("P");
	grRate_2->Draw("P");
	grRate_3->Draw("P");
	grRate_4->Draw("P");
	TLegend *leg = new TLegend(0.9,0.72,0.99,0.99);
	leg->AddEntry(grRate_0,"4.0 #sigma");
	leg->AddEntry(grRate_1,"3.8 #sigma");
	leg->AddEntry(grRate_2,"3.5 #sigma");
	leg->AddEntry(grRate_3,"3.0 #sigma");
	leg->AddEntry(grRate_4,"2.5 #sigma");
	leg->Draw();	
	cRate->SaveAs(Form("%s/RateVsCell_%sCut.pdf",gSystem->Getenv("AD_AC227_PLOTS"),pdfVar.c_str()));

	TCanvas *chRate = new TCanvas("chRate","Rate Per Cell",1);
	hRelRate_0->Draw("HIST");
	hRelRate_1->Draw("HIST&SAME");
	hRelRate_2->Draw("HIST&SAME");
	hRelRate_3->Draw("HIST&SAME");
	hRelRate_4->Draw("HIST&SAME");
	leg = new TLegend(0.9,0.72,0.99,0.99);
	leg->AddEntry(hRelRate_0,"4.0 #sigma","l");
	leg->AddEntry(hRelRate_1,"3.8 #sigma","l");
	leg->AddEntry(hRelRate_2,"3.5 #sigma","l");
	leg->AddEntry(hRelRate_3,"3.0 #sigma","l");
	leg->AddEntry(hRelRate_4,"2.5 #sigma","l");
	leg->Draw();	
	chRate->SaveAs(Form("%s/RateVsCellHist_%sCut.pdf",gSystem->Getenv("AD_AC227_PLOTS"),pdfVar.c_str()));

	TCanvas *chResid = new TCanvas("chResid","Resid Per Cell",1);
	gPad->SetRightMargin(0.1);
	hResid_1->Draw("HIST");
	hResid_2->Draw("HIST&SAME");
	hResid_3->Draw("HIST&SAME");
	hResid_4->Draw("HIST&SAME");
	hRelResid_1->Draw("HIST&SAME");
	hRelResid_2->Draw("HIST&SAME");
	hRelResid_3->Draw("HIST&SAME");
	hRelResid_4->Draw("HIST&SAME");
	leg = new TLegend(0.9,0.72,0.99,0.99);
	leg->AddEntry(hResid_1,"3.8 #sigma","l");
	leg->AddEntry(hResid_2,"3.5 #sigma","l");
	leg->AddEntry(hResid_3,"3.0 #sigma","l");
	leg->AddEntry(hResid_4,"2.5 #sigma","l");
	leg->Draw();	
	chResid->SaveAs(Form("%s/ResidVsCell_%sCut.pdf",gSystem->Getenv("AD_AC227_PLOTS"),pdfVar.c_str()));

	TCanvas *cChi2 = new TCanvas("cChi2","chi2",1);
	//grChi2->Draw("AP");
	grRelChi2->Draw("AP");
	cChi2->SaveAs(Form("%s/Chi2_%sCut.pdf",gSystem->Getenv("AD_AC227_PLOTS"),pdfVar.c_str()));

	TCanvas *cP0 = new TCanvas("cP0","p0",1);
	//grP0->Draw("AP");
	grRelP0->Draw("AP");
	cP0->SaveAs(Form("%s/P0_%sCut.pdf",gSystem->Getenv("AD_AC227_PLOTS"),pdfVar.c_str()));

	TCanvas *c = new TCanvas("c","c",1);
	grEff_0->SetMarkerColor(1);
	grEff_1->SetMarkerColor(2);
	grEff_2->SetMarkerColor(6);
	grEff_3->SetMarkerColor(4);
	grEff_4->SetMarkerColor(8);

	grEff_0->Draw("AP");
	grEff_1->Draw("P");
	grEff_2->Draw("P");
	grEff_3->Draw("P");
	grEff_4->Draw("P");


//	TFile *f = new TFile(Form("%s/CutStudiesResults.root",gSystem->Getenv("AD_AC227ANALYSIS_RESULTS")),"UPDATE");
//	grRelP0->Write(Form("%s",pdfVar.c_str()));
//	f->Close();

	return 0;
}
