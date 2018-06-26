#include <iostream>
using namespace std;



#define MULTIBIN 1.0

//filetype = 0, without fitting
//filetype = 1, with fitting

//type = 0, PDF
//type = 1, PVALUE
//type = 2, CDF
//.x ~/grid_scripts2/ts/plot_histos3.C (0, "PVALUE","", 1, kBlue, 2, "PVALUE", "", 1, kRed, 2, "PVALUE", "", 1, kBlack, 2, "PVALUE", "", 1, kBlack, 2)

void plot_histos3(int filetype = 0, TString type1 = "PDF", TString filenameinput1 = "", int style1 = 1, int color1 = kBlack, int width1 = 1, TString type2 = "PDF", TString filenameinput2 = "", int style2 = 1, int color2 = kBlack, int width2 = 1, TString type3 = "PDF", TString filenameinput3 = "", int style3 = 1, int color3 = kBlack, int width3 = 1, TString type4 = "PDF", TString filenameinput4 = "", int style4 = 1, int color4 = kBlack, int width4 = 1) {
	
	
	
	gStyle->SetPalette(1);gStyle->SetPadColor(0);gStyle->SetFrameFillColor(0);gStyle->SetCanvasColor(0);gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(0000000);gStyle->SetOptFit(0000);gStyle->SetPalette(1);gStyle->SetLabelFont(42, "xyz");gStyle->SetTitleFont(42, "xyz");
	/*
	TCanvas c1("Distribution", "Distribution", 1280,1024); c1.Divide(2,1);c1.cd(1);gPad->SetLogy();c1.cd(2);gPad->SetLogy();
	*/
	Double_t dimWindowX = 1400;
	Double_t dimWindowY = 900;
	
	TString filenameinput[4];
	int style[4];
	int color[4];
	int width[4];
	TString type[4];
	TH1D* h1[4];
	TH1D* hf[4];
	
	Double_t eta = 0.5;
	Double_t cumTSdisbinx, max_cumTSdisbinx ;
	max_cumTSdisbinx = 36;
	cumTSdisbinx = max_cumTSdisbinx*MULTIBIN;
	
	TString nameh, namef;
	filenameinput[0] = filenameinput1;
	style[0] = style1;
	color[0] = color1;
	width[0] = width1;
	type[0] = type1;
	nameh = "h0"; nameh += type[0];
	namef = "f0"; namef += type[0];
	h1[0] = new TH1D(nameh, nameh, cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[0] = new TH1D(namef, namef, cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[0]->SetLineWidth(width[0]);
	hf[0]->SetLineWidth(2);
	h1[0]->SetLineColor(color[0]);
	hf[0]->SetLineColor(kBlack);
	h1[0]->SetLineStyle(style[0]);
	
	if(type[0] == "PVALUE") {
		h1[0]->GetXaxis()->SetTitle("h");
		h1[0]->GetYaxis()->SetTitle("P(TS >= h)");
	} 
	if(type[0] == "PDF") {
		h1[0]->GetXaxis()->SetTitle("TS");
		h1[0]->GetYaxis()->SetTitle("counts (normalized)");	
	}
	if(type[0] == "CDF") {
		h1[0]->GetXaxis()->SetTitle("h");
		h1[0]->GetYaxis()->SetTitle("P(TS <= h)");	
	}
	h1[0]->GetYaxis()->SetTitleSize(0.03);
	h1[0]->GetYaxis()->SetTitleOffset(1.60);
	
	
	filenameinput[1] = filenameinput2;
	style[1] = style2;
	color[1] = color2;
	width[1] = width2;
	type[1] = type2;
	nameh = "h1"; nameh += type[1];
	namef = "f1"; namef += type[1];
	h1[1] = new TH1D(nameh, nameh, cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[1] = new TH1D(namef, namef, cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[1]->SetLineWidth(width[1]);
	hf[1]->SetLineWidth(2);
	h1[1]->SetLineColor(color[1]);
	hf[1]->SetLineColor(kBlack);
	h1[1]->SetLineStyle(style[1]);
	
	filenameinput[2] = filenameinput3;
	style[2] = style3;
	color[2] = color3;
	width[2] = width3;
	type[2] = type3;
	nameh = "h2"; nameh += type[2];
	namef = "f2"; namef += type[2];
	h1[2] = new TH1D(nameh, nameh, cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[2] = new TH1D(namef, namef, cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[2]->SetLineWidth(width[2]);
	hf[2]->SetLineWidth(2);
	h1[2]->SetLineColor(color[2]);
	hf[2]->SetLineColor(kBlack);
	h1[2]->SetLineStyle(style[2]);
	
	filenameinput[3] = filenameinput4;
	style[3] = style4;
	color[3] = color4;
	width[3] = width4;
	type[3] = type4;
	nameh = "h3"; nameh += type[3];
	namef = "f3"; namef += type[3];
	h1[3] = new TH1D(nameh, nameh,  cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[3] = new TH1D(namef, namef,  cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[3]->SetLineWidth(width[3]);
	hf[3]->SetLineWidth(2);
	h1[3]->SetLineColor(color[3]);
	hf[3]->SetLineColor(kBlack);
	h1[3]->SetLineStyle(style[3]);
	

	cout << "Number of bins " << cumTSdisbinx << " from 0 to " << max_cumTSdisbinx << endl;
	
	TH1D* chi2_1 = new TH1D("\\chi^{2}_{1}", "\\chi^{2}_{1}", cumTSdisbinx , 0, max_cumTSdisbinx );
	TString histname;
	histname += "1/2 \\chi^{2}_{1}";
	TH1D* chi2_1_05 = new TH1D(histname, histname, cumTSdisbinx , 0, max_cumTSdisbinx );
	TString histname33;
	histname33 += "1/2 \\chi^{2}_{3}";
	TH1D* chi2_3_05 = new TH1D(histname33, histname33, cumTSdisbinx , 0, max_cumTSdisbinx );
	
	TH1D* integral_chi2_1 = new TH1D("\\int \\chi^{2}_{1}", "int \\chi^{2}_{1}", cumTSdisbinx , 0, max_cumTSdisbinx );
	TString histname2, histname3;
	histname2 += "1/2";
	histname3 += "1/2";
	histname2 += " \\int \\chi^{2}_{1}";
	histname3 += " \\int \\chi^{2}_{3}";
	TH1D* integral_chi2_1_05 = new TH1D(histname2,histname2, cumTSdisbinx , 0, max_cumTSdisbinx );
	TH1D* integral_chi2_3_05 = new TH1D(histname3,histname3, cumTSdisbinx , 0, max_cumTSdisbinx );
	

	
	Int_t counts = 0;


	Long64_t nlines;
	TTree* T;
	
	//TCanvas* c1_h1 = new TCanvas("Distribution", "Distribution", dimWindowX, dimWindowY); gPad->SetLogy(); 	gPad->SetTitle();

	
	Float_t BIN, BINCENTER, DPDF, EPDF, DCDF, ECDF, D, E, FPDF, FPVALUE, F, DPVALUE, EPVALUE;
	for (int i=0; i<4; i++) {
		if (filenameinput[i] == "") 
			continue;
	
		T = new TTree("DATA", "") ;
		
		if(filetype == 0) {
			nlines = T->ReadFile(filenameinput[i], "BIN:DPDF:EPDF:DPVALUE:EPVALUE:DCDF:ECDF");
		
			T->SetBranchAddress("BIN", &BIN);
			T->SetBranchAddress("DPDF", &DPDF);
			T->SetBranchAddress("EPDF", &EPDF);
			T->SetBranchAddress("DPVALUE", &DPVALUE);
			T->SetBranchAddress("EPVALUE", &EPVALUE);
			T->SetBranchAddress("DCDF", &DCDF);
			T->SetBranchAddress("ECDF", &ECDF);
		
		} else {
			nlines = T->ReadFile(filenameinput[i], "BIN:DPDF:EPDF:DPVALUE:EPVALUE:DCDF:ECDF:FPDF:FPVALUE");
			
			T->SetBranchAddress("BIN", &BIN);
			T->SetBranchAddress("DPDF", &DPDF);
			T->SetBranchAddress("EPDF", &EPDF);
			T->SetBranchAddress("DPVALUE", &DPVALUE);
			T->SetBranchAddress("EPVALUE", &EPVALUE);
			T->SetBranchAddress("DCDF", &DCDF);
			T->SetBranchAddress("ECDF", &ECDF);
			T->SetBranchAddress("FPDF", &FPDF);
			T->SetBranchAddress("FPVALUE", &FPVALUE);
		}
		cout << "N lines: " << nlines << endl;
		
		if(nlines == 0) {
			cout << "End of the procedure" << endl;
			return;
		}
		
		TString clname;
		clname += (int)(gRandom->Uniform()*1000);
		TH1F* h1cl = h1[i]->Clone(clname);
		h1cl->SetName(clname);
		h1cl->SetTitle(clname);

		for(Long64_t j = 0; j<nlines; j++) {
			T->GetEntry(j);
			if (type[i] == "PDF") {
				D = DPDF;
				E = EPDF;
			} 
			if (type[i] == "CDF") {
				D = DCDF;
				E = ECDF;
			}
			if (type[i] == "PVALUE") {
				D = DPVALUE;
				E = EPVALUE;
			}
			BIN = BIN+1;
			h1[i]->SetBinContent(BIN, D);
			h1cl->SetBinContent(BIN, D);
			h1[i]->SetBinError(BIN, E);
			cout << BIN << " " << D << " " << E;
			if (filetype==1) {
				if (type[i] == "PDF") {
					F = FPDF;
				}
				if (type[1] == "PVALUE") {
					F = FPVALUE;
				}
				cout << " " << F;
				hf[i]->SetBinContent(BIN, F);
			}
			cout << endl;
		}
		TString tit;
		tit += type[i];
		tit += " ";
		tit += filenameinput[i];
		h1[i]->SetTitle(tit);
		
		if(i==0) {
			h1[i]->Draw("E");
			h1cl->Draw("SAMEL");
			hf[i]->Draw("SAMEL");
		} else {
			h1[i]->Draw("SAMEE");
			h1cl->Draw("SAMEL");
			hf[i]->Draw("SAMEL");
		}
	}	
	
	
	
	//fitting
	
	//fitting2
	Double_t dimbin = max_cumTSdisbinx / (double) cumTSdisbinx;
	TF1* f1sf = new TF1("fchi2", "[0] * exp(-x/2.0) / (sqrt(2*pi * x)) ", dimbin, max_cumTSdisbinx);
	TF1* f2sfN= new TF1("fchi2N", " [0] * exp(-x/2.0) * pow(x, [1]/2.0 - 1) / ( pow(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", dimbin, max_cumTSdisbinx);
	TF1* f2sf1= new TF1("f2sf1", " [0] * exp(-x/2.0) * pow(x, [1]/2.0 - 1) / ( pow(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", dimbin, max_cumTSdisbinx);
	TF1* f2sf1_05= new TF1("f2sf1_05", " [0] * exp(-x/2.0) * pow(x, [1]/2.0 - 1) / ( pow(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", dimbin, max_cumTSdisbinx);
	TF1* f2sf3_05= new TF1("f2sf3_05", " [0] * exp(-x/2.0) * pow(x, [1]/2.0 - 1) / ( pow(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", dimbin, max_cumTSdisbinx);
	TF1* f2sf = new TF1("delta", "[0]", 0, dimbin);
	TF1* f3sf = new TF1("sumchi2deltaNO", " delta + fchi2", 0, max_cumTSdisbinx);
	TF1* f4sf = new TF1("sumchi2NdeltaNO", "delta + fchi2N", 0, max_cumTSdisbinx);
	
	//reference
	f2sf3_05->SetParLimits(0, 0.5, 0.5);
	f2sf3_05->SetParameter(0, 0.5);
	f2sf3_05->SetParLimits(1, 3, 3);
	f2sf3_05->SetParameter(1, 3);
	
	//reference
	f2sf1_05->SetParLimits(0, 0.5, 0.5);
	f2sf1_05->SetParameter(0, 0.5);
	f2sf1_05->SetParLimits(1, 1, 1);
	f2sf1_05->SetParameter(1, 1);
	//reference
	f2sf1->SetParLimits(0, 1, 1);
	f2sf1->SetParameter(0, 1);
	f2sf1->SetParLimits(1, 1, 1);
	f2sf1->SetParameter(1, 1);
	
	
	
	Double_t evalfirst1 = -1, evalfirst2 = -1;
	
	
	for(int i=1; i<cumTSdisbinx; i++) {
		Double_t eval, eval2;
		Double_t x;
		//histogram of chi2/1
		x = chi2_1->GetBinCenter(i);
		eval = f2sf1->Eval(x);
		if(evalfirst1 == -1)
			evalfirst1 = eval;
		chi2_1->SetBinContent(i, eval);
		
		//histogram of chi2_1/2
		eval2 = f2sf1_05->Eval(x);
		if(evalfirst2 == -1)
			evalfirst2 = eval2;
		chi2_1_05->SetBinContent(i, eval2);
		
		//histogram of chi2_3/2
		eval2 = f2sf3_05->Eval(x);
		if(evalfirst2 == -1)
			evalfirst2 = eval2;
		chi2_3_05->SetBinContent(i, eval2);
	}
	for(int i=0; i<cumTSdisbinx; i++) {
		Double_t integral, integral2, integral3;
		
		int jjj=i+1;
		
		//histogram of integral of chi2_1
		integral = f2sf1->Integral(i, cumTSdisbinx);
		integral_chi2_1->SetBinContent(jjj, integral);
		//histogram of integral of chi2_1/2
		integral2 = f2sf1_05->Integral(i, cumTSdisbinx);
		integral_chi2_1_05->SetBinContent(jjj, integral2);
		//histogram of integral of chi2_3/2
		integral2 = f2sf3_05->Integral(i, cumTSdisbinx);
		integral_chi2_3_05->SetBinContent(jjj, integral2);
		
	}	

	if(type[0] == "PVALUE") {
		integral_chi2_1->SetLineColor(kGreen+2);
		integral_chi2_1->SetLineWidth(1);
		integral_chi2_1->SetLineStyle(2);
		integral_chi2_1->Draw("SAMEL");
		
		integral_chi2_1_05->SetLineColor(kRed);
		integral_chi2_1_05->SetLineWidth(1);
		integral_chi2_1_05->SetLineStyle(3);
		integral_chi2_1_05->Draw("SAMEL");
		
		
		integral_chi2_3_05->SetLineColor(kCyan+1);
		integral_chi2_3_05->SetLineWidth(1);
		integral_chi2_3_05->SetLineStyle(5);
		integral_chi2_3_05->Draw("SAMEL");
	} 
	if(type[0] == "PDF") {

		chi2_1->SetLineColor(kGreen+2);
		chi2_1->SetLineWidth(1);
		chi2_1->SetLineStyle(2);
		chi2_1->Draw("SAMEL");
		chi2_1_05->SetLineColor(kRed);
		chi2_1_05->SetLineWidth(1);
		chi2_1_05->SetLineStyle(3);
		chi2_1_05->Draw("SAMEL");
		chi2_3_05->SetLineColor(kCyan+3);
		chi2_3_05->SetLineWidth(1);
		chi2_3_05->SetLineStyle(5);
		chi2_3_05->Draw("SAMEL");
	}
}