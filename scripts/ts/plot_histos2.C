#include <iostream>
using namespace std;



#define MULTIBIN 1.0

//0 delta
//1 eta1
//2 N1
//3 eta2
//4 N2
//5 translation2
//6 eta3
//7 N3
//8 translation3
//9 eta4
//10 N4
//11 translation4
Double_t traslchi2N1N2sumN3N4deltafuncV2(Double_t *x, Double_t *par) {
	
	Double_t f=0;
	Double_t trasl2 = par[5]; 
	Double_t x1 = x[0]-trasl2;
	Double_t trasl3 = par[8]; 
	Double_t x2 = x[0]-trasl3;
	Double_t trasl4 = par[11]; 
	
	Double_t f = (x[0]<=trasl2) ? 
	par[0]*(x[0]<=1.0/MULTIBIN?1:0)   + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) :  
	(x[0]>trasl2 && x[0] <= trasl3) ? 
	par[3] * exp(-x1/2.0) * pow(x1, par[4]/2.0 - 1) / ( pow(2, par[4]/2.0 ) * TMath::Gamma(par[4]/2.0) ) : 
	(x[0]>trasl3 && x[0] <= trasl4) ?
	par[6] * exp(-x[0]/2.0) * pow(x[0], par[7]/2.0 - 1) / ( pow(2, par[7]/2.0 ) * TMath::Gamma(par[7]/2.0) ) :
	par[9] * exp(-x[0]/2.0) * pow(x[0], par[10]/2.0 - 1) / ( pow(2, par[10]/2.0 ) * TMath::Gamma(par[10]/2.0) );
	
	
	
	
	//((x[0]>par[2]/MULTIBIN)?1:0)
	return f;
	
}


//filetype = 0, without fitting
//filetype = 1, with fitting

//type = 0, PDF
//type = 1, PVALUE
//type = 2, CDF

void plot_histos2(int filetype = 0, TString type = "PDF", TString filenameinput1 = "", int style1 = 1, int color1 = kBlack, int width1 = 1, TString filenameinput2 = "", int style2 = 1, int color2 = kBlack, int width2 = 1, TString filenameinput3 = "", int style3 = 1, int color3 = kBlack, int width3 = 1, TString filenameinput4 = "", int style4 = 1, int color4 = kBlack, int width4 = 1, TString filenameinput5 = "", int style5 = 1, int color5 = kBlack, int width5 = 1, TString filenameinput6 = "", int style6 = 1, int color6 = kBlack, int width6 = 1, TString filenameinput7 = "", int style7 = 1, int color7 = kBlack, int width7 = 1, TString filenameinput8 = "", int style8 = 1, int color8 = kBlack, int width8 = 1) {
	Double_t eta = 0.5;
	Double_t cumTSdisbinx, max_cumTSdisbinx ;
	max_cumTSdisbinx = 36;
	cumTSdisbinx = max_cumTSdisbinx*MULTIBIN;

	//FITTING
	int fittingfunction = 0;
	//M95+ICL=9**(23)*********************************************************************************************************************
	//M95*******************************************************************************************************************************
	//M95*******************************************************************************************************************************
	TF1* f4sf_traslN1N2sumlN3N4v2 = new TF1("traslchi2N1N2sumN3N4deltafuncV2", traslchi2N1N2sumN3N4deltafuncV2, 0, max_cumTSdisbinx,12);
	f4sf_traslN1N2sumlN3N4v2->SetParName(0, "delta");
	f4sf_traslN1N2sumlN3N4v2->SetParName(1, "eta1");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(1, 0, 1);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 1, 0.5);
	f4sf_traslN1N2sumlN3N4v2->SetParName(2, "N1");
	Int_t  f4sf_test_N1=1;
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(2, 0, 5);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 2, f4sf_test_N1);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //!!!!!!!!!!!!
	f4sf_traslN1N2sumlN3N4v2->SetParName(3, "eta2");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(3, 0, 1);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 3, 0.5);
	f4sf_traslN1N2sumlN3N4v2->SetParName(4, "N2");
	Int_t  f4sf_test_N2=5;
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(4, 0, 5);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 4, f4sf_test_N2);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_traslN1N2sumlN3N4v2->SetParName(5, "translation2");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(5, 0, 36);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	f4sf_traslN1N2sumlN3N4v2->SetParName(6, "eta3");
	double eta3 = 1.43127e-02;
	//f4sf_traslN1N2sumlN3N4v2->SetParLimits(6, eta3, eta3);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 6, eta3);
	f4sf_traslN1N2sumlN3N4v2->SetParName(7, "N3");
	Int_t  f4sf_test_N3=5;
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(7, 1, 10);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(7, f4sf_test_N3);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(7, f4sf_test_N3, f4sf_test_N3); //!!!!!!!!!!!!
	Int_t  f4sf_trasl_N1N2v2_trasl = 9;
	f4sf_traslN1N2sumlN3N4v2->SetParName(8, "translation3");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(8, 0, 36);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	
	double eta4 = 3;
	f4sf_traslN1N2sumlN3N4v2->SetParName(9, "eta4");
	//f4sf_traslN1N2sumlN3N4v2->SetParLimits(9, eta4, eta4);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 9, eta4);
	f4sf_traslN1N2sumlN3N4v2->SetParName(10, "N4");
	Int_t  f4sf_test_N4=1;
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(10, 0, 5);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(10, f4sf_test_N4);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(10, f4sf_test_N4, f4sf_test_N4); //!!!!!!!!!!!!
	Int_t  f4sf_trasl_N1N2v2_trasl = 14;
	f4sf_traslN1N2sumlN3N4v2->SetParName(11, "translation3");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(11, 0, 36);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(11, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(11, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	
	f4sf_traslN1N2sumlN3N4v2->SetLineColor(kBlack);
		
	//ENDFITTING
	
		gStyle->SetPalette(1);gStyle->SetPadColor(0);gStyle->SetFrameFillColor(0);gStyle->SetCanvasColor(0);gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(0000000);gStyle->SetOptFit(0000);gStyle->SetPalette(1);gStyle->SetLabelFont(42, "xyz");gStyle->SetTitleFont(42, "xyz");
	/*
	TCanvas c1("Distribution", "Distribution", 1280,1024); c1.Divide(2,1);c1.cd(1);gPad->SetLogy();c1.cd(2);gPad->SetLogy();
	*/
	TCanvas* c1 = new TCanvas("Distribution", "Distribution", 1280,1024);
	Double_t dimWindowX = 1400;
	Double_t dimWindowY = 900;
	
	TString filenameinput[8];
	int style[8];
	int color[8];
	int width[8];
	TH1D* h1[8];
	TH1D* hf[8];
	
	
	TString nameh, namef;
	filenameinput[0] = filenameinput1;
	style[0] = style1;
	color[0] = color1;
	width[0] = width1;
	nameh = "h0"; nameh += type;
	namef = "f0"; namef += type;
	h1[0] = new TH1D(nameh, nameh, cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[0] = new TH1D(namef, namef, cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[0]->SetLineWidth(width[0]);
	hf[0]->SetLineWidth(2);
	h1[0]->SetLineColor(color[0]);
	hf[0]->SetLineColor(kBlack);
	h1[0]->SetLineStyle(style[0]);
	
	if(type == "PVALUE") {
		h1[0]->GetXaxis()->SetTitle("h");
		h1[0]->GetYaxis()->SetTitle("P(TS >= h)");
	} 
	if(type == "PDF") {
		h1[0]->GetXaxis()->SetTitle("TS");
		h1[0]->GetYaxis()->SetTitle("counts (normalized)");	
	}
	if(type == "CDF") {
		h1[0]->GetXaxis()->SetTitle("h");
		h1[0]->GetYaxis()->SetTitle("P(TS <= h)");	
	}
	h1[0]->GetYaxis()->SetTitleSize(0.03);
	h1[0]->GetYaxis()->SetTitleOffset(1.60);
	
	
	filenameinput[1] = filenameinput2;
	style[1] = style2;
	color[1] = color2;
	width[1] = width2;
	nameh = "h1"; nameh += type;
	namef = "f1"; namef += type;
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
	nameh = "h2"; nameh += type;
	namef = "f2"; namef += type;
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
	nameh = "h3"; nameh += type;
	namef = "f3"; namef += type;
	h1[3] = new TH1D(nameh, nameh,  cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[3] = new TH1D(namef, namef,  cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[3]->SetLineWidth(width[3]);
	hf[3]->SetLineWidth(2);
	h1[3]->SetLineColor(color[3]);
	hf[3]->SetLineColor(kBlack);
	h1[3]->SetLineStyle(style[3]);
	
	filenameinput[4] = filenameinput5;
	style[4] = style5;
	color[4] = color5;
	width[4] = width5;
	nameh = "h4"; nameh += type;
	namef = "f4"; namef += type;
	h1[4] = new TH1D(nameh, nameh,  cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[4] = new TH1D(namef, namef,  cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[4]->SetLineWidth(width[4]);
	hf[4]->SetLineWidth(2);
	h1[4]->SetLineColor(color[4]);
	hf[4]->SetLineColor(kBlack);
	h1[4]->SetLineStyle(style[4]);
	
	filenameinput[5] = filenameinput6;
	style[5] = style6;
	color[5] = color6;
	width[5] = width6;
	nameh = "h5"; nameh += type;
	namef = "f5"; namef += type;
	h1[5] = new TH1D(nameh, nameh,  cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[5] = new TH1D(namef, namef,  cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[5]->SetLineWidth(width[5]);
	hf[5]->SetLineWidth(2);
	h1[5]->SetLineColor(color[5]);
	hf[5]->SetLineColor(kBlack);
	h1[5]->SetLineStyle(style[5]);
	

	filenameinput[6] = filenameinput7;
	style[6] = style7;
	color[6] = color7;
	width[6] = width7;
	nameh = "h6"; nameh += type;
	namef = "f6"; namef += type;
	h1[6] = new TH1D(nameh, nameh,  cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[6] = new TH1D(namef, namef,  cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[6]->SetLineWidth(width[6]);
	hf[6]->SetLineWidth(2);
	h1[6]->SetLineColor(color[6]);
	hf[6]->SetLineColor(kBlack);
	h1[6]->SetLineStyle(style[6]);
	
	filenameinput[7] = filenameinput8;
	style[7] = style8;
	color[7] = color8;
	width[7] = width8;
	nameh = "h6"; nameh += type;
	namef = "f6"; namef += type;
	h1[7] = new TH1D(nameh, nameh,  cumTSdisbinx, 0, max_cumTSdisbinx );
	hf[7] = new TH1D(namef, namef,  cumTSdisbinx, 0, max_cumTSdisbinx );
	h1[7]->SetLineWidth(width[7]);
	hf[7]->SetLineWidth(2);
	h1[7]->SetLineColor(color[7]);
	hf[7]->SetLineColor(kBlack);
	h1[7]->SetLineStyle(style[7]);
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
	for (int i=0; i<8; i++) {
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
			if (type == "PDF") {
				D = DPDF;
				E = EPDF;
			} 
			if (type == "CDF") {
				D = DCDF;
				E = ECDF;
			}
			if (type == "PVALUE") {
				D = DPVALUE;
				E = EPVALUE;
			}
			BIN = BIN+1;
			h1[i]->SetBinContent(BIN, D);
			h1cl->SetBinContent(BIN, D);
			h1[i]->SetBinError(BIN, E);
			cout << BIN << " " << D << " " << E;
			if (filetype==1) {
				if (type == "PDF") {
					F = FPDF;
				}
				if (type == "PVALUE") {
					F = FPVALUE;
				}
				cout << " " << F;
				hf[i]->SetBinContent(BIN, F);
			}
			cout << endl;
		}
		TString tit;
		tit += type;
		tit += " ";
		tit += filenameinput[i];
		h1[i]->SetTitle(tit);
		
		if(i==0) {
			
			switch(fittingfunction) {
					
				case 23:
					cout << "####################### traslchi2N1N2sumN3N4deltafuncV2 (M95 with ICL) " << endl;
					h1[i]->Fit("traslchi2N1N2sumN3N4deltafuncV2", "+");
					bestfit = f4sf_traslN1N2sumlN3N4v2;
					double params[12];
					params[0] = bestfit->GetParameter(0);
					params[1] = bestfit->GetParameter(1);
					params[2] = bestfit->GetParameter(2);
					params[3] = bestfit->GetParameter(3);
					params[4] = bestfit->GetParameter(4);
					params[5] = bestfit->GetParameter(5);
					params[6] = bestfit->GetParameter(6);
					params[7] = bestfit->GetParameter(7);
					params[8] = bestfit->GetParameter(8);
					params[9] = bestfit->GetParameter(9);
					params[10] = bestfit->GetParameter(10);
					params[11] = bestfit->GetParameter(11);
					break;
				default:
					break;
			}
			
			
			
		}
		
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

	if(type == "PVALUE") {
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
	if(type == "PDF") {

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