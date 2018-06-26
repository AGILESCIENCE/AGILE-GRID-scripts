#include <iostream>
using namespace std;



#define MULTIBIN 1.0



//funzione per fare il fitting
Double_t sumchi2deltafunc(Double_t *x, Double_t *par)
{
	Double_t f = par[0]*(x[0]<=1.0/MULTIBIN?1:0) + par[1]*exp(-x[0]/2.0) / (sqrt(2*TMath::Pi() * x[0]));
	return f;
}

Double_t sumchi2Ndeltafunc(Double_t *x, Double_t *par)
{
	Double_t f = par[0]*(x[0]<=1.0/MULTIBIN?1:0)   + par[1] * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) ;
	return f;
}

Double_t sumchi2Ndeltafunc_V2(Double_t *x, Double_t *par)
{
	Double_t f = par[0]*((x[0]<=1.0/MULTIBIN)?1:0) + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) ;
	return f;
}

Double_t sumchi2test(Double_t *x, Double_t *par)
{
	Double_t f=0;
	Double_t trasl = 6; 
	Double_t x1 = x[0]-trasl;

	Double_t f = (x[0]<=trasl) ? par[0]*(x[0]<=1.0/MULTIBIN?1:0)   + par[1] * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) : par[3] * exp(-x1/2.0) * pow(x1, par[4]/2.0 - 1) / ( pow(2, par[4]/2.0 ) * TMath::Gamma(par[4]/2.0) );

	//Double_t f2 = (x[0]>trasl) ? par[3] * exp(-x1/2.0) * pow(x1, par[4]/2.0 - 1) / ( pow(2, par[4]/2.0 ) * TMath::Gamma(par[4]/2.0) ) : 0 ;

	
	//((x[0]>par[2]/MULTIBIN)?1:0)
	return f;
}

//0 delta
//1 eta1
//2 N1
//3 eta2
//4 N2
//5 translation2
Double_t traslchi2N1N2deltafuncV2(Double_t *x, Double_t *par) {
	
	Double_t f=0;
	Double_t trasl = par[5]; 
	Double_t x1 = x[0]-trasl;
	
	Double_t f = (x[0]<=trasl) ? par[0]*(x[0]<=1.0/MULTIBIN?1:0)   + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) : par[3] * exp(-x1/2.0) * pow(x1, par[4]/2.0 - 1) / ( pow(2, par[4]/2.0 ) * TMath::Gamma(par[4]/2.0) );
	
	
	
	//((x[0]>par[2]/MULTIBIN)?1:0)
	return f;
	
}

//0 delta
//1 eta1
//2 N1
//3 eta2
//4 N2
//5 translation2
//6 eta3
//7 N3
//8 translation3
Double_t traslchi2N1N2N3deltafuncV2(Double_t *x, Double_t *par) {
	
	Double_t f=0;
	Double_t trasl2 = par[5]; 
	Double_t x1 = x[0]-trasl2;
	Double_t trasl3 = par[8]; 
	Double_t x2 = x[0]-trasl3;
	
	Double_t f = (x[0]<=trasl2) ? par[0]*(x[0]<=1.0/MULTIBIN?1:0)   + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) :  (x[0]>trasl2 && x[0] <= trasl3) ? par[3] * exp(-x1/2.0) * pow(x1, par[4]/2.0 - 1) / ( pow(2, par[4]/2.0 ) * TMath::Gamma(par[4]/2.0) ) : par[6] * exp(-x2/2.0) * pow(x2, par[7]/2.0 - 1) / ( pow(2, par[7]/2.0 ) * TMath::Gamma(par[7]/2.0) );
	
	
	
	//((x[0]>par[2]/MULTIBIN)?1:0)
	return f;
	
}

//0 delta
//1 eta1
//2 N1
//3 eta2
//4 N2
//5 translation2
//6 eta3
//7 N3
//8 translation3
Double_t traslchi2N1N2sumN3deltafuncV2(Double_t *x, Double_t *par) {
	
	Double_t f=0;
	Double_t trasl2 = par[5]; 
	Double_t x1 = x[0]-trasl2;
	Double_t trasl3 = par[8]; 
	Double_t x2 = x[0]-trasl3;
	
	Double_t f = (x[0]<=trasl2) ? 
	par[0]*(x[0]<=1.0/MULTIBIN?1:0)   + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) :  
	(x[0]>trasl2 && x[0] <= trasl3) ? 
	par[3] * exp(-x1/2.0) * pow(x1, par[4]/2.0 - 1) / ( pow(2, par[4]/2.0 ) * TMath::Gamma(par[4]/2.0) ) : 
	par[6] * exp(-x[0]/2.0) * pow(x[0], par[7]/2.0 - 1) / ( pow(2, par[7]/2.0 ) * TMath::Gamma(par[7]/2.0) );
	
	
	
	//((x[0]>par[2]/MULTIBIN)?1:0)
	return f;
	
}


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


//0 delta
//1 eta1
//2 N1
//3 eta2
//4 N2
//5 translation2 (T1=t_icl)
//6 eta3
//7 N3
//8 translation3 (T2)
//9 eta4
//10 N4
//11 translation4 (T3)
Double_t traslchi2N1N2sumN3N4deltafuncV3(Double_t *x, Double_t *par) {
	
	Double_t f=0;
	Double_t trasl2 = par[5]; 
	Double_t x1 = x[0]-trasl2;
	Double_t trasl3 = par[8];
	Double_t x2 = x[0];
	Double_t trasl4 = par[11]; 
	Double_t x3 = x[0]; 
	
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

//0 delta
//1 eta1
//2 N1
//3 eta2
//4 N2
//5 translation2
//6 eta3
//7 N3
//8 translation3
Double_t chi2delta_sumN1taslN2sumN3funcV2(Double_t *x, Double_t *par) {
	
	Double_t f=0;
	Double_t trasl2 = par[5]; 
	Double_t x1 = x[0] -trasl2;
	Double_t trasl3 = par[8]; 
	Double_t x2 = x[0];
	
	Double_t f = (x[0]<=trasl2) ? par[0]*(x[0]<=1.0/MULTIBIN?1:0)   + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) :  (x[0]>trasl2 && x[0] <= trasl3) ? par[3] * exp(-x1/2.0) * pow(x1, par[4]/2.0 - 1) / ( pow(2, par[4]/2.0 ) * TMath::Gamma(par[4]/2.0) ) : par[6] * exp(-x2/2.0) * pow(x2, par[7]/2.0 - 1) / ( pow(2, par[7]/2.0 ) * TMath::Gamma(par[7]/2.0) );
	
	
	
	//((x[0]>par[2]/MULTIBIN)?1:0)
	return f;
	
}


//delta + two-component fitting function
Double_t sumchi2N1N2deltafunc(Double_t *x, Double_t *par)
{
	
	/*
	 par[0] = normalization of delta function
	 par[1] = locCL TS threshold
	 par[2] = normalization of below threshold chi^2
	 par[3] = dof of below threshold chi^2
	 par[4] = normalization of above threshold chi^2
	 par[5] = dof of above threshold chi^2
	 */
	Double_t firstbin = 1.0/MULTIBIN;
	Double_t f;
	
	f = par[0] * ((x[0]<=firstbin)?1:0); 
	f += ((x[0]<=par[1]/MULTIBIN)?1:0)? (par[2] * exp(-x[0]/2.0) * pow(x[0], par[3]/2.0 - 1) / ( pow(2, par[3]/2.0 ) * TMath::Gamma(par[3]/2.0) ) ) : (par[4] * exp(-x[0]/2.0) * pow(x[0], par[5]/2.0 - 1) / ( pow(2, par[5]/2.0 ) * TMath::Gamma(par[5]/2.0) ) ) ;

	return f;
}

Double_t sumchi2N1N2deltafuncV2(Double_t *x, Double_t *par)
{
	
	/*
	 par[0] = normalization of delta function
	 par[1] = locCL TS threshold
	 par[2] = normalization of below threshold chi^2
	 par[3] = dof of below threshold chi^2
	 par[4] = normalization of above threshold chi^2
	 par[5] = dof of above threshold chi^2
	 */
	Double_t firstbin = 1.0/MULTIBIN;
	Double_t x1;
	x1 = x[0]-firstbin;
	Double_t x2;
	x2 = x[0]-par[1];
	Double_t f;
	x1 = x[0];
	x2 = x[0];
	
	f = par[0] * ((x[0]<=firstbin)?1:0); 
	
	f += ((x[0]>firstbin && x[0] <= par[1]/MULTIBIN)?1:0) * (par[2] * exp(-x1/2.0) * pow(x1, par[3]/2.0 - 1) / ( pow(2, par[3]/2.0 ) * TMath::Gamma(par[3]/2.0) ) ); 
	f += ((x[0] >  par[1]/MULTIBIN)?1:0) * (par[4] * exp(-x2/2.0) * pow(x2, par[5]/2.0 - 1) / ( pow(2, par[5]/2.0 ) * TMath::Gamma(par[5]/2.0) ) );
	
	
	
	return f;
}

Double_t sumchi2N1N2N3deltafuncV2(Double_t *x, Double_t *par)
{
	
	/*
	 par[0] = normalization of delta function
	 par[1] = threshold 1
	 par[2] = normalization of below threshold chi^2
	 par[3] = dof of below threshold chi^2
	 par[4] = normalization of above threshold chi^2
	 par[5] = dof of above threshold chi^2
	 par[6] = threshodl 2
	 par[7] = normalization of above threshold chi^2
	 par[8] = dof of above threshold chi^2
	 */
	Double_t firstbin = 1.0/MULTIBIN;
	Double_t x1;
	x1 = x[0]-firstbin;
	Double_t x2, x3;
	x2 = x[0]-par[1];
	Double_t f;
	x1 = x[0];
	x2 = x[0];
	x3 = x[0];
	
	f = par[0] * ((x[0]<=firstbin)?1:0); 
	
	f += ((x[0]>firstbin && x[0] <= par[1]/MULTIBIN)?1:0) * (par[2] * exp(-x1/2.0) * pow(x1, par[3]/2.0 - 1) / ( pow(2, par[3]/2.0 ) * TMath::Gamma(par[3]/2.0) ) ); 
	f += ((x[0] >  par[1]/MULTIBIN && x[0] <=  par[6]/MULTIBIN)?1:0) * (par[4] * exp(-x2/2.0) * pow(x2, par[5]/2.0 - 1) / ( pow(2, par[5]/2.0 ) * TMath::Gamma(par[5]/2.0) ) );
	f += ((x[0] >  par[6]/MULTIBIN)?1:0) * (par[7] * exp(-x3/2.0) * pow(x3, par[8]/2.0 - 1) / ( pow(2, par[8]/2.0 ) * TMath::Gamma(par[8]/2.0) ) );
	
	
	
	return f;
}


void fit_histo(TString filenameinput, double hhhh, int fittingfunction, int n, TString titleExternal = "") {
	
	double eta = 0.5;
	gStyle->SetPalette(1);
	gStyle->SetPadColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetOptStat(0000000);
	gStyle->SetOptFit(0000);
	gStyle->SetPalette(1);
	gStyle->SetLabelFont(42, "xyz");
	gStyle->SetTitleFont(42, "xyz");
	Double_t dimWindowX = 1400;
	Double_t dimWindowY = 900;


	
	TString titTS, titcumTS;

	titTS= "TS distribution " + titleExternal;
	titcumTS = "P(TS >= h) ";

	TH1D *h1, *hh1, *hh2, *hh3, *hh4;

	TH1D *h4;//se un bin vale 0 significa che non è stato analizzato. Se vale un numero piccolo prossimo a zero
	TH1D* chi2_1, *chi2_1_05, *chi2_3_05;
	TH1D* integral_chi2_1, *integral_chi2_1_05, *integral_chi2_3_05;

	Double_t cumTSdisbinx, max_cumTSdisbinx ;

	max_cumTSdisbinx = 36;
	cumTSdisbinx = max_cumTSdisbinx*MULTIBIN;
		
	h1 = new TH1D("TS2", titTS, cumTSdisbinx, 0, max_cumTSdisbinx );
	h4 = new TH1D("cumsqrTS4", titcumTS + titleExternal, cumTSdisbinx , 0, max_cumTSdisbinx );
	hPVALUE = new TH1D("PVALUE", "PVALUE", cumTSdisbinx , 0, max_cumTSdisbinx );
	hCDF = new TH1D("CDF", "CDF", cumTSdisbinx , 0, max_cumTSdisbinx );
	chi2_1 = new TH1D("\\chi^{2}_{1}", "\\chi^{2}_{1}", cumTSdisbinx , 0, max_cumTSdisbinx );
	TString histname;
	histname += eta;
	histname += "*\\chi^{2}_{1}";
	chi2_1_05 = new TH1D(histname, histname, cumTSdisbinx , 0, max_cumTSdisbinx );
	TString histname33;
	histname33 += eta;
	histname33 += "*\\chi^{2}_{3}";
	chi2_3_05 = new TH1D(histname33, histname33, cumTSdisbinx , 0, max_cumTSdisbinx );
	
	integral_chi2_1 = new TH1D("int \\chi^{2}_{1}", "int \\chi^{2}_{1}", cumTSdisbinx , 0, max_cumTSdisbinx );
	TString histname2, histname3;
	histname2 += eta;
	histname3 += eta;
	//histname2 = histname2.Chop();
	histname2 += "* int \\chi^{2}_{1}";
	histname3 += "* int \\chi^{2}_{3}";
	integral_chi2_1_05 = new TH1D(histname2,histname2, cumTSdisbinx , 0, max_cumTSdisbinx );
	integral_chi2_3_05 = new TH1D(histname3,histname3, cumTSdisbinx , 0, max_cumTSdisbinx );
	TH1D* integral_bestfit = new TH1D("Best fit","Best fit", cumTSdisbinx , 0, max_cumTSdisbinx );
	

	
	TRandom m;
	Int_t counts = 0;
	TF1* f1 = new TF1("f_chi2", "TMath::Exp(-x/2.0) / (TMath::Sqrt(TMath::TwoPi() * x)) ", 0, max_cumTSdisbinx);
	TF1* f2 = new TF1("eta * f_chi2", "[0] * TMath::Exp(-x/2.0) / ( TMath::Sqrt(TMath::TwoPi() * x)) ", 0, max_cumTSdisbinx);
	f2->SetParameter(0, eta);
	TF1* f1_der = new TF1("f_chi2_der", "TMath::Exp(-x/2.0) / (TMath::Sqrt(TMath::TwoPi() * x)) * (1 + 1 / x)", 0, max_cumTSdisbinx);
	TF1* f2_der = new TF1("eta * f_chi2_der", "[0] * TMath::Exp(-x/2.0) / (TMath::Sqrt(TMath::TwoPi() * x)) * (1 + 1 / x)", 0, max_cumTSdisbinx);
	f2_der->SetParameter(0, eta);
	TF1* f3 = new TF1("1/sqrt(N)", "1.0 / (TMath::Sqrt(x)) ", 0, max_cumTSdisbinx); 
	
	Int_t counts = 0;
	
	
	Long64_t nlines;
	TTree* T;
	
	//TCanvas* c1_h1 = new TCanvas("Distribution", "Distribution", dimWindowX, dimWindowY); gPad->SetLogy(); 	gPad->SetTitle();
	
	
	Float_t BIN, BINCENTER, DPDF, EPDF, DPVALUE, EPVALUE, D, E, FPDF, FPVALUE, F, DPVALUE, EPVALUE, DCDF, ECDF;
	for (int i=0; i<4; i++) {
		if (filenameinput == "") 
			continue;
		
		T = new TTree("DATA", "") ;
		
		nlines = T->ReadFile(filenameinput, "BIN:DPDF:EPDF:DPVALUE:EPVALUE:DCDF:ECDF");
		
		T->SetBranchAddress("BIN", &BIN);
		T->SetBranchAddress("DPDF", &DPDF);
		T->SetBranchAddress("EPDF", &EPDF);
		T->SetBranchAddress("DPVALUE", &DPVALUE);
		T->SetBranchAddress("EPVALUE", &EPVALUE);
		T->SetBranchAddress("DCDF", &DCDF);
		T->SetBranchAddress("ECDF", &ECDF);
		
		cout << "N lines: " << nlines << endl;
		
		if(nlines == 0) {
			cout << "End of the procedure" << endl;
			return;
		}
		
		TString clname;
		clname += (int)(gRandom->Uniform()*1000);
		TH1F* h1cl = h1->Clone(clname);
		h1cl->SetName(clname);
		h1cl->SetTitle(clname);
		
		clname = "";
		clname += (int)(gRandom->Uniform()*1000);
		TH1F* h4cl = h4->Clone(clname);
		h4cl->SetName(clname);
		h4cl->SetTitle(clname);
		
		for(Long64_t j = 0; j<nlines; j++) {
			T->GetEntry(j);
			BIN=BIN+1;
			h1->SetBinContent(BIN, DPDF);
			h1cl->SetBinContent(BIN, DPDF);
			h1->SetBinError(BIN, EPDF);
			h4->SetBinContent(BIN, DPVALUE);
			hPVALUE->SetBinContent(BIN, DPVALUE);
			hPVALUE->SetBinError(BIN, EPVALUE);
			h4->SetBinError(BIN, EPVALUE);
			hCDF->SetBinContent(BIN, DCDF);
			hCDF->SetBinError(BIN, ECDF);
			
		}
		
	}	
	
	
	
	
	
	for(int i=0; i<cumTSdisbinx; i++) {
		int jjj=i+1;
		
		cout << jjj <<  "\t";
		cout << h1->GetBinContent(jjj) <<  "\t" << h1->GetBinError(jjj) << "\t";
		cout << hPVALUE->GetBinContent(jjj) <<  "\t" << hPVALUE->GetBinError(jjj) << endl;
	}	
	
	
	
	
	
	//fitting
	TF1* f2fitter = new TF1("fchi2_1fit", "[0] * TMath::Exp(-x/2.0) / ( TMath::Sqrt(TMath::TwoPi() * x)) ", 0, max_cumTSdisbinx);
	f2fitter->SetParName(0, "eta");
	f2fitter->SetParLimits(0, 0, 1);
	f2fitter->SetParameter(0, 0.5);
	TF1* f3fitter = new TF1("fchi2_Nfit", " [0] * TMath::Exp(-x/2.0) * TMath::Power(x, [1]/2.0 - 1) / ( TMath::Power(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", 0, max_cumTSdisbinx);
	f3fitter->SetParName(0, "eta");
	f3fitter->SetParLimits(0, 0, 1);
	f3fitter->SetParameter(0, 0.5);
	f3fitter->SetParName(1, "N");
	f3fitter->SetParLimits(1, 1, 10);
	f3fitter->SetParameter(1, 2);
	//TFormula* ff1 = new TFormula("fkingfunctionsum", "(1. - 1./[1]) * TMath::Power(1. +  ((x/[0])**2.)/(2.0*[1]), -[1]) + (1. - 1./[3]) * TMath::Power(1. +  ((x/[2])**2.)/(2.0*[3]), -[3])");
	TFormula* ff2 = new TFormula("fkingfunction", "(1. - 1./[1]) * TMath::Power(1. +  ((x/[0])**2.)/(2.0*[1]), -[1])");
	TF1* f3fitter = new TF1("fchi2_Nfit_kf", "(1. - 1./[3]) * TMath::Power(1. +  ((x/[2])**2.)/(2.0*[3]), -[3]) + [0] * TMath::Exp(-x/2.0) * TMath::Power(x, [1]/2.0 - 1) / ( TMath::Power(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", 0, max_cumTSdisbinx);
	TF1* f3fitter = new TF1("fchi2_Nfit_gaus1", "( 1.0 / TMath::Sqrt(TMath::TwoPi() * TMath::Power([2],2) ) ) *  TMath::Exp(-x**2 / (2.0 * TMath::Power([2],2) )) + [0] * TMath::Exp(-x/2.0) * TMath::Power(x, [1]/2.0 - 1) / ( TMath::Power(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", 0, max_cumTSdisbinx);
	TF1* f3fitter = new TF1("fchi2_Nfit_gaus2", " gaus + [3] * TMath::Exp(-x/2.0) * TMath::Power(x, [4]/2.0 - 1) / ( TMath::Power(2, [4]/2.0 ) * TMath::Gamma([4]/2.0) ) ", 0, max_cumTSdisbinx);
	TF1* f3fitter = new TF1("fdouble_chi2_Nfit", " [0] * TMath::Exp(-x/2.0) * TMath::Power(x, [1]/2.0 - 1) / ( TMath::Power(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) + [2] * TMath::Exp(-x/2.0) * TMath::Power(x, [3]/2.0 - 1) / ( TMath::Power(2, [3]/2.0 ) * TMath::Gamma([3]/2.0) ) ", 0, max_cumTSdisbinx);
	
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
	
	
	
	//fitting con range
	TF1* f3sf = new TF1("sumchi2delta", sumchi2deltafunc, 0, max_cumTSdisbinx,2);
	f3sf->SetParName(0, "delta");
	f3sf->SetParName(1, "eta");
	f3sf->SetParLimits(1, 0, 1);
	f3sf->SetParameter( 1, 0.5);
	f3sf->SetLineColor(kBlue);
	TF1* f4sf = new TF1("sumchi2(N)delta", sumchi2Ndeltafunc, 0, max_cumTSdisbinx,3);
	f4sf->SetParName(0, "delta");
	f4sf->SetParName(1, "eta");
	f4sf->SetParLimits(1, 0, 1);
	f4sf->SetParameter( 1, 0.5);
	f4sf->SetParName(2, "N");
	f4sf->SetParLimits(2, 0, 10);
	f4sf->SetParameter( 2, 1);
	f4sf->SetLineColor(kGreen);
	TF1* f4sf1 = new TF1("sumchi2(1)delta", sumchi2Ndeltafunc, 0, max_cumTSdisbinx,3);
	f4sf1->SetParName(0, "delta");
	f4sf1->SetParName(1, "eta");
	f4sf1->SetParLimits(1, 0, 1);
	f4sf1->SetParameter( 1, 0.5);
	f4sf1->SetParName(2, "N");
	f4sf1->SetParLimits(2, 1, 1);
	f4sf1->SetParameter( 2, 1);
	f4sf1->SetLineColor(kBlack);
	TF1* f4sf2 = new TF1("sumchi2(2)delta", sumchi2Ndeltafunc, 0, max_cumTSdisbinx,3);
	f4sf2->SetParName(0, "delta");
	f4sf2->SetParName(1, "eta");
	f4sf2->SetParLimits(1, 0, 1);
	f4sf2->SetParameter( 1, 0.5);
	f4sf2->SetParName(2, "N");
	f4sf2->SetParLimits(2, 2, 2);
	f4sf2->SetParameter( 2, 2);
	f4sf2->SetLineColor(kRed);
	TF1* f4sf3 = new TF1("sumchi2(3)delta", sumchi2Ndeltafunc, 0, max_cumTSdisbinx,3);
	f4sf3->SetParName(0, "delta");
	f4sf3->SetParName(1, "eta");
	f4sf3->SetParLimits(1, 0, 1);
	f4sf3->SetParameter( 1, 0.5);
	f4sf3->SetParName(2, "N");
	f4sf3->SetParLimits(2, 3, 3);
	f4sf3->SetParameter( 2, 3);
	f4sf3->SetLineColor(kGreen);
	
	TF1* f4sfv2 = new TF1("sumchi2(N)deltaV2", sumchi2Ndeltafunc_V2, 0, max_cumTSdisbinx,3);
	f4sfv2->SetParName(0, "delta");
	f4sfv2->SetParName(1, "eta");
	f4sfv2->SetParLimits(1, 0, 1);
	f4sfv2->SetParameter( 1, 0.5);
	f4sfv2->SetParName(2, "N");
	f4sfv2->SetParLimits(2, 0, 10);
	f4sfv2->SetParameter( 2, 1);
	f4sfv2->SetLineColor(kBlack);
	TF1* f4sf1v2 = new TF1("sumchi2(1)deltaV2", sumchi2Ndeltafunc_V2, 0, max_cumTSdisbinx,3);
	f4sf1v2->SetParName(0, "delta");
	f4sf1v2->SetParName(1, "eta");
	f4sf1v2->SetParLimits(1, 0, 1);
	f4sf1v2->SetParameter( 1, 0.5);
	f4sf1v2->SetParName(2, "N");
	f4sf1v2->SetParLimits(2, 1, 1);
	f4sf1v2->SetParameter( 2, 1);
	f4sf1v2->SetLineColor(kBlack);
	TF1* f4sf2v2 = new TF1("sumchi2(2)deltaV2", sumchi2Ndeltafunc_V2, 0, max_cumTSdisbinx,3);
	f4sf2v2->SetParName(0, "delta");
	f4sf2v2->SetParName(1, "eta");
	f4sf2v2->SetParLimits(1, 0, 1);
	f4sf2v2->SetParameter( 1, 0.5);
	f4sf2v2->SetParName(2, "N");
	f4sf2v2->SetParLimits(2, 2, 2);
	f4sf2v2->SetParameter( 2, 2);
	f4sf2v2->SetLineColor(kRed);
	TF1* f4sf3v2 = new TF1("sumchi2(3)deltaV2", sumchi2Ndeltafunc_V2, 0, max_cumTSdisbinx,3);
	f4sf3v2->SetParName(0, "delta");
	f4sf3v2->SetParName(1, "eta");
	f4sf3v2->SetParLimits(1, 0, 1);
	f4sf3v2->SetParameter( 1, 0.5);
	f4sf3v2->SetParName(2, "N");
	f4sf3v2->SetParLimits(2, 3, 3);
	f4sf3v2->SetParameter( 2, 3);
	f4sf3v2->SetLineColor(kBlack);
	TF1* f4sf4v2 = new TF1("sumchi2(4)deltaV2", sumchi2Ndeltafunc_V2, 0, max_cumTSdisbinx,3);
	f4sf4v2->SetParName(0, "delta");
	f4sf4v2->SetParName(1, "eta");
	f4sf4v2->SetParLimits(1, 0, 1);
	f4sf4v2->SetParameter( 1, 0.5);
	f4sf4v2->SetParName(2, "N");
	f4sf4v2->SetParLimits(2, 4, 4);
	f4sf4v2->SetParameter( 2, 4);
	f4sf4v2->SetLineColor(kBlack);
	
	
	TF1* f4sf_trasl_N1N2v2 = new TF1("traslchi2N1N2deltafuncV2", traslchi2N1N2deltafuncV2, 0, max_cumTSdisbinx,6);
	f4sf_trasl_N1N2v2->SetParName(0, "delta");
	f4sf_trasl_N1N2v2->SetParName(1, "eta1");
	f4sf_trasl_N1N2v2->SetParLimits(1, 0, 1);
	f4sf_trasl_N1N2v2->SetParameter( 1, 0.5);
	f4sf_trasl_N1N2v2->SetParName(2, "N1");
	Int_t f4sf_test_N1=1; //fix N1
	f4sf_trasl_N1N2v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //fix it!
	f4sf_trasl_N1N2v2->SetParameter( 2, f4sf_test_N1);
	f4sf_trasl_N1N2v2->SetParName(3, "eta2");
	f4sf_trasl_N1N2v2->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2v2->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2v2->SetParName(4, "N2");
	Int_t f4sf_test_N2=5; //fix N2
	f4sf_trasl_N1N2v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2); //fix it!
	f4sf_trasl_N1N2v2->SetParameter( 4, f4sf_test_N2);
	Int_t f4sf_trasl_N1N2v2_trasl = 6; //fix translation
	f4sf_trasl_N1N2v2->SetParName(5, "translation2");
	f4sf_trasl_N1N2v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //fix it!
	f4sf_trasl_N1N2v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2v2->SetLineColor(kBlack);
	
	TF1* f4sf_trasl_N1N2v2_F = new TF1("traslchi2N1N2deltafuncV2_ALLFREE", traslchi2N1N2deltafuncV2, 0, max_cumTSdisbinx,6);
	f4sf_trasl_N1N2v2_F->SetParName(0, "delta");
	f4sf_trasl_N1N2v2_F->SetParName(1, "eta1");
	f4sf_trasl_N1N2v2_F->SetParLimits(1, 0, 1);
	f4sf_trasl_N1N2v2_F->SetParameter( 1, 0.5);
	f4sf_trasl_N1N2v2_F->SetParName(2, "N1");
	f4sf_test_N1=1; //starting value
	f4sf_trasl_N1N2v2_F->SetParLimits(2, 0, 5);
	f4sf_trasl_N1N2v2_F->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //fix it!
	f4sf_trasl_N1N2v2_F->SetParameter( 2, f4sf_test_N1);
	f4sf_trasl_N1N2v2_F->SetParName(3, "eta2");
	f4sf_trasl_N1N2v2_F->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2v2_F->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2v2_F->SetParName(4, "N2");
	f4sf_test_N2=3; //starting value
	f4sf_trasl_N1N2v2_F->SetParLimits(4, 0, 10);
	f4sf_trasl_N1N2v2_F->SetParameter( 4, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 6; //starting value
	f4sf_trasl_N1N2v2_F->SetParName(5, "translation2");
	f4sf_trasl_N1N2v2_F->SetParLimits(5, 0, 36);
	f4sf_trasl_N1N2v2_F->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //fix it!
	f4sf_trasl_N1N2v2_F->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2v2_F->SetLineColor(kBlack);
	
	
	TF1* f4sf_trasl_N1N2N3v2_F = new TF1("traslchi2N1N2N3deltafuncV2_ALLFREE", traslchi2N1N2N3deltafuncV2, 0, max_cumTSdisbinx,9);
	f4sf_trasl_N1N2N3v2_F->SetParName(0, "delta");
	f4sf_trasl_N1N2N3v2_F->SetParName(1, "eta1");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(1, 0, 1);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 1, 0.5);
	f4sf_trasl_N1N2N3v2_F->SetParName(2, "N1");
	f4sf_test_N1=1;
	f4sf_trasl_N1N2N3v2_F->SetParLimits(2, 0, 5);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 2, f4sf_test_N1);
	//f4sf_trasl_N1N2N3v2_F->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //!!!!!!!!!!!!
	f4sf_trasl_N1N2N3v2_F->SetParName(3, "eta2");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2N3v2_F->SetParName(4, "N2");
	f4sf_test_N2=5;
	f4sf_trasl_N1N2N3v2_F->SetParLimits(4, 0, 5);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 4, f4sf_test_N2);
	//f4sf_trasl_N1N2N3v2_F->SetParLimits(4, f4sf_test_N2, f4sf_test_N2); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_trasl_N1N2N3v2_F->SetParName(5, "translation2");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(5, 0, 36);
	f4sf_trasl_N1N2N3v2_F->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2_F->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	f4sf_trasl_N1N2N3v2_F->SetParName(6, "eta3");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(6, 0, 1);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 6, 0.5);
	f4sf_trasl_N1N2N3v2_F->SetParName(7, "N3");
	f4sf_test_N2=3;
	f4sf_trasl_N1N2N3v2_F->SetParLimits(7, 0, 5);
	f4sf_trasl_N1N2N3v2_F->SetParameter(7, f4sf_test_N2);
	//f4sf_trasl_N1N2N3v2_F->SetParLimits(7, f4sf_test_N2, f4sf_test_N2); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 12;
	f4sf_trasl_N1N2N3v2_F->SetParName(8, "translation3");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(8, 0, 36);
	f4sf_trasl_N1N2N3v2_F->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2_F->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	
	f4sf_trasl_N1N2N3v2_F->SetLineColor(kBlack);

	
	TF1* f4sf_trasl_N1N2N3v2 = new TF1("traslchi2N1N2N3deltafuncV2", traslchi2N1N2N3deltafuncV2, 0, max_cumTSdisbinx,9);
	f4sf_trasl_N1N2N3v2->SetParName(0, "delta");
	f4sf_trasl_N1N2N3v2->SetParName(1, "eta1");
	f4sf_trasl_N1N2N3v2->SetParLimits(1, 0, 1);
	f4sf_trasl_N1N2N3v2->SetParameter( 1, 0.5);
	f4sf_trasl_N1N2N3v2->SetParName(2, "N1");
	f4sf_test_N1=1;
	f4sf_trasl_N1N2N3v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1);
	f4sf_trasl_N1N2N3v2->SetParameter( 2, f4sf_test_N1);
	f4sf_trasl_N1N2N3v2->SetParName(3, "eta2");
	f4sf_trasl_N1N2N3v2->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2N3v2->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2N3v2->SetParName(4, "N2");
	f4sf_test_N2=5;
	f4sf_trasl_N1N2N3v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2);
	f4sf_trasl_N1N2N3v2->SetParameter( 4, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_trasl_N1N2N3v2->SetParName(5, "translation2");
	f4sf_trasl_N1N2N3v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2->SetParName(6, "eta3");
	f4sf_trasl_N1N2N3v2->SetParLimits(6, 0, 1);
	f4sf_trasl_N1N2N3v2->SetParameter( 6, 0.5);
	f4sf_trasl_N1N2N3v2->SetParName(7, "N3");
	f4sf_test_N2=3;
	f4sf_trasl_N1N2N3v2->SetParLimits(7, f4sf_test_N2, f4sf_test_N2);
	f4sf_trasl_N1N2N3v2->SetParameter(7, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 9;
	f4sf_trasl_N1N2N3v2->SetParName(8, "translation3");
	f4sf_trasl_N1N2N3v2->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	
	f4sf_trasl_N1N2N3v2->SetLineColor(kBlack);
	
	
	
	TF1* f4sf_sumN1N2traslN3v2 = new TF1("chi2delta_sumN1taslN2sumN3funcV2", chi2delta_sumN1taslN2sumN3funcV2, 0, max_cumTSdisbinx,9);
	f4sf_sumN1N2traslN3v2->SetParName(0, "delta");
	f4sf_sumN1N2traslN3v2->SetParName(1, "eta1");
	f4sf_sumN1N2traslN3v2->SetParLimits(1, 0, 1);
	f4sf_sumN1N2traslN3v2->SetParameter( 1, 0.5);
	f4sf_sumN1N2traslN3v2->SetParName(2, "N1");
	f4sf_test_N1=1;
	f4sf_sumN1N2traslN3v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1);
	f4sf_sumN1N2traslN3v2->SetParameter( 2, f4sf_test_N1);
	f4sf_sumN1N2traslN3v2->SetParName(3, "eta2");
	f4sf_sumN1N2traslN3v2->SetParLimits(3, 0, 1);
	f4sf_sumN1N2traslN3v2->SetParameter( 3, 0.5);
	f4sf_sumN1N2traslN3v2->SetParName(4, "N2");
	f4sf_test_N2=5;
	f4sf_sumN1N2traslN3v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2);
	f4sf_sumN1N2traslN3v2->SetParameter( 4, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_sumN1N2traslN3v2->SetParName(5, "translation2");
	f4sf_sumN1N2traslN3v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_sumN1N2traslN3v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_sumN1N2traslN3v2->SetParName(6, "eta3");
	f4sf_sumN1N2traslN3v2->SetParLimits(6, 0, 1);
	f4sf_sumN1N2traslN3v2->SetParameter( 6, 0.5);
	f4sf_sumN1N2traslN3v2->SetParName(7, "N3");
	f4sf_test_N2=3;
	//f4sf_sumN1N2traslN3v2->SetParLimits(7, f4sf_test_N2, f4sf_test_N2);
	f4sf_sumN1N2traslN3v2->SetParameter(7, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 14;
	f4sf_sumN1N2traslN3v2->SetParName(8, "translation3");
	f4sf_sumN1N2traslN3v2->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_sumN1N2traslN3v2->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	
	f4sf_sumN1N2traslN3v2->SetLineColor(kBlack);
	
	
	/*
	 par[0] = normalization of delta function
	 par[1] = locCL TS threshold
	 par[2] = normalization of below threshold chi^2
	 par[3] = dof of below threshold chi^2
	 par[4] = normalization of above threshold chi^2
	 par[5] = dof of above threshold chi^2
	 */
	TF1* f4sfN1N2 = new TF1("sumchi2(N1)(N2)delta", sumchi2N1N2deltafunc, 0, max_cumTSdisbinx,6);
	f4sfN1N2->SetParName(0, "delta");
	f4sfN1N2->SetParLimits(0, 0.0, 1);
	f4sfN1N2->SetParameter( 0, 0.7);
	int constloclc = 4;
	f4sfN1N2->SetParName(1, "loccl");
	f4sfN1N2->SetParLimits(1, constloclc,constloclc);
	f4sfN1N2->SetParameter( 1, constloclc);
	f4sfN1N2->SetParName(2, "eta1");
	f4sfN1N2->SetParLimits(2, 0, 1);
	f4sfN1N2->SetParameter( 2, 0.5);
	Int_t f4sfN1N2_N1 = 3;
	f4sfN1N2->SetParName(3, "N1");
	f4sfN1N2->SetParLimits(3, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2->SetParameter(3, f4sfN1N2_N1);
	f4sfN1N2->SetParName(4, "eta2");
	f4sfN1N2->SetParLimits(4, 0, 1);
	f4sfN1N2->SetParameter( 4, 0.5);
	Int_t f4sfN1N2_N2 = 3;
	f4sfN1N2->SetParName(5, "N2");
	f4sfN1N2->SetParLimits(5, f4sfN1N2_N2, f4sfN1N2_N2);
	f4sfN1N2->SetParameter(5, f4sfN1N2_N2);
	f4sfN1N2->SetLineColor(kBlack);

	TF1* f4sfN1N2v2 = new TF1("sumchi2(N1)(N2)deltaV2", sumchi2N1N2deltafuncV2, 0, max_cumTSdisbinx,6);
	f4sfN1N2v2->SetParName(0, "delta");
	f4sfN1N2v2->SetParLimits(0, 0.0, 1);
	//f4sfN1N2->SetParameter( 0, 0.7);
	constloclc=8;
	f4sfN1N2v2->SetParName(1, "loccl");
	f4sfN1N2v2->SetParLimits(1, constloclc,constloclc);
	f4sfN1N2v2->SetParameter( 1, constloclc);
	
	f4sfN1N2v2->SetParName(2, "eta1");
	f4sfN1N2v2->SetParLimits(2, 0, 1);
	f4sfN1N2v2->SetParameter( 2, 0.5);
	
	f4sfN1N2_N1 = 3;
	f4sfN1N2v2->SetParName(3, "N1");
	//f4sfN1N2v2->SetParLimits(3, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2v2->SetParameter(3, f4sfN1N2_N1);
	
	f4sfN1N2v2->SetParName(4, "eta2");
	f4sfN1N2v2->SetParLimits(4, 0, 10);
	f4sfN1N2v2->SetParameter( 4, 0.5);
	
	f4sfN1N2_N2 = 2;
	f4sfN1N2v2->SetParName(5, "N2");
	//f4sfN1N2v2->SetParLimits(5, f4sfN1N2_N2, f4sfN1N2_N2);
	f4sfN1N2v2->SetParameter(5, f4sfN1N2_N2);
	
	f4sfN1N2v2->SetLineColor(kBlack);
	
	
	/*
	 par[0] = normalization of delta function
	 par[1] = threshold 1
	 par[2] = normalization of below threshold chi^2
	 par[3] = dof of below threshold chi^2
	 par[4] = normalization of above threshold chi^2
	 par[5] = dof of above threshold chi^2
	 par[6] = threshodl 2
	 par[7] = normalization of above threshold chi^2
	 par[8] = dof of above threshold chi^2
	 */
	TF1* f4sfN1N2N3v2 = new TF1("sumchi2(N1)(N2)(N3)deltaV2", sumchi2N1N2N3deltafuncV2, 0, max_cumTSdisbinx,9);
	f4sfN1N2N3v2->SetParName(0, "delta");
	f4sfN1N2N3v2->SetParLimits(0, 0.0, 1);
	f4sfN1N2N3v2->SetParameter( 0, 0.7);
	constloclc=4;
	f4sfN1N2N3v2->SetParName(1, "thr1");
	f4sfN1N2N3v2->SetParLimits(1, constloclc,constloclc);
	f4sfN1N2N3v2->SetParameter( 1, constloclc);
	
	f4sfN1N2N3v2->SetParName(2, "eta1");
	f4sfN1N2N3v2->SetParLimits(2, 0, 1);
	f4sfN1N2N3v2->SetParameter( 2, 0.5);
	
	f4sfN1N2_N1 = 1;
	f4sfN1N2N3v2->SetParName(3, "N1");
	f4sfN1N2N3v2->SetParLimits(3, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2N3v2->SetParameter(3, f4sfN1N2_N1);
	
	f4sfN1N2N3v2->SetParName(4, "eta2");
	f4sfN1N2N3v2->SetParLimits(4, 0, 1);
	f4sfN1N2N3v2->SetParameter( 4, 0.5);
	
	f4sfN1N2_N2 = 3;
	f4sfN1N2N3v2->SetParName(5, "N2");
	f4sfN1N2N3v2->SetParLimits(5, f4sfN1N2_N2, f4sfN1N2_N2);
	f4sfN1N2N3v2->SetParameter(5, f4sfN1N2_N2);
	
	constloclc=9;
	f4sfN1N2N3v2->SetParName(6, "thr2");
	f4sfN1N2N3v2->SetParLimits(6, constloclc,constloclc);
	f4sfN1N2N3v2->SetParameter(6, constloclc);
	
	
	f4sfN1N2N3v2->SetParName(7, "eta3");
	f4sfN1N2N3v2->SetParLimits(7, 0, 200);
	f4sfN1N2N3v2->SetParameter(7, 0.5);
	
	f4sfN1N2_N1 = 1;
	f4sfN1N2N3v2->SetParName(8, "N3");
	//f4sfN1N2N3v2->SetParLimits(8, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2N3v2->SetParameter(8, f4sfN1N2_N1);
	
	f4sfN1N2N3v2->SetLineColor(kBlack);
	
	double f4sf_test_N3 = 0;
	
	TF1* f4sf_traslN1N2sumlN3v2 = new TF1("traslchi2N1N2sumN3deltafuncV2", traslchi2N1N2sumN3deltafuncV2, 0, max_cumTSdisbinx,9);
	f4sf_traslN1N2sumlN3v2->SetParName(0, "delta");
	f4sf_traslN1N2sumlN3v2->SetParName(1, "eta1");
	f4sf_traslN1N2sumlN3v2->SetParLimits(1, 0, 1);
	f4sf_traslN1N2sumlN3v2->SetParameter( 1, 0.5);
	f4sf_traslN1N2sumlN3v2->SetParName(2, "N1");
	f4sf_test_N1=1;
	f4sf_traslN1N2sumlN3v2->SetParLimits(2, 0, 5);
	f4sf_traslN1N2sumlN3v2->SetParameter( 2, f4sf_test_N1);
	f4sf_traslN1N2sumlN3v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //!!!!!!!!!!!!
	f4sf_traslN1N2sumlN3v2->SetParName(3, "eta2");
	f4sf_traslN1N2sumlN3v2->SetParLimits(3, 0, 1);
	f4sf_traslN1N2sumlN3v2->SetParameter( 3, 0.5);
	f4sf_traslN1N2sumlN3v2->SetParName(4, "N2");
	f4sf_test_N2=5;
	f4sf_traslN1N2sumlN3v2->SetParLimits(4, 0, 5);
	f4sf_traslN1N2sumlN3v2->SetParameter( 4, f4sf_test_N2);
	f4sf_traslN1N2sumlN3v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 5;
	f4sf_traslN1N2sumlN3v2->SetParName(5, "translation2");
	f4sf_traslN1N2sumlN3v2->SetParLimits(5, 0, 36);
	f4sf_traslN1N2sumlN3v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	f4sf_traslN1N2sumlN3v2->SetParName(6, "eta3");
	f4sf_traslN1N2sumlN3v2->SetParLimits(6, 0, 1);
	f4sf_traslN1N2sumlN3v2->SetParameter( 6, 0.5);
	f4sf_traslN1N2sumlN3v2->SetParName(7, "N3");
	f4sf_test_N3=3;
	f4sf_traslN1N2sumlN3v2->SetParLimits(7, 0, 5);
	f4sf_traslN1N2sumlN3v2->SetParameter(7, f4sf_test_N3);
	f4sf_traslN1N2sumlN3v2->SetParLimits(7, f4sf_test_N3, f4sf_test_N3); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 8;
	f4sf_traslN1N2sumlN3v2->SetParName(8, "translation3");
	f4sf_traslN1N2sumlN3v2->SetParLimits(8, 0, 36);
	f4sf_traslN1N2sumlN3v2->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3v2->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	
	f4sf_traslN1N2sumlN3v2->SetLineColor(kBlack);
	
	double f4sf_test_N4 = 0;
	
	//M95+ICL=9**(23)*********************************************************************************************************************
	//M95*******************************************************************************************************************************
	//M95*******************************************************************************************************************************
	TF1* f4sf_traslN1N2sumlN3N4v2 = new TF1("traslchi2N1N2sumN3N4deltafuncV2", traslchi2N1N2sumN3N4deltafuncV2, 0, max_cumTSdisbinx,12);
	f4sf_traslN1N2sumlN3N4v2->SetParName(0, "delta");
	f4sf_traslN1N2sumlN3N4v2->SetParName(1, "eta1");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(1, 0, 1);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 1, 0.5);
	f4sf_traslN1N2sumlN3N4v2->SetParName(2, "N1");
	f4sf_test_N1=1;
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(2, 0, 5);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 2, f4sf_test_N1);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //!!!!!!!!!!!!
	f4sf_traslN1N2sumlN3N4v2->SetParName(3, "eta2");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(3, 0, 1);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 3, 0.5);
	f4sf_traslN1N2sumlN3N4v2->SetParName(4, "N2");
	f4sf_test_N2=5;
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
	f4sf_test_N3=5;
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(7, 1, 10);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(7, f4sf_test_N3);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(7, f4sf_test_N3, f4sf_test_N3); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 9;
	f4sf_traslN1N2sumlN3N4v2->SetParName(8, "translation3");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(8, 0, 36);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	double eta4 = 3;
	f4sf_traslN1N2sumlN3N4v2->SetParName(9, "eta4");
	//f4sf_traslN1N2sumlN3N4v2->SetParLimits(9, eta4, eta4);
	f4sf_traslN1N2sumlN3N4v2->SetParameter( 9, eta4);
	f4sf_traslN1N2sumlN3N4v2->SetParName(10, "N4");
	f4sf_test_N4=2;
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(10, 0, 5);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(10, f4sf_test_N4);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(10, f4sf_test_N4, f4sf_test_N4); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 14;
	f4sf_traslN1N2sumlN3N4v2->SetParName(11, "translation3");
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(11, 0, 36);
	f4sf_traslN1N2sumlN3N4v2->SetParameter(11, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v2->SetParLimits(11, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	
	f4sf_traslN1N2sumlN3N4v2->SetLineColor(kBlack);
	
	//M95+ICL=6**(24)*********************************************************************************************************************
	//M95*******************************************************************************************************************************
	//M95*******************************************************************************************************************************
	TF1* f4sf_traslN1N2sumlN3N4v3 = new TF1("traslchi2N1N2sumN3N4deltafuncV3", traslchi2N1N2sumN3N4deltafuncV3, 0, max_cumTSdisbinx, 12);
	f4sf_traslN1N2sumlN3N4v3->SetParName(0, "delta");
	f4sf_traslN1N2sumlN3N4v3->SetParName(1, "eta1");
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(1, 0, 1);
	f4sf_traslN1N2sumlN3N4v3->SetParameter( 1, 0.5);
	f4sf_traslN1N2sumlN3N4v3->SetParName(2, "N1");
	f4sf_test_N1=1;
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(2, 0, 5);
	f4sf_traslN1N2sumlN3N4v3->SetParameter( 2, f4sf_test_N1);
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //!!!!!!!!!!!!
	f4sf_traslN1N2sumlN3N4v3->SetParName(3, "eta2");
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(3, 0, 1);
	f4sf_traslN1N2sumlN3N4v3->SetParameter( 3, 0.5);
	f4sf_traslN1N2sumlN3N4v3->SetParName(4, "N2");
	f4sf_test_N2=5;
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(4, 0, 5);
	f4sf_traslN1N2sumlN3N4v3->SetParameter( 4, f4sf_test_N2);
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(4, f4sf_test_N2, f4sf_test_N2); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_traslN1N2sumlN3N4v3->SetParName(5, "translation2");
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(5, 0, 36);
	f4sf_traslN1N2sumlN3N4v3->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	f4sf_traslN1N2sumlN3N4v3->SetParName(6, "eta3");
	
	f4sf_traslN1N2sumlN3N4v3->SetParameter( 6, eta3);
	f4sf_traslN1N2sumlN3N4v3->SetParName(7, "N3");
	f4sf_test_N3=4;
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(7, 1, 5);
	f4sf_traslN1N2sumlN3N4v3->SetParameter(7, f4sf_test_N3);
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(7, f4sf_test_N3, f4sf_test_N3); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 9;
	f4sf_traslN1N2sumlN3N4v3->SetParName(8, "translation3");
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(8, 0, 36);
	f4sf_traslN1N2sumlN3N4v3->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	
	f4sf_traslN1N2sumlN3N4v3->SetParName(9, "eta4");
	f4sf_traslN1N2sumlN3N4v3->SetParameter( 9, eta4);
	f4sf_traslN1N2sumlN3N4v3->SetParName(10, "N4");
	f4sf_test_N4=3;
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(10, 0, 5);
	f4sf_traslN1N2sumlN3N4v3->SetParameter(10, f4sf_test_N4);
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(10, f4sf_test_N4, f4sf_test_N4); //!!!!!!!!!!!!
	f4sf_trasl_N1N2v2_trasl = 14;
	f4sf_traslN1N2sumlN3N4v3->SetParName(11, "translation3");
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(11, 0, 36);
	f4sf_traslN1N2sumlN3N4v3->SetParameter(11, f4sf_trasl_N1N2v2_trasl);
	f4sf_traslN1N2sumlN3N4v3->SetParLimits(11, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!!!!!!!!
	
	f4sf_traslN1N2sumlN3N4v3->SetLineColor(kBlack);
	
	
	
	/*h1->Fit("fchi2_1fit", "+");
	h1->Fit("fchi2_Nfit", "+");
	h1->Fit("sumchi2delta", "+");
	h1->Fit("sumchi2Ndelta", "+");
	h1->Fit("delta", "+", 0, 0.1);*/
//	h1->Fit("fchi2_Nfit_kf", "+");
//	h1->Fit("fchi2_Nfit_gaus1", "+");
//	h1->Fit("fchi2_Nfit_gaus2", "+");
//	h1->Fit("fdouble_chi2_Nfit", "+");
	//h1->Fit("fkingfunctionsum", "+");
	//h1->Fit("fkingfunction", "+");
	h1->SaveAs("TS2.root");
	//**************************
	TF1* bestfit;
	
	switch(fittingfunction) {
			
		case 1:
	
			cout << "####################### fchi2_1fit " << endl;
			//h1->Fit("fchi2_1fit", "+", "", 0, max_cumTSdisbinx);
			break;
		case 2:
			cout << "####################### fchi2_Nfit " << endl;
			//h1->Fit("fchi2_Nfit", "+", "", 0, max_cumTSdisbinx);
			break;
		case 3:
			cout << "####################### delta " << endl;
			//h1->Fit("delta", "+", "", 0, 0.1);
			break;
		case 4:
			cout << "####################### sumchi2delta " << endl;
			//h1->Fit("sumchi2delta", "+");
			break;
		case 5:
			cout << "####################### sumchi2(N)delta " << endl;
			//h1->Fit("sumchi2(N)delta", "+");
			break;
		case 6:
			cout << "####################### sumchi2(1)delta " << endl;
			h1->Fit("sumchi2(1)delta", "+");
			bestfit = f4sf1;
			break;
		case 7:
			cout << "####################### sumchi2(2)delta " << endl;
			//h1->Fit("sumchi2(2)delta", "+");
			break;
		case 8:
			cout << "####################### sumchi2(3)delta " << endl;
			//h1->Fit("sumchi2(3)delta", "+");
		
			break;
		case 9:	
			cout << "####################### sumchi2(N)deltaV2 - Equation (3) of A&A 540 A79, 2012 free" << endl;
			h1->Fit("sumchi2(N)deltaV2", "+");
			bestfit = f4sfv2;
			break;
		case 10:	
			cout << "####################### sumchi2(1)deltaV2 - Equation (3) of A&A 540 A79, 2012 N=1" << endl;
			h1->Fit("sumchi2(1)deltaV2", "+");
			bestfit = f4sf1v2;
			break;
		case 11:
			cout << "####################### sumchi2(2)deltaV2 - Equation (3) of A&A 540 A79, 2012 N=2" << endl;
			h1->Fit("sumchi2(2)deltaV2", "+");
			bestfit = f4sf2v2;
			break;
		case 12:
			cout << "####################### sumchi2(3)deltaV2 - Equation (3) of A&A 540 A79, 2012 N=3" << endl;
			h1->Fit("sumchi2(3)deltaV2", "+");
			bestfit = f4sf3v2;
			break;
		case 13:
			cout << "####################### sumchi2(4)deltaV2 - Equation (3) of A&A 540 A79, 2012 N=4" << endl;
			h1->Fit("sumchi2(4)deltaV2", "+");
			bestfit = f4sf4v2;
			break;
		case 14:
			cout << "####################### sumchi2(N1)(N2)delta " << endl;
			h1->Fit("sumchi2(N1)(N2)delta", "+");
			bestfit = f4sfN1N2;
			break;
		case 15:
			cout << "####################### sumchi2(N1)(N2)deltaV2 " << endl;
			h1->Fit("sumchi2(N1)(N2)deltaV2", "+");
			bestfit = f4sfN1N2v2;
			break;
		case 16:
			cout << "####################### sumchi2(N1)(N2)(N3)deltaV2 " << endl;
			h1->Fit("sumchi2(N1)(N2)(N3)deltaV2", "+");
			bestfit = f4sfN1N2N3v2;
			break;
		case 17:
			cout << "####################### traslchi2N1N2deltafuncV2 - Equation (4) of A&A 540 A79, 2012" << endl;
			h1->Fit("traslchi2N1N2deltafuncV2", "+");
			bestfit = f4sf_trasl_N1N2v2;
			break;
		case 18:
			cout << "####################### traslchi2N1N2deltafuncV2_ALLFREE - Equation (4) of A&A 540 A79, 2012 but with ALL paramters free" << endl;
			h1->Fit("traslchi2N1N2deltafuncV2_ALLFREE", "+");
			bestfit = f4sf_trasl_N1N2v2_F;
			break;
		case 19:
			cout << "####################### traslchi2N1N2N3deltafuncV2" << endl;
			h1->Fit("traslchi2N1N2N3deltafuncV2", "+");
			bestfit = f4sf_trasl_N1N2N3v2_F;		
			break;
		case 20:
			cout << "####################### traslchi2N1N2N3deltafuncV2_ALLFREE " << endl;
			h1->Fit("traslchi2N1N2N3deltafuncV2_ALLFREE", "+");
			bestfit = f4sf_trasl_N1N2N3v2;	
			break;
		case 21:
			cout << "####################### chi2delta_sumN1taslN2sumN3funcV2 " << endl;
			h1->Fit("chi2delta_sumN1taslN2sumN3funcV2", "+");
			bestfit = f4sf_sumN1N2traslN3v2;
			break;
		case 22:
			cout << "####################### traslchi2N1N2sumN3deltafuncV2 " << endl;
			h1->Fit("traslchi2N1N2sumN3deltafuncV2", "+");
			bestfit = f4sf_traslN1N2sumlN3v2;
			break;
		case 23:
			cout << "####################### traslchi2N1N2sumN3N4deltafuncV2 (M95 with ICL) - sliding window" << endl;
			h1->Fit("traslchi2N1N2sumN3N4deltafuncV2", "+");
			bestfit = f4sf_traslN1N2sumlN3N4v2;
			break;
		case 24:
			cout << "####################### traslchi2N1N2sumN3N4deltafuncV3 (M95 without ICL) " << endl;
			h1->Fit("traslchi2N1N2sumN3N4deltafuncV3", "+");
			bestfit = f4sf_traslN1N2sumlN3N4v3;
			break;
	
	}
	
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
	//params[12] = bestfit->GetParameter(12);
	/*cout << "Fitting with sumchi2(1)delta" << endl;
	cout << "delta " << params[0] << endl;
	cout << "eta " << params[1] << endl;
	cout << "N " << params[2] << endl;
	cout << "chisquare " << f3sf->GetChisquare() << endl;
	cout << "Fitting with sumchi2Ndelta" << endl;
	cout << "delta " << f4sf->GetParameter(0) << endl;
	cout << "eta " << f4sf->GetParameter(1) << endl;
	cout << "N " << f4sf->GetParameter(2) << endl;
	cout << "chisquare " << f4sf->GetChisquare() << endl;
	*/
	
	Double_t evalfirst1 = -1, evalfirst2 = -1;
	
	
	for(int i=1; i<cumTSdisbinx; i++) {
		Double_t eval, eval2;
		Double_t x;
		//histogram of chi2/1
		x = chi2_1->GetBinCenter(i);
		eval = f1->Eval(x);
		if(evalfirst1 == -1)
			evalfirst1 = eval;
		chi2_1->SetBinContent(i, eval);
		
		//histogram of chi2_1/2
		eval2 = f2->Eval(x);
		if(evalfirst2 == -1)
			evalfirst2 = eval2;
		chi2_1_05->SetBinContent(i, eval2);
		
		//histogram of chi2_3/2
		eval2 = f2sf3_05->Eval(x);
		if(evalfirst2 == -1)
			evalfirst2 = eval2;
		chi2_3_05->SetBinContent(i, eval2);
		
		//histogram of derivativo of chi2/1
		//x = derivative_chi2_1->GetBinCenter(i);
		//eval = f1_der->Eval(x);
		//derivative_chi2_1->SetBinContent(i, eval);

		//histogram of derivative of chi2/2
		//x = derivative_chi2_1_05->GetBinCenter(i);
		//eval = f2_der->Eval(x);
		//derivative_chi2_1_05->SetBinContent(i, eval);

	}	
	
	//cout << "Print bestfit" << endl;
	//bestfit->Print();
	
	


	TCanvas* c1_h1 = new TCanvas("TS Distribution", "TS Distribution", dimWindowX, dimWindowY);
	/*
	DUE
	c1_h1->Divide(2,1);
	c1_h1->cd(1);*/
	TH1D* h1_clone = h1->Clone("NH1");
// 	c1_h1->SetLogx();
	gPad->SetLogy();
	h1->GetXaxis()->SetTitle("TS");
	h1->GetYaxis()->SetTitle("counts (normalized)");
//  	h1->Scale(1.0/h1->GetBinContent(h1->GetMaximumBin()) * evalfirst2);
	
// 	h1->Scale(eta * MULTIBIN);
	cout << "INTEGRAL " << h1->Integral(1, cumTSdisbinx-1) << endl;
// 	h1->GetXaxis()->SetRangeUser(0, 16);
	h1->SetLineColor(kBlue);
	h1->SetLineWidth(2);
	h1->GetYaxis()->SetTitleSize(0.03);
	h1->GetYaxis()->SetTitleOffset(1.80);
	h1->Draw("L");
	f1->SetLineColor(kGreen);
// 	f1->Draw("SAME");
	f2->SetLineColor(kRed);
// 	f2->Draw("SAME");
// 	f1_der->Draw("SAME");
// 	f2_der->Draw("SAME");
// 	c1_h1->SetGridx();
// 	c1_h1->SetGridy();
	
	c1_h1->cd(1);
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
	
	double p;
	cout << "P(TS=h) and P(TS>=h) -> p-value" << endl;
	cout << "Best fit: " << bestfit->Eval(hhhh) << endl;
	p= bestfit->Integral(hhhh, cumTSdisbinx, params);
	cout << "p-value bestfit = " << p  << endl;


	/*

	Double_t gmb;
	//TCanvas* c2 = new TCanvas("Cumulative TS distribution", "Cumulative TS distribution", dimWindowX, dimWindowY);
	c1_h1->cd(2);
	//cout << "chi2_3 reference: " << endl;
	//f2sfN->Print();
	cout << filenameinput << ".distf" << endl;
	cout << "(TS) (PDF) (PDF ERROR) (PVALUE) (PVALUE ERROR) (CDF) (CDF ERROR) (BEST FIT) (BEST FIT ERROR)" << endl;
	for(int i=0; i<cumTSdisbinx; i++) {
		Double_t integral, integral2, integral3;
		
		int jjj=i+1;
		integral = h1->Integral(i, cumTSdisbinx);
		h4->SetBinContent(i, integral);
		
		//histogram of integral of chi2_1
		integral = f2sf1->Integral(i, cumTSdisbinx);
		integral_chi2_1->SetBinContent(jjj, integral);
		//histogram of integral of chi2_1/2
		integral2 = f2sf1_05->Integral(i, cumTSdisbinx);
		integral_chi2_1_05->SetBinContent(jjj, integral2);
		//histogram of integral of chi2_3/2
		integral2 = f2sf3_05->Integral(i, cumTSdisbinx);
		integral_chi2_3_05->SetBinContent(jjj, integral2);
		
		//AB AB AB da verificare i ranges ################################################################
		integral3 = bestfit->Integral(i, cumTSdisbinx, params);
		integral_bestfit->SetBinContent(jjj, integral3);
		
		
		cout << jjj << "\t" << h1->GetBinContent(jjj) <<  "\t" << h1->GetBinError(jjj) << "\t";
		cout << hPVALUE->GetBinContent(jjj) <<  "\t" << hPVALUE->GetBinError(jjj) << "\t";
		cout << hCDF->GetBinContent(jjj) <<  "\t" << hCDF->GetBinError(jjj);
		cout << "\t" << bestfit->Eval(jjj) << "\t0\n";
		//cout << integral_bestfit->GetBinCenter(jjj) << "\t";
		//cout << integral_bestfit->GetBinContent(jjj) <<  endl;
	}
	//cout << filenameinput << ".fit" << endl;
	//for(int i=0; i<cumTSdisbinx; i++) {
	//	int jjj=i+1;
	//	cout << jjj;		
	//	cout << "\t" << bestfit->Eval(jjj) << "\t0\t";
	//	cout << integral_bestfit->GetBinContent(jjj) <<  "\t0" << endl;
	//}

	double p;
	cout << "P(TS=h) and P(TS>=h) -> p-value" << endl;
	cout << bestfit->Eval(hhhh) << endl;
	p= bestfit->Integral(hhhh, cumTSdisbinx, params);
	cout << "p-value bestfit = " << p  << endl;
	
	//****************************
	cout << "--- Pre trial" << endl;
	double h = hhhh;
	printf("Probability of TS>%d: %6.12f (", h, p);
	cout << p << ") for single map " << endl;
	cout << "It correspond to " << TMath::ErfInverse(1-p)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << p << endl;
	Double_t postp = 1 - (1 - TMath::Power(1-p, n));
	
	printf("--- Post trial for %d maps  %6.12e\n", n, postp);
	
	cout << "It correspond to " << TMath::ErfInverse(postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << 1-postp << endl;
	
	double start = 0;
	double end = 100;
	TF1* f2sf1   = new TF1("f2sf1",    " [0] * exp(-x/2.0) * pow(x, [1]/2.0 - 1) / ( pow(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", start, end);
	TF1* f2sf1_05= new TF1("f2sf1_05", " [0] * exp(-x/2.0) * pow(x, [1]/2.0 - 1) / ( pow(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", start, end);
	TF1* f2sf3_05= new TF1("f2sf3_05", " [0] * exp(-x/2.0) * pow(x, [1]/2.0 - 1) / ( pow(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", start, end);
	
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
	
	Double_t ppp = 0;
	ppp=f2sf1->Integral(h, end); postp = 1 - (1 - TMath::Power(1-ppp, n));
	cout << "chi^2_1    \t" << f2sf1->Eval(h) << "\t " <<    ppp    << "\t" << TMath::ErfInverse(1-ppp)*TMath::Sqrt(2) << "\t" << postp << endl;
	ppp=f2sf1_05->Integral(h, end); postp = 1 - (1 - TMath::Power(1-ppp, n));
	cout << "1/2 chi^2_1\t" << f2sf1_05->Eval(h) << "\t " << ppp << "\t" << TMath::ErfInverse(1-ppp)*TMath::Sqrt(2) << "\t" << postp<< endl;
	ppp=f2sf3_05->Integral(h, end); postp = 1 - (1 - TMath::Power(1-ppp, n));
	cout << "1/2 chi^2_3\t" << f2sf3_05->Eval(h) << "\t " << ppp << "\t" << TMath::ErfInverse(1-ppp)*TMath::Sqrt(2) << "\t" << postp<< endl;

	//*****************************

	gmb = h4->GetMaximumBin();
//  	h4->Scale(1.0/h4->GetBinContent(gmb) * evalfirst2);
	h4->GetXaxis()->SetTitle("h");
	h4->GetYaxis()->SetTitle("P(TS >= h)");
	h4->SetLineWidth(2);
	h4->SetLineColor(4);
	
	hPVALUE->SetLineWidth(2);
	hPVALUE->SetLineColor(4);
	

	//AB: rimuovi questo se non vuoi il disegno dell'istogramma
	h4->GetYaxis()->SetTitleSize(0.03);
	h4->GetYaxis()->SetTitleOffset(1.80);
	h4->Draw("L");
	
	
	hPVALUE->Draw("SAMEE");
	
	gPad->SetLogy();

 	gmb = integral_chi2_1->GetMaximumBin();
//     	integral_chi2_1->Scale(1.0/integral_chi2_1->GetBinContent(gmb) );
	gmb = integral_chi2_1_05->GetMaximumBin();
//   	integral_chi2_1_05->Scale(1.0/integral_chi2_1_05->GetBinContent(gmb) * evalfirst2);
	integral_chi2_1->SetLineColor(kGreen+2);
	integral_chi2_1->SetLineWidth(1);
	integral_chi2_1->SetLineStyle(2);
	integral_chi2_1->Draw("SAMEL");
	
	integral_chi2_1_05->SetLineColor(kRed);
	integral_chi2_1_05->SetLineWidth(1);
	integral_chi2_1_05->SetLineStyle(3);
	integral_chi2_1_05->Draw("SAMEL");
	
	integral_bestfit->SetLineColor(kBlack);
	integral_bestfit->SetLineWidth(2);
	integral_bestfit->Draw("SAMEL");
	
	integral_chi2_3_05->SetLineColor(kCyan+3);
	integral_chi2_3_05->SetLineWidth(1);
	integral_chi2_3_05->SetLineStyle(5);
	integral_chi2_3_05->Draw("SAMEL");
	



	///////////////
	gmb = chi2_1_05->GetMaximumBin();


	*/

}


