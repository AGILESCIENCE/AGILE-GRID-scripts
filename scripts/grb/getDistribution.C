#include <iostream>
using namespace std;

#define MULTIBIN 1.0
#define MAXBIN 36

void getDistribution(TString filename) {
	gStyle->SetPalette(1);gStyle->SetOptStat(0000000);gStyle->SetOptFit(0000);
	
	Double_t dimWindowX = 1400;
	Double_t dimWindowY = 900;

	TTree* T;
	T = new TTree("DATA", "");
	nlines = T->ReadFile(filename, "STSTART/F:STSTOP/F:SCOUNTS/F:SEXP/F:SRATE/F:B1TSTART/F:B1TSTOP/F:B1COUNTS/F:B1EXP/F:B2TSTART/F:B2TSTOP/F:B2COUNTS/F:B2EXP/F:ALPHA/F:OFF/C:BCOUNTS/F:BEXP/F:BRATE/F:SIGMA/F");
	//415132181.0 415132211.0 2 67.69 0.0295443215 415131681.0 415132181.0 5 449.38 415132211.0 415132811.0 5 318.30 0.09  off 10 767.68 0.0130262827 0.9634100942
	
	cout << nlines << endl;
	Float_t SIGMA, SEXP, BEXP;
	//Int_t
	
	T->SetBranchAddress("SIGMA", &SIGMA);
	T->SetBranchAddress("SEXP", &SEXP);
	T->SetBranchAddress("BEXP", &BEXP);

	Int_t max_cumTSdisbinx = MAXBIN;
	Int_t cumTSdisbinx = max_cumTSdisbinx*MULTIBIN;
	
	TH1D* h1 = new TH1D("TS2", "SIGMA", cumTSdisbinx, 0, max_cumTSdisbinx );
	TH1D* hPVALUE = new TH1D("PVALUE", "PVALUE", cumTSdisbinx , 0, max_cumTSdisbinx ); //P(X>=x)
	TH1D* hCDF = new TH1D("CDF", "CDF", cumTSdisbinx , 0, max_cumTSdisbinx ); //P(X<=x)

	for(Long64_t j = 0; j<nlines; j++) {
        T->GetEntry(j);
        if(SEXP > 0 && BEXP > 0) {
        	//cout << SIGMA << endl;
        	h1->Fill(SIGMA);
        }
    }
    
    TH1F* h1cl = h1->Clone("TS2c");
	double scalefactor = h1->Integral();
	
	//calculate integral of TS
	cout << "integral" << endl;
	for(Long64_t i=0; i<cumTSdisbinx; i++) {
		Double_t integral;
		integral = h1->Integral(i, cumTSdisbinx);
		
		hPVALUE->SetBinContent(i, integral);
	}
	for(Long64_t i=1; i<cumTSdisbinx; i++) {
		Double_t integral;
		integral = h1->Integral(0, i);
		hCDF->SetBinContent(i+1, integral);
	}
	
	TH1F* hPVALUEcl = hPVALUE->Clone("PVALUEcl");
	double scalefactor_integral = hPVALUE->GetBinContent(1);
	//double scalefactor_integral = hPVALUE->Integral(0, cumTSdisbinx);
	hPVALUE->Scale(1.0/scalefactor_integral);
	TH1F* hCDFcl = hCDF->Clone("CDFcl");
	double scalefactor_integral2 = hCDF->GetBinContent(cumTSdisbinx-1);
	hCDF->Scale(1.0/scalefactor_integral2);
	
    h1->Scale(1.0/scalefactor);
	for(int ii1=1; ii1<=h1->GetNbinsX(); ii1++) {
		double ct2 = h1cl->GetBinContent(ii1);
		h1->SetBinError(ii1, TMath::Sqrt(ct2) * 1.0/scalefactor);
		double ct3 = hPVALUEcl->GetBinContent(ii1);
		hPVALUE->SetBinError(ii1, TMath::Sqrt(ct3) * 1.0/scalefactor);
		double ct4 = hCDFcl->GetBinContent(ii1);
		hCDF->SetBinError(ii1, TMath::Sqrt(ct4) * 1.0/scalefactor);
	}

	
    
	cout << "DISTRIBUTION ************************************" << endl;
	cout << "*************************************************" << endl;
	cout << "bin pdf pdf_error p-value p-value_error cdf cdf_error" << endl;
	for(Long64_t i=0; i<cumTSdisbinx; i++) {
		int jjj=i+1;
		
		cout << jjj-1 <<  "\t";
		cout << h1->GetBinContent(jjj) <<  "\t" << h1->GetBinError(jjj) << "\t";
		cout << hPVALUE->GetBinContent(jjj) <<  "\t" << hPVALUE->GetBinError(jjj) << "\t";
		cout << hCDF->GetBinContent(jjj) <<  "\t" << hCDF->GetBinError(jjj) << endl;
	}	
	cout << "*************************************************" << endl;

	hPVALUE->SetLineWidth(2);
	hPVALUE->SetLineColor(4);
	hCDF->SetLineWidth(2);
	hCDF->SetLineColor(kGreen+4);

	TCanvas* c1_h1 = new TCanvas("TS Distribution", "TS Distribution", dimWindowX, dimWindowY);
	gPad->SetLogy();
	c1_h1->Divide(2,1);
	c1_h1->cd(1);
	
	h1->SetLineColor(kBlue);
	h1->SetLineWidth(2);
	h1->GetYaxis()->SetTitleSize(0.03);
	h1->GetYaxis()->SetTitleOffset(1.80);
	h1->Draw("L");	
	c1_h1->cd(2);
	hPVALUE->Draw("SAMEE");
}
