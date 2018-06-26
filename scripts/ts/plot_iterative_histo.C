#include <iostream>
using namespace std;


/*
Parameters:
filename: name of the file that contain the result of the alike scan
catalog: name of a source catalog. This catalog is a text file with l b separated by a tab character
useTS: use TS instead of sqrt(TS)
titleExternal: an additional string for the titles of the histograms
outfileprefix: prefix for the output files

Formato file di input
oFile << srcData.label << "\t";
    oFile << srcData.srcL << "\t";
    oFile << srcData.srcB << "\t";
    oFile << srcData.TS << "\t";
    oFile << srcData.flux << "\t";
    oFile << srcData.index << "\t";
    oFile << srcData.fixflag << "\t";
    oFile << srcData.minTS << "\t";
    oFile << srcData.gal << "\t";
    oFile << srcData.iso << "\t";
    oFile <<  endl;
*/

/*
NOTA SUI VALORI DI TS NEL FILE TXT O NEL FITS
Se un bin vale TS=0 significa che non è stato analizzato. 
Se vale un numero piccolo prossimo a zero allora è stato analizzato e lo possiamo considerare 0. Lo contiamo nei counts < 0
Se vale un numero > limit allora è stato analizzato. Lo contiamo nei counts > 0
*/

#define MULTIBIN 1.0
#define MAXBIN 36

Double_t distance(double ll, double bl, double lf, double bf) {
	double d1 = bl - bf;
	double d2 = ll - lf;
	
	double bl1 = TMath::TwoPi() / 4.0 - (bl * TMath::TwoPi()  / 360.0);
	double bf1 = TMath::TwoPi() / 4.0 - (bf * TMath::TwoPi()  / 360.0);
	double m4 = TMath::Cos(bl1) * TMath::Cos(bf1)  + TMath::Sin(bl1) * TMath::Sin(bf1) * TMath::Cos(d2 * TMath::TwoPi()  / 360.0);                    
	double d4 = TMath::ACos(m4) *  360.0 / TMath::TwoPi();
	return d4;
	
}

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
Double_t sumchi2N1N2translN3deltafuncV2(Double_t *x, Double_t *par) {
	
	Double_t f=0;
	Double_t trasl2 = par[5]; 
	Double_t x1 = x[0];
	Double_t trasl3 = par[8]; 
	Double_t x2 = x[0]-trasl3;
	
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

//inputtype=5 output restituito da AG_multisim4
TTree* plot_iterative_histo(TString filenameinput, int inputtype=1, double enabledistsel = 0, Double_t eta = 0.5, Double_t limitFLUX = 0, Bool_t useTS = true,  TString titleExternal = "", TString outfileprefix = "") {
	
	//parametri di analisi da modificare
	//centro della mappa (caso empty field)
	//double centerl = 0;
	//double centerb = -90;
	double centerl = 79;
	double centerb = 0;
	//valore massimo di gal
	double galmax = 9999;
	
	cout << "Analysis parameters" << endl;
	cout << "center (l,b) of input map " << centerl << " " << centerb << " with max distance from center of " << enabledistsel << endl;
	cout << "GAL parameters selection. Max gal value " << galmax << endl;
	cout << "Eta vaule for plot: " << eta << endl;
	cout << "Limit flux selection: " << limitFLUX << endl;
	
	
	
	Bool_t saveRootFile = false;
// 	Double_t limitTS = 10e-3; //10e-4 - limite con cui definire lo 0
		
	//gStyle->SetPalette(1);gStyle->SetPadColor(0);gStyle->SetFrameFillColor(0);gStyle->SetCanvasColor(0);gStyle->SetFrameBorderMode(0);gStyle->SetOptStat(0000000);gStyle->SetOptFit(0000);gStyle->SetPalette(1);gStyle->SetLabelFont(42, "xyz");gStyle->SetTitleFont(42, "xyz");
	
	gStyle->SetPalette(1);gStyle->SetOptStat(0000000);gStyle->SetOptFit(0000);
	
	Double_t dimWindowX = 1400;
	Double_t dimWindowY = 900;


	
	TString titTS, titcumTS;
	if(useTS) {
		titTS= "TS distribution " + titleExternal;
		titcumTS = "P(TS >= h) ";
	} else {
		titTS = "sqr(TS) distribution " + titleExternal;
		titcumTS = "Cumulative sqr(TS) distribution ";
	}
	TH1D *h1, *hh1, *hh2, *hh3, *hh4;

	TH1D *h4;//se un bin vale 0 significa che non è stato analizzato. Se vale un numero piccolo prossimo a zero
	TH1D* chi2_1, *chi2_1_05, *chi2_3_05;
	TH1D* integral_chi2_1, *integral_chi2_1_05, *integral_chi2_3_05;
// 	Int_t nbinsall = 80;
	Double_t cumTSdisbinx, max_cumTSdisbinx ;
	if(useTS) {
		max_cumTSdisbinx = MAXBIN;
		cumTSdisbinx = max_cumTSdisbinx*MULTIBIN;
		
	} else {
		max_cumTSdisbinx = MAXBIN;
		cumTSdisbinx = max_cumTSdisbinx;
	}
	cout << "Number of bins " << cumTSdisbinx << " from 0 to " << max_cumTSdisbinx << endl;
	if(useTS)
		h1 = new TH1D("TS2", titTS, cumTSdisbinx, 0, max_cumTSdisbinx );
	else
		h1 = new TH1D("sqrTS2", titTS, cumTSdisbinx, 0, max_cumTSdisbinx );
// 	h3 = new TH1D("cumsqrTS2", titcumTS, cumTSdisbinx , 0, max_cumTSdisbinx );
	h4 = new TH1D("cumsqrTS4", titcumTS + titleExternal, cumTSdisbinx , 0, max_cumTSdisbinx );
	hPVALUE = new TH1D("PVALUE", "PVALUE", cumTSdisbinx , 0, max_cumTSdisbinx ); //P(X>=x)
	hCDF = new TH1D("CDF", "CDF", cumTSdisbinx , 0, max_cumTSdisbinx ); //P(X<=x)
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
	TH1D* derivative_chi2_1 = new TH1D("D(CHI2)", "D(CHI2)", cumTSdisbinx , 0, max_cumTSdisbinx );
	TH1D* derivative_chi2_1_05 = new TH1D("eta * D(CHI2)", "eta * D(CHI2)", cumTSdisbinx , 0, max_cumTSdisbinx );
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
	Int_t numberofTSgood = 0;
	Int_t numberofFLUXgood = 0;
	Float_t SOURCE, L, B, TS, FLUX, FLUXERR, FLUXUL, GAL, ISO, R, EXP, SI, MINTS, FIXFLAG , CTS, CTSERR, CTSUL, TOTEXP, TOTNCOUNTS, FCN0, FCN1, EDM0, EDM1, ITER0, ITER1;
	L = B = TS = FLUX = FLUXERR = FLUXUL = GAL = ISO = R = EXP = SI = MINTS = FIXFLAG = CTS = CTSERR = CTSUL = TOTEXP = TOTNCOUNTS = FCN0 = FCN1 = EDM0 =  EDM1 = ITER0 =  ITER1 = 1;

	
	
	Long64_t nlines;
	TTree* T;
	
	TFile* fin = new TFile();
	int resfile = fin->Open( filenameinput, "read");
	cout << resfile << endl;
	
	if(resfile == 0) {
		T = new TTree("DATA", "") ;
		if(inputtype == 1) {
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS:GAL:ISO:R");
			cout << "******************************************* " << endl;
			cout << "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS:GAL:ISO:R" << endl;
			cout << "******************************************* " << endl;
		}
		if(inputtype == 2) {
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:EXP");
			cout << "******************************************* " << endl;
			cout << "SOURCE:L:B:TS:EXP" << endl;
			cout << "******************************************* " << endl;
		}
		//nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS:GAL:ISO");
		//nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS");
		if(inputtype == 3) {
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS");
			cout << "******************************************* " << endl;
			cout << "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS" << endl;
			cout << "******************************************* " << endl;
		}
		if(inputtype == 4) {
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:SPECTRALINDEX:FIXFLAG:MINTS:GAL:ISO:R");
			cout << "******************************************* " << endl;
			cout << "SOURCE:L:B:TS:FLUX:SPECTRALINDEX:FIXFLAG:MINTS:GAL:ISO:R" << endl;
			cout << "******************************************* " << endl;
		}
		//output restituito da AG_multisim4
		if(inputtype == 5) {
			/*
			After R:
			double exposure = GetTotalExposure(i);
			srcout << " " << exposure*m_sources[i].GetFlux();
			srcout << " " << exposure*m_sources[i].GetFluxerr();
			srcout << " " << exposure*m_sources[i].GetFluxul();
			srcout << " " << exposure;
			srcout << " " << m_fitInfo[i].counts;

			srcout << " " << m_fitInfo[i].fcn0;
			srcout << " " << m_fitInfo[i].fcn1;
			srcout << " " << m_fitInfo[i].edm0;
			srcout << " " << m_fitInfo[i].edm1;
			srcout << " " << m_fitInfo[i].iter0;
			srcout << " " << m_fitInfo[i].iter1 << endl;
			*/
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS:GAL:ISO:R:CTS:CTSERR:CTSUL:TOTEXP:TOTNCOUNTS:FCN0:FCN1:EDM0:EDM1:ITER0:ITER1");
			cout << "******************************************* " << endl;
			cout << "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS:GAL:ISO:R:CTS:CTSERR:CTSUL:TOTEXP:TOTNCOUNTS:FCN0:FCN1:EDM0:EDM1:ITER0:ITER1" << endl;
			cout << "******************************************* " << endl;
		}
		if(inputtype == 6) {
			//NUMIT L B TS FLUX FLUXERR FLUXUL SPECTRAL_INDEX FIXFLAG MINTS R EXP CTS CTSERROR CTSUL TOTEXP TOTNCOUNTS FCN0 FCN1 EDM0 EDM1 ITER0 ITER1 GAL ISO
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:FLUXERR:FLUXUL:INDEX:FIXFLAG:MINTS:R:EXP:CTS:CTSERR:CTSUL:TOTEXP:TOTNCOUNTS:FCN0:FCN1:EDM0:EDM1:ITER0:ITER1:GAL:ISO");
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:FLUXERR:FLUXUL:INDEX:FIXFLAG:MINTS:R:EXP:CTS:CTSERR:CTSUL:TOTEXP:TOTNCOUNTS:FCN0:FCN1:EDM0:EDM1:ITER0:ITER1");
			cout << "******************************************* " << endl;
			cout << "SOURCE:L:B:TS:FLUX:FLUXERR:FLUXUL:INDEX:FIXFLAG:MINTS:R:EXP:CTS:CTSERR:CTSUL:TOTEXP:TOTNCOUNTS:FCN0:FCN1:EDM0:EDM1:ITER0:ITER1:GAL:ISO" << endl;
			cout << "******************************************* " << endl;
		}
	} else {
		T = DATA;
		cout << T << endl;
		nlines = T->GetEntry();
	}
	
	cout << "N lines: " << nlines << endl;
	
	TFile* fout = 0;
	/*if(resfile == 0) {
		cout << "Save the tree" << endl;
		TString fnout = filenameinput;
		fnout += ".root";
		TFile* fout = new TFile(fnout, "recreate");
		T->Write();
	}*/
	
	if(nlines == 0) {
		cout << "End of the procedure" << endl;
		return;
	}
	if(inputtype == 1) {
		T->SetBranchAddress("SOURCE", &SOURCE);
		T->SetBranchAddress("L", &L);
		T->SetBranchAddress("B", &B);
		T->SetBranchAddress("TS", &TS);
		T->SetBranchAddress("FLUX", &FLUX);
		T->SetBranchAddress("GAL", &GAL);
		T->SetBranchAddress("ISO", &ISO);
		T->SetBranchAddress("R", &R);
	}
	if(inputtype == 2) {
		T->SetBranchAddress("SOURCE", &SOURCE);
		T->SetBranchAddress("L", &L);
		T->SetBranchAddress("B", &B);
		T->SetBranchAddress("TS", &TS);
		T->SetBranchAddress("EXP", &EXP);
	}
	if(inputtype == 3) {
		T->SetBranchAddress("SOURCE", &SOURCE);
		T->SetBranchAddress("L", &L);
		T->SetBranchAddress("B", &B);
		T->SetBranchAddress("TS", &TS);
		T->SetBranchAddress("FLUX", &FLUX);
	}
	if(inputtype == 4) {
		T->SetBranchAddress("SOURCE", &SOURCE);
		T->SetBranchAddress("L", &L);
		T->SetBranchAddress("B", &B);
		T->SetBranchAddress("TS", &TS);
		T->SetBranchAddress("FLUX", &FLUX);
		T->SetBranchAddress("SPECTRALINDEX", &SI);
		T->SetBranchAddress("FIXFLAG", &FIXFLAG);
		T->SetBranchAddress("MINTS", &MINTS);
		T->SetBranchAddress("GAL", &GAL);
		T->SetBranchAddress("ISO", &ISO);
		T->SetBranchAddress("R", &R);
	}
	if(inputtype == 5 || inputtype == 6) {
		T->SetBranchAddress("SOURCE", &SOURCE);
		T->SetBranchAddress("L", &L);
		T->SetBranchAddress("B", &B);
		T->SetBranchAddress("TS", &TS);
		T->SetBranchAddress("FLUX", &FLUX);
		T->SetBranchAddress("R", &R);
		T->SetBranchAddress("CTS", &CTS);
		T->SetBranchAddress("TOTNCOUNTS", &TOTNCOUNTS);
		T->SetBranchAddress("FCN0", &FCN0);
		T->SetBranchAddress("FCN1", &FCN1);
		T->SetBranchAddress("EDM0", &EDM0);
		T->SetBranchAddress("EDM1", &EDM1);
		T->SetBranchAddress("ITER0", &ITER0);
		T->SetBranchAddress("ITER1", &ITER1);
		
		//T->SetBranchAddress("GAL", &GAL);
		//T->SetBranchAddress("ISO", &ISO);
	}
	Int_t index = 0;
	
	hh_FLUX = new TH1D("H_FLUX", "H_FLUX", 100, 0, 1000 );
	hh_ISO = new TH1D("H_ISO", "H_ISO", 200, 0, 50 );
	hh_GAL = new TH1D("H_GAL", "H_GAL", 200, 0, 5 );
	
	//hh_FCN0_L = new TH1D("FCN0_L", "FCN0_L");
	//hh_FCN0_H = new TH1D("FCN0_H", "FCN0_H");
	/*TCanvas* cc = new TCanvas("Fit parameters");
	gPad->Divide(2,2);
	gPad->cd(1);
	T->Draw("FCN0 >> hfcn");
	gPad->cd(2);
	T->Draw("FCN0 >> hfcn20", "TS>20");
	*/
	
	TH2D* hhISOGAL = new TH2D("HGALISO", "HGALISO",160, 0, 20, 400, 0, 50);
	hhISOGAL->GetXaxis()->SetTitle("GAL");
	hhISOGAL->GetYaxis()->SetTitle("ISO");
	
	if(enabledistsel > 0)
		cout << "!Distance selection: enabled=" << enabledistsel <<  " - l: " << centerl << " b:" << centerb << " distmax " << enabledistsel << endl;
	//cout << "!GAL selection enabled" << endl;
	long nentries = 0;
	//nlines = 100;
	
	TH2D* hlb = new TH2D("lb", "lb", 200, centerl-10, centerl+10, 200, centerb-10, centerb+10);
	hlb->GetXaxis()->SetTitle("l");
	hlb->GetYaxis()->SetTitle("b");
	
	TH2D* htsr = new TH2D("tsr", "tsr", cumTSdisbinx , 0, max_cumTSdisbinx, 30, 0, 3);
	htsr->GetXaxis()->SetTitle("ts");
	htsr->GetYaxis()->SetTitle("r");
	TH2D* htsd = new TH2D("tsd", "tsd", cumTSdisbinx , 0, max_cumTSdisbinx, 60, 0, 6);
	htsd->GetXaxis()->SetTitle("ts");
	htsd->GetYaxis()->SetTitle("d");
	TH2D* hrd = new TH2D("dr", "dr", 30 , 0, 3, 30, 0, 3);
	hrd->GetXaxis()->SetTitle("d");
	hrd->GetYaxis()->SetTitle("r");
	

	/*TH1D* hl = new TH1D("l", "l", 200, 79.9031-10, 79.9031+10);
	hl->GetXaxis()->SetTitle("l");
	TH1D* hb = new TH1D("b", "b", 200, 0.558037-10, 0.558037+10);
	hl->GetXaxis()->SetTitle("b");
	*/
	
	TH1D* hl = new TH1D("l", "l", 20, centerl-2, centerl+2);
	hl->GetXaxis()->SetTitle("l");
	TH1D* hb = new TH1D("b", "b", 20, centerb-2, centerb+2);
	hl->GetXaxis()->SetTitle("b");
	
	TH1D* rsel = new TH1D("rsel", "rsel", 2, 0,2);
	rsel->GetXaxis()->SetTitle("(r+0.1)<dist");
	
	TH2D* htsflux = new TH2D("tsflux", "tsflux", 36, 0, 36, 100, 0, 1000);
	htsflux->GetXaxis()->SetTitle("TS");
	htsflux->GetYaxis()->SetTitle("Flux");
	int numHTS = 0;
	
	int countsremove[10]; 
	for(int cri=0; cri<10; cri++) countsremove[cri]=0;
	
	int discarded = 0;
	for(Long64_t j = 0; j<nlines; j++) {
		T->GetEntry(j);
		
		
		if(EDM0 != 0.5) {
			discarded++;
			cout << "* " << SOURCE << " " << L << " " << B  << " " << (FLUX) << " " << TS << " " << EDM0 << " " << EDM1 << endl;
			continue;
		}
		
		
		if(inputtype == 2) {
			TS = TS * TS; //BE CAREFULL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(EXP == 0) continue;
		}
		
			
		//if(TS < 12)	continue;
		if(TS > 0.0) {
			numberofTSgood++;
		}
		if(FLUX >= limitFLUX) {
			numberofFLUXgood++;
			
		}
		
		//fill the histogram
		
		//if(GAL <= galmax) {
		
			double dist = distance(L, B, centerl, centerb);
			if(enabledistsel > 0) {
				
				if(dist < enabledistsel) {
					if(TS>25) {
						//continue;
						cout << SOURCE << " " << L << " " << B << " " << dist << " " << (FLUX) << " " << TS << " " << GAL << " " << ISO << endl;
						numHTS++;
					}
					//h1->Fill((useTS?TS:TMath::Sqrt(TS)));
					h1->Fill(2*(FCN0-FCN1));
					
					hh_FLUX->Fill(FLUX / 1e-08);
					hh_GAL->Fill(GAL);
					hh_ISO->Fill(ISO);
					hhISOGAL->Fill(GAL, ISO);
					hlb->Fill(L, B);
					hl->Fill(L);
					hb->Fill(B);
					htsflux->Fill((useTS?TS:TMath::Sqrt(TS)), FLUX / 1e-08);
					if(R >0) {
						htsr->Fill(TS, R);
						hrd->Fill(dist, R);
						if( (R+0.1)>dist )
							rsel->Fill(1); //associ ok
						else {
							rsel->Fill(0);
						}

					}
					htsd->Fill(TS, dist);
					
					nentries++;
				}
			} else {
				if(TS>25) {
					//continue;
					cout << SOURCE << " " << L << " " << B << " NODIST " << FLUX << " " << TS << " " << GAL << " " << ISO << endl;
					numHTS++;
				}
				/*if(TS>30 && TS<=31 && countsremove[0] <=14) {
					cout << "TS>30 && TS<=31 " << endl;
					countsremove[0]++;
					continue;
				}
				if(TS>31 && TS<=32 && countsremove[1] <=1) {
					cout << "TS>31 && TS<=32 " << endl;
					countsremove[1]++;
					continue;
				}
				if(TS>32 && TS<=33 && countsremove[1] <=2) {
					cout << "TS>32 && TS<=33 " << endl;
					countsremove[1]++;
					continue;
				}
				if(TS>34 && TS<=34 && countsremove[1] <=2) {
					cout << "TS>33 && TS<=34 " << endl;
					countsremove[1]++;
					continue;
				}
				if(TS>34 && TS<=35 && countsremove[2] <=1) {
					cout << "TS>34 && TS<=35 " << endl;
					countsremove[2]++;
					continue;
				}
				if(TS>35 && TS<=36 && countsremove[3] <=1) {
					cout << "TS>35 && TS<=36 " << endl;
					countsremove[3]++;
					continue;
				}
				*/
				h1->Fill((useTS?TS:TMath::Sqrt(TS)));
				
				hh_FLUX->Fill(FLUX / 1e-08);
				hh_GAL->Fill(GAL);
				hh_ISO->Fill(ISO);
				hhISOGAL->Fill(GAL, ISO);
				hlb->Fill(L, B);
				hl->Fill(L);
				hb->Fill(B);
				htsflux->Fill((useTS?TS:TMath::Sqrt(TS)), FLUX / 1e-08);
				if(R >0) {
					htsr->Fill(TS, R);
					hrd->Fill(dist, R);
					if( (R+0.1)>dist )
						rsel->Fill(1); //ass ok
					else {
						rsel->Fill(0);
					}
				}
				htsd->Fill(TS, dist);
				nentries++;
			}
		//}
	}
	
	
	
	cout << "nentries: " << nentries << endl;
	cout << "numHTS: " << numHTS << endl;
	TH1F* h1cl = h1->Clone("TS2c");
	
	new TCanvas("TRIALS");
	h1cl->Draw("");
	new TCanvas();
	
	
	double scalefactor = h1->Integral();
	
	/*
	cout << "TS distribution value" << endl;
	for(int ii1=0; ii1<=h1->GetNbinsX(); ii1++) {
		cout << ii1 << " " << h1->GetBinCenter(ii1) << " " << h1->GetBinContent(ii1) << endl;
	}
	*/
	
	//calculate integral of TS
	cout << "integral" << endl;
	for(Long64_t i=0; i<cumTSdisbinx; i++) {
		Double_t integral;
		integral = h1->Integral(i, cumTSdisbinx);
		
		hPVALUE->SetBinContent(i, integral);
		//hPVALUE->SetBinError(i, TMath::Sqrt(integral) * 1.0/scalefactor);
	//	cout << i << " " << integral << " " << TMath::Sqrt(integral) * 1.0/scalefactor << endl;
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
	
	//calcolo errore nei bin
	h1->Scale(1.0/scalefactor);
	for(int ii1=1; ii1<=h1->GetNbinsX(); ii1++) {
		double ct2 = h1cl->GetBinContent(ii1);
		h1->SetBinError(ii1, TMath::Sqrt(ct2) * 1.0/scalefactor);
		double ct3 = hPVALUEcl->GetBinContent(ii1);
		hPVALUE->SetBinError(ii1, TMath::Sqrt(ct3) * 1.0/scalefactor);
		double ct4 = hCDFcl->GetBinContent(ii1);
		hCDF->SetBinError(ii1, TMath::Sqrt(ct4) * 1.0/scalefactor);
	}
	//TH1D* hhhhisto = hPVALUE;
	//TH1D* hhhhisto = h1;
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
	
	cout << "MAP FITS: Number of FLUX > " << limitFLUX << "  are " << numberofFLUXgood << endl;
	//****************
	//**************** 
	//****************
	
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
	Int_t f4sf_test_N1=1;
	//f4sf_trasl_N1N2v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1);
	f4sf_trasl_N1N2v2->SetParameter( 2, f4sf_test_N1);
	f4sf_trasl_N1N2v2->SetParName(3, "eta2");
	f4sf_trasl_N1N2v2->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2v2->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2v2->SetParName(4, "N2");
	Int_t f4sf_test_N2=3;
	//f4sf_trasl_N1N2v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2);
	f4sf_trasl_N1N2v2->SetParameter( 4, f4sf_test_N2);
	Int_t f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_trasl_N1N2v2->SetParName(5, "translation2");
	f4sf_trasl_N1N2v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2v2->SetLineColor(kBlack);
	
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	TF1* f4sf_trasl_N1N2v2_F = new TF1("traslchi2N1N2deltafuncV2_ALLFREE", traslchi2N1N2deltafuncV2, 0, max_cumTSdisbinx,6);
	f4sf_trasl_N1N2v2_F->SetParName(0, "delta");
	f4sf_trasl_N1N2v2_F->SetParName(1, "eta1");
	f4sf_trasl_N1N2v2_F->SetParLimits(1, 0, 1);
	f4sf_trasl_N1N2v2_F->SetParameter( 1, 0.5);
	f4sf_trasl_N1N2v2_F->SetParName(2, "N1");
	f4sf_test_N1=1;
	f4sf_trasl_N1N2v2_F->SetParLimits(2, 0, 5);
	//f4sf_trasl_N1N2v2_F->SetParLimits(2, f4sf_test_N1, f4sf_test_N1); //!!!!!
	f4sf_trasl_N1N2v2_F->SetParameter( 2, f4sf_test_N1);
	f4sf_trasl_N1N2v2_F->SetParName(3, "eta2");
	f4sf_trasl_N1N2v2_F->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2v2_F->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2v2_F->SetParName(4, "N2");
	f4sf_test_N2=5;
	f4sf_trasl_N1N2v2_F->SetParLimits(4, 0, 5);
	//f4sf_trasl_N1N2v2_F->SetParLimits(4, f4sf_test_N2, f4sf_test_N2); //!!!!!!
	f4sf_trasl_N1N2v2_F->SetParameter( 4, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_trasl_N1N2v2_F->SetParName(5, "translation2");
	f4sf_trasl_N1N2v2_F->SetParLimits(5, 0, 36);
	f4sf_trasl_N1N2v2_F->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl); //!!!!!
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
	f4sf_trasl_N1N2N3v2_F->SetParName(3, "eta2");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2N3v2_F->SetParName(4, "N2");
	f4sf_test_N2=4;
	f4sf_trasl_N1N2N3v2_F->SetParLimits(4, 0, 5);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 4, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_trasl_N1N2N3v2_F->SetParName(5, "translation2");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(5, 0, 36);
	f4sf_trasl_N1N2N3v2_F->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2_F->SetParName(6, "eta3");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(6, 0, 1);
	f4sf_trasl_N1N2N3v2_F->SetParameter( 6, 0.5);
	f4sf_trasl_N1N2N3v2_F->SetParName(7, "N3");
	f4sf_test_N2=4;
	f4sf_trasl_N1N2N3v2_F->SetParLimits(7, 0, 5);
	f4sf_trasl_N1N2N3v2_F->SetParameter(7, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 6;
	f4sf_trasl_N1N2N3v2_F->SetParName(8, "translation3");
	f4sf_trasl_N1N2N3v2_F->SetParLimits(8, 0, 36);
	f4sf_trasl_N1N2N3v2_F->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	
	f4sf_trasl_N1N2N3v2_F->SetLineColor(kBlack);

	
	TF1* f4sf_trasl_N1N2N3v2 = new TF1("traslchi2N1N2N3deltafuncV2", traslchi2N1N2N3deltafuncV2, 0, max_cumTSdisbinx,9);
	f4sf_trasl_N1N2N3v2->SetParName(0, "delta");
	f4sf_trasl_N1N2N3v2->SetParName(1, "eta1");
	f4sf_trasl_N1N2N3v2->SetParLimits(1, 0, 1);
	f4sf_trasl_N1N2N3v2->SetParameter( 1, 0.5);
	f4sf_trasl_N1N2N3v2->SetParName(2, "N1");
	f4sf_test_N1=3;
	f4sf_trasl_N1N2N3v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1);
	f4sf_trasl_N1N2N3v2->SetParameter( 2, f4sf_test_N1);
	f4sf_trasl_N1N2N3v2->SetParName(3, "eta2");
	f4sf_trasl_N1N2N3v2->SetParLimits(3, 0, 1);
	f4sf_trasl_N1N2N3v2->SetParameter( 3, 0.5);
	f4sf_trasl_N1N2N3v2->SetParName(4, "N2");
	f4sf_test_N2=3;
	f4sf_trasl_N1N2N3v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2);
	f4sf_trasl_N1N2N3v2->SetParameter( 4, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 4;
	f4sf_trasl_N1N2N3v2->SetParName(5, "translation2");
	f4sf_trasl_N1N2N3v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2->SetParName(6, "eta3");
	f4sf_trasl_N1N2N3v2->SetParLimits(6, 0, 1);
	f4sf_trasl_N1N2N3v2->SetParameter( 6, 0.5);
	f4sf_trasl_N1N2N3v2->SetParName(7, "N3");
	f4sf_test_N2=1;
	f4sf_trasl_N1N2N3v2->SetParLimits(7, f4sf_test_N2, f4sf_test_N2);
	f4sf_trasl_N1N2N3v2->SetParameter(7, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 9;
	f4sf_trasl_N1N2N3v2->SetParName(8, "translation3");
	f4sf_trasl_N1N2N3v2->SetParLimits(8, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_trasl_N1N2N3v2->SetParameter(8, f4sf_trasl_N1N2v2_trasl);
	
	f4sf_trasl_N1N2N3v2->SetLineColor(kBlack);
	
	
	
	TF1* f4sf_sumN1N2traslN3v2 = new TF1("sumchi2N1N2translN3deltafuncV2", sumchi2N1N2translN3deltafuncV2, 0, max_cumTSdisbinx,9);
	f4sf_sumN1N2traslN3v2->SetParName(0, "delta");
	f4sf_sumN1N2traslN3v2->SetParName(1, "eta1");
	f4sf_sumN1N2traslN3v2->SetParLimits(1, 0, 1);
	f4sf_sumN1N2traslN3v2->SetParameter( 1, 0.5);
	f4sf_sumN1N2traslN3v2->SetParName(2, "N1");
	f4sf_test_N1=3;
	f4sf_sumN1N2traslN3v2->SetParLimits(2, f4sf_test_N1, f4sf_test_N1);
	f4sf_sumN1N2traslN3v2->SetParameter( 2, f4sf_test_N1);
	f4sf_sumN1N2traslN3v2->SetParName(3, "eta2");
	f4sf_sumN1N2traslN3v2->SetParLimits(3, 0, 1);
	f4sf_sumN1N2traslN3v2->SetParameter( 3, 0.5);
	f4sf_sumN1N2traslN3v2->SetParName(4, "N2");
	f4sf_test_N2=4;
	f4sf_sumN1N2traslN3v2->SetParLimits(4, f4sf_test_N2, f4sf_test_N2);
	f4sf_sumN1N2traslN3v2->SetParameter( 4, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 4;
	f4sf_sumN1N2traslN3v2->SetParName(5, "translation2");
	f4sf_sumN1N2traslN3v2->SetParLimits(5, f4sf_trasl_N1N2v2_trasl, f4sf_trasl_N1N2v2_trasl);
	f4sf_sumN1N2traslN3v2->SetParameter(5, f4sf_trasl_N1N2v2_trasl);
	f4sf_sumN1N2traslN3v2->SetParName(6, "eta3");
	f4sf_sumN1N2traslN3v2->SetParLimits(6, 0, 1);
	f4sf_sumN1N2traslN3v2->SetParameter( 6, 0.5);
	f4sf_sumN1N2traslN3v2->SetParName(7, "N3");
	f4sf_test_N2=1;
	f4sf_sumN1N2traslN3v2->SetParLimits(7, f4sf_test_N2, f4sf_test_N2);
	f4sf_sumN1N2traslN3v2->SetParameter(7, f4sf_test_N2);
	f4sf_trasl_N1N2v2_trasl = 9;
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
	Int_t f4sfN1N2_N1 = 1;
	f4sfN1N2->SetParName(3, "N1");
	f4sfN1N2->SetParLimits(3, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2->SetParameter(3, f4sfN1N2_N1);
	f4sfN1N2->SetParName(4, "eta2");
	f4sfN1N2->SetParLimits(4, 0, 1);
	f4sfN1N2->SetParameter( 4, 0.5);
	Int_t f4sfN1N2_N2 = 1;
	f4sfN1N2->SetParName(5, "N2");
	f4sfN1N2->SetParLimits(5, f4sfN1N2_N2, f4sfN1N2_N2);
	f4sfN1N2->SetParameter(5, f4sfN1N2_N2);
	f4sfN1N2->SetLineColor(kBlack);

	TF1* f4sfN1N2v2 = new TF1("sumchi2(N1)(N2)deltaV2", sumchi2N1N2deltafuncV2, 0, max_cumTSdisbinx,6);
	f4sfN1N2v2->SetParName(0, "delta");
	f4sfN1N2v2->SetParLimits(0, 0.0, 1);
	//f4sfN1N2->SetParameter( 0, 0.7);
	constloclc=4;
	f4sfN1N2v2->SetParName(1, "loccl");
	f4sfN1N2v2->SetParLimits(1, constloclc,constloclc);
	f4sfN1N2v2->SetParameter( 1, constloclc);
	
	f4sfN1N2v2->SetParName(2, "eta1");
	f4sfN1N2v2->SetParLimits(2, 0, 1);
	f4sfN1N2v2->SetParameter( 2, 0.5);
	
	f4sfN1N2_N1 = 3;
	f4sfN1N2v2->SetParName(3, "N1");
	f4sfN1N2v2->SetParLimits(3, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2v2->SetParameter(3, f4sfN1N2_N1);
	
	f4sfN1N2v2->SetParName(4, "eta2");
	f4sfN1N2v2->SetParLimits(4, 0, 1);
	f4sfN1N2v2->SetParameter( 4, 0.5);
	
	f4sfN1N2_N2 = 3;
	f4sfN1N2v2->SetParName(5, "N2");
	f4sfN1N2v2->SetParLimits(5, f4sfN1N2_N2, f4sfN1N2_N2);
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
	
	f4sfN1N2_N1 = 3;
	f4sfN1N2N3v2->SetParName(3, "N1");
	f4sfN1N2N3v2->SetParLimits(3, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2N3v2->SetParameter(3, f4sfN1N2_N1);
	
	f4sfN1N2N3v2->SetParName(4, "eta2");
	f4sfN1N2N3v2->SetParLimits(4, 0, 1);
	f4sfN1N2N3v2->SetParameter( 4, 0.5);
	
	f4sfN1N2_N2 = 4;
	f4sfN1N2N3v2->SetParName(5, "N2");
	f4sfN1N2N3v2->SetParLimits(5, f4sfN1N2_N2, f4sfN1N2_N2);
	f4sfN1N2N3v2->SetParameter(5, f4sfN1N2_N2);
	
	constloclc=9;
	f4sfN1N2N3v2->SetParName(6, "thr2");
	f4sfN1N2N3v2->SetParLimits(6, constloclc,constloclc);
	f4sfN1N2N3v2->SetParameter(6, constloclc);
	
	
	f4sfN1N2N3v2->SetParName(7, "eta3");
	f4sfN1N2N3v2->SetParLimits(7, 0, 1);
	f4sfN1N2N3v2->SetParameter(7, 0.5);
	
	f4sfN1N2_N1 = 1;
	f4sfN1N2N3v2->SetParName(8, "N3");
	f4sfN1N2N3v2->SetParLimits(8, f4sfN1N2_N1, f4sfN1N2_N1);
	f4sfN1N2N3v2->SetParameter(8, f4sfN1N2_N1);
	
	f4sfN1N2N3v2->SetLineColor(kBlack);
	
	
	
	
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
	//h1->SaveAs("TS2.root");
	//**************************
	TF1* bestfit;
	
	
	cout << "####################### fchi2_1fit " << endl;
	//h1->Fit("fchi2_1fit", "+", "", 0, max_cumTSdisbinx);
	cout << "####################### fchi2_Nfit " << endl;
	//h1->Fit("fchi2_Nfit", "+", "", 0, max_cumTSdisbinx);
	cout << "####################### delta " << endl;
	//h1->Fit("delta", "+", "", 0, 0.1);
	
	
	cout << "####################### sumchi2delta " << endl;
	//h1->Fit("sumchi2delta", "+");
	cout << "####################### sumchi2(N)delta " << endl;
	//h1->Fit("sumchi2(N)delta", "+");
	cout << "####################### sumchi2(1)delta " << endl;
	//h1->Fit("sumchi2(1)delta", "+");
	cout << "####################### sumchi2(2)delta " << endl;
	//h1->Fit("sumchi2(2)delta", "+");
	cout << "####################### sumchi2(3)delta " << endl;
	//h1->Fit("sumchi2(3)delta", "+");
	
	cout << "####################### sumchi2(N)deltaV2 " << endl;
	//h1->Fit("sumchi2(N)deltaV2", "+");
	//bestfit = f4sfv2;
	
	cout << "####################### sumchi2(1)deltaV2 " << endl;
	//h1->Fit("sumchi2(1)deltaV2", "+");
	//bestfit = f4sf1v2;
	
	cout << "####################### sumchi2(2)deltaV2 " << endl;
	//h1->Fit("sumchi2(2)deltaV2", "+");
	//bestfit = f4sf2v2;
	
	cout << "####################### sumchi2(3)deltaV2 " << endl;
	//h1->Fit("sumchi2(3)deltaV2", "+");
	//bestfit = f4sf3v2;
	
	cout << "####################### sumchi2(4)deltaV2 " << endl;
	//h1->Fit("sumchi2(4)deltaV2", "+");
	//bestfit = f4sf4v2;
	
	cout << "####################### sumchi2(N1)(N2)delta " << endl;
	//h1->Fit("sumchi2(N1)(N2)delta", "+");
	//bestfit = f4sfN1N2;
	
	cout << "####################### sumchi2(N1)(N2)deltaV2 " << endl;
	//h1->Fit("sumchi2(N1)(N2)deltaV2", "+");
	//bestfit = f4sfN1N2v2;
	
	cout << "####################### sumchi2(N1)(N2)(N3)deltaV2 " << endl;
	//h1->Fit("sumchi2(N1)(N2)(N3)deltaV2", "+");
	//bestfit = f4sfN1N2N3v2;
	
	//loccl=95
	cout << "####################### traslchi2N1N2deltafuncV2 " << endl;
	h1->Fit("traslchi2N1N2deltafuncV2", "+");
	bestfit = f4sf_trasl_N1N2v2;
	
	cout << "####################### traslchi2N1N2deltafuncV2_ALLFREE " << endl;
	//h1->Fit("traslchi2N1N2deltafuncV2_ALLFREE", "+");
	//bestfit = f4sf_trasl_N1N2v2_F;
	
	cout << "####################### traslchi2N1N2N3deltafuncV2" << endl;
	//h1->Fit("traslchi2N1N2N3deltafuncV2", "+");
	//bestfit = f4sf_trasl_N1N2N3v2_F;		
	
	cout << "####################### traslchi2N1N2N3deltafuncV2_ALLFREE " << endl;
	//h1->Fit("traslchi2N1N2N3deltafuncV2_ALLFREE", "+");
	//bestfit = f4sf_trasl_N1N2N3v2;	
	
	cout << "####################### sumchi2N1N2translN3deltafuncV2 " << endl;
	//h1->Fit("sumchi2N1N2translN3deltafuncV2", "+");
	//bestfit = f4sf_sumN1N2traslN3v2;
	
	
	double params[10];
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
	
	
	for(Long64_t i=1; i<cumTSdisbinx; i++) {
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
	c1_h1->Divide(2,1);
	c1_h1->cd(1);
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

	

	//****************
	//**************** calcolo dell'errore tra le distribuzione TS reale e quella teorica
	//****************
	/*
	TCanvas* c1_error = new TCanvas("TS Distribution error", "TS Distribution error", dimWindowX, dimWindowY);
	TH1D* TSerror;
	Int_t bbbin = 0;
	if(useTS)
		TSerror = new TH1D("TS ERROR", "TS ERROR", bbbin=(40*MULTIBIN/4), 0, 10*MULTIBIN );
	else
		TSerror = new TH1D("TS ERROR", "TS ERROR", bbbin=(40/4), 0, 10 );
	for(Long64_t i=2; i<bbbin; i++) {
		Double_t a1 = bestfit->Eval(h1->GetBinCenter(i));
		//cout << h1->GetBinCenter(i) << " " << a1 << endl;
		Double_t a2 = h1->GetBinContent(i);
		TSerror->SetBinContent(i, TMath::Abs(a2-a1));
		
	}
	TSerror->Draw();
	TSerror->Fit("pol0");
	TF1 *f1_error = TSerror->GetFunction("pol0");
	Double_t Errors = f1_error->GetParameter(0);
	*/

	Double_t gmb;
	//TCanvas* c2 = new TCanvas("Cumulative TS distribution", "Cumulative TS distribution", dimWindowX, dimWindowY);
	c1_h1->cd(2);
	//cout << "chi2_3 reference: " << endl;
	//f2sfN->Print();
	cout << "(TS) (bin center) (best fit integral) (P(TS>=h) integral) (chi2_1 int) (chi2_1/2 int) (chi2_3/2 int)" << endl;
	for(Long64_t i=0; i<cumTSdisbinx; i++) {
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
		
		cout << i << " " << integral_bestfit->GetBinCenter(jjj) << " ";
		cout << integral_bestfit->GetBinContent(jjj) <<  " " << h4->GetBinContent(i) << " ";
		cout << integral_chi2_1->GetBinContent(jjj) << " " <<  integral_chi2_1_05->GetBinContent(jjj) << " " << integral_chi2_3_05->GetBinContent(jjj) <<  endl;
	}

	

	gmb = h4->GetMaximumBin();
//  	h4->Scale(1.0/h4->GetBinContent(gmb) * evalfirst2);
	h4->GetXaxis()->SetTitle("h");
	h4->GetYaxis()->SetTitle("P(TS >= h)");
	h4->SetLineWidth(2);
	h4->SetLineColor(4);
	
	hPVALUE->SetLineWidth(2);
	hPVALUE->SetLineColor(4);
	hCDF->SetLineWidth(2);
	hCDF->SetLineColor(kGreen+4);

	//AB: rimuovi questo se non vuoi il disegno dell'istogramma
	h4->GetYaxis()->SetTitleSize(0.03);
	h4->GetYaxis()->SetTitleOffset(1.80);
	h4->Draw("L");
	
	
	hPVALUE->Draw("SAMEE");
	//hCDF->Draw("SAMEE");
	
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
	

// 	f1->Draw("same");
//     	f2->Draw("SAME");
// 	f1_der->Draw("SAME");
// 	f2_der->Draw("SAME");

	///////////////
	gmb = chi2_1_05->GetMaximumBin();
// 	chi2_1_05->Scale(1.0/chi2_1_05->GetBinContent(gmb) * evalfirst2);

//  	chi2_1->Draw("SAMEL");
 	//chi2_1_05->Draw("SAMEL");
	//gmb = derivative_chi2_1_05->GetMaximumBin();
	//derivative_chi2_1_05->Scale(1.0/derivative_chi2_1_05->GetBinContent(gmb) * evalfirst2);
// 	derivative_chi2_1->SetLineColor(3);
// 	derivative_chi2_1_05->SetLineColor(3);
// 	derivative_chi2_1->Draw("SAME");
// 	derivative_chi2_1_05->Draw("SAMEL");
	

	//c2->BuildLegend();

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
	//chi2_1_05->Draw("SAMEL");
//   	integral_chi2_1->Draw("SAMEL");
// 	integral_chi2_1_05->Draw("SAMEL");
// 	derivative_chi2_1->Draw("SAME");
// 	derivative_chi2_1_05->Draw("SAMEL");
	
	//c1_h1->BuildLegend();

	TCanvas* c1_galiso = new TCanvas("GALISO", "GALISO", dimWindowX, dimWindowY);
	c1_galiso->Divide(4,1);
	c1_galiso->cd(1);
	c1_galiso->SetLogy();
	hh_FLUX->Draw("");
	c1_galiso->cd(2);
	hh_GAL->Draw("");
	c1_galiso->cd(3);
	hh_ISO->Draw("");
	c1_galiso->cd(4);
	hhISOGAL->Draw("COL");
	
	//source coordinates **********************************
	
	TCanvas* canlb = new TCanvas("LB", "LB", 1024, 768);
	canlb->Divide(2,1);
	canlb->cd(1);
	gPad->SetLogz();
	hlb->Draw("HISTCOLZ");
	canlb->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	gPad->SetLogz();
	htsflux->Draw("HISTCOLZ");
	
	TCanvas* canlb2 = new TCanvas("LB2", "LB2", 1024, 768);
	canlb2->Divide(2,1);
	canlb2->cd(1);
	gPad->SetLogy();
	hl->Draw("HISTCOLZ");
	canlb2->cd(2);
	gPad->SetLogy();
	hb->Draw("HISTCOLZ");
	
	
	TCanvas* cantsr = new TCanvas("TSR", "TSR", 1024, 768);
	cantsr->Divide(3,1);
	cantsr->cd(1);
	htsr->Draw("HISTCOLZ");
	cantsr->cd(2);
	htsd->Draw("HISTCOLZ");
	cantsr->cd(3);
	hrd->Draw("HISTCOLZ");
	/*cantsr->cd(4);
	rsel->Draw("HIST");
	*/
	if(outfileprefix != "") {
		if(filename != "") {
			c1_map->cd();
			c1_map->SaveAs(outfileprefix + "_TSmap.jpg");
			if(saveRootFile) c1_map->SaveAs(outfileprefix + "_TSmap.root");
		}
		c1_h1->cd();
		c1_h1->SaveAs(outfileprefix + "_tsdist.jpg");
		if(saveRootFile) c1_h1->SaveAs(outfileprefix + "_tsdist.root");
		c2->cd();
		c2->SaveAs(outfileprefix + "_cumtsdist.jpg");
		if(saveRootFile) c2->SaveAs(outfileprefix + "_cumtsdist.root");

		c1_error->cd();
// 		c1_error->SaveAs(outfileprefix + "_tsdisterror.jpg");
// 		if(saveRootFile) c1_error->SaveAs(outfileprefix + "_tsdisterror.root");
	}
	
	if(fout) fout->Close();
	
	return T;

}