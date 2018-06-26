#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>       // std::string
#include "TGraphAsymmErrors.h"

using namespace std;

TVectorF x(100);
TVectorF y(100);
TVectorF exl(100);
TVectorF exh(100);
TVectorF eyl(100);
TVectorF eyh(100);
Double_t* edges;
float spectral_index = 0;

//funzione per fare il fitting
/*
<source name="PowerLaw_source" type="PointSource">
<!-- point source units are cm^-2 s^-1 MeV^-1 -->
<spectrum type="PowerLaw">
<parameter free="1" max="1000.0" min="0.001" name="Prefactor" scale="1e-09" value="1"/>
<parameter free="1" max="-1.0" min="-5." name="Index" scale="1.0" value="-2.1"/>
<parameter free="0" max="2000.0" min="30.0" name="Scale" scale="1.0" value="100.0"/>
</spectrum>
<spatialModel type="SkyDirFunction">
<parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="83.45"/>
<parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="21.72"/>
</spatialModel>
</source>
*/
Double_t powerlaw(Double_t *x, Double_t *par)
{
		//Prefactor = $N_0$ par[0]
		//Index = $\gamma$ par[1]
		//Scale = $E_0$ par[2]
		//dN / dE = N_0 * ( E / E_0 ) ^ \gamma
		Double_t gamma = -par[1]; //ATTENZIONE!!!!!!!!!
        Double_t f = par[0] * pow((x[0]/par[2]), gamma);
        return f;
}	

/*
<source name="LogParabola_source" type="PointSource">
<!-- point source units are cm^-2 s^-1 MeV^-1 -->
<spectrum type="LogParabola">
<parameter free="1" max="1000.0" min="0.001" name="norm" scale="1e-9" value="1"/>
<parameter free="1" max="10" min="0" name="alpha" scale="1.0" value="1"/>
<parameter free="1" max="1e4" min="20" name="Eb" scale="1" value="300."/>
<parameter free="1" max="10" min="0" name="beta" scale="1.0" value="2"/>
</spectrum>
<spatialModel type="SkyDirFunction">
<parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="83.45"/>
<parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="21.72"/>
</spatialModel>
</source>
*/
Double_t logparabola(Double_t *x, Double_t *par)
{
		//norm = $N_0$ par[0]
		//alpha = par[1]
		//beta par[2]
		//Eb = $E_0$ par[3]
		Double_t alpha = par[1]; //ATTENZIONE
        Double_t f = par[0] * pow((x[0]/par[3]), -(alpha + par[2] * TMath::Log(x[0] / par[3])));
        return f;
}

/*
<source name="PLExpCutoff_source" type="PointSource">
<!-- point source units are cm^-2 s^-1 MeV^-1 -->
<spectrum type="PLSuperExpCutoff">
<parameter free="1" max="1000" min="1e-05" name="Prefactor" scale="1e-07" value="1"/>
<parameter free="1" max="0" min="-5" name="Index1" scale="1" value="-1.7"/>
<parameter free="0" max="1000" min="50" name="Scale" scale="1" value="200"/>
<parameter free="1" max="30000" min="500" name="Cutoff" scale="1" value="3000"/>
</spectrum>
<spatialModel type="SkyDirFunction">
<parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="83.45"/>
<parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="21.72"/>
</spatialModel>
</source>
*/
Double_t plexpcutoff(Double_t *x, Double_t *par)
{
		//norm = $N_0$ par[0]
		//gamma = par[1]
		//E0 par[2]
		//Ec = par[3] //cutoff energy
		Double_t gamma = -par[1]; //ATTENZIONE!!!!!!!!!
        Double_t f = par[0] * pow((x[0]/par[2]), gamma) * exp(- (x[0] - par[2]) / par[3] );
        return f;
}

/*
<source name="PLSuperExpCutoff_source" type="PointSource">
<!-- point source units are cm^-2 s^-1 MeV^-1 -->
<spectrum type="PLSuperExpCutoff">
<parameter free="1" max="1000" min="1e-05" name="Prefactor" scale="1e-07" value="1"/>
<parameter free="1" max="0" min="-5" name="Index1" scale="1" value="-1.7"/>
<parameter free="0" max="1000" min="50" name="Scale" scale="1" value="200"/>
<parameter free="1" max="30000" min="500" name="Cutoff" scale="1" value="3000"/>
<parameter free="1" max="5" min="0" name="Index2" scale="1" value="1.5"/>
</spectrum>
<spatialModel type="SkyDirFunction">
<parameter free="0" max="360." min="-360." name="RA" scale="1.0" value="83.45"/>
<parameter free="0" max="90." min="-90." name="DEC" scale="1.0" value="21.72"/>
</spatialModel>
</source>
*/
Double_t plsuperexpcutoff(Double_t *x, Double_t *par)
{
		//norm = $N_0$ par[0]
		//gamma1 = par[1]
		//E0 par[2]
		//Ec = par[3] //cutoff energy
		//gamma2 = par[4]
		Double_t gamma1 = -par[1]; //ATTENZIONE!!!!!!!!!
		Double_t gamma2 = par[4];
        Double_t f = par[0] * pow((x[0]/par[2]), gamma1) * exp(- pow(x[0] / par[3], gamma2) );
        return f;
}

//(TMath::Log10(400)-TMath::Log10(200))/2 + TMath::Log10(200)

int loadAGILESpectra(TString filename, int um) {
	//(0)sqrtts (1)flux[ph/cm2/s] (2)flux_err (3)erg[erg/cm2/s] (4)Erg_err (5)Emin (6)Emax (7)E_log_center (8)exp (9)flux_ul (10)spectral_index (11)spectral_index_err (12)id_detection
	
	double correctionfactor = 1.0;
	double corrfact[100];
	for(int i=0; i<100; i++ ) corrfact[i] = 1;
	/*
	corrfact[0] = 1.79421;
	corrfact[1] = 2.06059;
	corrfact[2] = 2.73078*0.75;
	corrfact[3] = 1.7181*0.75;
	corrfact[4] = 0.08104;
	*/

	int i = 0;
	float sqrtts, flux, flux_err, erg, erg_err, emin, emax, e_log_center, exp, flux_ul, erglogul;
	
	TTree* t = new TTree("DATA", "");
	//METTERE IL DATA TYPE IN MODO ESPLICITO:
	Long64_t nlines = t->ReadFile(filename, "SQRTTS/F:FLUX/F:FLUXE/F:ERG/F:ERGE/F:EMIN/F:EMAX/F:ELOGC/F:EXP/F:FLUXUL/F:SI/F:SIE/F:ID/F:ERGLOGUL/F");
	edges = new Double_t[nlines + 1];
	cout << nlines << endl;


	t->SetBranchAddress("ERG", &erg);
	t->SetBranchAddress("ERGLOGUL", &erglogul);
	t->SetBranchAddress("ERGE", &erg_err);
	t->SetBranchAddress("EMIN", &emin);
	t->SetBranchAddress("EMAX", &emax);
	t->SetBranchAddress("ELOGC", &e_log_center);
	t->SetBranchAddress("SI", &spectral_index);
	t->SetBranchAddress("FLUX", &flux);
	t->SetBranchAddress("FLUXE", &flux_err);
	t->SetBranchAddress("FLUXUL", &flux_ul);
	t->SetBranchAddress("SQRTTS", &sqrtts);
	Long64_t j = 0;
	for(j = 0; j<nlines; j++) {
    	t->GetEntry(j);
		
    	//use this if you are working with erg
    	x[j] = e_log_center;
    	exl[j] = e_log_center - emin;
    	exh[j] = emax - e_log_center;
		
    	
    	//use this if are working with ph/cm2/s
    	/*
    	x[j] = emin + (emax-emin)/2.0;
    	exl[j] = (emax-emin)/2.0;
    	exh[j] = (emax-emin)/2.0;
		*/
		double fp, fpe, ful;
		if(um == 0) {
			fp = flux;
			fpe = flux_err;
			ful = flux_ul;
		} else {
			fp = erg;
			fpe = erg_err;
			ful = erglogul;
		}
		if(sqrtts >= 3) {
    		y[j] = fp * corrfact[j];
    		eyl[j] = eyh[j] = fpe * corrfact[j];
  		} else {
  			cout << j << " considering flux ul " << endl;
  			y[j] = ful * corrfact[j];
    		eyl[j]  = corrfact[j] * ful / 2.0;
			eyh[j] = 0;
  		}
  		/*
    	if(j == 3)
    		y[j] *= 1.4;
    	if(j == 4)
    		y[j] *= 1.8;
    	*/
    	//edges[j] = emin + (emax-emin)/2 ;
    	//edges[j] = e_log_center;
  
    	edges[j] = emin;
  	
		cout << j << ": sqrtts: " << sqrtts << " energy " << x[j] << " [" << emin << "-" << emax << "] ph " << y[j] << " +/- " << " (-" << eyl[j] << ", " << eyh[j] << ") cor factor " << corrfact[j] << endl;
	}
	edges[nlines] = emax;
	cout << "dim " << x.GetNrows() << endl;
	x.ResizeTo(j);
	y.ResizeTo(j);
	exl.ResizeTo(j);
	exh.ResizeTo(j);
	eyl.ResizeTo(j);
	eyh.ResizeTo(j);
	return x.GetNrows();
}

int loadFERMISpectra(TString filename, int um) {
	//(0)sqrtts (1)flux[ph/cm2/s] (2)erg[erg/cm2/s] (3)Emin (4)Emax (5)E_log_center
	

	int i = 0;
	float sqrtts, flux, flux_err, erg, erg_err, emin, emax, e_log_center, exp, flux_ul;
	
	TTree* t = new TTree("DATA", "");
	//METTERE IL DATA TYPE IN MODO ESPLICITO:
	Long64_t nlines = t->ReadFile(filename, "SQRTTS/F:FLUX/F:ERG/F:EMIN/F:EMAX/F:ELOGC/F:SI/F");
	edges = new Double_t[nlines + 1];
	cout << nlines << endl;

	t->SetBranchAddress("SQRTTS", &sqrtts);
	t->SetBranchAddress("FLUX", &flux);
	t->SetBranchAddress("ERG", &erg);
	t->SetBranchAddress("EMIN", &emin);
	t->SetBranchAddress("EMAX", &emax);
	t->SetBranchAddress("ELOGC", &e_log_center);
//	t->SetBranchAddress("SI", &spectral_index);
	Long64_t j = 0;
	for(j = 0; j<nlines; j++) {
    	t->GetEntry(j);
    	
    	//use this if you are working with erg
    	x[j] = e_log_center;
    	exl[j] = e_log_center - emin;
    	exh[j] = emax - e_log_center;
    	
    	
    	/*
    	//use this if are working with ph/cm2/s
    	x[j] = emin + (emax-emin)/2.0;
    	exl[j] = (emax-emin)/2.0;
    	exh[j] = (emax-emin)/2.0;
    	*/
		if(um == 0) {
    		y[j] = flux; //erg
			eyl[j] = eyh[j] = y[j] / 100.;
		}
		if(um == 1) {
			y[j] = erg;
			eyl[j] = eyh[j] = 0;
		}
		
    	
    	//edges[j] = emin + (emax-emin)/2 ;
    	//edges[j] = e_log_center;
    	edges[j] = emin;
    	
    	//eyl[j] = eyh[j] = flux / 10.;
    	
    	//cout << x[j] << " " << exl[j] << " " << exh[j] << " ph " << y[j] << " " << " " << eyl[j] << " " << eyh[j] << endl;
		cout << j << ": sqrtts: " << sqrtts << " energy " << x[j] << " [" << emin << "-" << emax << "] ph " << y[j] << " +/- " << " (-" << eyl[j] << ", " << eyh[j] << ") " << endl;
	}
	edges[nlines] = emax;
	cout << "Edges" << endl;
	for(int i=0; i<nlines; i++) {
		cout << i << ": " << edges[i] << " " << edges[i+1] << " center: " << x[i] << endl;
		
	}
	cout << "dim " << x.GetNrows() << endl;
	x.ResizeTo(j);
	y.ResizeTo(j);
	exl.ResizeTo(j);
	exh.ResizeTo(j);
	eyl.ResizeTo(j);
	eyh.ResizeTo(j);
	return x.GetNrows();
}

void drawSpectra() {

}

//experiment AGILE=0, Fermi=1
//overlap 0=do not overlap, 1=overlap
//um 0=ph cm-2 s-1 1=erg
//TCanvas* c1 = new TCanvas; c1->Divide(3,1);
void fitSpectra(string filename, int experiment, TCanvas* c1, int overlap, int color=kBlue, int um = 0) {

	if(c1 == 0) {
		c1 = new TCanvas();
		c1->Divide(3,1);
	}
	c1->cd(1);
	gPad->SetLogx();
   	gPad->SetLogy();
   	c1->cd(2);
   	gPad->SetLogx();
   	gPad->SetLogy();
	c1->cd(3);
	gPad->SetLogx();
	gPad->SetLogy();

	int n=0;
	if(experiment == 0)
		n = loadAGILESpectra(filename, um);
	if(experiment == 1)
		n = loadFERMISpectra(filename, um);
		
	double maxenergy = x[n-1] + exh[n-1];
   	double minenergy = x[0];// - exl[0];
   	cout << "min energy: " << minenergy << endl;
   	cout << "max energy: " << maxenergy << endl;
	TString tit = filename;
	
	//////////////graphs
	TGraphAsymmErrors* gr = new TGraphAsymmErrors(x,y,exl,exh,eyl,eyh);
	//TGraph* gr = new TGraph(x, y); 
   	
	//////////////histos
	/*
	TH1D* gr = new TH1D("spectra", tit, n, edges);
	for(int i=0; i<n+1; i++)
		cout << edges[i] << endl;
	for(int i=1; i<=n; i++) {
		gr->SetBinContent(i, y[i-1]);
		cout << "integral [" << edges[i-1] << ", " << edges[i] << "] " << ": " << gr->Integral(i, i) << endl;
		//correctoin
		
		if(i == -1) {
			Double_t integpre = gr->Integral(i, i);
			gr->SetBinContent(i, y[i-1] * 2.5);
			Double_t integpost = gr->Integral(i, i);
			Double_t diffint = integpost - integpre;
			cout << "integral B [" << edges[i-1] << ", " << edges[i] << "] " << ": " << gr->Integral(i, i) << endl;
			Double_t intprebinold = gr->Integral(i-1, i-1);
			intprebinold -= diffint;
			//Double_t newv = intprebinold  / (edges[i-2] - edges[i-1]);
			//gr->SetBinContent(i-1, newv);
		}
	}
	*/
   	gr->SetTitle(tit);
   	gr->SetMarkerColor(4);
   	gr->SetMarkerStyle(21);
   	c1->cd(1);
   	if(overlap == 1) {
   		gr->Draw("LP*"); //also AL for histos
		gr->GetXaxis()->SetTitle("MeV");
		if(um==0)
			gr->GetYaxis()->SetTitle("#ph/cm^{2}/s");
		else
			gr->GetYaxis()->SetTitle("#erg/cm^{2}/s");
   		gr->SetLineColor(color);
		gr->SetLineWidth(2);
   	} else {
   		gr->Draw("ALP*");
		gr->GetXaxis()->SetTitle("MeV");
		if(um==0)
			gr->GetYaxis()->SetTitle("ph/cm^{2}/s");
		else
			gr->GetYaxis()->SetTitle("#erg/cm^{2}/s");
   		gr->SetLineColor(color);
		gr->SetLineWidth(2);
   	}
	
	////////////////////////////////////////////////////////////////////////////////
   	//Prefactor = $N_0$ par[0]
		//Index = $\gamma$ par[1]
		//Scale = $E_0$ par[2]
   	TF1* plaw = new TF1("powerlaw", powerlaw, minenergy, maxenergy, 3);
    plaw->SetParName(0, "N0");
    plaw->SetParName(1, "gamma");
    plaw->SetParName(2, "E0");
    
    //<parameter free="1" max="1000.0" min="0.001" name="Prefactor" scale="1e-09" value="1"/>
   // plaw->SetParLimits(0, 1e-12, 1e-5);
    plaw->SetParameter(0, 1e-9);
    //plaw->FixParameter(0, 1e-9);
    
    //<parameter free="1" max="-1.0" min="-5." name="Index" scale="1.0" value="-2.1"/>
    plaw->SetParLimits(1, 0.1, 5);
    plaw->SetParameter(1, 2.1);
    //plaw->FixParameter(1, 2.1);
    
    //<parameter free="0" max="2000.0" min="30.0" name="Scale" scale="1.0" value="100.0"/>
    //plaw->SetParLimits(2, 30, 200000);
    plaw->SetParameter(2, 100);
    //plaw->FixParameter(2, 100);
    
    if(overlap == 1) 
    	plaw->SetLineColor(color);
    else
    	plaw->SetLineColor(color);
    
    ///////////////////////////////////////////////////////////////////////////////
    TF1* logp = new TF1("logparabola", logparabola, minenergy, maxenergy, 4);
    logp->SetParName(0, "norm");
    logp->SetParName(1, "alpha");
    logp->SetParName(2, "beta");
    logp->SetParName(3, "eb");
    
    //N0
    //<parameter free="1" max="1000.0" min="0.001" name="norm" scale="1e-9" value="1"/>
    logp->SetParLimits(0, 1e-12, 1e-6);
    logp->SetParameter(0, 1e-9);
    //logp->FixParameter(0, 1e-9);
    
    //alpha
    //<parameter free="1" max="10" min="0" name="alpha" scale="1.0" value="1"/>
    logp->SetParLimits(1, 0, 10);
    logp->SetParameter(1, 1);
    //logp->FixParameter(1, 2.44);
    
    //beta
    //<parameter free="1" max="10" min="0" name="beta" scale="1.0" value="2"/>
    logp->SetParLimits(2, 0.2, 10);
    logp->SetParameter(2, 2);
    //logp->FixParameter(2, 0.25); //0.0386
    
    //EB
    //<parameter free="1" max="1e4" min="20" name="Eb" scale="1" value="300."/>
    logp->SetParLimits(3, 20, 1e4);
    logp->SetParameter(3, 300);
    //logp->FixParameter(3, 1894.78);

    if(overlap == 1) 
    	logp->SetLineColor(color);
    else
    	logp->SetLineColor(color);
    
    //////////////////////////////////////////////////////////////////////////////
	TF1* ple = new TF1("plexpcutoff", plexpcutoff, minenergy, maxenergy, 4);
    ple->SetParName(0, "norm");
    ple->SetParName(1, "gamma");
    ple->SetParName(2, "E0");
    ple->SetParName(3, "Ec");	//cutoff
	
	//<parameter free="1" max="1000" min="1e-05" name="Prefactor" scale="1e-07" value="1"/>
	ple->SetParLimits(0, 1e-12, 1e-4);	
	ple->SetParameter(0, 1e-07);
    //ple->FixParameter(0, 1e-07);	
    
    //<parameter free="1" max="0" min="-5" name="Index1" scale="1" value="-1.7"/>
    ple->SetParLimits(1, 0, 5);
    ple->SetParameter(1, 0.5);
    //ple->FixParameter(1, 1.45);
    
    //<parameter free="0" max="1000" min="50" name="Scale" scale="1" value="200"/>
    ple->SetParLimits(2, 50, 10000);
    ple->SetParameter(2, 200);
    //ple->FixParameter(2, 1000);
    
    //<parameter free="1" max="30000" min="500" name="Cutoff" scale="1" value="3000"/>
    ple->SetParLimits(3, 500, 30000);
    ple->SetParameter(3, 3000);
    //ple->FixParameter(3,2100); //500, 30000
    
    if(overlap == 1) 
    	ple->SetLineColor(color);
    else
    	ple->SetLineColor(color);
    
    ///////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

	TF1* psle = new TF1("plsuperexpcutoff", plexpcutoff, minenergy, maxenergy, 5);
    psle->SetParName(0, "Prefactor");
    psle->SetParName(1, "Index1-spectral_index");
    psle->SetParName(2, "Scale");
    psle->SetParName(3, "Ecutoff");	//cutoff
	psle->SetParName(4, "Index2-exp_index");
	
	//<parameter free="1" max="1000" min="1e-05" name="Prefactor" scale="1e-07" value="1"/>
	psle->SetParLimits(0, 1e-12, 1e-4);	
	psle->SetParameter(0, 1e-07);
    //psle->FixParameter(0, 7.55801e-06);	
    
    //<parameter free="1" max="0" min="-5" name="Index1" scale="1" value="-1.7"/>
    //3FGL Spectral_Index
    psle->SetParLimits(1, 0, 5);
    psle->SetParameter(1, 1.7);
    psle->FixParameter(1, 1.0295);
    
    //<parameter free="0" max="1000" min="50" name="Scale" scale="1" value="200"/>
    psle->SetParLimits(2, 50, 1000);
    psle->SetParameter(2, 200);
    //psle->FixParameter(2, 72);
    
    //<parameter free="1" max="30000" min="500" name="Cutoff" scale="1" value="3000"/>
    //Cutoff 3FGL
    psle->SetParLimits(3, 500, 30000);
    psle->SetParameter(3, 3000);
    //psle->FixParameter(3, 3700); 
    
    //<parameter free="1" max="5" min="0" name="Index2" scale="1" value="1.5"/>
    //Exp_Index 
    psle->SetParLimits(4, 0, 5);
    psle->SetParameter(4, 1.5);
    //psle->FixParameter(4, 5);
    
    if(overlap == 1) 
    	psle->SetLineColor(color);
    else
    	psle->SetLineColor(color);
    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    
    //plaw -> PowerLaw
    //logp -> LogParabole
    //ple -> PLExpCutoff
    //psle -> PLSuperExpCutoff
    TF1* fitfun = plaw;

    /*

* Minuit (library libMinuit). Old version of Minuit, based on the TMinuit class. The list of possible algorithms are:
Migrad (default one)
Simplex
Minimize (it is a combination of Migrad and Simplex)
MigradImproved
Scan
Seek

* Minuit2 (library libMinuit2). New C++ version of Minuit. The list of possible algorithm is :
Migrad (default)
Simplex
Minimize
Scan

*Fumili . This is the same algorithm of TFumili, but implemented in the Minuit2 library.

* GSLMultiMin (library libMathMore). Minimizer based on the Multidimensional Minimization routines of the Gnu Scientific Library (GSL). The list of available algorithms is
BFGS2 (default) : second version of the vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm;
BFGS : old version of the vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm;
ConjugateFR : Fletcher-Reeves conjugate gradient algorithm;
ConjugatePR : Polak-Ribiere conjugate gradient algorithm;
SteepestDescent: steepest descent algorithm;

* GSLMultiFit (library libMathMore). Minimizer based on the Non-Linear Least-Square routines of GSL. This minimizer can be used only for least-square fits.

* GSLSimAn (library libMathMore). Minimizer based on simulated annealing.

* Genetic (library libGenetic). Genetic minimizer based on an algorithm implemented in the TMVA package.
*/

/*
Each minimizer can be configured using the ROOT::Math::MinimizerOptions class. The list of possible option that can be set are:

* Minimizer type (MinimizerOptions::SetMinimizerType(const char *)) .
* Minimizer algorithm (MinimizerOptions::SetMinimizerAlgorithm(const char *)).
* Print Level (MinimizerOptions::SetPrintLevel(int )) to set the verbose printing level (default is 0).
* Tolerance (MinimizerOptions::SetTolerance(double )) tolerance used to control the iterations.
* Maximum number of function calls (MinimizerOptions::SetMaxFunctionCalls(int )).
* Maximum number of iterations (MinimizerOptions::SetMaxIterations(int )). Note that this is not used by Minuit
FCN Upper value for Error Definition (MinimizerOptions::SetMaxIterations(int )). Value in the minimization function used to compute the parameter errors. The default is to get the uncertainties at the 68% CL is a value of 1 for a chi-squared function minimization and 0.5 for a log-likelihood function.
* Strategy (MinimizerOptions::SetStrategy(int )), minimization strategy used. For each minimization strategy Minuit uses different configuration parameters (e.g. different requirements in computing derivatives, computing full Hessian (strategy = 2) or an approximate version. The default is a value of 1. In this case the full Hessian matrix is computed only after the minimization.
* Precision (MinimizerOptions::SetTolerance(double )). Precision value in the evaluation of the minimization function. Default is numerical double precision.
*/
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
	//cout << ROOT::Math::MinimizerOptions::DefaultTolerance() << endl;
	//cout << ROOT::Math::MinimizerOptions::DefaultPrecision() << endl;
	ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.000000000001);
	ROOT::Math::MinimizerOptions::SetDefaultPrecision(0.0000000001);
	//ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(2 );
	ROOT::Math::MinimizerOptions opt;
	opt.Print();

    
    //OK gr->Fit(fitfun, "");
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
	//gr->Fit(fitfun);
   
   // gr->Fit(fitfun);
    //gr->Fit(fitfun);
    //gr->Fit(fitfun);
    
	/*
	 
	if(overlap > 0)
    	c1->cd(overlap);
	*/
	
	/*
	 OK
    fitfun->Draw("");
    gr->Draw("LP");
	*/
   // gr->Fit(ple);
    //f3sf->Draw("L");

}
