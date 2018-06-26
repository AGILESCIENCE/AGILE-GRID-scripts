

//h -> is the is the minimum level of TS chosen to define a detection;
//n-> number of maps

#define MULTIBIN 10.0

Double_t sumchi2Ndeltafunc(Double_t *x, Double_t *par)
{
	//Double_t f = par[0]*(x[0]<=1.0/MULTIBIN?1:0) + par[1]*exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) ;
	Double_t f = par[0]*((x[0]<=1.0/MULTIBIN)?1:0) + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) ;

	return f;
}

Double_t SMProbP(Double_t h, Double_t n = 1, Double_t delta, Double_t eta = 0.5, Int_t  dof=1, TString histofilename = "") {
	return 0; //non usare questa funzione
	cout << histofilename << endl;
	Double_t p;
	if(histofilename == "") {
		cout << "calculation with chi_2" << endl;
		Int_t max_cumTSdisbinx = 1000;
		//TF1* f2 = new TF1("f_chi2/2", "[0] * TMath::Exp(-x/2.0) / ( TMath::Sqrt(TMath::TwoPi() * x)) ", 0, max_cumTSdisbinx);
		//f2->SetParameter(0, eta);

		
		TF1* f4sf1 = new TF1("sumchi2(N)delta", sumchi2Ndeltafunc, 0, max_cumTSdisbinx,3);
		f4sf1->SetParName(0, "delta");
		f4sf1->SetParameter( 1, delta);
		f4sf1->SetParName(1, "eta");
		f4sf1->SetParameter( 1, eta);
		f4sf1->SetParName(2, "N");
		f4sf1->SetParameter( 2, dof);
		
		TF1* f3 = new TF1("fchi2_Nfit", " [0] * TMath::Exp(-x/2.0) * TMath::Power(x, [1]/2.0 - 1) / ( TMath::Power(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", 0, max_cumTSdisbinx);

		
		
		f3->SetParameter(0, eta);
		
		f3->SetParameter(1, dof);
		//Double_t p = f2->Integral(h, max_cumTSdisbinx); //p(h) = P(TS >= h)
		//p = f3->Integral(h, max_cumTSdisbinx); //p(h) = P
	
		if(delta == 0)
			p = f3->Integral(h, max_cumTSdisbinx); //p(h) = P
		else {
			p = f4sf1->Integral(h, max_cumTSdisbinx); //p(h) = P
			
		}

	} else {
		cout << "calculation with histogram" << endl;
		Int_t max_cumTSdisbinx = 36;
		TFile f(histofilename);
		
		p = TS2->Integral(h, max_cumTSdisbinx);
	}
	
	//p = p1;
	cout << "--- Pre trial" << endl;
	printf("Probability of TS>%d: %6.12f (", h, p);
	cout << p << ") for single map " << endl;
	cout << "It correspond to " << TMath::ErfInverse(1-p)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << p << endl;
	Double_t postp = 1 - (1 - TMath::Power(1-p, n));
	
	printf("--- Post trial for %d maps  %6.12e\n", n, postp);

	cout << "It correspond to " << TMath::ErfInverse(postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << 1-postp << endl;
	std::ofstream asciiFile("gooddetectionp");

	asciiFile << 1-p << endl;
	asciiFile.close();

	//return (1-postp) * 100;
}