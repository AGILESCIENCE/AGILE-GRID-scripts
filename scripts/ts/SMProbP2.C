
//p p-value
//h -> is the is the minimum level of TS chosen to define a detection;
//n-> number of maps


Double_t SMProbP2(Double_t p, Double_t h, Double_t n) {
	return 0; //non usare questa funzione
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
	ppp=f2sf1->Integral(h, end); postp = 1 - (1 - TMath::Power(ppp, n));
	cout << "function  \t eval    \t integral \t signficance \t pvalue " << endl;
	cout << "chi^2_1    \t" << f2sf1->Eval(h) << "\t " <<    ppp    << "\t" << TMath::ErfInverse(1-ppp)*TMath::Sqrt(2) << "\t" << postp << endl;
	ppp=f2sf1_05->Integral(h, end); postp = 1 - (1 - TMath::Power(1-ppp, n));
	cout << "1/2 chi^2_1\t" << f2sf1_05->Eval(h) << "\t " << ppp << "\t" << TMath::ErfInverse(1-ppp)*TMath::Sqrt(2) << "\t" << postp<< endl;
	ppp=f2sf3_05->Integral(h, end); postp = 1 - (1 - TMath::Power(1-ppp, n));
	cout << "1/2 chi^2_3\t" << f2sf3_05->Eval(h) << "\t " << ppp << "\t" << TMath::ErfInverse(1-ppp)*TMath::Sqrt(2) << "\t" << postp<< endl;

}