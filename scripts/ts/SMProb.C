
//n-> number of maps
//m-> numder of bins for each map (or number of trials for each map)
//p -> is the p-value
//h -> is the is the minimum level of TS chosen to define a detection;
//k -> is the number of trials in the given position with TS ≥ h (a detection). This is the result of our experiment.


Double_t SMProb(Double_t n, Double_t m, Double_t p, Double_t h, Double_t k , Double_t drawgraph = false, Double_t eta = 0.1, int dof=1	) {
	Int_t max_cumTSdisbinx = 1000;
		//TF1* f2 = new TF1("f_chi2/2", "[0] * TMath::Exp(-x/2.0) / ( TMath::Sqrt(TMath::TwoPi() * x)) ", 0, max_cumTSdisbinx);
	//f2->SetParameter(0, eta);
/*	TF1* f3 = new TF1("f_chi2_N/2", " [0] * TMath::Exp(-x/2.0) * TMath::Power(x, [1]/2.0 - 1) / ( TMath::Power(2, [1]/2.0 ) * TMath::Gamma([1]/2.0) ) ", 0, max_cumTSdisbinx);
	f3->SetParameter(0, eta);
	f3->SetParameter(1, dof);
	
	//Double_t p = f2->Integral(h, max_cumTSdisbinx); //p(h) = P(TS >= h)
	Double_t p = f3->Integral(h, max_cumTSdisbinx); //p(h) = P
*/	
	double resultexit;
	cout << "---------- Single map" << endl;
	cout << "--- Pre trial" << endl;
	printf("Probability of TS>%d: %6.12f (", h, p);
	cout << p << ") for single map " << endl;
	cout << "It correspond to " << TMath::ErfInverse(1-p)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << p << endl;
	
	double postp = 1 - TMath::Power(1-p, n);
	if(postp == 0)
		postp = p * n;
	
	printf("--- Post trial for %d maps  %6.12e\n", n, postp);

	resultexit=TMath::ErfInverse(postp)*TMath::Sqrt(2);
	cout << "It correspond to " << ROOT::Math::gaussian_quantile_c(postp,1.0) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfInverse(1-postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfcInverse(postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << 1-postp << endl;
	
	Double_t sum=0;
	for(int j=0; j<=k-1; j++)
		sum  += ( TMath::Binomial(n, j) * TMath::Power(p, j) * TMath::Power(1-p, n-j) );

// 	Double_t sum2=0;
// 	for(int j=k; j<=n; j++)
// 		sum2 += ( TMath::Binomial(n, j) * TMath::Power(p, j) * TMath::Power(1-p, n-j) );

// 	Double_t sum3=0;
// 	for(int j=0; j<=n; j++)
// 		sum3 += ( TMath::Binomial(n, j) * TMath::Power(p, j) * TMath::Power(1-p, n-j) );


// 	printf("%6.9f\n", sum);
// 
// 	printf("%6.9f\n", 1 - TMath::BinomialI(p, n, k));
// 
// 	printf("%6.9f\n", 1 - sum2);

// 	printf("%6.9f\n", sum3);

	
	Double_t resultsum = 1 - TMath::Power(sum, m);
	
	cout << "\n---------- " << n << " maps" << endl;
	
	cout << "Probability of wrong " << k << " (or more) detection in "<< n <<" maps of " << m << " bins : \n"   ;
	cout << "Pretrial significance" << endl;
	printf("1) sum      %6.16e\n", resultsum);
	//cout << "It correspond to " << TMath::ErfInverse(1-resultsum)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	cout << "It correspond to " <<  ROOT::Math::gaussian_quantile_c(resultsum,1.0) << " (normal) sigma \n" << endl;
	Double_t result = 1 - TMath::Power(1 - TMath::BinomialI(p, n, k), m);
	
 	printf("2) binomial %6.16e\n\n", result);
	
	cout << "\n---------- " << endl;
	cout << "Probability of wrong " << k << " (exactly) detection in "<< n <<" maps of 1 bins : \n"   ;
	Double_t resultK =  TMath::Binomial(n, j) * TMath::Power(p, j) * TMath::Power(1-p, n-j) ;
	printf("1) p      %6.16e\n", resultK);
	//resultexit=TMath::ErfInverse(1-resultK)*TMath::Sqrt(2);
	//cout << "It correspond to " << TMath::ErfInverse(1-resultK)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	cout << "It correspond to " <<  ROOT::Math::gaussian_quantile_c(resultK,1.0) << " (normal) sigma \n" << endl;
	
 	/*cout << "Posttrial significance" << endl;
 	
 	Double_t postp2 = 1 - (1 - TMath::Power(1-resultsum, n));
 	printf("3) post sum  %6.16e\n", postp2);
 	cout << "It correspond to " << TMath::ErfInverse(postp2)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	cout << endl;
	*/
	
	std::ofstream asciiFile("gooddetection");
	asciiFile << 1-result << endl;
// 	asciiFile << 1-p << endl;
	asciiFile.close();

// 	printf("%6.9e\n", 1 - TMath::Power(1 - TMath::BinomialI(p, n, k), m));
	if(drawgraph == false)
		return resultexit;

	cout << "draw graph" << endl;
	Int_t np = 720;
	Double_t startx=1;
	Double_t* x = new Double_t[np-startx];
	Double_t* y = new Double_t[np-startx];
	
	for(int i=startx; i<np; i++) {
		x[i-startx] = i;

		//n
		Double_t res = SMProb(i, m, p, h, k, false);
		//m
// 		Double_t res = SMProb(n, startx, k, h, false);
// 		startx += 1000;
		//k
// 		Double_t res = SMProb(n, m, startx, h, false);
// 		startx += 1;

		y[i-startx] = res;
		
	}
	TGraph* gg = new TGraph(np-startx, x, y);
	//gg->SetLineColor(8);
	gg->Draw("AL");
	
	for(int i=startx; i<np; i++) {
		cout << x[i-startx] << " " << y[i-startx] << endl;
	}
	return result;
}