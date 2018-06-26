
//n-> number of maps
//p -> is the p-value
//k -> is the number of trials in the given position with TS ≥ h (a detection). This is the result of our experiment.


Double_t SMProbWithPValue(Double_t n, Double_t p, Double_t k , Double_t drawgraph = false, Double_t eta = 0.1, int dof=1	) {
	Int_t max_cumTSdisbinx = 1000;
	double resultexit;
	cout << "#########################################################" << endl;
	cout << "---------- Single map" << endl;
	cout << "--- Pre trial: ";
	cout << p << " for single map " << endl;
	cout << "It correspond to " << TMath::ErfInverse(1-p)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << p << endl;
	
	double postp = 1 - TMath::Power(1-p, n);
	if(postp == 0)
		postp = p * n;
	
	printf("--- Post trial %6.12e for %i maps\n", postp, n);

	resultexit=TMath::ErfInverse(postp)*TMath::Sqrt(2);
	cout << "It correspond to " << ROOT::Math::gaussian_quantile_c(postp,1.0) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfInverse(1-postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfcInverse(postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	cout << "p value: " << 1-postp << endl;
	
	Double_t sum=0;
	for(int j=0; j<=k-1; j++)
		sum  += ( TMath::Binomial(n, j) * TMath::Power(p, j) * TMath::Power(1-p, n-j) );

	
	Double_t resultsum = 1 - TMath::Power(sum, 1);
	
	cout << "\n---------- " << n << " maps" << endl;
	
	cout << "Probability of wrong " << k << " (or more) detection in "<< n <<" maps: \n"   ;
	cout << "Pretrial significance" << endl;
	printf("1) sum      %6.16e\n", resultsum);
	//cout << "It correspond to " << TMath::ErfInverse(1-resultsum)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	cout << "It correspond to " <<  ROOT::Math::gaussian_quantile_c(resultsum,1.0) << " (normal) sigma \n" << endl;
	Double_t result = 1 - TMath::Power(1 - TMath::BinomialI(p, n, k), 1);
	
 	printf("2) binomial %6.16e\n\n", result);
	
	cout << "\n---------- " << endl;
	cout << "Probability of wrong " << k << " (exactly) detection in "<< n <<" maps : \n"   ;
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
