#include <stdio.h>
#include <math.h>
//Versione finale dello script. Questo fa tutto. 28 Novembre 2011
//Riferimento per il paper sul TS
/*

/*
 hp
 1) chi2_1 * delta
 2) chi2_3 * delta
 3) ICL C1 Cygnus region FF=3, t_ICL=9, minTS=4
 4) ICL C1 Empty Galactic region FF=3, t_ICL=6, minTS=4
 5) FF=1 Galactic region
 6) ALL Cygnus Galactic region FF=3 minTS=4
 7) ALL Empty Galactic region FF=3 minTS=4
 8) IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4 NOSOURCE
 9) IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4 SOURCE
 10) Equation 3 - A&A 540 A79, 2012
 11) Equation 4 - A&A 540 A79, 2012
 */

//N-> number of maps
//M-> numder of bins for each map (or number of trials for each map)
//TS or H -> is the is the minimum level of TS chosen to define a detection;
//K -> is the number of trials in the given position with TS ≥ h (a detection). This is the result of our experiment.


void getPValueAndPostTrial(double TS = 25, int N = 180, Double_t M = 1, int K = 10, int hp = 1, double delta=0.5, int N1 = 1, double eta1=0.45, int N2 = 2, double eta2=4.8e-03 ) {
	
	
		double pre = 0;
		
		if(hp == 1)
			pre = ROOT::Math::chisquared_cdf_c(TS, 1)  * delta;
		if(hp == 2)
			pre = ROOT::Math::chisquared_cdf_c(TS, 3) * delta;
		
		//ICL C1 (t_ICL = 9)
		if(hp == 3) {
			double loccl = 6;
			double ticl = 9;
			if(TS<1)
				pre = 8.47994e-01;
			if(TS>=1 && TS < loccl)
				pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 4.99384e-01;
			if(TS>=loccl && TS < ticl)
				pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 6.71985e-03 ;
			if(TS>=ticl && TS < 14)
				pre = ROOT::Math::chisquared_cdf_c(TS, 5) * 1.82840e-02; 
			if(TS>=14)
				pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 1.32454e+00; 
		}
	
		//ICL C1 Empty Galacti region
		if(hp == 4) {
			double loccl = 6;
			double ticl = 9;
			if(TS<1)
				pre = 8.94878e-01;
			if(TS>=1 && TS < loccl)
				pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 3.46840e-01;
			if(TS>=loccl && TS < ticl)
				pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 4.50982e-03 ;
			if(TS>=ticl && TS < 14)
				pre = ROOT::Math::chisquared_cdf_c(TS, 5) * 1.27508e-02; 
			if(TS>=14)
				pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 9.17136e-01;
		}
	
		//FF=1 Galactic region
		if(hp == 5) {
			if(TS<1)
				pre = 6.50849e-01;
				if(TS>=1 )
					pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 4.67050e-01;
		}
	
	//ALL Cygnus Galacti region
	if(hp == 6) {
		double loccl = 6;
		if(TS<1)
			pre = 8.56219e-01;
		if(TS>=1 && TS < loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 4.67802e-01;
		if(TS>=loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 6.17182e-03 ;
	}
	
	//ALL Empty Galacti region
	if(hp == 7) {
		double loccl = 6;
		if(TS<1)
			pre = 8.93596e-01;
		if(TS>=1 && TS < loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 3.46840e-01;
		if(TS>=loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 3.96381e-03 ;
	}
		
		
	//IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4 NOSOURCE
		if(hp == 8) {
			double loccl = 6;
			double ticl = 9;
			if(TS<1)
				pre = 8.6e-01;
			if(TS>=1 && TS < loccl)
				pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 0.47;
			if(TS>=loccl && TS < ticl)
				pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 6e-03 ;
			if(TS>=ticl && TS < 14)
				pre = ROOT::Math::chisquared_cdf_c(TS, 5) * 1.22e-02; 
			if(TS>=14)
				pre = ROOT::Math::chisquared_cdf_c(TS, 3) * 0.053;
		}
	
	//IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4 SOURCE
	if(hp == 9) {
		double loccl = 6;
		double ticl = 9;
		if(TS<1)
			pre = 7.98e-01;
		if(TS>=1 && TS < loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 6.56748e-01;
		if(TS>=loccl && TS < ticl)
			pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 1.25973e-02 ;
		if(TS>=ticl && TS < 14)
			pre = ROOT::Math::chisquared_cdf_c(TS, 5) * 3.38648e-02; 
		if(TS>=14)
			pre = ROOT::Math::chisquared_cdf_c(TS, 2) * 6.01559e-01;
	}
	
	//Equation 3 - A&A 540 A79, 2012
	if(hp == 10) {
		if(TS<1)
			pre = delta1;
			if(TS>=1 )
				pre = ROOT::Math::chisquared_cdf_c(TS, N1) * eta1;
	}
		
	//Equation 4 - A&A 540 A79, 2012
	if(hp == 11) {
		double loccl = 6;
		if(TS<1)
			pre = delta;
		if(TS>=1 && TS < loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS, N1) * eta1;
		if(TS>=loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS - loccl, N2) * eta2 ;
	}
	
	long double post = 1 - TMath::Power(1.-pre, N);
	if(post == 0) {
		cout << "Warning: Post trial calculated with the approxiimation" << endl;
		post = pre * N;
	}
	cout << "TS = " << TS << endl;
	cout << "p-value pre-trial  = " << pre << endl;
	cout << "p-value post-trial = " << post << endl;
	
	long double sigma_pre = TMath::ErfInverse(1-pre)*TMath::Sqrt(2);
	long double sigma_post = TMath::ErfInverse(1-post)*TMath::Sqrt(2);
	//long double sigma_pre = ROOT::Math::gaussian_quantile_c(pre,1.0);
	//long double sigma_post = ROOT::Math::gaussian_quantile_c(post,1.0);
	
	cout << "sigma pre-trial  = " << sigma_pre  << " " << endl;
	cout << "sigma pre-trial (ErfInverse) = " << TMath::ErfInverse(1.0-pre)*TMath::Sqrt(2)  << endl;
	cout << "sigma pre-trial (normal_quantile_c) = " << ROOT::Math::normal_quantile_c(pre,1.0)  << endl;
	cout << "sigma post-trial = " << sigma_post  << endl;
	
	cout << "---------------------------------------------------" << endl;
	
	long double resultexit;
	cout << "---------- Single map" << endl;
	cout << "--- Pre trial" << endl;
	printf("Probability of TS>%d: %3.20Lf (", TS, pre);
	cout << pre << ") for single map " << endl;
	cout << "p-value: " << pre << endl;
	cout << "It correspond to " << TMath::ErfInverse(1-pre)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	
	
	long double postp = 1 - TMath::Power(1-pre, N);
	if(postp == 0)
		postp = pre * N;
	
	printf("--- Post trial for %d maps  %6.12e\n", N, postp);
	
	resultexit=TMath::ErfInverse(postp)*TMath::Sqrt(2);
	cout << "p-value: " << 1-postp << endl;
	cout << "It correspond to " << TMath::ErfInverse(1-post)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfInverse(1-postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfcInverse(postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	
	
	long double sumA=0;
	
	/*for(int j=0; j<=K-1; j++) {
		long double pr2 = ( TMath::Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) );
		sumA  += pr2;
		printf("%i      %1.60Lf\n", j, sumA);
		//printf("%E \n", pr2);
		
	}
	
	long double resultsumA = 1. - TMath::Power(sumA, M);
	printf("resultsumA     %4.30Lf\n", resultsumA);*/
	/*
	sumA = 0;
	for(int j=K-1; j>=0; j--) {
		long double pr2 = ( TMath::Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) );
		sumA  += pr2;
		//printf("%i      %1.60Lf	%1.60Lf\n", j, sumA, pr2);
		
	}
	
	long double resultsumA = 1. - TMath::Power(sumA, M);
	printf("resultsumA     %4.30Lf\n", resultsumA);
	*/
	long double sum=0.;
	for(int j=K-1; j>=0; j--){
		long double bin = TMath::Binomial(N, j);
		long double p1 = TMath::Power(pre, j);
		//printf("p1      %1.30Lf\n", p1);
		//long double p2 = (1. - pre)**(N-j);
		long double p2 = TMath::Power(1.-pre, N-j);
		//printf("p2      %1.30Lf\n", p2);
		sum = sum + (bin * p1 * p2);
		long double pippo;
		pippo = bin * p1 * p2;
		//printf("A%2i %1.50Lf\n", j, sum);
		//printf("B%2i %1.50Lf\n", j, pippo);
	}
	//printf("--> %1.50Lf\n", sum);
	
	
	long double resultsum = 1. - TMath::Power(sum, 1);
	//long double resultsum = 1. - (sum**M);
	printf("resultsum     %4.20Lf\n", resultsum);
	cout << "resultsum=" << resultsum << endl;
	//
	
	
	cout << "\n---------- " << N << " maps" << endl;
	
	cout << "Probability of wrong " << K << " (or more) detection in "<< N <<" maps of " << M << " bins : \n"   ;
	cout << "Pretrial significance" << endl;
	printf("1) p-value sum      %6.16e\n", resultsum);
	//NO cout << "It correspond to " << TMath::ErfInverse(1-resultsum)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	cout << "It correspond to " <<  TMath::ErfInverse(1-resultsum)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	
	
	long double result = 1 - TMath::Power(1 - TMath::BinomialI(pre , N, K), M);
	printf("2) p-value binomial %6.16e\n", result);
	cout << result << endl;
	cout << "It correspond to " << (long double) TMath::ErfInverse(1-result)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	
	cout << "\n---------- " << endl;
	cout << "Probability of wrong " << K << " (exactly) detection in "<< N <<" maps of 1 bins : \n"   ;
	long double resultK =  TMath::Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) ;
	printf("1) p-value      %6.16e\n", resultK);
	cout << resultK << endl;
	//NO resultexit=TMath::ErfInverse(1-resultK)*TMath::Sqrt(2);
	//NO cout << "It correspond to " << TMath::ErfInverse(1-resultK)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	cout << "It correspond to " <<  TMath::ErfInverse(1-resultK)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	

	
}