#include <stdio.h>
#include <math.h>
//Versione finale dello script. Questo fa tutto.
//Riferimento per il paper sul TS
/*

/*
 hp
 1) chi2_1 * delta
 3) chi2_3 * delta
 4) ICL C1 Cygnus region FF=3, t_ICL=9, minTS=4
 5) ICL C1 Empty Galactic region FF=3, t_ICL=6, minTS=4
 6) FF=1 Galactic region
 7) ALL Cygnus Galactic region FF=3 minTS=4
 8) ALL Empty Galactic region FF=3 minTS=4
 9) IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4 NOSOURCE
 10) IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4 SOURCE
 11) Equation 3 - A&A 540 A79, 2012
 12) Equation 4 - A&A 540 A79, 2012
 */

//N-> number of maps
//TS or H -> is the is the minimum level of TS chosen to define a detection;
//K -> is the number of trials in the given position with TS ≥ h (a detection). This is the result of our experiment.
//r = 1 -> hp = 1
//r = 3 -> hp = 3
//M-> numder of bins for each map (or number of trials for each map)

void repeatedPostTrial(int N = 180, int K = 10, int hp = 1, double TS = 25, double delta=1, int M = 1, int N1 = 1, double eta1=0.45, int N2 = 2, double eta2=4.8e-03 ) {
	
	
	Double_t pre = 0;
		
	cout << "\n---------- " << endl;
	if(hp == 1) {
		pre = ROOT::Math::chisquared_cdf_c(TS, 1)  * delta;
		cout << "\chi^2_1 * " << delta << endl;
	}
	if(hp == 3) {
		pre = ROOT::Math::chisquared_cdf_c(TS, 3) * delta;
		cout << "\chi^2_3 * " << delta << endl;
	}
		
		//ICL C1 (t_ICL = 9)
		if(hp == 4) {
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
		if(hp == 5) {
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
		if(hp == 6) {
			if(TS<1)
				pre = 6.50849e-01;
				if(TS>=1 )
					pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 4.67050e-01;
		}
	
	//ALL Cygnus Galacti region
	if(hp == 7) {
		double loccl = 6;
		if(TS<1)
			pre = 8.56219e-01;
		if(TS>=1 && TS < loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 4.67802e-01;
		if(TS>=loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 6.17182e-03 ;
	}
	
	//ALL Empty Galacti region
	if(hp == 8) {
		double loccl = 6;
		if(TS<1)
			pre = 8.93596e-01;
		if(TS>=1 && TS < loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS, 1) * 3.46840e-01;
		if(TS>=loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS - loccl, 5) * 3.96381e-03 ;
	}
		
		
	//IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4 NOSOURCE
		if(hp == 9) {
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
	if(hp == 10) {
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
	if(hp == 11) {
		if(TS<1)
			pre = delta1;
			if(TS>=1 )
				pre = ROOT::Math::chisquared_cdf_c(TS, N1) * eta1;
	}
		
	//Equation 4 - A&A 540 A79, 2012
	if(hp == 12) {
		double loccl = 6;
		if(TS<1)
			pre = delta;
		if(TS>=1 && TS < loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS, N1) * eta1;
		if(TS>=loccl)
			pre = ROOT::Math::chisquared_cdf_c(TS - loccl, N2) * eta2 ;
	}
	
	Double_t post = 1.0 - TMath::Power(1. - pre, N);
	if(post == 0) {
		cout << "Warning: Post trial calculated with the approxiimation" << endl;
		post = pre * N;
	}
	cout << "h = TS = " << TS << endl;
	printf("p  = %E\n", pre);
	
	cout << endl << "N  = " << N << endl;
	//cout << "p = " << pre << endl;
	
	if(hp == 1) cout << "r = 1 (\delta = " << delta << ")" << endl;
	if(hp == 2) cout << "r = 3 (\delta = " << delta << ")" << endl;
	cout << "M  = " << M << endl;
	cout << "This is the M of the A&A paper" << endl;
	//cout << "p-value pre-trial  = " << pre << endl;
	printf("p-value pre-trial  = %E\n", pre);
	//cout << "p-value post-trial = " << post << endl;
	printf("p-value post-trial = 1 - (1 - p)^N = %E\n", post);
	
	//double sigma_pre = TMath::ErfcInverse(pre)*TMath::Sqrt(2);
	//double sigma_post = TMath::ErfcInverse(post)*TMath::Sqrt(2);
	//long double sigma_pre = ROOT::Math::gaussian_quantile_c(pre,1.0);
	//long double sigma_post = ROOT::Math::gaussian_quantile_c(post,1.0);
	
	//cout << "sigma pre-trial  = " << sigma_pre  << endl;
	cout << "sigma pre-trial  (ErfInverse) = " << TMath::ErfInverse(1.0-pre)*TMath::Sqrt(2)  << endl;
	//cout << "sigma pre-trial (normal_quantile_c) = " << ROOT::Math::normal_quantile_c(pre,1.0)  << endl;
	cout << "sigma post-trial (ErfInverse) = " << TMath::ErfInverse(1.0-post)*TMath::Sqrt(2)  << endl;
	//cout << "sigma post-trial (normal_quantile_c) = " << sigma_post  << endl;
	
	/*
	cout << "---------------------------------------------------" << endl;
	
	long double resultexit;
	cout << "---------- Single map" << endl;
	cout << "--- Pre trial" << endl;
	printf("Probability of TS>%d: %3.15f ", TS, pre);
	cout <<  " for single map " << endl;
	cout << "p-value: " << pre << endl;
	//cout << "It correspond to " << ROOT::Math::gaussian_quantile_c(pre,1.0) << " (normal) sigma" << endl;
	cout << "It correspond to " << TMath::ErfInverse(1.0-pre)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	
	long double postp = 1 - TMath::Power(1-pre, N);
	if(postp == 0)
		postp = pre * N;
	
	printf("--- Post trial for %d maps  %6.12e\n", N, postp);
	
	resultexit=TMath::ErfInverse(postp)*TMath::Sqrt(2);
	cout << "p-value: " << 1-postp << endl;
	//cout << "It correspond to " << ROOT::Math::gaussian_quantile_c(postp,1.0) << " (normal) sigma" << endl;
	cout << "It correspond to " << TMath::ErfInverse(1-postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfcInverse(postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	*/
	//cout << "---------------" << endl;
	
	Double_t sumA = 0.0;
	
	/*for(int j=0; j<=K-1; j++) {
		long double pr2 = ( TMath::Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) );
		sumA  += pr2;
		printf("%i      %1.60Lf\n", j, sumA);
		//printf("%E \n", pr2);
		
	}
	
	long double resultsumA = 1. - TMath::Power(sumA, M);
	printf("resultsumA     %4.30Lf\n", resultsumA);*/
	
	sumA = 0.0;
	for(Int_t j=K-1; j>=0; j--) {
		Double_t pr2 = ( TMath::Binomial(N, j) * TMath::Power(pre, (Int_t) j) * TMath::Power(1.0-pre, (Int_t) N-j) );
		sumA  += pr2;
		//printf("%i      %1.60Lf	%1.60Lf\n", j, sumA, pr2);
	}
	//cout << sumA << endl;
	Double_t resultsumA = 1. - TMath::Power(sumA, (Int_t) M);
	//printf("resultsumA     %4.30Lf\n", resultsumA);
	
	Double_t sum=0.;
	for(int j=K-1; j>=0; j--){
		Double_t bin = TMath::Binomial(N, j);
		Double_t p1 = TMath::Power(pre, j);
		//printf("p1      %1.30Lf\n", p1);
		//long double p2 = (1. - pre)**(N-j);
		Double_t p2 = TMath::Power(1.-pre, N-j);
		//printf("p2      %1.30Lf\n", p2);
		sum = sum + (bin * p1 * p2);
		//printf("A%2i %1.50Lf\n", j, sum);
		//printf("B%2i %1.50Lf\n", j, pippo);
	}
	//printf("--> %1.50Lf\n", sum);
	
	
	double resultsum = 1. - TMath::Power(sum, 1);
	//long double resultsum = 1. - (sum**M);
	
	//cout << "resultsum=" << resultsum << endl;
	//
	
	
	cout << "\n---------- " << N << " maps" << endl;
	cout << "Probability of wrong " << K << " (or more) detection in "<< N <<" maps with M=1 trial for each map:\n";
	
	
	cout << "N = " << N << endl;
	cout << "k = " << K << endl;
	cout << "G = 1" << endl;
	if(hp == 1) cout << "r = 1 (\delta = " << delta << ")" << endl;
	if(hp == 2) cout << "r = 3 (\delta = " << delta << ")" << endl;
	//cout << "p = " << pre << endl;
	printf("p = %E\n\n", pre);
	//cout << "Pretrial significance" << endl;
	//printf("1) p-value sum      %6.16e\n", resultsum);
	printf("P(N,X≥k)=1-\sum_{j=0}^{k-1} [ (N j) p^j (1-p)^{N-j} ] = %E\n", resultsum);
	cout << "It correspond to " <<	TMath::ErfInverse(1.0-resultsum)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
				
	//cout << "It correspond to " <<  ROOT::Math::gaussian_quantile_c(resultsum,1.0) << " (normal) sigma \n" << endl;
	
	cout << "\nProbability of wrong " << K << " (or more) detection in "<< N <<" maps with M=" << M << " trials for each map: \n";

	cout << "N = " << N << endl;
	cout << "k = " << K << endl;
	cout << "G = 1" << endl;
	if(hp == 1) cout << "r = 1 (\delta = " << delta << ")" << endl;
	if(hp == 2) cout << "r = 3 (\delta = " << delta << ")" << endl;
	printf("p = %E\n", pre);
	cout << "M = " << M << endl;
	
	double result = 1 - TMath::Power(1 - TMath::BinomialI(pre , N, K), M);
	cout << "\nMethod 1:"<<endl;
	printf("P(N,X≥k)=1-( \sum_{j=0}^{k-1} [ (N j) p^j (1-p)^{N-j} ] )^M = %.6e\n", result);
	cout << "It correspond to " <<	TMath::ErfInverse(1.0-result)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	
	double result2 = 1 - TMath::Power(sum, M);
	cout << "Method 2:"<<endl;
	printf("P(N,X≥k)=1-( \sum_{j=0}^{k-1} [ (N j) p^j (1-p)^{N-j} ] )^M = %.6e\n", result2);
	cout << "It correspond to " <<	TMath::ErfInverse(1.0-result2)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	//cout << "It correspond to " << (long double) ROOT::Math::gaussian_quantile_c(result,1.0) << " (normal) sigma \n" << endl;
	
	/*
	cout << "\n---------- " << endl;
	cout << "Probability of wrong " << K << " (exactly) detection in "<< N <<" maps of 1 bins : \n"   ;
	long double resultK =  TMath::Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) ;
	printf("1) p-value      %6.16e\n", resultK);
	cout << resultK << endl;
	//NO resultexit=TMath::ErfInverse(1-resultK)*TMath::Sqrt(2);
	cout << "It correspond to " << TMath::ErfInverse(1-resultK)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	//cout << "It correspond to " <<  ROOT::Math::gaussian_quantile_c(resultK,1.0) << " (normal) sigma \n" << endl;
	*/

	
}
