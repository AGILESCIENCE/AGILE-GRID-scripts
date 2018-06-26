
#include <stdio.h>
#include <math.h>
#include <TMath.h>
#include <Math/QuantFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>
#include <TObject.h>
#include <iostream>
#include <gmp.h>

using namespace std;
using namespace ROOT::Math;

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
 */

//N-> number of maps
//M-> numder of bins for each map (or number of trials for each map)
//TS or H -> is the is the minimum level of TS chosen to define a detection;
//K -> is the number of trials in the given position with TS ≥ h (a detection). This is the result of our experiment.


long double Binomial(Int_t n,Int_t k)
{
   // Calculate the binomial coefficient n over k.

   if (k==0 || n==k) return 1;
   if (n<=0 || k<0 || n<k) return 0;

   Int_t k1=TMath::Min(k,n-k);
   Int_t k2=n-k1;
   long double fact=k2+1;
   for (long double i=k1;i>1.;--i)
      fact *= (k2+i)/i;
   return fact;
}



void getPValueAndPostTrial(double TS = 25, int N = 180, Double_t M = 1, int K = 10, int hp = 1, double delta=0.5 ) {
	
	
		double pre = 0;
		
		if(hp == 1)
			pre = ROOT::Math::chisquared_cdf_c(TS, 1.0)  * delta;
		if(hp == 2)
			pre = ROOT::Math::chisquared_cdf_c(TS, 3.0) * delta;
		
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
		
		
	//IGRJ17354-3255 region ICL C1 FF=3, t_ICL=9, minTS=4
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
	
	long double post = 1 - TMath::Power(1.-pre, N);
	if(post == 0) {
		cout << "Warning: Post trial calculated with the approxiimation" << endl;
		post = pre * N;
	}
	cout << "TS = " << TS << endl;
	cout << "p-value pre-trial  = " << pre << endl;
	cout << "p-value post-trial = " << post << endl;
	
	//double sigma_pre = TMath::ErfcInverse(pre)*TMath::Sqrt(2);
	//double sigma_post = TMath::ErfcInverse(post)*TMath::Sqrt(2);
	long double sigma_pre = ROOT::Math::gaussian_quantile_c(pre,1.0);
	long double sigma_post = ROOT::Math::gaussian_quantile_c(post,1.0);
	
	cout << "sigma pre-trial  = " << sigma_pre  << endl;
	cout << "sigma post-trial = " << sigma_post  << endl;
	
	cout << "---------------------------------------------------" << endl;
	
	long double resultexit;
	cout << "---------- Single map" << endl;
	cout << "--- Pre trial" << endl;
	printf("Probability of TS>%d: %3.20Lf (", TS, pre);
	cout << pre << ") for single map " << endl;
	cout << "p-value: " << pre << endl;
	cout << "It correspond to " << ROOT::Math::gaussian_quantile_c(pre,1.0) << " (normal) sigma" << endl;
	
	
	long double postp = 1 - TMath::Power(1-pre, N);
	if(postp == 0)
		postp = pre * N;
	
	printf("--- Post trial for %d maps  %6.12e\n", N, postp);
	
	resultexit=TMath::ErfInverse(postp)*TMath::Sqrt(2);
	cout << "p-value: " << 1-postp << endl;
	cout << "It correspond to " << ROOT::Math::gaussian_quantile_c(postp,1.0) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfInverse(1-postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	//cout << "It correspond to " << TMath::ErfcInverse(postp)*TMath::Sqrt(2) << " (normal) sigma" << endl;
	
	
	long double sumA=0;
	
	for(int j=0; j<=K-1; j++) {
		long double pr2 = ( Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) );
		sumA  += pr2;
		//printf("%i      %1.60Lf\n", j, pr2);
		//printf("%E \n", pr2);
		
	}
	printf("sumA     %4.30Lf\n", sumA);
	long double resultsumA = 1. - TMath::Power(sumA, M);
	//printf("resultsumA     %4.30Lf\n", resultsumA);
	
	sumA = 0;
	for(int j=K-1; j>=0; j--) {
		long double pr2 = ( Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) );
		sumA  = sumA +  pr2;
		printf("%i      %1.60Lf	\n", j, sumA);
		
	}
	printf("sumA     %4.30Lf\n", sumA);
	resultsumA = 1. - TMath::Power(sumA, M);
	//printf("resultsumA     %4.30Lf\n", resultsumA);
	
	/*sumA = 0;
	for(int j=K-1; j>=0; j--) {
		long double pr2 = ( Binomial(N, j) * pow((long double)pre, (int)j) * pow((long double)1-pre, (int)N-j) );
		sumA  = sumA +  pr2;
		//printf("%i      %1.60Lf	%1.60Lf\n", j, sumA, pr2);
		
	}
	printf("sumA     %4.30Lf\n", sumA);
	resultsumA = 1. - TMath::Power(sumA, M);*/
	
	mpf_t mpf_sumA;
	mpf_t mpf_tmp;
	mpf_init2 (mpf_sumA, 256); 
	mpf_init2 (mpf_tmp, 256); 
	int   nprint = 60;
	for(int j=K-1; j>=0; j--) {
		long double pr2 = ( Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) );
		mpf_set_d(mpf_tmp, pr2);
		mpf_add(mpf_sumA, mpf_sumA, mpf_tmp);
		
     	
     	gmp_printf ("fixed point mpf %.*Ff with %d digits\n", nprint, mpf_sumA, nprint);
	}
	
	mpf_t mpf_resultsumA;
	mpf_init2 (mpf_resultsumA, 256); 
	mpf_pow_ui(mpf_resultsumA, mpf_sumA, M);
	mpf_ui_sub(mpf_resultsumA, 1, mpf_resultsumA);
	gmp_printf ("resultSUMA gmp %.*Ff with %d digits\n", nprint, mpf_resultsumA, nprint);
	
	double resultsum = mpf_get_d(mpf_resultsumA);
	printf("resultsum      %1.60f\n", resultsum);
	
	
	/*long double sum=0.;
	for(int j=0; j<=K-1; j++){
		long double bin = Binomial(N, j);
		long double p1 = pow((long double)pre, (long double)j);
		//printf("p1      %1.30Lf\n", p1);
		//long double p2 = (1. - pre)**(N-j);
		long double p2 = pow((long double)1.-pre, (long double)N-j);
		//printf("p2      %1.30Lf\n", p2);
		
		long double mult;
		
		mult = bin * p1 * p2;
		
		sum = sum + mult;
		
		//printf("A%2i %1.50Lf\n", j, sum);
		//printf("%1.80Lf\n", mult);
	}
	//printf("--> %1.80Lf\n", sum);
	
	
	long double resultsum = 1. - pow((long double)sum, (long double)1);
	//long double resultsum = 1. - (sum**M);
	printf("resultsum     %4.80Lf\n", resultsum);
	cout << "resultsum=" << resultsum << endl;
	//
	
	*/
	cout << "\n---------- " << N << " maps" << endl;
	
	cout << "Probability of wrong " << K << " (or more) detection in "<< N <<" maps of " << M << " bins : \n";
	cout << "Pretrial significance" << endl;
	printf("1) p-value sum      %6.16e\n", resultsum);
	//NO cout << "It correspond to " << TMath::ErfInverse(1-resultsum)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	printf("It correspond to %1.20f (normal) sigma \n\n", ROOT::Math::gaussian_quantile_c(resultsum,1.0));
	
	
	/*long double result = 1 - TMath::Power(1 - TMath::BinomialI(pre , N, K), M);
	printf("2) p-value binomial %6.16e\n", result);
	cout << result << endl;
	printf("It correspond to %1.60f (normal) sigma \n", ROOT::Math::gaussian_quantile_c(result,1.0));*/
	
	cout << "\n---------- " << endl;
	cout << "Probability of wrong " << K << " (exactly) detection in "<< N <<" maps of 1 bins : \n"   ;
	///////long double resultK =  TMath::Binomial(N, j) * TMath::Power(pre, j) * TMath::Power(1-pre, N-j) ;
	///////printf("1) p-value      %6.16e\n", resultK);
	///////cout << resultK << endl;
	//NO resultexit=TMath::ErfInverse(1-resultK)*TMath::Sqrt(2);
	//NO cout << "It correspond to " << TMath::ErfInverse(1-resultK)*TMath::Sqrt(2) << " (normal) sigma \n" << endl;
	//////cout << "It correspond to " <<  ROOT::Math::gaussian_quantile_c(resultK,1.0) << " (normal) sigma \n" << endl;
	

	
}

int main(int argc, char** argv) {
	double TS = atof(argv[1]);
	int N  = atoi(argv[2]);
	double M = atof(argv[3]);
	int K = atoi(argv[4]);
	int hp = atoi(argv[5]);
	
	cout << "----" << endl;
	long double pippo;
	cout << sizeof(pippo);
	cout << "-----\n" << endl;
	
	mpf_t x, y;
       mpf_init (x);           /* use default precision */
       mpf_init2 (y, 256);     /* precision at least 256 bits */
    cout << sizeof(y) << endl;
	
	getPValueAndPostTrial(TS, N, M, K, hp);
	
}
