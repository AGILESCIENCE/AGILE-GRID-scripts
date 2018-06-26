

/*
 hp
 1) chi2_1 * delta
 2) chi2_3 * delta
 3) ICL C1 Cygnus region FF=3, t_ICL=9, minTS=4
 4) ICL C1 Empty Galactic region FF=3, t_ICL=6, minTS=4
 5) FF=1 Galactic region
 6) ALL Cygnus Galactic region FF=3 minTS=4
 7) ALL Empty Galactic region FF=3 minTS=4
 */

void getTS(double pvalue = 1e-2, int N = 180, int hp = 1, double delta=0.5) {
	double TS = 1.0;
	double step = 0.0001;
	double pre = 10;
	cout << hp << " " << TS << " " << step << " " << pre << " " << pvalue << endl;
	
	while(pre > pvalue) {
		TS += step;
		
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
		
		
		//cout << pre << "\t";
		//cout << TS << " " << pre << " " << pvalue << endl;
	}
	
	double post = 1 - TMath::Power(1-pre, N);
	if(post == 0)
		post = pre * N;
	
	cout << "TS = " << TS << endl;
	cout << "p-value pre-trial  = " << pre << endl;
	cout << "p-value post-trial = " << post << endl;
	
	//double sigma_pre = TMath::ErfcInverse(pre)*TMath::Sqrt(2);
	//double sigma_post = TMath::ErfcInverse(post)*TMath::Sqrt(2);
	double sigma_pre = ROOT::Math::gaussian_quantile_c(pre,1.0);
	double sigma_post = ROOT::Math::gaussian_quantile_c(post,1.0);
	
	cout << "sigma pre-trial  = " << sigma_pre  << endl;
	cout << "sigma post-trial = " << sigma_post  << endl;

	
}