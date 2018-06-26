



void getSigma(double sigma = 5) {
	
	double step = 0.00001;
	double pre = 0.1;
	double sigmac = 0; 
	
	cout << sigmac << " " << sigma << endl;
	while(sigmac <= sigma) {
		pre -= step;
		sigmac = ROOT::Math::gaussian_quantile_c(pre,1.0);
		//cout << sigmac << endl;
		
	}
	
	
	cout << "p-value pre-trial  = " << pre << endl;
	
	
	double sigma_pre = ROOT::Math::gaussian_quantile_c(pre,1.0);
	
	cout << "sigma pre-trial  = " << sigma_pre  << endl;
	
	
}