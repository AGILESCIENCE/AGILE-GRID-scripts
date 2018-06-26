#define MULTIBIN 1.0
Double_t sumchi2Ndeltafunc(Double_t *x, Double_t *par)
{
	Double_t f = par[0]*((x[0]<=1.0/MULTIBIN)?1:0) + par[1] * ((x[0]>1.0/MULTIBIN)?1:0) * exp(-x[0]/2.0) * pow(x[0], par[2]/2.0 - 1) / ( pow(2, par[2]/2.0 ) * TMath::Gamma(par[2]/2.0) ) ;
	return f;
}

void m(Double_t x, Double_t *par) {
	int max_cumTSdisbinx = 36;
	TF1* f4sf = new TF1("sumchi2(N)delta", sumchi2Ndeltafunc, 0, max_cumTSdisbinx,3);
	f4sf->SetParName(0, "delta");
	f4sf->SetParameter(0, par[0]);
	f4sf->SetParName(1, "eta");
	f4sf->SetParLimits(1, 0, 1);
	f4sf->SetParameter( 1, par[1]);
	f4sf->SetParName(2, "N");
	f4sf->SetParLimits(2, 0, 10);
	f4sf->SetParameter( 2, par[2]);
	f4sf->SetLineColor(kGreen);

	double y = f4sf->Eval(x);
	cout << y << endl;

	double integral = f4sf->Integral(0, 35);
	cout << integral << endl;
}
