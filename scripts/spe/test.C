void test() {
	int m_mapCount = 5;
	double edges2[6] = {100, 300, 1000, 3000, 10000, 50000};
	TH1D gr2("spectra", "spectra", m_mapCount, edges2);
	for (int m=0; m<=m_mapCount; ++m)
		cout << edges2[m] << endl;
	gr2.SetBinContent(1, 0.671431);
	gr2.SetBinContent(2, 1.12314);
	gr2.SetBinContent(3, 0.94071);
	gr2.SetBinContent(4, 0.655919);
	gr2.SetBinContent(5, 0.751574);
	for (int m=0; m<m_mapCount; ++m) {
		cout << "val " << m+1 << " " << gr2.GetBinContent(m+1) << endl;
	}
	
	//TF1 f3("PL", "[0] * x^(-[1])", 0, m_mapCount);
	TF1 f3("PL", "[0] ", 0, m_mapCount);
	
	gr2.Fit(&f3, "");

	for(int i=0; i<m_mapCount; i++)
		cout << edges2[i] + (edges2[i+1]-edges2[i])/2 << " " << f3(edges2[i] + (edges2[i+1]-edges2[i])/2) << endl;
}
