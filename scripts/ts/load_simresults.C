
Double_t distance(double ll, double bl, double lf, double bf) {
	double d1 = bl - bf;
	double d2 = ll - lf;
	
	double bl1 = TMath::TwoPi() / 4.0 - (bl * TMath::TwoPi()  / 360.0);
	double bf1 = TMath::TwoPi() / 4.0 - (bf * TMath::TwoPi()  / 360.0);
	double m4 = TMath::Cos(bl1) * TMath::Cos(bf1)  + TMath::Sin(bl1) * TMath::Sin(bf1) * TMath::Cos(d2 * TMath::TwoPi()  / 360.0);                    
	double d4 = TMath::ACos(m4) *  360.0 / TMath::TwoPi();
	return d4;
	
}

void load_simresults(TString filenameinput, int inputtype=5, double enabledistsel = 0) {
	Float_t SOURCE, L, B, TS, FLUX, FLUXERR, FLUXUL, R, EXP, SI, MINTS, FIXFLAG , CTS, CTSERR, CTSUL, TOTEXP, TOTNCOUNTS, FCN0, FCN1, EDM0, EDM1, ITER0, ITER1;
	
	Float_t GAL, ISO;
	
	double centerl = 79;
	double centerb = 0;
	
	L = B = TS = FLUX = FLUXERR = FLUXUL = R = EXP = SI = MINTS = FIXFLAG = CTS = CTSERR = CTSUL = TOTEXP = TOTNCOUNTS = FCN0 = FCN1 = EDM0 =  EDM1 = ITER0 =  ITER1 = 1;
	
	GAL=ISO=-1;
	
	Long64_t nlines;
	TTree* T;
	T = new TTree("DATA", "") ;
	
		if(inputtype == 5) {
			/*
			After R:
			double exposure = GetTotalExposure(i);
			srcout << " " << exposure*m_sources[i].GetFlux();
			srcout << " " << exposure*m_sources[i].GetFluxerr();
			srcout << " " << exposure*m_sources[i].GetFluxul();
			srcout << " " << exposure;
			srcout << " " << m_fitInfo[i].counts;

			srcout << " " << m_fitInfo[i].fcn0;
			srcout << " " << m_fitInfo[i].fcn1;
			srcout << " " << m_fitInfo[i].edm0;
			srcout << " " << m_fitInfo[i].edm1;
			srcout << " " << m_fitInfo[i].iter0;
			srcout << " " << m_fitInfo[i].iter1 << endl;
			*/
			nlines = T->ReadFile(filenameinput, "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS:GAL:ISO:R:CTS:CTSERR:CTSUL:TOTEXP:TOTNCOUNTS:FCN0:FCN1:EDM0:EDM1:ITER0:ITER1");
			cout << "******************************************* " << endl;
			cout << "SOURCE:L:B:TS:FLUX:INDEX:FIXFLAG:MINTS:GAL:ISO:R:CTS:CTSERR:CTSUL:TOTEXP:TOTNCOUNTS:FCN0:FCN1:EDM0:EDM1:ITER0:ITER1" << endl;
			cout << "******************************************* " << endl;
		}
		cout << "N lines: " << nlines << endl;
		
		if(inputtype == 5) {
			T->SetBranchAddress("SOURCE", &SOURCE);
			T->SetBranchAddress("L", &L);
			T->SetBranchAddress("B", &B);
			T->SetBranchAddress("TS", &TS);
			T->SetBranchAddress("FLUX", &FLUX);
			T->SetBranchAddress("R", &R);
			T->SetBranchAddress("CTS", &CTS);
			T->SetBranchAddress("TOTNCOUNTS", &TOTNCOUNTS);
			T->SetBranchAddress("FCN0", &FCN0);
			T->SetBranchAddress("FCN1", &FCN1);
			T->SetBranchAddress("EDM0", &EDM0);
			T->SetBranchAddress("EDM1", &EDM1);
			T->SetBranchAddress("ITER0", &ITER0);
			T->SetBranchAddress("ITER1", &ITER1);
		
			//T->SetBranchAddress("GAL", &GAL);
			//T->SetBranchAddress("ISO", &ISO);
		}
	
		TCanvas* c1 = new TCanvas;
		c1->Divide(3,3);
		
		//FCN1 = 2 ln L1
		//FCN0 = 2 ln L0
		//TS = TS = -2(ln L0 - ln L1) = 2 ln L1 - 2 ln L0 = FCN1 - FCN0
		
		TH2D* h_tsdist = new TH2D("h_tsdist", "h_tsdist", 400, 0, 400, 100, 0, 10);
		
		TH1D* h_fcn0 = new TH1D("h_fcn0", "h_fcn0", 40000, 0, 40000);
		TH1D* h_fcn1 = new TH1D("h_fcn1", "h_fcn1", 40000, 0, 40000);
		TH1D* h_fcndiff = new TH1D("h_fcndiff", "h_fcndiff", 400, 0, 400);
		TH1D* h_ts = new TH1D("h_ts", "h_ts", 400, 0, 400);
		TH1D* h_ts_thr = new TH1D("h_ts_thr", "h_ts_thr", 400, 0, 400);
		TH1D* h_fcn0_thr = new TH1D("h_fcn0_thr", "h_fcn0_thr", 40000, 0, 40000);
		TH1D* h_fcn0_thr = new TH1D("h_fcn1_thr", "h_fcn1_thr", 40000, 0, 40000);
		TH1D* h_fcndiff_thr = new TH1D("h_fcndiff_thr", "h_fcndiff_thr", 400, 0, 400);
		TH1D* h_ts_thr = new TH1D("h_ts_thr", "h_ts_thr", 400, 0, 400);
		
		TH1D* h_edm0 = new TH1D("h_emd0", "h_edm0", 10000, 0, 100);
		TH1D* h_edm0_thr = new TH1D("h_emd0_thr", "h_edm0_thr", 10000, 0, 100);
		
		TH1D* h_edm1 = new TH1D("h_emd1", "h_edm1", 10000, 0, 0.02);
		TH1D* h_edm1_thr = new TH1D("h_emd1_thr", "h_edm1_thr", 10000, 0, 0.02);
		
		TH1D* h_iter0 = new TH1D("h_iter0", "h_iter0", 1000, 0, 1000);
		TH1D* h_iter0_thr = new TH1D("h_iter0_thr", "h_iter0_thr", 1000, 0, 1000);
		
		TH1D* h_iter1 = new TH1D("h_iter1", "h_iter1", 10000, 0, 10000);
		TH1D* h_iter1_thr = new TH1D("h_iter1_thr", "h_iter1_thr", 10000, 0, 10000);
		
		int discarded = 0;
		for(Long64_t j = 0; j<nlines; j++) {
			T->GetEntry(j);
			
			if(EDM0 != 0.5) {
				cout << "* " << SOURCE << " " << L << " " << B  << " " << (FLUX) << " " << TS << " " << EDM0 << " " << EDM1 << endl;
				discarded++;
				continue;
			}
				
			double dist = distance(L, B, centerl, centerb);

			h_tsdist->Fill(TS, dist);
			h_fcn0->Fill(FCN0);
			h_fcn1->Fill(FCN1);
			h_fcndiff->Fill(2*(FCN0-FCN1));
			h_ts->Fill(TS);
			h_edm0->Fill(EDM0);
			h_edm1->Fill(EDM1);
			h_iter0->Fill(ITER0);
			h_iter1->Fill(ITER1);
			if(dist > 1) {
				h_fcn0_thr->Fill(FCN0);
				h_fcn1_thr->Fill(FCN1);
				h_fcndiff_thr->Fill(2*(FCN0-FCN1));
				h_ts_thr->Fill(TS);
				h_edm0_thr->Fill(EDM0);
				h_edm1_thr->Fill(EDM1);
				h_iter0_thr->Fill(ITER0);
				h_iter1_thr->Fill(ITER1);
				h_ts->Fill(TS);
			}
			
		}
		
		cout << "discarded: " << discarded << endl;
		
		double scf = h_fcn0->Integral();
		h_fcn0->Scale(1.0/scf);
		h_fcn0_thr->Scale(1.0/scf);
		scf = h_fcn1->Integral();
		h_fcn1->Scale(1.0/scf);
		h_fcn1_thr->Scale(1.0/scf);
		scf = h_fcndiff->Integral();
		h_fcndiff->Scale(1.0/scf);
		h_fcndiff_thr->Scale(1.0/scf);
		scf = h_ts->Integral();
		h_ts->Scale(1.0/scf);
		h_ts_thr->Scale(1.0/scf);
		
		scf = h_edm0->Integral();
		h_edm0->Scale(1.0/scf);
		h_edm0_thr->Scale(1.0/scf);
		scf = h_edm1->Integral();
		h_edm1->Scale(1.0/scf);
		h_edm1_thr->Scale(1.0/scf);
		
		scf = h_iter0->Integral();
		h_iter0->Scale(1.0/scf);
		h_iter0_thr->Scale(1.0/scf);
		scf = h_iter1->Integral();
		h_iter1->Scale(1.0/scf);
		h_iter1_thr->Scale(1.0/scf);
		
		c1->cd(1);
		gPad->SetLogy();
		h_fcn0->Draw();
		h_fcn0_thr->SetLineColor(kRed);
		h_fcn0_thr->Draw("SAME");
		c1->cd(2);
		gPad->SetLogy();
		h_fcn1->Draw();
		h_fcn1_thr->SetLineColor(kRed);
		h_fcn1_thr->Draw("SAME");
		c1->cd(3);
		gPad->SetLogy();
		h_fcndiff->Draw();
		h_ts->SetLineColor(kRed);
		h_ts->Draw("SAME");
		c1->cd(4);
		gPad->SetLogy();
		T->Draw("EDM0");
		h_edm0->Draw();
		h_edm0_thr->SetLineColor(kRed);
		h_edm0_thr->Draw("SAME");
		
		c1->cd(5);
		gPad->SetLogy();
		T->Draw("EDM1");
		h_edm1->Draw();
		h_edm1_thr->SetLineColor(kRed);
		h_edm1_thr->Draw("SAME");
		c1->cd(6);
		gPad->SetLogy();
		h_iter0->Draw();
		h_iter0_thr->SetLineColor(kRed);
		h_iter0_thr->Draw("SAME");
		c1->cd(7);
		gPad->SetLogy();
		h_iter1->Draw();
		h_iter1_thr->SetLineColor(kRed);
		h_iter1_thr->Draw("SAME");
		c1->cd(8);
		gPad->SetLogy();
		T->Draw("FCN1", "TS>9");
		
		c1->cd(9);
		gPad->SetLogy();
		gPad->SetLogx();
		h_tsdist->Draw("BOXCOL");
		
		
}