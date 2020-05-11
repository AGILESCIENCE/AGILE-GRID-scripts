
void FermiSensitivity(const TString FileNameSens= "sens.txt", TString drawString="APS", Int_t colorMap = 0) {
	Int_t len = 0;

	TString FileNameSensTemp = FileNameSens;
	Int_t indp = FileNameSensTemp.Index("[.]txt",&len,0);

	ifstream inputFile ;
	inputFile.open(FileNameSens.Data()) ;
	if(!inputFile.is_open())
		return;

	const Int_t LineLen = 15120;
	char line[LineLen];

	const Double_t XScaleFactor = 1.05;
	const Double_t YScaleFactor = 1.1;

	const Int_t LenCurve = 100;
	Double_t xmin = 0;
	Double_t xmax =0 ;
	Double_t x2[LenCurve];
	Double_t y2[LenCurve];
	Double_t ex2[LenCurve];
	Double_t ey2[LenCurve];
	Double_t x3[LenCurve];
	Double_t y3[LenCurve];
	Double_t ex3[LenCurve];
	Double_t ey3[LenCurve];

	Double_t StartTime = 0;
	Double_t EndTime = 0;
	Double_t Time = 0;
	Double_t Time2 = 0;
	Int_t nl = -1;
	Double_t sens4;
	Double_t sens5;
	Double_t maxflux= -1000;
	Double_t minflux = 1000;
	Int_t n2 = 0;
	Int_t n3 = 0;
	Double_t DeltaTime = 4;

	char slash;
	while( !(inputFile.eof()) ) {
		inputFile.getline(line,LineLen);

		if( strlen(line) <= 0)
			break;
		nl++;
		TString SLine(line);
		TObjArray* tokenobj = SLine.Tokenize(" ") ;


		TObjString* stt = (TObjString *) (tokenobj->At(0)) ;
		Time = (stt->GetString()).Atof();
		stt = (TObjString *) (tokenobj->At(1)) ;
		Time2 = (stt->GetString()).Atof();
		
		DeltaTime = Time2 -Time;
		
		if(nl == 1) {
			x2[0] = Time - DeltaTime*0.5;
			xmin = Time - DeltaTime*2;
		}
		xmax = Time2 + DeltaTime;

		TObjString* stf = (TObjString *) (tokenobj->At(2)) ;
		sens4  = (stf->GetString()).Atof();

		TObjString* stfe = (TObjString *) (tokenobj->At(3)) ;
		sens5  = (stfe->GetString()).Atof();
		
		cout << setprecision(15)  << Time << " " << Time2 << " " << DeltaTime << " " << sens4 << " " << sens5 << endl;

		x2[n2] = Time + DeltaTime*0.5;
		ex2[n2] = DeltaTime*0.5;
		y2[n2] = sens4;
		ey2[n2] = 0;
		n2++;

		x3[n3] = Time + DeltaTime*0.5;
		ex3[n3] = DeltaTime*0.5;
		y3[n3] = sens5;
		ey3[n3] = 0;
		n3++;

		if(nl==0) {
			StartTime = Time;
		}

		minflux = TMath::Min(minflux,sens4);
		maxflux = TMath::Max(maxflux,sens5);
	}
	EndTime = Time + DeltaTime;
	inputFile.close();

	TGraphErrors *g3 = new TGraphErrors(n3,x3,y3,ex3,ey3);
	TString title = "";
	g3->SetTitle(title.Data());
	if(colorMap) {
		g3->SetMarkerColor(colorMap);
		g3->SetLineColor(colorMap);
	} else {
		g3->SetMarkerColor(kBlue);
		g3->SetLineColor(kBlue);
	}
	//g3->SetMarkerStyle(20);

	TGraphErrors *g2 = new TGraphErrors(n2,x2,y2,ex2,ey2);
	if(colorMap) {
		g2->SetMarkerColor(colorMap);
		g2->SetLineColor(colorMap);
	}
	else {
		g2->SetMarkerColor(kBlack);
		g2->SetLineColor(kBlack);
	}
	//g2->SetMarkerStyle(20);

	g3->Draw(drawString);
	
	Double_t MinHisto = (minflux > 0) ? minflux/YScaleFactor : minflux*YScaleFactor ;
	g3->GetHistogram()->SetMinimum(minflux - (maxflux - minflux)*(YScaleFactor-1));
	g3->GetHistogram()->SetMaximum(maxflux + (maxflux - minflux)*(YScaleFactor-1));
	TAxis *xaxis = g3->GetXaxis();
	xaxis->SetTitle("MJD [days]");
	xaxis->SetTitleSize(0.050);
	xaxis->SetTitleOffset(0.9);
	TAxis *yaxis = g3->GetYaxis();
	TString YAxisTitle = TString::Format("Flux [10^{-8}ph cm^{-2}s^{-1}]");
	yaxis->SetTitle(YAxisTitle.Data());
	yaxis->SetTitleSize(0.050);
	yaxis->SetTitleOffset(0.9);

	xaxis->SetLimits(xmin,xmax);

	
	g3->Draw("PS");
	g2->Draw("PS");

	return;
	
	
}

//colorMap=0 -> differnt color based on TS
//colotMap > 0 use the same color for all
void AgileDrawLLightCurve(const TString FileNameLightCurve= "analysisres.resfull", TString drawString="APS", Int_t colorMap = 0, Int_t markerStyle=20, Double_t flux_ave=0, Double_t flux_ave_err=0)
{
	Int_t len = 0;

	TString FileNameLightCurveTemp = FileNameLightCurve;
	Int_t indp = FileNameLightCurveTemp.Index("[.]resfull",&len,0);

	TString sourcename;

	ifstream inputFile ;
	inputFile.open(FileNameLightCurve.Data()) ;
	if(!inputFile.is_open())
		return;

	const Int_t LineLen = 15120;
	char line[LineLen];

	Double_t UPArrowLength = 25.;
	const Double_t XScaleFactor = 1.05;
	const Double_t YScaleFactor = 1.1;

	const Int_t LenCurve = 100;
	Double_t xmin = 0;
	Double_t xmax =0 ;
	Double_t x2[LenCurve];
	Double_t y2[LenCurve];
	Double_t ex2[LenCurve];
	Double_t ey2[LenCurve];

	Double_t x3[LenCurve];
	Double_t y3[LenCurve];
	Double_t ex3[LenCurve];
	Double_t ey3[LenCurve];

	Double_t x4[LenCurve];
	Double_t y4[LenCurve];
	Double_t ex4[LenCurve];
	Double_t ey4[LenCurve];

	Double_t x5[LenCurve];
	Double_t y5[LenCurve];
	Double_t ex5[LenCurve];
	Double_t ey5[LenCurve];

	Double_t StartTime = 0;
	Double_t EndTime = 0;
	Double_t Time = 0;
	Double_t Time2 = 0;
	Int_t nl = -1;
	Double_t TS[1000];
	Double_t temp;
	Double_t temp2;
	Double_t flux;
	Double_t fluxerr;
	Double_t cts;
	Double_t ctserr;
	Double_t uflux;
	Double_t maxflux = -1000;
	Double_t minflux = 1000;
	Int_t n2 = 0;
	Int_t n3 = 0;
	Int_t n4 = 0;
	Int_t n5 = 0;
	Double_t DeltaTime = 4;

	char slash;
	while( !(inputFile.eof()) ) {
		inputFile.getline(line,LineLen);

		if( strlen(line) <= 0)
			break;
		nl++;
		TString SLine(line);
		TObjArray* tokenobj = SLine.Tokenize(" ") ;

		TObjString* ss = (TObjString *) (tokenobj->At(1)) ;
		sourcename  = (ss->GetString());
		cout << sourcename << endl;

		TObjString* sts = (TObjString *) (tokenobj->At(2)) ;
		TS[nl]  = (sts->GetString()).Atof();

		/*
		TObjString* stt = (TObjString *) (tokenobj->At(6)) ;
		Time = (stt->GetString()).Atof();
		Int_t lent = 0;
		Int_t indt = (stt->GetString()).Index("/[0-9.]*",&lent,0)+1;
		lent--;
		Time2 = TString((stt->GetString())(indt,lent)).Atof();
		*/

		TObjString* stt = (TObjString *) (tokenobj->At(41)) ;
		Time = (stt->GetString()).Atof();
		stt = (TObjString *) (tokenobj->At(42)) ;
		Time2 = (stt->GetString()).Atof();
		
		DeltaTime = Time2 -Time;
		
		if(nl == 1) {
			x2[0] = Time - DeltaTime*0.5;
			x3[0] = Time - DeltaTime*0.5;
			x4[0] = Time - DeltaTime*0.5;
			x5[0] = Time - DeltaTime*0.5;
			xmin = Time - DeltaTime*2;
		}
		xmax = Time2 + DeltaTime;

		TObjString* stf = (TObjString *) (tokenobj->At(21)) ;
		flux  = (stf->GetString()).Atof();
		flux *= 1e8	;

		TObjString* stfe = (TObjString *) (tokenobj->At(22)) ;
		fluxerr  = (stfe->GetString()).Atof();
		fluxerr *= 1e8;

		TObjString* stfu = (TObjString *) (tokenobj->At(25)) ;
		uflux  = (stfu->GetString()).Atof();
		uflux *= 1e8;

		TObjString* stcts = (TObjString *) (tokenobj->At(15)) ;
		cts  = (stcts->GetString()).Atof();

		TObjString* stctserr = (TObjString *) (tokenobj->At(16)) ;
		ctserr  = (stctserr->GetString()).Atof();
		
		cout << TS[nl] << " " << setprecision(15)  << Time << " " << Time2 << " " << DeltaTime << " " << flux << " " << fluxerr << " " << uflux << endl;

		if(TS[nl] < 3) {
			x2[n2] = Time + DeltaTime*0.5;
			ex2[n2] = DeltaTime*0.5;
			y2[n2] = uflux;
			ey2[n2] = 0;
			n2++;
		}
		else if(TS[nl] >= 3 && TS[nl] < 4) {
			x3[n3] = Time + DeltaTime*0.5;
			ex3[n3] = DeltaTime*0.5;
			y3[n3] = flux;
			ey3[n3] = fluxerr;
			n3++;
		}
		else if(TS[nl] > 4 && TS[nl] < 5) {
			x4[n4] = Time + DeltaTime*0.5;
			ex4[n4] = DeltaTime*0.5;
			y4[n4] = flux;
			ey4[n4] = fluxerr;
			n4++;
		}
		else if(TS[nl] > 5 ) {
			x5[n5] = Time + DeltaTime*0.5;
			ex5[n5] = DeltaTime*0.5;
			y5[n5] = flux;
			ey5[n5] = fluxerr;
			n5++;
		}
		if(nl==0) {
			StartTime = Time;
		}
		Double_t MinPoint = (TS[nl] > 3) ? flux-fluxerr : uflux - UPArrowLength ;
		minflux = TMath::Min(minflux,MinPoint);
		maxflux = TMath::Max(maxflux,flux+fluxerr);
	}
	EndTime = Time + DeltaTime;
	minflux = TMath::Min(minflux,flux_ave-flux_ave_err);
	maxflux = TMath::Max(maxflux,flux_ave+flux_ave_err);
	UPArrowLength = (maxflux - minflux)/6.;
	inputFile.close();

	TGraphErrors *g3 = new TGraphErrors(n3,x3,y3,ex3,ey3);
	//TString title = "Light curve of " + sourcename;
	TString title = "";
	g3->SetTitle(title.Data());
	if(colorMap) {
		g3->SetMarkerColor(colorMap);
		g3->SetLineColor(colorMap);
	} else {
		g3->SetMarkerColor(kBlue);
		g3->SetLineColor(kBlue);
	}
	g3->SetMarkerStyle(markerStyle);

	TGraphErrors *g2 = new TGraphErrors(n2,x2,y2,ex2,ey2);
	if(colorMap) {
		g2->SetMarkerColor(colorMap);
		g2->SetLineColor(colorMap);
	}
	else {
		g2->SetMarkerColor(kBlack);
		g2->SetLineColor(kBlack);
	}
	g2->SetMarkerStyle(markerStyle);

	TGraphErrors *g4 = new TGraphErrors(n4,x4,y4,ex4,ey4);
	if(colorMap) {
		g4->SetMarkerColor(colorMap);
		g4->SetLineColor(colorMap);
	}
	else {
		g4->SetMarkerColor(kOrange+7);
		g4->SetLineColor(kOrange+7);		
	}
	g4->SetMarkerStyle(markerStyle);

	TGraphErrors *g5 = new TGraphErrors(n5,x5,y5,ex5,ey5);
	if(colorMap) {
		g5->SetMarkerColor(colorMap);
		g5->SetLineColor(colorMap);
	}
	else {
		g5->SetMarkerColor(kRed);
		g5->SetLineColor(kRed);
	}
	g5->SetMarkerStyle(markerStyle);

	g3->Draw(drawString);
	Double_t MinHisto = (minflux > 0) ? minflux/YScaleFactor : minflux*YScaleFactor ;
	g3->GetHistogram()->SetMinimum(minflux - (maxflux - minflux)*(YScaleFactor-1));
	g3->GetHistogram()->SetMaximum(maxflux + (maxflux - minflux)*(YScaleFactor-1));
	TAxis *xaxis = g3->GetXaxis();
	xaxis->SetTitle("MJD [days]");
	xaxis->SetTitleSize(0.050);
	xaxis->SetTitleOffset(0.9);
	TAxis *yaxis = g3->GetYaxis();
	TString YAxisTitle = TString::Format("Flux [10^{-8}ph cm^{-2}s^{-1}]");
	yaxis->SetTitle(YAxisTitle.Data());
	yaxis->SetTitleSize(0.050);
	yaxis->SetTitleOffset(0.9);

	/*
	Double_t xmin = TMath::Min(x3[0]-ex3[0],x4[0]-ex4[0]);
	xmin = TMath::Min(xmin,x5[0]-ex5[0]);

	Double_t xmax = TMath::Max(x3[n3-1]+ex3[0],x4[n4-1]+ex4[0]);
	xmax = TMath::Max(xmax,x5[n5-1]+ex5[0]);
	xmax = TMath::Max(xmax,x2[n5-1]+ex2[0]);
	Double_t xmin_temp = xmin - (xmax - xmin)*(XScaleFactor-1) ;
	Double_t xmax_temp = xmax + (xmax - xmin)*(XScaleFactor-1) ;
	xmin = xmin_temp ;
	xmax = xmax_temp ;
	*/
	xaxis->SetLimits(xmin,xmax);

	if(flux_ave > 0) {
		TLine* fluxupdown[2] ;
		for (Int_t i=0; i<2; i++) {
			Double_t yline = flux_ave+(2*i-1)*flux_ave_err;
			fluxupdown[i] = new TLine(xmin,yline,xmax,yline);
			fluxupdown[i]->SetLineColor(kGreen+1);
			fluxupdown[i]->SetLineStyle(9);
			fluxupdown[i]->SetLineWidth(3);
			fluxupdown[i]->Draw();
		}
	}
	
	g3->Draw("PS");
	g2->Draw("PS");
	g4->Draw("PS");
	g5->Draw("PS");

	Double_t y2l;
	for (Int_t i=0; i<n2; i++) {
		if(ey2[i] == 0) {
			y2l = y2[i]-UPArrowLength;
			TArrow* a = new TArrow(x2[i],y2[i],x2[i],y2l,0.02,">");
			if(colorMap)
				a->SetLineColor(colorMap);
			else
				a->SetLineColor(kBlack);
			a->Draw();
		}
	}

	

	return;
}

void AgileCatalog2GLLightCurve(const TString FileNameLightCurve= "analysisres.resfull", Double_t flux_ave=0, Double_t flux_ave_err=0) {
	TString sourcename = "IGRJ17354-3255";
	const TString CanvasTitle = "Light curve of "+sourcename;
	TCanvas* c1 = new TCanvas("c1",CanvasTitle.Data(),200,10,700,500);
	c1->SetGrid();
	
	AgileDrawLLightCurve("analysisres14400.resfull", "APS", kRed, 20, flux_ave, flux_ave_err);
	AgileDrawLLightCurve("analysisres21600.resfull", "PS", kBlue);
	FermiSensitivity("fermisens14400.txt", "PS", kBlack);
	
	TString Filename = "lc_"+sourcename;
	Filename += ".png";
	c1->SaveAs(Filename.Data());
}
