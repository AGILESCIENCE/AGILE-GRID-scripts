void DrawGraphArrow(Int_t n, Color_t, const Double_t *, const Double_t *, const Double_t *,
     const Double_t *, const Double_t *, const Double_t *);

TString title ;

void Fig10AgileCatalog2AGLSpectra(Int_t NFig = 0)
{
   const Int_t NFigure = 7;
   const Int_t NPointAGILE = 7;
   const Double_t xAGILE[NPointAGILE] = {38.7,70.7,173.2,547.7,1732.0,5477.2,22360.6};
   Double_t yAGILE[NFigure][NPointAGILE] = {{3.99e-6,4.38e-6,4.98e-6,3.11e-6,1.10e-6,2.18e-7,1.19e-9},
                                                  {1.77e-6,7.81e-7,1.52e-6,1.23e-6,3.61e-7,1.19e-7,2.58e-9},
                                                  {9.49e-7,7.47e-7,6.77e-7,3.72e-7,1.03e-7,2.14e-8,1.68e-10},
                                                  {8.29e-7,5.35e-7,7.27e-7,5.04e-7,2.38e-7,4.74e-8,6.72e-10},
                                                  {7.09e-7,2.21e-7,2.66e-7,2.61e-7,1.50e-7,3.32e-8,1.73e-10},
                                                  {0.00000383605,0.00000218654,0.00000134517,4.25609e-7,2.96438e-8,3.48939e-9,1.36713e-8},
                                                  {0.00000113037,7.44895e-7,5.21792e-7,1.64637e-7,2.65843e-8,8.32692e-9,3.0905e-9}
                                                  };
   Double_t exlAGILE[NPointAGILE] ;
   Double_t exhAGILE[NPointAGILE] ;
   const Double_t eylAGILE[NFigure][NPointAGILE] = {{3.62e-7,1.56e-7,8.83e-8,4.83e-8,3.93e-8,2.83e-8,1.30e-10},
                                                   {0.00e-7,1.40e-7,7.11e-8,3.89e-8,2.80e-8,2.99e-8,9.41e-10},
                                                   {0.00e-7,7.48e-8,3.45e-8,1.46e-8,9.99e-9,8.66e-9,2.46e-10},
                                                   {0.00e-7,1.07e-7,4.80e-8,2.14e-8,2.06e-8,1.61e-8,2.20e-10},
                                                   {0.00e-7,7.19e-8,2.56e-8,1.27e-8,1.39e-8,1.17e-8,4.74e-11},
                                                   {2.79602e-7,1.03665e-7,4.66108e-8,1.91148e-8,4.12405e-9,1.4307e-9,0.00e-7},
                                                   {2.01993e-7,5.98697e-8,1.97651e-8,6.41949e-9,3.05404e-9,3.3232e-9,1.04981e-9}
                                                   };
   Double_t eyhAGILE[NFigure][NPointAGILE] ;
   const Double_t xboundary[NPointAGILE+1] = {30.,50.,100.,300.,1000.,3000.,10000.,50000.};

   const Int_t NPointFERMI = 7;
   const Double_t yFERMI[NFigure][NPointFERMI] = {{5.37e-6,3.22e-6,1.07e-6,2.16e-7,1.14e-8},
                                                 {1.78e-6,1.37e-6,5.72e-7,1.17e-7,3.28e-9},
                                                 {6.93e-7,3.70e-7,1.02e-7,1.52e-8,4.35e-10},
                                                 {7.71e-7,4.44e-7,1.52e-7,3.41e-8,2.51e-9},
                                                 {2.70e-7,2.11e-7,8.34e-8,1.52e-8,3.26e-10},
                                                 {1.6458885e-06,5.4454074e-07,1.2553377e-07,2.5455243e-08,2.2631892e-09},
                                                 {1.7635266e-06,5.0496607e-07,8.9191275e-08,1.4981598e-08,1.0909111e-09}
                                                 };
   const Double_t eylFERMI[NFigure][NPointFERMI] = {{1.37e-9,2.29e-9,2.52e-9,1.32e-9,1.05e-10},
                                                   {4.65e-9,1.01e-9,1.37e-9,6.95e-10,2.88e-11},
                                                   {1.72e-10,2.60e-10,2.38e-10,8.63e-11,4.03e-12},
                                                   {1.96e-10,3.16e-10,3.61e-10,2.13e-10,2.32e-11},
                                                   {7.06e-11,1.55e-10,1.98e-10,8.88e-11,3.02e-12},
                                                   {4.65e-9,1.01e-9,1.37e-9,6.95e-10,0.00e-7},
                                                   {4.65e-9,1.01e-9,1.37e-9,6.95e-10,0.00e-7}
                                                   };
   Double_t eyhFERMI[NFigure][NPointFERMI] ;

   printf("NFig %d \n", NFig);

   NFig--;
   if(NFig<0 || NFig >= NFigure) {
     printf("NFig out of range\n");
     return -1;
   }

   TString sourcename ;
   switch(NFig) {
      case 0:
        sourcename = "Vela";
        break;
      case 1:
        sourcename = "Geminga";
        break;
      case 2:
        sourcename = "J2021+4029";
        break;
      case 3:
        sourcename = "J1710-4429";
        break;
      case 4:
        sourcename = "J1836+5924";
        break;
      case 5:
          sourcename = "Crab";
          break;
      case 6:
          sourcename = "3C 454.3";
          break;
   }
   title = sourcename + " spectrum";

   TCanvas* c1 = new TCanvas("c1",title.Data(),200,10,700,500);
   c1->SetGrid();
   c1->SetLogx();
   c1->SetLogy();

   for(Int_t nf = 0; nf < NFigure; nf++) {
   for(Int_t n = 0; n < NPointAGILE; n++) {
     eyhAGILE[nf][n] = eylAGILE[nf][n];
     exlAGILE[n] = xAGILE[n] - xboundary[n];
     exhAGILE[n] = xboundary[n+1] - xAGILE[n] ;
   }
   for(Int_t n = 0; n < NPointFERMI; n++) {
     eyhFERMI[nf][n] = eylFERMI[nf][n];
   }
   }

   //remove edpcorrection=0.75
   Bool_t corr = true;
   if(corr) {
     yAGILE[NFig][4] = yAGILE[NFig][4] / 0.75;
     yAGILE[NFig][5] = yAGILE[NFig][5] / 0.75;
     yAGILE[NFig][6] = yAGILE[NFig][6] / 0.75;
   }
   const Int_t FAOffset = 2;


   DrawGraphArrow(NPointAGILE-1,kBlack, xAGILE,yAGILE[NFig],exlAGILE,exhAGILE,eylAGILE[NFig],eyhAGILE[NFig]);
   DrawGraphArrow(NPointFERMI-1,kRed,   xAGILE+FAOffset,yFERMI[NFig],exlAGILE+FAOffset,
        exhAGILE+FAOffset,eylFERMI[NFig],eyhFERMI[NFig]);


   TString Filename = "Spectrum"+sourcename+"Fig11";
   TString letter = TString::Format("%c",'a'+NFig);
   Filename += letter+".png";
   c1->SaveAs(Filename.Data());
//   c1->Save();

  return;
}

void DrawGraphArrow(Int_t n, Color_t color, const Double_t *x, const Double_t *y, const Double_t *exl, const Double_t *exh, const Double_t *eyl, const Double_t *eyh)
{
   TGraphAsymmErrors *g = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
   g->SetTitle(title.Data());
   g->SetMarkerColor(color);
   g->SetMarkerStyle(20);
   if(color == kBlack)
     g->Draw("APS");
   else
     g->Draw("PS");
   TAxis *xaxis = g->GetXaxis();
   xaxis->SetTitle("Energy [MeV]");
   xaxis->SetTitleSize(0.050);
   xaxis->SetTitleOffset(0.9);

   TAxis *yaxis = g->GetYaxis();
   TString YAxisTitle = TString::Format("Flux[ph cm^{-2}s^{-1}]");
   yaxis->SetTitle(YAxisTitle.Data());
   yaxis->SetTitleSize(0.050);
   yaxis->SetTitleOffset(0.9);

   if(color == kBlack) {
   TArrow *a;
   Double_t yl;
   for (Int_t i=0; i<n; i++) {
      if(eyl[i] == 0) {
        yl = TMath::Log10(y[i])-0.5;
        yl = TMath::Power(10,yl);
        a = new TArrow(x[i],y[i],x[i],yl,0.02,">");
        a->Draw();
      }
   }
   }

   return;
}
