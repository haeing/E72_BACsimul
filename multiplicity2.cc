void multiplicity2(){

  gStyle->SetOptStat(0);

  TFile *file =  new TFile("build/version2.root","read");
  TTree* data = (TTree*)file->Get("tree");

  Int_t nhMppc;
  Double_t x[1000];
  Double_t y[1000];
  Double_t z[1000];

  Double_t evtposx;
  Double_t evtposy;
  

  data->SetBranchAddress("nhMppc",&nhMppc);
  data->SetBranchAddress("evtposx",&evtposx);
  data->SetBranchAddress("evtposy",&evtposy);
  data->SetBranchAddress("mppcposx",x);
  data->SetBranchAddress("mppcposy",y);
  data->SetBranchAddress("mppcposz",z);

  Double_t total = data->GetEntries();
  
  Int_t numx = 4;
  Int_t numy = 4;
  Int_t numz = 2;
  Int_t multi[4][4][2];
  Int_t cell_num[4][4][2];

  Double_t one = 6;
  Int_t thre = 8;
  Int_t result;
  Double_t effi_num[5];

  Int_t num_check=10;
  Int_t multi_effi[num_check];
  Int_t effi_thre[num_check];
  for(int i=0;i<num_check;i++){
    multi_effi[i]=0;
    effi_thre[i]=0;
  }

  TH1D* hist_multi = new TH1D("hist_multi","hist_multi",32,0,32);
  TH1D* hist1 = new TH1D("hist1","hist1",25,0,25);

  TH2D* evt = new TH2D("evt","evt",50,-40,40,50,-40,40);

  TGraph *effi = new TGraph();
  //TGraph *effi_tight = new TGraph();

  Int_t multi_thre = 0;
  for(int n=0;n<total;n++){
    data->GetEntry(n);
    for(int nx=0;nx<numx;nx++){
      for(int ny=0;ny<numy;ny++){
	for(int nz=0;nz<numz;nz++){
	  multi[nx][ny][nz]=0;
	  cell_num[nx][ny][nz]=0;
	}
      }
    }
  
  result=0;
  multi_thre = 0;
  for(int i=0;i<num_check;i++)multi_effi[i]=0;
    //result_tight=0;
    //if(evtposx>-40&&evtposx<40&&evtposy>-40&&evtposy<40)total_tight++;

  for(int i=0;i<nhMppc;i++){
    for(int nx=0;nx<numx;nx++){
      for(int ny=0;ny<numy;ny++){
	for(int nz=0;nz<numz;nz++){
	  if(-12+one*nx<=x[i]&&-12+one*(nx+1)>x[i]&&-12+one*ny<=y[i]&&-12+one*(ny+1)>y[i]&&-200+nz*200<=z[i]&&-200+(nz+1)*200>z[i]){
	      multi[nx][ny][nz]=1;
	      cell_num[nx][ny][nz]++;
	  }
	}
      }
    }
  }
        
  
  for(int nx=0;nx<numx;nx++){
    for(int ny=0;ny<numy;ny++){
      for(int nz=0;nz<numz;nz++){
	result+=multi[nx][ny][nz];
	if(cell_num[nx][ny][nz]>0)hist1->Fill(cell_num[nx][ny][nz]);
	if(cell_num[nx][ny][nz]>thre)multi_thre++;
	for(int l=1;l<num_check;l++){
	    if(cell_num[nx][ny][nz]>l)multi_effi[l-1]++;
	  }
	
	
      }
    }
  }
  hist_multi->Fill(multi_thre);
  
  //hist_multi->Fill(result);

  

  
    
  if(multi_thre>0)evt->Fill(evtposx,evtposy);
  for(int f=0;f<num_check;f++){
    if(multi_effi[f]>0)effi_thre[f]++;
  }
  multi_thre=0;
  }
    /*
    for(int l=0;l<5;l++){
      thre=l;
      if(result>=thre)effi_num[l]++;
      if(evtposx>-40&&evtposx<40&&evtposy>-40&&evtposy<40){
	if(result_tight>=thre)effi_num_tight[l]++;
      }

    }
  }
    

  Double_t divide;
  for(int i=0;i<5;i++){
    effi->AddPoint(i,effi_num[i]/total);
    effi_tight->AddPoint(i,effi_num_tight[i]/total_tight);
  }
    

  effi->SetMarkerStyle(21);
  effi_tight->SetMarkerStyle(21);
  effi_tight->SetMarkerColor(kRed);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(effi);
  mg->Add(effi_tight);
  mg->SetTitle("Efficiency;Thre. of multiplicity;efficiency");

  TLegend *l1 = new TLegend(0.8,0.5,0.48,0.6);
  l1->AddEntry(effi,"whole range");
  l1->AddEntry(effi_tight,"8 cm * 8 cm");
    */

  TGraph *graph_multi = new TGraph();
  for(int i=0;i<num_check;i++){
    std::cout<<effi_thre[i]<<std::endl;
    graph_multi->AddPoint(1+i,effi_thre[i]/total);
  }
  graph_multi->SetMarkerStyle(21);
  graph_multi->SetTitle("Efficiency;threshold of # of photons;efficiency");

  
  hist_multi->SetTitle("Whole range;multiplicity;n");
  //hist_multi_tight->SetTitle("8 cm * 8 cm;multiplicity;n");
  //hist1->SetTitle("# of photons at each cell;# of photons;n");
  TCanvas *c1 = new TCanvas("c1","c1",800,650);
  //c1->Divide(2);
  c1->cd(1);
  hist_multi->Draw();
  //c1->cd(2);
  //hist_multi_tight->Draw();
  
  
  TCanvas *c2 = new TCanvas("c2","c2",800,650);
  c2->cd();
  hist1->Draw();

  evt->SetTitle("Position of the beam starting point;x [mm];y [mm]");
  TCanvas *c3 = new TCanvas("c3","c3",800,650);
  c3->cd();
  evt->Draw("colz");
  /*
  TCanvas *c4 = new TCanvas("c4","c4",800,650);
  c4->cd();
  mg->Draw("AP");
  l1->Draw();
  */

  TCanvas *c5 = new TCanvas("c5","c5",800,650);
  c5->cd();
  graph_multi->Draw("AP");
}
