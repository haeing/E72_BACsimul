void multiplicity3(){

  gStyle->SetOptStat(0);

  int N = 3;
  TFile *file[N];
  file[0] =  new TFile("build/v1_110.root","read");
  file[1] =  new TFile("build/v2_110.root","read");
  file[2] =  new TFile("build/v4_110.root","read");
  TTree *data[N];
  for(int i=0;i<N;i++){
    data[i] = (TTree*)file[i]->Get("tree");
  }
    Int_t nhMppc[N];
    Double_t x[N][1000];
    Double_t y[N][1000];
    Double_t z[N][1000];
    Int_t copynum[N][1000];

    Double_t evtposx[N];
    Double_t evtposy[N];




  for(int i=0;i<N;i++){
    data[i]->SetBranchAddress("nhMppc",&nhMppc[i]);
    data[i]->SetBranchAddress("evtposx",&evtposx[i]);
    data[i]->SetBranchAddress("evtposy",&evtposy[i]);
    data[i]->SetBranchAddress("mppcposx",x[i]);
    data[i]->SetBranchAddress("mppcposy",y[i]);
    data[i]->SetBranchAddress("mppcposz",z[i]);
    data[i]->SetBranchAddress("mppcnum",copynum[i]);
  }


  //Condition
  
  Int_t numx = 1;
  Int_t numy = 1;
  Int_t numz = 5;


  Double_t onex = 24/numx;
  Double_t oney = 24/numy;
  Int_t thre = 10;

  Int_t mul_th = 1;

  Int_t num_check=10;
  

  //Parameters
  Int_t multi[N][numx][numy][numz];
  Int_t cell_num[N][numx][numy][numz];
  Int_t numpercell[N];
  
  Int_t multi_effi[N][num_check];
  Int_t effi_thre[N][num_check];
  for(int i=0;i<num_check;i++){
    for(int j=0;j<N;j++){
      multi_effi[j][i]=0;
      effi_thre[j][i]=0;
    }
  }
  
  Double_t total = data[0]->GetEntries();

  
  

  
  
  Double_t effi_num[5];

  
  
  

  //TH1D* hist_multi = new TH1D("hist_multi","hist_multi",32,0,32);
  TH1D* hist_numpercell[N];
  TGraph *effi[N];
  TH2D* evt[N][num_check];
  for(int i=0;i<N;i++){
    hist_numpercell[i] = new TH1D(Form("numpercell%d",i),Form("numpercell%d",i),50,0,50);
    effi[i]= new TGraph();
    for(int j=0;j<num_check;j++){
      evt[i][j] = new TH2D(Form("evt%d%d",i,j),Form("evt%d%d",i,j),200,-65,65,200,-65,65);
      evt[i][j]->SetTitle(Form("Threshold of # of photons : %d;x [mm];y [mm]",j+1));
    }
  }

  //TH2D* evt = new TH2D("evt","evt",200,-65,65,200,-65,65);

  
  

  //Int_t multi_thre = 0;

  for(int v=0;v<N;v++){
    for(int n=0;n<total;n++){
      data[v]->GetEntry(n);
      for(int nx=0;nx<numx;nx++){
	for(int ny=0;ny<numy;ny++){
	  for(int nz=0;nz<numz;nz++){
	    multi[v][nx][ny][nz]=0;
	    cell_num[v][nx][ny][nz]=0;
	  }
	}
      }
      
      numpercell[v]=0;
      for(int i=0;i<num_check;i++){
	multi_effi[v][i]=0;
      }

      for(int i=0;i<nhMppc[v];i++){
	for(int nx=0;nx<numx;nx++){
	  for(int ny=0;ny<numy;ny++){
	    for(int nz=0;nz<numz;nz++){
	  
	      if(-12+onex*nx<=x[v][i]&&-12+onex*(nx+1)>x[v][i]&&-12+oney*ny<=y[v][i]&&-12+oney*(ny+1)>y[v][i]&&nz+1<=copynum[v][i]&&nz+2>copynum[v][i]){
		multi[v][nx][ny][nz]=1;
		cell_num[v][nx][ny][nz]++;
	      }
	    }
	  }
	}

      }
        


      for(int nx=0;nx<numx;nx++){
	for(int ny=0;ny<numy;ny++){
	  for(int nz=0;nz<numz;nz++){
	      
	    numpercell[v]+=multi[v][nx][ny][nz];
	    if(cell_num[v][nx][ny][nz]>0)hist_numpercell[v]->Fill(cell_num[v][nx][ny][nz]);
	    for(int l=1;l<num_check+1;l++)if(cell_num[v][nx][ny][nz]>l)multi_effi[v][l-1]++;
	  }
	}
      }
      for(int f=0;f<num_check;f++){
	if(multi_effi[v][f]>mul_th)effi_thre[v][f]++;
	if(multi_effi[v][f]<=mul_th)evt[v][f]->Fill(evtposx[v],evtposy[v]);
      }


    }
  }
    

  
  TGraph *graph_multi[N];
  for(int i=0;i<N;i++){
    graph_multi[i]= new TGraph();
    graph_multi[i]->SetMarkerStyle(4);
    graph_multi[i]->SetMarkerColor(1+i);
    graph_multi[i]->SetTitle(Form("Efficiency version %d;threshold of # of photons;efficiency",i+1));
    hist_numpercell[i]->SetTitle("# of photons at each cell;# of photons;n");
  }
  for(int i=0;i<num_check;i++){
    for(int j=0;j<N;j++){
      graph_multi[j]->AddPoint(1+i,effi_thre[j][i]/total);
    }
  }
  TMultiGraph *mg_multi = new TMultiGraph();
  for(int i=0;i<N;i++)mg_multi->Add(graph_multi[i]);
  mg_multi->SetTitle("Efficiency;threshold of # of photons;n");
  
  

  
  //hist_multi->SetTitle("Whole range;multiplicity;n");
  //hist_multi_tight->SetTitle("8 cm * 8 cm;multiplicity;n");

  hist_numpercell[0]->SetLineColor(kBlack);
  hist_numpercell[1]->SetLineColor(kBlue);
  hist_numpercell[2]->SetLineColor(kRed);
  
  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<N;i++){
    le->AddEntry(hist_numpercell[i],Form("version %d",i+1));
  }

  TLegend *le1 = new TLegend(0.8,0.5,0.48,0.6);
  for(int i=0;i<N;i++){
    le1->AddEntry(graph_multi[i],Form("version %d",i+1));
  }
  
  
  TCanvas *c1 = new TCanvas("c1","c1",800,650);
  c1->cd();
  
  hist_numpercell[0]->Draw();
  hist_numpercell[1]->Draw("same");
  hist_numpercell[2]->Draw("same");
  le->Draw();
  
  
  TCanvas *c2 = new TCanvas("c2","c2",800,650);
  c2->cd();
  mg_multi->Draw("AP");
  le1->Draw();

  TCanvas *c3[N];
  for(int i=0;i<N;i++){
    c3[i] = new TCanvas(Form("c3%d",i),Form("version %d",i+1),800,650);
    c3[i]->Divide(num_check/2,2);
  }
  
  for(int i=0;i<num_check/2;i++){
    for(int j=0;j<N;j++){
      c3[j]->cd(i+1);
      evt[j][i]->Draw("colz");
      c3[j]->cd(num_check/2+i+1);
      evt[j][num_check/2+i]->Draw("colz");
    }
  }


 
  
}
