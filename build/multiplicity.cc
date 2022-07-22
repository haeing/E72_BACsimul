void multiplicity(){

  TFile *file = new TFile("reflect_140_lg_50_mid_20_beam_z2mm.root","read");
  TTree *data = (TTree*)file->Get("data");
  Int_t num;
  Double_t posx[50];
  Double_t posy[50];
  Double_t posz[50];

  data->SetBranchAddress("nhMppc",&num);
  data->SetBranchAddress("mppcposx",posx);
  data->SetBranchAddress("mppcposy",posy);
  data->SetBranchAddress("mppcposz",posz);

  Double_t total = data->GetEntries();
  Int_t multi = 0;

  TH1D* hist = new TH1D("hist","hist",50,0,50);
  for(int n=0;n<total;n++){
    data->GetEntry(n);
    for(int i=0;i<num;i++){
      for(int j=i
    }
    
      
    
}
