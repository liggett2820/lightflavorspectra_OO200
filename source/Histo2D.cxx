#ifndef _HISTO_2D_CXX_
#define _HISTO_2D_CXX_

// Ported verbatim from lightflavorspectra_etof/source/Histo2D.cxx -- standalone
// variable-Y-binning histogram utility, no physics-config dependence, so no changes
// needed.

#include "../headers/Histo2D.h"
#include "TCanvas.h"

ClassImp(Histo2D);

Histo2D::Histo2D() : TNamed("name", "title"){
  m_num_bins_x = 0;
  m_bin_edges_x.clear();
  m_weak_construction = true;
  m_histogram_drawing_base = nullptr;

  vector<TH1D*> temp = {};
  m_vertical_histogram_projections = temp;
  vector<TH2D*> temp2D = {};
  m_vertical_2D_drawing_projections = temp2D;
}

//________________________________________________________________________________________________________________________________________________


Histo2D::Histo2D(string a_name, string a_title, unsigned int a_num_x_bins, const double* a_x_bin_edges, vector<TH1D*> a_y_projections)
   : TNamed(a_name.c_str(), a_title.c_str()){

  m_num_bins_x = a_num_x_bins;
  m_bin_edges_x.clear();
  m_bin_edges_x.push_back(a_x_bin_edges[0]);
  for(unsigned int index = 1; index <= a_num_x_bins; index++){
    m_bin_edges_x.push_back(a_x_bin_edges[index]);
    if(a_x_bin_edges[index] <= a_x_bin_edges[index-1]){
      cout << "ERROR: Histo2D::Histo2D  You supplied bin edges that have no width or are not in order! " << a_x_bin_edges[index-1] << "  ->  " << a_x_bin_edges[index] << endl;
    }
  }
  if(a_y_projections.size() > 0){
    m_vertical_histogram_projections = a_y_projections;
    m_weak_construction = false;
  }else{
    m_weak_construction = true;
    m_vertical_histogram_projections.resize(m_num_bins_x, nullptr);
  }
  m_histogram_drawing_base = nullptr;
}

//________________________________________________________________________________________________________________________________________________

Histo2D::Histo2D(string a_name, string a_title, unsigned int a_num_x_bins, double a_x_min, double a_x_max, vector<TH1D*> a_y_projections)
  : TNamed(a_name.c_str(), a_title.c_str()){
  vector<double> bin_edges_x = HistogramUtilities::makeEvenBinEdges(a_num_x_bins,a_x_min,a_x_max);
  //Histo2D(a_name,a_title,a_num_x_bins, &bin_edges_x[0],a_y_projections);
  SetName(a_name.c_str());
  SetTitle(a_title.c_str());
  m_num_bins_x = a_num_x_bins;
  m_bin_edges_x = bin_edges_x;
  if (a_y_projections.size() > 0) {
    m_vertical_histogram_projections = a_y_projections;
    m_weak_construction = false;
  } else {
    m_weak_construction = true;
    m_vertical_histogram_projections.resize(m_num_bins_x, nullptr);
  }
  m_histogram_drawing_base = nullptr;
}

//________________________________________________________________________________________________________________________________________________


Histo2D::Histo2D(TH2D* a_histo): TNamed(a_histo->GetName(), a_histo->GetTitle()){
  if(!a_histo){
    cout << "ERROR: Histo2D::Histo2D(TH2D* a_histo) : Handed Null Pointer!" << endl;
    return;
  }
  m_weak_construction = false;
  SetName(a_histo->GetName());
  SetTitle(a_histo->GetTitle());
  m_num_bins_x = a_histo->GetNbinsX();
  m_bin_edges_x = HistogramUtilities::getBinEdges_X(a_histo);
  m_vertical_histogram_projections.clear();
  for(int index_x = 0; index_x < m_num_bins_x; index_x++){
    m_vertical_histogram_projections.push_back(a_histo->ProjectionY(Form("%s_xIndex%d",GetName(),index_x),index_x+1,index_x+1));
  }

  m_histogram_drawing_base = nullptr;
}

//________________________________________________________________________________________________________________________________________________

Histo2D::Histo2D(const Histo2D& other): TNamed(other) {
    copy(other);
}

//________________________________________________________________________________________________________________________________________________


Histo2D& Histo2D::operator=(const Histo2D& other) {
    if (this != &other) {
        if (m_histogram_drawing_base) delete m_histogram_drawing_base;
        for (auto& proj : m_vertical_histogram_projections) {
            delete proj;
        }
        for (auto& proj : m_vertical_2D_drawing_projections) {
            delete proj;
        }
        copy(other);
    }
    return *this;
}

//________________________________________________________________________________________________________________________________________________


Histo2D::~Histo2D(){
    if (m_histogram_drawing_base) delete m_histogram_drawing_base;
    for (auto& proj : m_vertical_histogram_projections) {
        delete proj;
    }
    for (auto& proj : m_vertical_2D_drawing_projections) {
        delete proj;
    }
}

//________________________________________________________________________________________________________________________________________________


void Histo2D::copy(const Histo2D& other) {
    m_weak_construction = other.m_weak_construction;
    SetName(other.GetName());
    SetTitle(other.GetTitle());
    //m_title = other.m_title;
    m_num_bins_x = other.m_num_bins_x;
    m_bin_edges_x = other.m_bin_edges_x;
    m_histogram_drawing_base = other.m_histogram_drawing_base ? new TH2D(*other.m_histogram_drawing_base) : nullptr;
    m_vertical_histogram_projections.clear();
    for (const auto& proj : other.m_vertical_histogram_projections) {
        m_vertical_histogram_projections.push_back(new TH1D(*proj));
    }
    m_vertical_2D_drawing_projections.clear();
    for (const auto& proj : other.m_vertical_2D_drawing_projections) {
        m_vertical_2D_drawing_projections.push_back(new TH2D(*proj));
    }
}

//________________________________________________________________________________________________________________________________________________


void Histo2D::PrintState() const {
    cout << "Histo2D State:" << endl;
    //cout << "Name: " << GetName() << endl;
    //cout << "Title: " << GetTitle() << endl;
    cout << "Number of X Bins: " << m_num_bins_x << endl;
    cout << "X Bin Edges: ";
    for (const auto& edge : m_bin_edges_x) {
        cout << edge << " ";
    }
    cout << endl;
    cout << "Weak Construction: " << (m_weak_construction ? "true" : "false") << endl;
    cout << "Number of Vertical Histogram Projections: " << m_vertical_histogram_projections.size() << endl;
    cout << "Histogram Drawing Base: " << (m_histogram_drawing_base ? "set" : "not set") << endl;
    cout << "Number of Vertical 2D Drawing Projections: " << m_vertical_2D_drawing_projections.size() << endl;
}

//________________________________________________________________________________________________________________________________________________


void Histo2D::SetYProjection(int a_index_x, TH1D* a_y_projection, bool a_copy){
  if(a_index_x < 0 || a_index_x >= m_num_bins_x){
    cout << "ERROR: Histo2D::SetYProjection(int a_index_x, TH1D* a_y_projection) : Index out of range" << endl;
    return;
  }
  if(a_index_x >= m_vertical_histogram_projections.size()){
    for(int index = m_vertical_histogram_projections.size(); index <= a_index_x; index++){
      m_vertical_histogram_projections.push_back(nullptr);
    }
  }
  // only delete if internal data are copies
  if(m_vertical_histogram_projections[a_index_x] && a_copy){
    delete m_vertical_histogram_projections[a_index_x];
  }
  if(a_copy){
    m_vertical_histogram_projections[a_index_x] = (TH1D*) a_y_projection->Clone();
  }else{
    m_vertical_histogram_projections[a_index_x] = a_y_projection;
  }
}



//________________________________________________________________________________________________________________________________________________

void Histo2D::SetYProjections(vector<TH1D*> a_y_projections, bool a_copy){
  if(a_y_projections.size() != m_num_bins_x){
    cout << "ERROR: Histo2D::SetYProjections(vector<TH1D*> a_y_projections) you supplied " << a_y_projections.size() << " histograms, but the internal structure required m_num_bins_x" << endl;  }else{
    m_weak_construction = false;
  }
  m_vertical_histogram_projections = a_y_projections;
}

//________________________________________________________________________________________________________________________________________________


void Histo2D::SetEmptyEvenYProjections(unsigned int a_num_y_bins, double a_y_min, double a_y_max){
  m_vertical_histogram_projections.clear();
  vector<double> bin_edges = HistogramUtilities::makeEvenBinEdges(a_num_y_bins,a_y_min,a_y_max);
  for(int index_x = 0; index_x < m_num_bins_x; index_x++){
    TH1D* base = new TH1D(Form("%s_xIndex%d",GetName(),index_x),Form("%s xIndex%d",GetTitle(),index_x),a_num_y_bins,&bin_edges[0]);
    m_vertical_histogram_projections.push_back(base);
  }
  m_weak_construction = false;
}

//________________________________________________________________________________________________________________________________________________


TH1D* Histo2D::GetInternalYProjection(int a_index_x){
  if(m_weak_construction){
    cout << "ERROR: set the y projections first!!!!" << endl;
    return NULL;
  }
  return m_vertical_histogram_projections[a_index_x];
}


//________________________________________________________________________________________________________________________________________________


void Histo2D::Draw(string a_option, bool a_low_content_range_is_zero, string a_xaxis_title, string a_yaxis_title){
  //if(m_vertical_2D_drawing_projections.size() == 0){
    make_vertical_drawing_projections();
  //}


  double min_content = 1E200;
  double max_content = -1E200;
  for(int index_x = 0; index_x < m_num_bins_x; index_x++){
    if(!m_vertical_2D_drawing_projections[index_x]) continue;
    for(int bin = 1; bin <= m_vertical_histogram_projections[index_x]->GetNbinsX(); bin++){
      double content = m_vertical_histogram_projections[index_x]->GetBinContent(bin);
      if(content < min_content) min_content = content;
      if(content > max_content) max_content = content;
    }
  }


  if(a_low_content_range_is_zero){
    m_histogram_drawing_base->SetMinimum(0);
    m_histogram_drawing_base->SetMaximum(1.05*max_content);

    for(int index_x = 0; index_x < m_num_bins_x; index_x++){
      if(!m_vertical_2D_drawing_projections[index_x]) continue;
      m_vertical_2D_drawing_projections[index_x]->SetMinimum(0);
      m_vertical_2D_drawing_projections[index_x]->SetMaximum(1.05*max_content);
    }
  }else{
    m_histogram_drawing_base->SetMinimum(min_content - 0.05*(max_content - min_content));
    m_histogram_drawing_base->SetMaximum(max_content + 0.05*(max_content - min_content));
    for(int index_x = 0; index_x < m_num_bins_x; index_x++){
      if(!m_vertical_2D_drawing_projections[index_x]) continue;
      m_vertical_2D_drawing_projections[index_x]->SetMinimum(min_content - 0.05*(max_content - min_content));
      m_vertical_2D_drawing_projections[index_x]->SetMaximum(max_content + 0.05*(max_content - min_content));
    }
  }
  m_histogram_drawing_base->SetStats(kFALSE);
  if(a_xaxis_title != "") m_histogram_drawing_base->GetXaxis()->SetTitle(a_xaxis_title.c_str());
  if(a_yaxis_title != "") m_histogram_drawing_base->GetYaxis()->SetTitle(a_yaxis_title.c_str());
  m_histogram_drawing_base->Draw(a_option.c_str());
  for(int index_x = 0; index_x < m_num_bins_x; index_x++){
    if(!m_vertical_2D_drawing_projections[index_x]) continue;
    m_vertical_2D_drawing_projections[index_x]->SetStats(kFALSE);
    m_vertical_2D_drawing_projections[index_x]->Draw(Form("%s SAME",a_option.c_str()));
  }



}

//________________________________________________________________________________________________________________________________________________


void Histo2D::SetXDrawRange(double a_min, double a_max){
  if(m_histogram_drawing_base) m_histogram_drawing_base->GetXaxis()->SetRangeUser(a_min,a_max);
}

//________________________________________________________________________________________________________________________________________________

void Histo2D::SetYDrawRange(double a_min, double a_max){
  if(m_histogram_drawing_base) m_histogram_drawing_base->GetYaxis()->SetRangeUser(a_min,a_max);
}

//________________________________________________________________________________________________________________________________________________

void Histo2D::SetZDrawRangeForLogScale(){
  if(!m_histogram_drawing_base) return;
  double min_nonzero_value = 1E200;
  double max_nonzero_value = 0.0;
  for(int index_x = 0; index_x < m_vertical_histogram_projections.size(); index_x++){
    if(m_vertical_histogram_projections[index_x]){
      for(int bin = 1; bin <= m_vertical_histogram_projections[index_x]->GetNbinsX(); bin++){
        double content = m_vertical_histogram_projections[index_x]->GetBinContent(bin);
        if(content < min_nonzero_value && content > 0) min_nonzero_value = content;
        if(content > max_nonzero_value && content > 0) max_nonzero_value = content;
      }
    }
  }
  m_histogram_drawing_base->SetMinimum(0.8*min_nonzero_value);
  m_histogram_drawing_base->SetMaximum(1.2*max_nonzero_value);
  for(int index_x = 0; index_x < m_vertical_2D_drawing_projections.size(); index_x++){
    if(m_vertical_2D_drawing_projections[index_x]){
      m_vertical_2D_drawing_projections[index_x]->SetMinimum(0.8*min_nonzero_value);
      m_vertical_2D_drawing_projections[index_x]->SetMaximum(1.2*max_nonzero_value);
    }
  }

}



//________________________________________________________________________________________________________________________________________________


void Histo2D::Fill(double a_x, double a_y, double a_weight){
  if(a_x < m_bin_edges_x[0]) return;
  int x_index = 0;
  while(a_x > m_bin_edges_x[x_index] && x_index < m_num_bins_x) x_index++;
  x_index--;
  if(x_index == m_num_bins_x) return;
  //cout << "x index found " << x_index << "  with value: " << a_y << "  and weight " << a_weight << endl;
  if(!m_vertical_histogram_projections[x_index]){
    cout << "ERROR: Histo2D::Fill()   no y projection found for x index " << x_index << endl;
    return;
  }
  m_vertical_histogram_projections[x_index]->Fill(a_y);//,a_weight);
  int bin = m_vertical_histogram_projections[x_index]->FindBin(a_y);
  //cout << "   at hor bin " << bin << " now with content " << m_vertical_histogram_projections[x_index]->GetBinContent(bin) << endl;
}


//________________________________________________________________________________________________________________________________________________


void Histo2D::make_vertical_drawing_projections(){
  //delete the old ones
  if(m_histogram_drawing_base) delete m_histogram_drawing_base;
  for(int index_x = 0; index_x < m_vertical_2D_drawing_projections.size(); index_x++){
    if(m_vertical_2D_drawing_projections[index_x]) delete m_vertical_2D_drawing_projections[index_x];
  }
  m_vertical_2D_drawing_projections.clear();


  //make the new ones
  double min_vert = 1E200;
  double max_vert = -1E200;
  for(int index_x = 0; index_x < m_num_bins_x; index_x++){
    //cout << "PTR: " << m_vertical_histogram_projections[index_x] << endl;
    if(!m_vertical_histogram_projections[index_x]){
      m_vertical_2D_drawing_projections.push_back(nullptr);
      continue;
    }
    vector<double> bin_edges_vert = HistogramUtilities::getBinEdges(m_vertical_histogram_projections[index_x]);
    if(bin_edges_vert[0] < min_vert)                       min_vert = bin_edges_vert[0];
    if(bin_edges_vert[bin_edges_vert.size()-1] > max_vert) max_vert = bin_edges_vert[bin_edges_vert.size()-1];
    double bin_edges_hor[2] = {m_bin_edges_x[index_x],m_bin_edges_x[index_x+1]};

    //cout << "j   from " << bin_edges_hor[0] << "   to   " << bin_edges_hor[1] << endl;
    TH2D* proj = new TH2D(Form("%s_xIndex%d_Draw",GetName(),index_x),"",1,&bin_edges_hor[0],((int) bin_edges_vert.size())-1,&bin_edges_vert[0]);
    for(int bin = 1; bin <= proj->GetNbinsY(); bin++){
      proj->SetBinContent(1,bin,m_vertical_histogram_projections[index_x]->GetBinContent(bin));
      proj->SetBinError(1,bin,m_vertical_histogram_projections[index_x]->GetBinError(bin));

    }
    proj->SetContour(99);
    m_vertical_2D_drawing_projections.push_back(proj);
  }
  //cout << " vert = " << min_vert << "  to  " << max_vert << endl;
  vector<double> bin_edges_drawing = HistogramUtilities::makeEvenBinEdges(30,min_vert,max_vert);
  m_histogram_drawing_base = new TH2D(Form("%s_drawing_base",GetName()),GetTitle(),m_num_bins_x,&m_bin_edges_x[0],30,&bin_edges_drawing[0]);

}

//________________________________________________________________________________________________________________________________________________

/*
void Histo2D::normalizeEachVerticalSlice(bool a_isDifferential){
  for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
    double normalizationFactor = 1.0;
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      if(a_isDifferential) normalizationFactor += a_histo->GetYaxis()->GetBinWidth(binY)*a_histo->GetBinContent(binX,binY);
      else                 normalizationFactor += a_histo->GetBinContent(binX,binY);
    }
    if(normalizationFactor == 0.0) continue;
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      a_histo->SetBinContent(binX,binY, a_histo->GetBinContent(binX,binY)/normalizationFactor);
      a_histo->SetBinError(binX,binY, a_histo->GetBinError(binX,binY)/normalizationFactor);
    }
  }
}*/






//________________________________________________________________________________________________________________________________________________




void test_Histo2D(){

  double x_bin_edges[4] = {-3,-2.5,1,4};
  vector<double> y_bin_edges_0 = {0,1,2,3,4,5,6,7,11,12,13,14,18,19,20};
  vector<double> y_bin_edges_1 = {0,1,4,6,7,12,13,14,15,16,17,18,19,20};
  vector<double> y_bin_edges_2 = {0,2,3,4,7,8,9,10,16,17,18,20};

  TH1D* hist_0 = new TH1D("","",((int)y_bin_edges_0.size()) - 1, &y_bin_edges_0[0]);
  TH1D* hist_1 = new TH1D("","",((int)y_bin_edges_1.size()) - 1, &y_bin_edges_1[0]);
  TH1D* hist_2 = new TH1D("","",((int)y_bin_edges_2.size()) - 1, &y_bin_edges_2[0]);

  vector<TH1D*> histos = {hist_0,hist_1,hist_2};
  cout << "Constructing" << endl;
  Histo2D* hist = new Histo2D("name", "title", 3, &x_bin_edges[0], histos);

  TRandom3* random_obj = new TRandom3(time(NULL));
  cout << "Filling" << endl;
  // loop for embedding sample (not the data distribution)
  for(int trial_index = 0; trial_index < 100000; trial_index++){
    double val = random_obj->Gaus(12,6);
    hist->Fill(-2.8,val,1.0);
    hist->Fill(0,val,5.0);
    hist->Fill(2,val,2.0);
  }

  TCanvas* canv = new TCanvas();
  cout << "Drawing" << endl;
  hist->Draw("colz");

  hist->PrintState();

  hist->SetName("name");
  hist->SetTitle("title");

  // Test writing the Histo2D to a file, then reading it back
  // Test writing the Histo2D to a file, then reading it back
  TFile* file = new TFile("test_Histo2D.root", "RECREATE");
  hist->Write();
  file->Close();

  TFile* file2 = new TFile("test_Histo2D.root", "READ");
  if (!file2 || file2->IsZombie()) {
      cout << "Error opening file" << endl;
      return;
  }

  file2->ls(); // List the contents of the file to check if the object is there


  Histo2D* hist2 = (Histo2D*)file2->Get("name");
  if (!hist2) {
      cout << "Error reading Histo2D from file" << endl;
      file2->Close();
      return;
  }
  hist2->PrintState();

  TCanvas* canv2 = new TCanvas();
  hist2->Draw("colz");
  //file2->Close();





}





#endif
