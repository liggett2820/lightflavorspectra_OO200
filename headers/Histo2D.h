#ifndef _HISTO_2D_
#define _HISTO_2D_

// Ported verbatim from lightflavorspectra_etof/headers/Histo2D.h -- a standalone
// variable-Y-binning histogram utility class with no dependence on particle
// species/energy/eventConfig, so nothing needed to change for the 200 GeV O+O pi/K/p
// restriction.

#include "TMemberInspector.h"
#include "TNamed.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TFile.h"
#include <vector>
#include <string>
#include <iostream>
#include <ctime>

#include "../headers/HistogramUtilities.h"
#include "../headers/PhysMath.h"

using namespace std;

// Histogram Class For Variable Bin Widths in Y Direction

class Histo2D : public TNamed{
public:
  // if you don't supply y projections here, you have to do it later!
  Histo2D(string a_name, string a_title, unsigned int a_num_x_bins, const double* a_x_bin_edges, vector<TH1D*> a_y_projections = {});
  Histo2D(string a_name, string a_title, unsigned int a_num_x_bins, double a_x_min, double a_x_max, vector<TH1D*> a_y_projections = {});
  Histo2D(TH2D* a_histo); // imports a TH2D to manipulate verticle slices independantly
  Histo2D(); // default constructor for ROOT I/O
  Histo2D(const Histo2D& other); // copy constructor
  Histo2D& operator=(const Histo2D& other); // assignment operator
  ~Histo2D();

  void SetYProjection(int a_index_x, TH1D* a_y_projection, bool a_copy = true); // this will use the same pointer (does not copy!)
  void SetYProjections(vector<TH1D*> a_y_projections, bool a_copy = true); // these will use the same pointers (does not copy!)
  void SetEmptyEvenYProjections(unsigned int a_num_y_bins, double a_y_min, double a_y_max);

  TH1D* GetInternalYProjection(int a_index_x);

  void Fill(double a_x, double a_y, double a_weight = 1.0);

  void Draw(string a_option = "", bool a_low_content_range_is_zero = true, string a_xaxis_title = "", string a_yaxis_title = "");
  void Draw(Option_t *option=""){Draw("colz",false);};
  void SetXDrawRange(double a_min, double a_max);
  void SetYDrawRange(double a_min, double a_max);
  void SetZDrawRangeForLogScale();

  void PrintState() const; // Function to print the state of the class


  //Operations to Change 1D projections
  //void normalizeEachVerticalSlice(bool a_isDifferential = false);



private:

  bool m_weak_construction; // doesn't contain y projections yet
  //string m_name;
  //string m_title;
  unsigned int    m_num_bins_x;
  vector<double>  m_bin_edges_x;
  vector<TH1D*>   m_vertical_histogram_projections;

  //vector<unsigned int>   m_num_bins_y;
  //vector<vector<double>> m_bin_edges_y;
  //vector<vector<double>> m_bin_content;
  //vector<vector<double>> m_bin_error_high;
  //vector<vector<double>> m_bin_error_low;
  TH2D* m_histogram_drawing_base; // does not have correct verticle bin structure (just for drawing)
  vector<TH2D*> m_vertical_2D_drawing_projections;

  void make_vertical_drawing_projections();
  void copy(const Histo2D& other); // helper function for copying

  ClassDef(Histo2D, 1); // ROOT class definition macro

};




// Test Functions
void test_Histo2D();


#endif
