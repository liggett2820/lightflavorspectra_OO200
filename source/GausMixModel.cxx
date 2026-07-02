#ifndef _GausMixModel_cxx_
#define _GausMixModel_cxx_

// Ported verbatim from lightflavorspectra_etof/source/GausMixModel.cxx -- standalone
// numerical utility, no physics-config dependence, so no changes needed.

#include "../headers/GausMixModel.h"


//_____________________________________________________________________________
GausMixModel::GausMixModel(){};


void GausMixModel::importDataAsHistogram(TH1D* a_histo){m_dataHistogram =(TH1D*) a_histo->Clone();};
void GausMixModel::importDataAs2DHistogram(TH2D* a_histo){m_dataHistogram2D = (TH2D*) a_histo->Clone();};
//void importDataAsArray(int a_numData, double* a_data);
void GausMixModel::setNumberOfGaussians(int a_numGaus){m_numGaus = a_numGaus;};

//_____________________________________________________________________________
//void GausMixModel::initializeParameters(int a_numParams, double* a_parameters){
//  m_parameters.clear();
//  for(int iii = 0; iii < a_numParams; iii++){
//    m_parameters.push_back(a_parameters[iii]);
//  }
//}

//_____________________________________________________________________________
bool GausMixModel::doKMeans(double* a_meansToFill){
  //Note: This putting the means at the edges of the histogram sometimes. If this happens, GMM
  //will just evenly divide the means on the range of the histogram
  srand(time(NULL));
  int nBins = m_dataHistogram->GetNbinsX();
  double lowHistoEdge = m_dataHistogram->GetBinLowEdge(1);
  double highHistoEdge = m_dataHistogram->GetBinLowEdge(nBins+1);

  int lowBin = 1;
  int highBin = nBins;
  //int lowBin = HistogramUtilities::firstBinWithContent(m_dataHistogram);
  //int highBin = HistogramUtilities::lastBinWithContent(m_dataHistogram);
  //double lowHistoEdge = m_dataHistogram->GetBinLowEdge(lowBin);
  //double highHistoEdge = m_dataHistogram->GetBinLowEdge(highBin+1);

  vector<double> means;
  //start the means off equally spaced
  for(int iii = 0; iii < m_numGaus; iii++){
    double step = (highHistoEdge - lowHistoEdge)/((double) (m_numGaus + 1));
    means.push_back(lowHistoEdge + step * ((double) ( iii + 1 )) );
  }

  //double meanTolerance = m_dataHistogram->GetBinWidth(1)/20.0;
  double smallestBinWidth = m_dataHistogram->GetBinWidth(1);
  for(int bin = 2; bin <= m_dataHistogram->GetNbinsX(); bin++){
    if(m_dataHistogram->GetBinWidth(bin) < smallestBinWidth) smallestBinWidth = m_dataHistogram->GetBinWidth(bin);
  }
  double meanTolerance = smallestBinWidth/20.0;


  vector<double> prevMeans(m_numGaus);
  vector< vector<double> > associationMatrix;
  for(int iii = 0; iii < nBins; iii++){
    vector<double> tempVect1(m_numGaus);
    associationMatrix.push_back(tempVect1);
  }

  for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    prevMeans[gausIndex] = lowHistoEdge - (highHistoEdge - lowHistoEdge);
  }

  bool meansOutOfTolerance = true;
  int trialMax = 200;
  int trialIndex = 0;
  while(meansOutOfTolerance){
    if(trialIndex > trialMax) break;

    //Calculate Assosiation Matrix
    for(int bin = lowBin; bin <= highBin; bin++){
      double smallestDistance = 1e200;
      int closestMean = -1;
      for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
        double distance = fabs(m_dataHistogram->GetBinCenter(bin)  //---------------------------------------- should this be bin + 1?
                               - means[gausIndex]);
        if(distance < smallestDistance){
          smallestDistance = distance;
          closestMean = gausIndex;
        }
        associationMatrix[bin-1][gausIndex] = 0.0; //for initialization
      }
      associationMatrix[bin-1][closestMean] = 1.0;
    }

    //See If Means are in all precise enough and copy mean array
    bool allInTolerance = true;
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      if(fabs(prevMeans[gausIndex] - means[gausIndex]) > meanTolerance){
        allInTolerance = false;
      }
      prevMeans[gausIndex] = means[gausIndex];

    }
    meansOutOfTolerance = !allInTolerance;



    //calculate new means
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      double newMeanNumerator = 0;
      double newMeanDenominator = 0;
      for(int bin = lowBin; bin <= highBin; bin++){
        newMeanNumerator += associationMatrix[bin-1][gausIndex]
                            * m_dataHistogram->GetBinCenter(bin)
                            * m_dataHistogram->GetBinContent(bin);
        newMeanDenominator += associationMatrix[bin-1][gausIndex]
                            * m_dataHistogram->GetBinContent(bin);
      }
      if(newMeanDenominator == 0){
        means[gausIndex] = m_dataHistogram->GetBinCenter(rand() % nBins + 1);
      }else{
        means[gausIndex] = newMeanNumerator/newMeanDenominator;
      }
    }
    trialIndex++;
  }

  //sort the means
  std::sort(means.begin(),means.end());
  bool anyInPeripheralBins = false;
  if(nBins > 7){
    double lowPeripheral = m_dataHistogram->GetBinCenter(3);
    double highPeripheral = m_dataHistogram->GetBinCenter(m_dataHistogram->GetNbinsX()-4);
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      a_meansToFill[gausIndex] = means[gausIndex];
      if(means[gausIndex] > highPeripheral || means[gausIndex] < lowPeripheral){
        anyInPeripheralBins = true;
      }
    }
  }
  //see if any hit the edge of histogram
  // will fail if in last 3 bins
  return !anyInPeripheralBins;
}

//_____________________________________________________________________________
//_____________________________________________________________________________


//######   KMEANS   2D    #########

bool GausMixModel::doKMeans2D(double* a_meansToFill){
  cout << "doKMeans2D does not function yet" << endl;
  /*
  srand(time(NULL));

  int nXBins = m_dataHistogram2D->GetNbinsX();
  int nYBins = m_dataHistogram2D->GetNbinsY();

  double lowHistoEdgeX = m_dataHistogram2D->GetXaxis()->GetBinLowEdge(1);
  double highHistoEdgeX = m_dataHistogram2D->GetXaxis()->GetBinLowEdge(nXBins);
  double histogramRangeX = highHistoEdgeX - lowHistoEdgeX;
  double lowHistoEdgeY = m_dataHistogram2D->GetYaxis()->GetBinLowEdge(1);
  double highHistoEdgeY = m_dataHistogram2D->GetYaxis()->GetBinLowEdge(nYBins);
  double histogramRangeY = highHistoEdgeY - lowHistoEdgeY;



  vector<double> means;
  //start the means off equally spaced on diagonal
  //from top right to bottom left
  std::cout << "Initial Means: ";
  double stepY = histogramRangeY/((double) (m_numGaus + 1));
  double stepX = histogramRangeX/((double) (m_numGaus + 1));
  for(int iii = 0; iii < m_numGaus; iii++){
    double meanX = lowHistoEdgeX + histogramRangeX - (iii + 1)*stepX;
    double meanY = lowHistoEdgeY + (iii + 1)*stepY;
    means.push_back(meanX);
    means.push_back(meanY);
    //std::cout << means[iii] << "  ";
  }
  std::cout << std::endl;

  double meanTolerance = m_dataHistogram->GetBinWidth(1)/20.0;
  vector<double> prevMeans(2*m_numGaus);

  vector< vector< vector<int> > > associationMatrix;//[nXBins][nYBins][m_numGaus];
  for(int iii = 0; iii < nXBins; iii++){
    vector< vector<int> > tempVect3;
    for(int jjj = 0; jjj < nYBins; jjj++){
      vector<int> tempVect4;
      for(int kkk = 0; kkk < m_numGaus; kkk++) tempVect4.push_back(0);
      tempVect3.push_back(tempVect4);
    }
    associationMatrix.push_back(tempVect3);
  }


  cout << "Setup Association Matrix" << endl;

  for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    prevMeans[2*gausIndex] = lowHistoEdgeX - histogramRangeX;
    prevMeans[2*gausIndex+1] = lowHistoEdgeY - histogramRangeY;
  }

  cout << "Starting Trial Loop of KMeans" << endl;
  bool meansOutOfTolerance = true;
  int trialMax = 200;
  int trialIndex = 0;
  while(meansOutOfTolerance){
    cout << "Trial Index: " << trialIndex << endl;
    if(trialIndex > trialMax) break;

    cout << "calculating association matrix" << endl;
    //Calculate Assosiation Matrix
    for(int binX = 1; binX <= nXBins; binX++){
      for(int binY = 1; binX <= nYBins; binY++){
        cout << "  bin X Y= " << binX << " " << binY << endl;
        double smallestDistance = 1e200;
        int closestMean = -1;
        for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
          double distance =
                            (m_dataHistogram2D->GetXaxis()->GetBinCenter(binX) - means[2*gausIndex])
                            *(m_dataHistogram2D->GetXaxis()->GetBinCenter(binX) - means[2*gausIndex])
                            +(m_dataHistogram2D->GetYaxis()->GetBinCenter(binY) - means[2*gausIndex+1])
                            *(m_dataHistogram2D->GetYaxis()->GetBinCenter(binY) - means[2*gausIndex+1]);

          if(distance < smallestDistance){
            smallestDistance = distance;
            closestMean = gausIndex;
          }
          associationMatrix[binX][binY][gausIndex] = 0; //for initialization
        }
        associationMatrix[binX][binY][closestMean] = 1;

    }
    cout << "finished association matrix" << endl;



    //See If Means are in all precise enough and copy mean array
    bool allInTolerance = true;
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      double distance =
                            (prevMeans[2*gausIndex] - means[2*gausIndex])
                            *(prevMeans[2*gausIndex] - means[2*gausIndex])
                            +(prevMeans[2*gausIndex+1] - means[2*gausIndex+1])
                            *(prevMeans[2*gausIndex+1] - means[2*gausIndex+1]);
      if(TMath::Sqrt(distance) > meanTolerance){
        allInTolerance = false;
      }
      prevMeans[2*gausIndex] = means[2*gausIndex];
      prevMeans[2*gausIndex+1] = means[2*gausIndex+1];

    }
    meansOutOfTolerance = !allInTolerance;

    / *
    stopped here 22invin023nvnna 23 23 23 2 3 23 2 3 23 avn oan 2323
    //calculate new means
    std::cout << "Current Means: ";
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      double newMeanNumerator1 = 0;
      double newMeanDenominator1 = 0;

      for(int binX = 1; binX <= nXBins; binX++){
        for(int binY = 1; binX <= nYBins; binY++){

          newMeanNumerator += associationMatrix[binX][binY][gausIndex]

                            * m_dataHistogram->GetBinCenter(binIndex+1)
                            * m_dataHistogram->GetBinContent(binIndex+1);
          newMeanDenominator += associationMatrix[binX][binY][gausIndex]
                            * m_dataHistogram->GetBinContent(binIndex+1);
        }
      }
      if(newMeanDenominator == 0){
        //std::cout << "\nInitialization was Bad: No Data Near Mean No"
        //          << gausIndex << std::endl;
        means[gausIndex] = m_dataHistogram->GetBinCenter(rand() % nBins + 1);
      }else{
        means[gausIndex] = newMeanNumerator/newMeanDenominator;
      }
      std::cout << means[gausIndex] << "  ";

    }
    std::cout << std::endl;
    trialIndex++;
  }

  //sort the means
  std::sort(means.begin(),means.end());
  //std::cout << "Final K-Means: ";
  bool anyInPeripheralBins = false;
  double lowPeripheral = m_dataHistogram->GetBinCenter(3);
  double highPeripheral = m_dataHistogram->GetBinCenter(m_dataHistogram->GetNbinsX()-4);
  for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    //std::cout << means[gausIndex] << "  ";
    a_meansToFill[gausIndex] = means[gausIndex];
    if(means[gausIndex] > highPeripheral || means[gausIndex] < lowPeripheral){
      anyInPeripheralBins = true;
    }
  }
  //std::cout << std::endl;

  //see if any hit the edge of histogram
  // will fail if in last 3 bins
  return !anyInPeripheralBins;
  */

  return false;

}



//_____________________________________________________________________________

void GausMixModel::doGMM(double* a_parametersToFill,
                        double* a_initialMeans,
                        double* a_initialAmps,
                        int a_maxTrials,
                        double a_tolerance,
                        double a_avgSigma,
                        double a_minBinContent){
  vector<double> initialSigmas;
  for(int iii = 0; iii < m_numGaus; iii++){
    initialSigmas.push_back(a_avgSigma);
  }
  doGMM(a_parametersToFill,a_initialMeans,a_initialAmps,&initialSigmas[0],a_maxTrials,a_tolerance,a_minBinContent);
}






//_____________________________________________________________________________

void GausMixModel::doGMM(double* a_parametersToFill,
                          double* a_initialMeans,
                          double* a_initialAmps,
                          double* a_initialSigmas,
                          int a_maxTrials,
                          double a_tolerance,
                          double a_minBinContent){
  //bool showCanvas = false;
  if(a_minBinContent > 0.0){
    HistogramUtilities::zeroBinsBelowThreshold(m_dataHistogram, a_minBinContent);
    a_minBinContent = 0;
  }
  int nBins = m_dataHistogram->GetNbinsX();
  if(nBins <= 0){
   cout << "ERROR: GausMixModel::doGMM number of bins handed to gaus mixture model was <= 0" << endl;
   return;
  }
  //int lowBinIndex = HistogramUtilities::firstBinWithContent(m_dataHistogram);
  //int highBinIndex = HistogramUtilities::lastBinWithContent(m_dataHistogram);
  //double lowHistoEdge = m_dataHistogram->GetBinLowEdge(lowBinIndex);
  //double highHistoEdge = m_dataHistogram->GetBinLowEdge(highBinIndex+1);
  //nBins = highBinIndex - lowBinIndex;


  double lowHistoEdge   = m_dataHistogram->GetBinLowEdge(1);
  double highHistoEdge  = m_dataHistogram->GetBinLowEdge(nBins+1);
  double histogramRange = highHistoEdge - lowHistoEdge;

  //SETUP THE GRAPHS TO MONITOR PARAMETERS WHILE FITTING
  std::vector< std::vector<TGraph*> > theGraphs;
  TCanvas* paramCanvas = NULL;
  /*if(showCanvas){
    paramCanvas = new TCanvas("Parameter Canvas");
    paramCanvas->Divide(3,m_numGaus);
    for(int iii = 0; iii < m_numGaus; iii++){
      //order = amp / mean / sig
      vector<TGraph*> tempGraphVect;
      tempGraphVect.push_back( new TGraph());
      tempGraphVect.push_back( new TGraph());
      tempGraphVect.push_back( new TGraph());
      theGraphs.push_back(tempGraphVect);
    }
  }*/


  vector<double> means(m_numGaus);
  //double means[m_numGaus];
  if(a_initialMeans){
    for(int iii = 0; iii < m_numGaus; iii++){
      means[iii] = a_initialMeans[iii];
    }
  }else{
    bool goodKMeans = doKMeans(&means[0]);
    if(!goodKMeans){
      cout << "KMeans Failed" << endl;
      for(int iii = 0; iii < m_numGaus; iii++){
        means[iii] = lowHistoEdge + histogramRange/((double) m_numGaus + 1)*((double) iii+1);
      }
    }
  }



  vector<double> localParams(3*m_numGaus);
  vector<double> converganceThresholds(3*m_numGaus);
  for(int iii = 0; iii < m_numGaus; iii++){
    //If the amps are specified, then use them
    if(a_initialAmps){
      localParams[3*iii] = a_initialAmps[iii];
    }else{
      localParams[3*iii] = 1.0/((double) m_numGaus);
    }

    localParams[3*iii + 1] = means[iii];

    double test_width = m_dataHistogram->GetBinWidth(m_dataHistogram->FindBin(means[iii]));


    //If the sigmas are specified, then use them
    if(a_initialSigmas){
      localParams[3*iii + 2] = a_initialSigmas[iii]*a_initialSigmas[iii];//a_avgSigma*a_avgSigma;
    }else{
      localParams[3*iii + 2] = test_width;
    }


    converganceThresholds[3*iii] = 0.001;
    converganceThresholds[3*iii+1] = test_width/50.0;
    converganceThresholds[3*iii+2] = test_width*test_width/50.0;
    //std::cout << "Mean " << iii << " initialized to "
    //     << localParams[3*iii + 1] << std::endl;
  }

  //std::cout << "Initial Parameters:" << std::endl;
  //for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
  //  std::cout << "   " << localParams[3*gausIndex+0] << "  "
  //                    << localParams[3*gausIndex+1] << "  "
  //                    << localParams[3*gausIndex+2] << std::endl;
  //}


  //Temporary Data Structures
  vector<double> numPerGaus(m_numGaus);
  vector< vector<double> > responsible;
  vector<double> prevParams(3*m_numGaus);

  for(int iii = 0; iii < nBins; iii++){
    vector<double> tempVect1(m_numGaus);
    responsible.push_back(tempVect1);
  }


  //Trial Loop
  bool converged = false;
  double prevLogLike = 0;
  int trialIndex = 0;
  for(trialIndex = 0; !converged; trialIndex++){
    if(trialIndex > a_maxTrials) break;
    for(int iii = 0; iii < 3*m_numGaus; iii++){
      prevParams[iii] = localParams[iii];
    }
    //Push to canvases
    /*if(showCanvas){
      for(int gausIndex4 = 0; gausIndex4 < m_numGaus; gausIndex4++){
        theGraphs[gausIndex4][0]->SetPoint(trialIndex,trialIndex, prevParams[3*gausIndex4]);
        theGraphs[gausIndex4][1]->SetPoint(trialIndex,trialIndex, prevParams[3*gausIndex4 + 1]);
        theGraphs[gausIndex4][2]->SetPoint(trialIndex,trialIndex, prevParams[3*gausIndex4 + 2]);

        if(trialIndex%25==0){
          paramCanvas->cd(3*gausIndex4 + 1);
          theGraphs[gausIndex4][0]->Draw();
          paramCanvas->cd(3*gausIndex4 + 2);
          theGraphs[gausIndex4][1]->Draw();
          paramCanvas->cd(3*gausIndex4 + 3);
          theGraphs[gausIndex4][2]->Draw();
          paramCanvas->Update();
          gSystem->ProcessEvents();
        }
      }
    }*/


    double numPerGausTotal = 0;
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      numPerGaus[gausIndex] = 0;
    }

    for(int binIndex = 0; binIndex < nBins; binIndex++){
      double x = m_dataHistogram->GetBinCenter(binIndex+1);
      double norm = 0.0;
      for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
        if(m_dataHistogram->GetBinContent(binIndex+1) < a_minBinContent){
          responsible[binIndex][gausIndex]=0;
        }else{
          responsible[binIndex][gausIndex]=localParams[3*gausIndex]
              *TMath::Exp(-((x - localParams[3*gausIndex + 1])
                           *(x - localParams[3*gausIndex + 1]))/
                           (2.0 * localParams[3*gausIndex + 2]))
              /TMath::Sqrt(2.0*3.14159265359*localParams[3*gausIndex + 2]);
          norm += responsible[binIndex][gausIndex];
        }

        // std::cout << "x=" << x << "  denom="
        //       << TMath::Sqrt(2.0*3.14159265359*localParams[3*gausIndex + 2])
        //      << std::endl;

        //std::cout << " resp_"<<binIndex<<"_"<<gausIndex<<"="
        //          << responsible[binIndex][gausIndex] << std::endl;
      }

      //normalize
      for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
        if(norm!=0){
          responsible[binIndex][gausIndex] = responsible[binIndex][gausIndex]/norm;
        }
        numPerGaus[gausIndex] += responsible[binIndex][gausIndex]
                                 *m_dataHistogram->GetBinContent(binIndex+1);
      }

    }

    //Now Compute New Parameters
    //double newParams[3*m_numGaus];
    //Initialize to Zero
    for(int iii = 0; iii < 3*m_numGaus; iii++){
      localParams[iii] = 0;
    }
    //Compute total number
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      numPerGausTotal += numPerGaus[gausIndex];
    }
    //std::cout << "Sum of Mixing Coefficients: " << numPerGausTotal << std::endl;



    //Calculate New Means And Mixing Coefficients
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      if(numPerGausTotal == 0) cout << "ERROR: GausMixModel numPerGausTotal is zero" << endl;
      localParams[3*gausIndex] = numPerGaus[gausIndex]/numPerGausTotal;
      if(localParams[3*gausIndex] == 0) cout << "ERROR: GausMixModel numPerGaus is zero" << endl;
      for(int binIndex = 0; binIndex < nBins; binIndex++){
        double x = m_dataHistogram->GetBinCenter(binIndex+1);
        localParams[3*gausIndex+1]+= x * responsible[binIndex][gausIndex]
                                  * m_dataHistogram->GetBinContent(binIndex+1);
      }
      localParams[3*gausIndex+1] = localParams[3*gausIndex+1]/numPerGaus[gausIndex];
    }

    //Calculate New Variances
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      for(int binIndex = 0; binIndex < nBins; binIndex++){
        double x = m_dataHistogram->GetBinCenter(binIndex+1);
        localParams[3*gausIndex+2]+=responsible[binIndex][gausIndex]
                                  *m_dataHistogram->GetBinContent(binIndex+1)
                                  *(x - localParams[3*gausIndex+1])
                                  *(x - localParams[3*gausIndex+1]);

      }
      localParams[3*gausIndex+2] = localParams[3*gausIndex+2]/numPerGaus[gausIndex];
    }

    //std::cout << "Parameters Trial " << trialIndex << ":" << std::endl;
    //for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    //  std::cout << "   " << localParams[3*gausIndex+0] << "  "
    //                    << localParams[3*gausIndex+1] << "  "
    //                    << localParams[3*gausIndex+2] << std::endl;
    //}

    //Test For Convergance
    //if(trialIndex > 4){
    //  converged = true;
    //  for(int paramIndex = 0; paramIndex < 3*m_numGaus; paramIndex++){
    //    double currentChange = fabs(localParams[paramIndex] - prevParams[paramIndex]);
    //    if(currentChange > converganceThresholds[paramIndex]) converged=false;
    //    //std::cout << "  Param: " << paramIndex << " currentChange: "
    //    //  << currentChange << " threshold: " << converganceThresholds[paramIndex] << std::endl;
    //  }
    //}


    //############# Calculate Log Likelihood #############
    double logLike = 0.0;
    for(int binIndex = 0; binIndex < nBins; binIndex++){
      double x = m_dataHistogram->GetXaxis()->GetBinCenter(binIndex+1);
      double prob = 0.0;
      for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
        if(m_dataHistogram->GetBinContent(binIndex+1) < a_minBinContent){
          continue;
        }else{

          double dx = x - localParams[3*gausIndex+1];
          double expFactor = -0.5 *dx*dx / localParams[3*gausIndex+2];

          prob += localParams[3*gausIndex]*TMath::Exp(expFactor)/TMath::Sqrt(2.0*3.14159265359*localParams[3*gausIndex + 2]);
          //cout << "Prob = " << prob << endl;
        }
      }
      if(prob > 0.0){
        logLike += TMath::Log(prob) * m_dataHistogram->GetBinContent(binIndex+1);
      }
    }// end X Bin Loop

    //cout << "LogLike: " << logLike << endl;
    //Test For Convergance
    if(trialIndex > 4){
      converged = true;
      double dev = fabs((logLike - prevLogLike)/prevLogLike);
      //cout << "Percent Deviation: " << dev << endl;
      if(dev > a_tolerance ){
        converged = false;
      }
    }

    prevLogLike = logLike;



  }//End of Trial Loop


  m_parameters.clear();
  for(int paramIndex = 0; paramIndex < 3*m_numGaus; paramIndex++){
    a_parametersToFill[paramIndex] = localParams[paramIndex];
    m_parameters.push_back(localParams[paramIndex]);
  }



  //chi squared convergence
  //double chidof = getChiSquared()/((double) m_dataHistogram->GetNbinsX());
  //cout << "Chi Squared per dof: " << chidof << endl;
  //if(chidof > 15){
  //  cout << "Failed" << endl;
  //}

  //std::cout << "Number of Trials: " << trialIndex << std::endl;
  //std::cout << "Final Parameters:" << std::endl;
  //std::cout << "      Weight      Mean    Variance" << std::endl;
  //for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
  //  std::cout << "     " << localParams[3*gausIndex+0] << "   "
  //                    << localParams[3*gausIndex+1] << "   "
  //                    << localParams[3*gausIndex+2] << std::endl;
  //}
  //  std::cout << std::endl << std::endl << std::endl;
}






//_____________________________________________________________________________
//_____________________________________________________________________________






void GausMixModel::doGMM2D(double* a_parametersToFill,
               double* a_initialMeans, double* a_initialAmps,
               double a_avgXSigma, double a_avgYSigma,
               int a_maxTrials, double a_tolerance, bool a_forceZeroCovariance){

  //cout << "Starting Gaussian Mixture Model in Two Dimensions" << endl;
  //cout << "Number of Gaussians: " << m_numGaus << endl;
  bool showCanvas = true;
  //double converganceThreshold = 1.0E-8;

  double minBinContent = 1;


  int nXBins = m_dataHistogram2D->GetNbinsX();
  int nYBins = m_dataHistogram2D->GetNbinsY();
  //int lowBinIndex = HistogramUtilities::firstBinWithContent(m_dataHistogram2D);
  //int highBinIndex = HistogramUtilities::lastBinWithContent(m_dataHistogram2D);
  //double lowHistoEdge = m_dataHistogram2D->GetBinLowEdge(lowBinIndex);
  //double highHistoEdge = m_dataHistogram2D->GetBinLowEdge(highBinIndex+1);
  //nBins = highBinIndex - lowBinIndex;


  double lowHistoEdgeX = m_dataHistogram2D->GetXaxis()->GetBinLowEdge(1);
  double highHistoEdgeX = m_dataHistogram2D->GetXaxis()->GetBinLowEdge(nXBins);
  double histogramRangeX = highHistoEdgeX - lowHistoEdgeX;
  double lowHistoEdgeY = m_dataHistogram2D->GetYaxis()->GetBinLowEdge(1);
  double highHistoEdgeY = m_dataHistogram2D->GetYaxis()->GetBinLowEdge(nYBins);
  double histogramRangeY = highHistoEdgeY - lowHistoEdgeY;


  //SETUP THE GRAPHS TO MONITOR PARAMETERS WHILE FITTING
  //vector< vector<TGraph*> > theGraphs;
  //TCanvas* paramCanvas = NULL;
  //if(showCanvas){
  //  paramCanvas = new TCanvas("Parameter Canvas");
  //  paramCanvas->Divide(3,m_numGaus);
  //  for(int iii = 0; iii < m_numGaus; iii++){
  //    //order = amp / mean / sig
  //    vector< TGraph* > tempGraphVect = { new TGraph(), new TGraph(), new TGraph()};
  //    theGraphs.push_back(tempGraphVect);
  //  }
  //}



  //INITIALIZE THE PARAMETERS
  // STRUCTURE [amp, meanX, meanY, cov11, cov12, cov22] for each gaus  so [gaus#][param]
  //vector< vector<double>> localParams;
  std::vector< std::vector<double> > localParams;//[m_numGaus][6];
  //cout << "initial local params: " << localParams.size()<< endl;
  for(int kkk = 0; kkk < m_numGaus; kkk++){
    //cout << kkk << " " << m_numGaus << " ";
    vector<double> tempVect51;
    for(int jjj = 0; jjj < 6; jjj++) tempVect51.push_back(0);
    //cout << tempVect51.size() << endl;
    localParams.push_back(tempVect51);
  }
  //cout << "local params: " << localParams.size() << " " << localParams[0].size() << endl;


  for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    //Amplitudes
    if(a_initialAmps){
      localParams[gausIndex][0] = a_initialAmps[gausIndex];
    }else{
      localParams[gausIndex][0] = 1.0/((double) m_numGaus);
    }

    //MEANS
    if(a_initialMeans){
      localParams[gausIndex][1] = a_initialMeans[2*gausIndex];
      localParams[gausIndex][2] = a_initialMeans[2*gausIndex+1];
    }else{
      localParams[gausIndex][1] = lowHistoEdgeX + histogramRangeX*((double) gausIndex +1)/((double) m_numGaus + 1);
      localParams[gausIndex][2] = lowHistoEdgeY + histogramRangeY/2.0;
    }

    //Covariances
    if(a_avgXSigma >= 0){
     localParams[gausIndex][3] = a_avgXSigma*a_avgXSigma;
    }else{
     localParams[gausIndex][3] = histogramRangeX*histogramRangeX/25.0;
    }
    if(a_avgYSigma >= 0){
     localParams[gausIndex][5] = a_avgYSigma*a_avgYSigma;
    }else{
     localParams[gausIndex][5] = histogramRangeY*histogramRangeY/25.0;
    }
    localParams[gausIndex][4] =0.0; //no starting correlation


  }

  //cout << "done intializing variables" << endl;


  //just means for convergence thresholds
  //double converganceThresholds[m_numGaus][2];
  //for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
  //  converganceThresholds[gausIndex][0] = 0.001;
  //  converganceThresholds[gausIndex][1] = 0.001;
  //}

  // PRINT INITIAL PARAMETERS
  //std::cout << endl << "Initial Parameters:" << std::endl;
  //for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
  //  for(int paramIndex = 0; paramIndex < 6; paramIndex++){
      //cout << "   "
      //     << localParams[gausIndex][paramIndex];
    //}
    //cout << endl;
  //}
  //cout << endl;
  //cout << "Done Initializing Parameters" << endl;


  //Temporary Data Structures
  vector<double> numPerGaus(m_numGaus);
  for(int jjj = 0; jjj < m_numGaus; jjj++) numPerGaus.push_back(0);

  vector< vector<double> > prevParams;//[m_numGaus][6];
  for(int iii = 0; iii < m_numGaus; iii++){
    vector<double> tempVect2;
    for(int jjj = 0; jjj < 6; jjj++) tempVect2.push_back(0);
    prevParams.push_back(tempVect2);
  }

  vector< vector< vector<double> > > responsible;//[nXBins][nYBins][m_numGaus];
  for(int iii = 0; iii < nXBins; iii++){
    vector< vector<double> > tempVect3;
    for(int jjj = 0; jjj < nYBins; jjj++){
      vector<double> tempVect4;
      for(int kkk = 0; kkk < m_numGaus; kkk++) tempVect4.push_back(0);
      tempVect3.push_back(tempVect4);
    }
    responsible.push_back(tempVect3);
  }





  double prevLogLike = 0;

  //############# Trial Loop #############
  bool converged = false;
  int trialIndex = 0;
  for(trialIndex = 0; !converged; trialIndex++){
    if(trialIndex > a_maxTrials) break;
    //store current parameters as past parameters
    //reset numPerGausTotal
    //cout << "Copying parameters to previous parameters" << endl;
    double numPerGausTotal = 0;
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      numPerGaus[gausIndex] = 0;
      for(int paramIndex = 0; paramIndex < 6; paramIndex++){
        prevParams[gausIndex][paramIndex] = localParams[gausIndex][paramIndex];
      }
    }

    //############# Calculate Prefactors and Determinants to speed up calculation #############
    //cout << "calculating prefactors" << endl;
    vector<double> prefactors;
    for(int jjj = 0; jjj < m_numGaus; jjj++) prefactors.push_back(0);
    vector<double> determinants;
    for(int jjj = 0; jjj < m_numGaus; jjj++) determinants.push_back(0);
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      determinants[gausIndex] = localParams[gausIndex][3]*localParams[gausIndex][5]
                               -localParams[gausIndex][4]*localParams[gausIndex][4];
      prefactors[gausIndex] = 0.5/(3.14159265359*TMath::Sqrt(determinants[gausIndex]));

      //cout << "det" << gausIndex << ": " << determinants[gausIndex] << endl;
      //cout << "pre" << gausIndex << ": " << prefactors[gausIndex] << endl;
    }
    //cout << endl;


    //########## Calculate Responsibilities #############
    //cout << "calculating responsibilities" << endl;
    for(int binXIndex = 0; binXIndex < nXBins; binXIndex++){
      double x = m_dataHistogram2D->GetXaxis()->GetBinCenter(binXIndex+1);

      for(int binYIndex = 0; binYIndex < nYBins; binYIndex++){
        double y = m_dataHistogram2D->GetYaxis()->GetBinCenter(binYIndex+1);

        double norm = 0.0;
        for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
          if(m_dataHistogram2D->GetBinContent(binXIndex+1,binYIndex+1) < minBinContent){
            responsible[binXIndex][binYIndex][gausIndex] = 0;
          }else{
            double tx = x - localParams[gausIndex][1];
            double ty = y - localParams[gausIndex][2];
            double expFactor = 0.5*(
                    2.0*localParams[gausIndex][4]*tx*ty
                    - localParams[gausIndex][3]*ty*ty
                    - localParams[gausIndex][5]*tx*tx
                 ) / determinants[gausIndex];

            responsible[binXIndex][binYIndex][gausIndex]
                = localParams[gausIndex][0] * prefactors[gausIndex]*TMath::Exp(expFactor);
            norm += responsible[binXIndex][binYIndex][gausIndex];
          }

        }

        //Normalize Responsibilities
        for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
          if(norm!=0){
            responsible[binXIndex][binYIndex][gausIndex] = responsible[binXIndex][binYIndex][gausIndex]/norm;
          }
          numPerGaus[gausIndex] += responsible[binXIndex][binYIndex][gausIndex]
                                   *m_dataHistogram2D->GetBinContent(binXIndex+1, binYIndex+1);
        }
      }//end Y Bin Loop
    }// end X Bin Loop

    //Now Compute New Parameters
    //cout << "Computing new parameters" << endl;
    //cout << "numGaus = " << m_numGaus << endl;
    //Initialize to Zero
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      for(int paramIndex = 0; paramIndex < 6; paramIndex++){
        localParams[gausIndex][paramIndex] = 0;
      }
    }
    //Compute total number
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      numPerGausTotal+= numPerGaus[gausIndex];
    }
    //std::cout << "Sum of Mixing Coefficients: " << numPerGausTotal << std::endl;



    //Calculate New Means And Mixing Coefficients
    //cout << "calculting means" << endl;
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      if(numPerGausTotal ==0) cout << "ERROR: GausMixModel numPerGausTotal is zero" << endl;
      localParams[gausIndex][0] = numPerGaus[gausIndex]/numPerGausTotal;
      if(localParams[gausIndex][0] == 0) cout << "   ERROR:   #############  numPerGaus is zero ###############" << endl;

      for(int binXIndex = 0; binXIndex < nXBins; binXIndex++){
        double x = m_dataHistogram2D->GetXaxis()->GetBinCenter(binXIndex+1);
        for(int binYIndex = 0; binYIndex < nYBins; binYIndex++){
          double y = m_dataHistogram2D->GetYaxis()->GetBinCenter(binYIndex+1);
          localParams[gausIndex][1] += x * responsible[binXIndex][binYIndex][gausIndex]
                                  *m_dataHistogram2D->GetBinContent(binXIndex+1, binYIndex+1);
          localParams[gausIndex][2] += y * responsible[binXIndex][binYIndex][gausIndex]
                                  *m_dataHistogram2D->GetBinContent(binXIndex+1, binYIndex+1);
        }//end y bin loop
      }
      localParams[gausIndex][1] = localParams[gausIndex][1]/numPerGaus[gausIndex];
      localParams[gausIndex][2] = localParams[gausIndex][2]/numPerGaus[gausIndex];
    }

    //Calculate New Variances
    //cout << "calculating variances" << endl;
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      for(int binXIndex = 0; binXIndex < nXBins; binXIndex++){
        double x = m_dataHistogram2D->GetXaxis()->GetBinCenter(binXIndex+1);

        for(int binYIndex = 0; binYIndex < nYBins; binYIndex++){
          double y = m_dataHistogram2D->GetYaxis()->GetBinCenter(binYIndex+1);
          double dx = x - localParams[gausIndex][1];
          double dy = y - localParams[gausIndex][2];
          localParams[gausIndex][3]+=responsible[binXIndex][binYIndex][gausIndex]
                    *m_dataHistogram2D->GetBinContent(binXIndex+1,binYIndex+1)
                    *dx*dx;
          if(!a_forceZeroCovariance){
            localParams[gausIndex][4]+=responsible[binXIndex][binYIndex][gausIndex]
                    *m_dataHistogram2D->GetBinContent(binXIndex+1,binYIndex+1)
                    *dx*dy;
          }
          localParams[gausIndex][5]+=responsible[binXIndex][binYIndex][gausIndex]
                    *m_dataHistogram2D->GetBinContent(binXIndex+1,binYIndex+1)
                    *dy*dy;

        }//end y bin loop
      }// end x bin loop
      //normalize
      localParams[gausIndex][3] = localParams[gausIndex][3]/numPerGaus[gausIndex];
      localParams[gausIndex][4] = localParams[gausIndex][4]/numPerGaus[gausIndex];
      localParams[gausIndex][5] = localParams[gausIndex][5]/numPerGaus[gausIndex];

    }

    //std::cout << "Parameters Trial " << trialIndex << ":" << std::endl;
    //for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    //  std::cout << "   " << localParams[3*gausIndex+0] << "  "
    //                    << localParams[3*gausIndex+1] << "  "
    //                    << localParams[3*gausIndex+2] << std::endl;
    //}



    //############# Calculate Prefactors and Determinants Again... #############
    //double prefactors[m_numGaus];
    //double determinants[m_numGaus];
    for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      determinants[gausIndex] = localParams[gausIndex][3]*localParams[gausIndex][5]
                               -localParams[gausIndex][4]*localParams[gausIndex][4];
      prefactors[gausIndex] = 0.5/(3.14159265359*TMath::Sqrt(determinants[gausIndex]));
    }


    //############# Calculate Log Likelihood #############
    double logLike = 0.0;
    for(int binXIndex = 0; binXIndex < nXBins; binXIndex++){
      double x = m_dataHistogram2D->GetXaxis()->GetBinCenter(binXIndex+1);

      for(int binYIndex = 0; binYIndex < nYBins; binYIndex++){
        double y = m_dataHistogram2D->GetYaxis()->GetBinCenter(binYIndex+1);

        double prob = 0.0;
        for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
          if(m_dataHistogram2D->GetBinContent(binXIndex+1,binYIndex+1) < minBinContent){
            continue;
          }else{
            double dx = x - localParams[gausIndex][1];
            double dy = y - localParams[gausIndex][2];
            double expFactor = 0.5*(
                    2.0*localParams[gausIndex][4]*dx*dy
                    - localParams[gausIndex][3]*dy*dy
                    - localParams[gausIndex][5]*dx*dx
                 ) / determinants[gausIndex];


            prob += localParams[gausIndex][0] * prefactors[gausIndex]*TMath::Exp(expFactor);
            //cout << " -p" << prob;
            //cout << "  e" << expFactor;
            //cout << "  pre" << prefactors[gausIndex];
            //cout << "  a" << localParams[gausIndex][0];
          }

        }
        if(prob > 0.0){
          logLike += TMath::Log(prob) * m_dataHistogram2D->GetBinContent(binXIndex+1, binYIndex+1);
        }


      }//end Y Bin Loop
    }// end X Bin Loop


    //print out the parameters
    //cout << "Trial " << trialIndex << ": " << "LogLike= " << logLike;
    //for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    //  cout << endl;
    //  for(int paramIndex = 0; paramIndex < 6; paramIndex++){
    //    cout << "  " << localParams[gausIndex][paramIndex];
    //  }
    //}
    //cout << endl;




    //Test For Convergance
    if(trialIndex > 4){
      converged = true;
      double dev = fabs((logLike - prevLogLike)/prevLogLike);
      //cout << "Percent Deviation: " << dev << endl;
      if(dev > a_tolerance ){
        converged = false;
      }
    }

    prevLogLike = logLike;



  }//End of Trial Loop

  //cout << "FINISHED TRIAL LOOP " << endl;
  m_parameters.clear();
  for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    for(int paramIndex = 0; paramIndex < 6; paramIndex++){
      a_parametersToFill[6*gausIndex + paramIndex] = localParams[gausIndex][paramIndex];
      m_parameters.push_back(localParams[gausIndex][paramIndex]);
    }
  }



  //chi squared convergence
  //double chidof = getChiSquared()/((double) m_dataHistogram2D->GetNbinsX());
  //cout << "Chi Squared per dof: " << chidof << endl;
  //if(chidof > 15){
  //  cout << "Failed" << endl;
  //}

}








TF1* GausMixModel::getFunction(){
     //double lowBound = m_dataHistogram->GetBinLowEdge(1);
     //double highBound = m_dataHistogram->GetBinLowEdge(m_dataHistogram->GetNbinsX()+1);
  int lowBinIndex = HistogramUtilities::firstBinWithContent(m_dataHistogram);
  int highBinIndex = HistogramUtilities::lastBinWithContent(m_dataHistogram);
  double lowHistoEdge = m_dataHistogram->GetBinLowEdge(lowBinIndex);
  double highHistoEdge = m_dataHistogram->GetBinLowEdge(highBinIndex+1);

  TF1* funct = new TF1("mixGausFunct",*this,lowHistoEdge,highHistoEdge,1);
  funct->SetNpx(2000);
  return funct;

}


TF1* GausMixModel::getSingleGaus(int a_gausIndex){
  double fraction = m_parameters[3*a_gausIndex];
  double mean     = m_parameters[3*a_gausIndex+1];
  double stdDev   = sqrt(m_parameters[3*a_gausIndex+2]);

  TF1* gausFunct = new TF1(Form("gausFunct_%02d",a_gausIndex),
    "[0]*0.3989422804 * TMath::Exp(-0.5*( x - [1] )*(x - [1])/([2]*[2])) /[2]",
    mean - 20*stdDev,mean + 20*stdDev);
  gausFunct->SetNpx(500);
  gausFunct->SetLineColor(kGreen);

  gausFunct->SetParameter(0,fraction*m_dataHistogram->Integral("width"));
  gausFunct->SetParameter(1,mean);
  gausFunct->SetParameter(2,stdDev);

  return gausFunct;
}



TH2D* GausMixModel::getTGraphFor2D(){
  TH2D* returnHisto = (TH2D*) m_dataHistogram2D->Clone();

  double integral = 0;
  for(int binX = 1; binX <= returnHisto->GetNbinsX(); binX++){
    for(int binY = 1; binY <= returnHisto->GetNbinsY(); binY++){
      integral += returnHisto->GetBinContent(binX,binY);
    }
  }

  vector<double> prefactors(m_numGaus);
  vector<double> determinants(m_numGaus);
  for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
      determinants[gausIndex] = m_parameters[6*gausIndex+3]*m_parameters[6*gausIndex+5]
                               -m_parameters[6*gausIndex+4]*m_parameters[6*gausIndex+4];
      prefactors[gausIndex] = 0.5/(3.14159265359*TMath::Sqrt(determinants[gausIndex]));
  }




  for(int binX = 0; binX <= returnHisto->GetNbinsX(); binX++){
   double x = returnHisto->GetXaxis()->GetBinCenter(binX+1);
    for(int binY = 0; binY <= returnHisto->GetNbinsY(); binY++){
      double y = returnHisto->GetYaxis()->GetBinCenter(binY+1);
      //double newContent = returnHisto->GetBinContent(binX,binY)/(integral
      //                    *returnHisto->GetXaxis()->GetBinWidth(binX)
      //                    *returnHisto->GetYaxis()->GetBinWidth(binY) );
      double newContent =0;
      for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
        //returnHisto->SetBinContent(binX,binY,newContent);
        double tx = x - m_parameters[6*gausIndex+1];
        double ty = y - m_parameters[6*gausIndex+2];
        double expFactor = 0.5*(
                    2.0*m_parameters[6*gausIndex+4]*tx*ty
                    - m_parameters[6*gausIndex+3]*ty*ty
                    - m_parameters[6*gausIndex+5]*tx*tx
                 ) / determinants[gausIndex];
        //newContent = newContent - m_parameters[6*gausIndex] * prefactors[gausIndex]*TMath::Exp(expFactor);
        newContent = newContent + m_parameters[6*gausIndex] * prefactors[gausIndex]*TMath::Exp(expFactor);
      }//end gaus loop

      returnHisto->SetBinContent(binX,binY,newContent);
    } // end y bin loop
  }//end x bin loop

  return returnHisto;
}






//_____________________________________________________________________________
double GausMixModel::operator() (double *x, double *p){
  double val = 0;
  for(int gausIndex = 0; gausIndex < m_numGaus; gausIndex++){
    val += m_parameters[3*gausIndex]
            *TMath::Exp(-((*x - m_parameters[3*gausIndex + 1])
                         *(*x - m_parameters[3*gausIndex + 1]))/
                         (2.0 * m_parameters[3*gausIndex + 2]))
            /TMath::Sqrt(2.0*3.14159265359*m_parameters[3*gausIndex + 2]);
  }
  val = val*HistogramUtilities::sumOfContent(m_dataHistogram) * m_dataHistogram->GetBinWidth(1);
  return val;
}

//_____________________________________________________________________________

double GausMixModel::getChiSquared(){
  double chiSqr = 0;
  TF1* function = getFunction();
  for(int iii = 1; iii <= m_dataHistogram->GetNbinsX(); iii++){
    double content = m_dataHistogram->GetBinContent(iii);
    double predict = function->Eval(m_dataHistogram->GetBinCenter(iii));
    if(content > 0.0){
      chiSqr+= (content - predict)*(content- predict)/content;
    }

  }
  return chiSqr;
}


//_____________________________________________________________________________
bool GausMixModel::isTooManyGaus(){
  double toleranceMean = 0.03; // how close means
  double toleranceVar = 0.05;
  bool isGood = true;
  for(int gausIndex1 = 0; gausIndex1 < m_numGaus; gausIndex1++){
    for(int gausIndex2 = 0; gausIndex2 < m_numGaus; gausIndex2++){
      double quantMean =  2.0*(m_parameters.at(gausIndex1*3 + 1) - m_parameters.at(gausIndex2*3 + 1))/
               (m_parameters.at(gausIndex1*3 + 1) + m_parameters.at(gausIndex2*3 + 1));
      double quantVar =  2.0*(m_parameters.at(gausIndex1*3 + 2) - m_parameters.at(gausIndex2*3 + 2))/
               (m_parameters.at(gausIndex1*3 + 2) + m_parameters.at(gausIndex2*3 + 2));

      if( fabs(quantMean) < toleranceMean || fabs(quantVar) < toleranceVar ) isGood = false;



    }
  }

  return isGood;


}


//_____________________________________________________________________________

double GausMixModel::getYield(int a_gausNumber){
 return m_parameters[3*a_gausNumber]*HistogramUtilities::sumOfContent(m_dataHistogram);
}




#endif
