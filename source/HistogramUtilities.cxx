#ifndef _HISTOGRAM_UTILITIES_CXX_
#define _HISTOGRAM_UTILITIES_CXX_

//#include "../headers/PhysMath.h"
#include "../headers/HistogramUtilities.h"

//##############  FUNCITON DEFINITIONS ################


//void HistogramUtilities::setUniqueNum(unsigned long a_num){g_histoUtil_uniqueIdentifier = a_num;}


void HistogramUtilities::setBinValues(TH1D* a_histo, int a_bin, double a_value, double a_error){
  if(!a_histo) return;
  a_histo->SetBinContent(a_bin,a_value);
  a_histo->SetBinError(a_bin,a_error);
}
void HistogramUtilities::setBinValues(TH2D* a_histo, int a_binX, int a_binY, double a_value, double a_error){
  if(!a_histo) return;
  a_histo->SetBinContent(a_binX,a_binY,a_value);
  a_histo->SetBinError(a_binX,a_binY,a_error);
}


// will return negative bin if empty;
int HistogramUtilities::firstBinWithContent(TH1D* a_histo){
  if(!a_histo) return -999;
  int numBins = a_histo->GetNbinsX();
  for(int bin = 1; bin <= numBins; bin++){
    if(a_histo->GetBinContent(bin) != 0.0) {
      return bin;
    }
  }
  return numBins+1;
}

int HistogramUtilities::lastBinWithContent(TH1D* a_histo){
  if(!a_histo) return -999;
  int numBins = a_histo->GetNbinsX();
  int returnBin = -1;
  for(int binIndex = numBins; binIndex > 0; binIndex--){
    if(a_histo->GetBinContent(binIndex) != 0.0) {
      returnBin = binIndex;
      break;
    }
  }
  return returnBin;
}
int HistogramUtilities::firstYBinWithContent(TH2D* a_histo, int a_binX){
  if(!a_histo) return -999;
  int numBins = a_histo->GetNbinsY();
  for(int binY = 1; binY <= numBins; binY++){
    double content = a_histo->GetBinContent(a_binX,binY);
    double error   = a_histo->GetBinError(a_binX,binY);
    if(std::isfinite(content) && ( content != 0.0 || error > 0.0) ){
      return binY;
    }
  }
  return numBins+1;
}

int HistogramUtilities::lastYBinWithContent(TH2D* a_histo, int a_binX){
  if(!a_histo) return -999;
  int numBins = a_histo->GetNbinsY();
  int returnBin = -1;
  for(int binY = numBins; binY > 0; binY--){
    double content = a_histo->GetBinContent(a_binX,binY);
    double error   = a_histo->GetBinError(a_binX,binY);
    if(std::isfinite(content) && ( content != 0.0 || error > 0.0) ){
      returnBin = binY;
      break;
    }
  }
  return returnBin;
}

int HistogramUtilities::numBinsWithContentAboveThreshold(TH1D* a_histo, double a_threshold){
  if(!a_histo) return -999;
  int bins = 0;
  for(int binIndex = 1; binIndex <= a_histo->GetNbinsX(); binIndex++){
    if(a_histo->GetBinContent(binIndex) > a_threshold){
      bins++;
    }
  }
  return bins;
}

bool HistogramUtilities::isEmpty(TH1D* a_histo){
  if(!a_histo) return true;
  for(int binIndex = 1; binIndex <= a_histo->GetNbinsX(); binIndex++){
    if(a_histo->GetBinContent(binIndex) != 0.0){
      return false;
    }
  }
  return true;
}

bool  HistogramUtilities::isEmpty2D(TH2D* a_histo){
  if(!a_histo) return true;
  for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      if(a_histo->GetBinContent(binX,binY) != 0.0){
        return false;
      }
    }
  }
  return true;
}

bool HistogramUtilities::isYSliceEmpty(TH2D* a_histo, int a_binX){
  if(!a_histo) return true;
  int binsWithContent = 0;
  for(int binIndex = 1; binIndex <= a_histo->GetNbinsY(); binIndex++){
    double content = a_histo->GetBinContent(a_binX,binIndex);
    if(std::isfinite(content) &&  content != 0.0){
      return false;
    }
  }
  return true;
}

bool HistogramUtilities::outsideHistogramRange(TH1D* a_histo, double a_value){
  if(!a_histo) return true;
  return (a_value < a_histo->GetBinLowEdge(1) || a_value > a_histo->GetBinLowEdge(a_histo->GetNbinsX()+1));
}


double HistogramUtilities::getLowEdge(TH1D* a_histo){
  if(!a_histo) return -999;
  return a_histo->GetXaxis()->GetBinLowEdge(1);
}

double HistogramUtilities::getHighEdge(TH1D* a_histo){
  if(!a_histo) return -999;
  return a_histo->GetXaxis()->GetBinLowEdge(a_histo->GetNbinsX()+1);
}


double HistogramUtilities::getBinContentAvg(TH1D* a_histo, int a_bin, int a_numBinsToAvg){
  if(!a_histo) return -999;
  double sum = a_histo->GetBinContent(a_bin);
  double numBinsUsed = 1.0;
  for(int binChange = 1; binChange <= a_numBinsToAvg; binChange++ ){
    sum += a_histo->GetBinContent(a_bin + binChange);
    sum += a_histo->GetBinContent(a_bin - binChange);
    numBinsUsed += 2.0;
  }
  return sum/numBinsUsed;
}


int HistogramUtilities::getMaxBinInRange(TH1D* a_histo, double a_lowRange, double a_highRange){
  if(!a_histo) return -999;
  int startBin = 1;
  if(a_lowRange != -999){
    startBin = a_histo->FindBin(a_lowRange);
    if(startBin <= 0) startBin = 1;
    if(startBin > a_histo->GetNbinsX()) startBin = a_histo->GetNbinsX();
  }


  int endBin = a_histo->GetNbinsX();
  if(a_lowRange != -999){
    endBin = a_histo->FindBin(a_highRange);
    if(startBin <= 0) endBin = 1;
    if(startBin > a_histo->GetNbinsX()) endBin = a_histo->GetNbinsX();
  }
  if(endBin < startBin){
    cout << " WARNING: HistogramUtilities::getMaxBinInRange the range supplied is out of order or off the range of the histogram. Range = [" << a_lowRange << ", " << a_highRange << "]" << endl;
    return 1;
  }

  int highestBin = startBin;
  double highestVal = a_histo->GetBinContent(highestBin);
  for(int bin = highestBin+1; bin <= endBin; bin++){
    if(a_histo->GetBinContent(bin) > highestVal){
      highestVal = a_histo->GetBinContent(bin);
      highestBin = bin;
    }
  }
  return highestBin;
}

double HistogramUtilities::getMaxInRange(TH1D* a_histo, double a_lowRange, double a_highRange){
  if(!a_histo) return -999;
  return a_histo->GetBinContent(getMaxBinInRange(a_histo,a_lowRange,a_highRange));
}

int HistogramUtilities::getMinBinInRange(TH1D* a_histo, double a_lowRange, double a_highRange){
  if(!a_histo) return -999;
  int lowestBin = -1;
  if(a_lowRange == -999) lowestBin = 1;
  else lowestBin = a_histo->FindBin(a_lowRange);
  double lowestVal = a_histo->GetBinContent(lowestBin);

  int endBin = -1;
  if(a_highRange == -999) endBin = a_histo->GetNbinsX();
  else endBin = a_histo->FindBin(a_highRange);

  for(int bin = lowestBin+1; bin <= endBin; bin++){
    if(a_histo->GetBinContent(bin) < lowestVal){
      lowestVal = a_histo->GetBinContent(bin);
      lowestBin = bin;
    }
  }
  return lowestBin;
}

double HistogramUtilities::getMinInRange(TH1D* a_histo, double a_lowRange, double a_highRange){
  if(!a_histo) return -999;
  return a_histo->GetBinContent(getMinBinInRange(a_histo,a_lowRange,a_highRange));
}



double HistogramUtilities::getSmallestBinWidth(TH1D* a_histo, int a_lowBin, int a_highBin){
  if(!a_histo) return -999;
  int lowBin = 1;
  int highBin = a_histo->GetNbinsX();
  if(a_lowBin > 0) lowBin = a_lowBin;
  if(a_highBin > 0) highBin = a_highBin;
  double smallestWidth = 1E99;
  for(int bin = lowBin; bin <= highBin; bin++){
    if(a_histo->GetBinWidth(bin) < smallestWidth){
      smallestWidth = a_histo->GetBinWidth(bin);
    }
  }
  return smallestWidth;
}



void HistogramUtilities::getPeakRange(TH1D* a_histo, int a_startBin, int &a_lowBinReturn, int &a_highBinReturn, double a_percentageRequired){
  if(!a_histo) return;
  double peakVal = a_histo->GetBinContent(a_startBin);
  double threshold = peakVal*a_percentageRequired;
  //track right
  for(int bin = a_startBin; bin <= a_histo->GetNbinsX(); bin++){
    if(a_histo->GetBinContent(bin) < threshold) break;
    a_highBinReturn = bin;
  }
  for(int bin = a_startBin; bin > 0; bin--){
    if(a_histo->GetBinContent(bin) < threshold) break;
    a_lowBinReturn = bin;
  }
  return;
}


void HistogramUtilities::printHistoInfo(TH1D* a_histo){
  if(!a_histo){
    cout << "HistogramUtilities::printHistoInfo() received an invalid pointer: " << a_histo << endl;
    return;
  }
  cout << "Histogram Name: " << a_histo->GetName() << endl;
  cout << "   Histogram Title: " << a_histo->GetTitle() << endl;
  cout << "   Histogram X axis: " << a_histo->GetXaxis()->GetTitle() << endl;
  cout << "   Histogram Y axis: " << a_histo->GetYaxis()->GetTitle() << endl;
  cout << "   Number of Entries: " << a_histo->GetEntries() << endl;
  cout << "   Number of Bins:  " << a_histo->GetNbinsX()<< endl;
  cout << "   Range of X-axis : [" << a_histo->GetBinLowEdge(1) << "," << a_histo->GetBinLowEdge(a_histo->GetNbinsX()+1) <<"]" << endl;

  cout << "   Underflow : Content=" << a_histo->GetBinContent(0)  << endl;
  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    cout << "   Bin : " << bin << "\tCenter : " << a_histo->GetBinCenter(bin) << "\tContent = " << a_histo->GetBinContent(bin) << "\tError = " << a_histo->GetBinError(bin)
         << "\tWidth = " << a_histo->GetBinWidth(bin) <<"\tLowEdge = " << a_histo->GetBinLowEdge(bin) << "\tHighEdge = " << a_histo->GetBinLowEdge(bin+1) << endl;
  }
  cout << "   Overflow : Content=" << a_histo->GetBinContent(a_histo->GetNbinsX()+1) << endl;

  return;

}




void HistogramUtilities::printHistoInfo(TH2D* a_histo){
  if(!a_histo){
    cout << "HistogramUtilities::printHistoInfo() received an invalid pointer: " << a_histo << endl;
    return;
  }
  cout << "Histogram Name: " << a_histo->GetName() << endl;
  cout << "   Histogram Title: " << a_histo->GetTitle() << endl;
  cout << "   Histogram X axis: " << a_histo->GetXaxis()->GetTitle() << endl;
  cout << "   Histogram Y axis: " << a_histo->GetYaxis()->GetTitle() << endl;
  cout << "   Number of Entries: " << a_histo->GetEntries() << endl;
  cout << "   Number of Bins X:  " << a_histo->GetNbinsX() << endl;
  cout << "   Range of X-axis : [" << a_histo->GetXaxis()->GetBinLowEdge(1) << "," << a_histo->GetXaxis()->GetBinLowEdge(a_histo->GetNbinsX()+1) <<"]" << endl;
  cout << "   Number of Bins Y:  " << a_histo->GetNbinsY() << endl;
  cout << "   Range of Y-axis : [" << a_histo->GetYaxis()->GetBinLowEdge(1) << "," << a_histo->GetYaxis()->GetBinLowEdge(a_histo->GetNbinsX()+1) <<"]" << endl;

  int numInfNanBins = 0;
  for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      double content = a_histo->GetBinContent(binX,binY);
      if(!std::isfinite(content)){
        cout << "   BinX : " << binX << " BinY : " << binY << "\tXCenter : " << a_histo->GetXaxis()->GetBinCenter(binX) << "\tYCenter : " << a_histo->GetYaxis()->GetBinCenter(binY) << "\tContent = " << a_histo->GetBinContent(binX,binY) << "\tError = " << a_histo->GetBinError(binX,binY) << endl;
      }
    }
  }

  return;

}






#ifdef _CPP11_
vector< vector< double > > HistogramUtilities::findPeaks(TH1D* a_histo){

  if(!a_histo){
    vector< vector< double > > blankVect;
    return blankVect;
  }

  vector< double > peakLoc;
  vector< double > peakBinRank;
  int numPeaks = 0;
  int numBins = a_histo->GetNbinsX();
  for(int binIndex = 1; binIndex <= numBins; binIndex++){
    peakLoc.push_back(a_histo->GetBinCenter(binIndex));
    peakBinRank.push_back(0);
    numPeaks += 1;
    double currentSum = 0.0;
    double numRangeBins = 0.0;
    for(int rangeIndex = 1; rangeIndex  <= numBins; rangeIndex++){
      int lowBin = binIndex - rangeIndex;
      int highBin = binIndex + rangeIndex;
      if(lowBin > 0 && lowBin <= numBins){
        currentSum += a_histo->GetBinContent(lowBin);
        numRangeBins += 1.0;
      }
      if(highBin > 0 && highBin <= numBins){
        currentSum += a_histo->GetBinContent(highBin);
        numRangeBins += 1.0;
      }

      if(currentSum < a_histo->GetBinContent(binIndex) * numRangeBins){
        peakBinRank[numPeaks - 1] = rangeIndex;
      }else{
        break;
      }
    }// range loop
  }// bin loop


  //sort the values to the best peaks first
  PhysMath::SimultaneousSortDoubleInt( peakBinRank, peakLoc);
  vector< vector< double > > returnVect = {peakLoc, peakBinRank};
  return returnVect;
}
#endif



void HistogramUtilities::interpolate(TH1D* a_histo, double a_x, double &a_y, double &a_yErr){
  if( !a_histo || a_histo->GetNbinsX() < 2 ) return;
  int binX = a_histo->FindBin(a_x);
  int firstBin  = 0;
  int secondBin = 0;
  if(binX <= 1 ){
    firstBin  = 1;
    secondBin = 2;
  }else{
    if( binX >= a_histo->GetNbinsX()){
      secondBin = a_histo->GetNbinsX();
      firstBin  = secondBin - 1;
    }else{
      if(a_x < a_histo->GetBinCenter(binX)){
        firstBin  = binX - 1;
        secondBin = binX;
      }else{
        firstBin  = binX;
        secondBin = binX + 1;
      }
    }
  }

  double x1    = a_histo->GetBinCenter(firstBin);
  double y1    = a_histo->GetBinContent(firstBin);
  double yErr1 = a_histo->GetBinError(firstBin);
  double x2    = a_histo->GetBinCenter(secondBin);
  double y2    = a_histo->GetBinContent(secondBin);
  double yErr2 = a_histo->GetBinError(secondBin);
  PhysMath::interpolate(a_x, x1, y1, yErr1, x2, y2, yErr2, a_y, a_yErr);
  //return (by1 - by2)/(bx1 - bx2)*(a_x - bx1) + by1; // slope*(x - x1) + y1
}

//_____________________________________________________________________________________________________________________________________

void HistogramUtilities::interpolate2D(TH2D* a_histo, int a_binX, double a_y, double &a_z, double &a_zError){
  if( !a_histo || a_histo->GetNbinsY() < 2 ) return;
  int binY = a_histo->GetYaxis()->FindBin(a_y);
  int firstBin  = 0;
  int secondBin = 0;
  if(binY <= 1 ){
    firstBin  = 1;
    secondBin = 2;
  }else{
    if( binY >= a_histo->GetNbinsY()){
      secondBin = a_histo->GetNbinsY();
      firstBin  = secondBin - 1;
    }else{
      if(a_y < a_histo->GetYaxis()->GetBinCenter(binY)){
        firstBin  = binY - 1;
        secondBin = binY;
      }else{
        firstBin  = binY;
        secondBin = binY + 1;
      }
    }
  }
  double y1    = a_histo->GetYaxis()->GetBinCenter(firstBin);
  double z1    = a_histo->GetBinContent(a_binX,firstBin);
  double zErr1 = a_histo->GetBinError(a_binX,firstBin);
  double y2    = a_histo->GetYaxis()->GetBinCenter(secondBin);
  double z2    = a_histo->GetBinContent(a_binX, secondBin);
  double zErr2 = a_histo->GetBinError(a_binX,secondBin);
  PhysMath::interpolate(a_y, y1, z1, zErr1, y2, z2, zErr2, a_z, a_zError);

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "PhysMath::interpolate2D " << a_histo << " XBin " << a_binX << endl;
    cout << "      y1  z1 zErr1    = " << y1 << " " << z1 << " " << zErr1 << endl;
    cout << "      y2  z2 zErr2    = " << y2 << " " << z2 << " " << zErr2 << endl;
    cout << "      y   z  zErr     = " << a_y << " " << a_z << " " << a_zError << endl;
  #endif
}


//_____________________________________________________________________________________________________________________________________


void HistogramUtilities::getGausData(double* a_returnedData, TH1D* a_histo, bool a_isDifferential, bool a_isntCounts, double a_lowBound, double a_highBound){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::getGausData(" << a_returnedData << ", " << a_histo << ", " << a_isDifferential << ", " << a_isntCounts << ", " << a_lowBound << ", " << a_highBound << ")" << endl;
  #endif
  if(!a_histo){
    cout << "HistogramUtilities::getGausData given a NULL histogram pointer..." << endl;
    return;
  }
  double lowBound  = getLowEdge(a_histo);
  double highBound = getHighEdge(a_histo);
  if(a_lowBound != -999)  lowBound = a_lowBound;
  if(a_highBound != -999) highBound = a_highBound;

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "  Call to HistogramUtilities::getGausData with lowBound: " << lowBound << " and highBound: " << highBound << endl;
    printHistoInfo(a_histo);
  #endif

  double weights  = 0.0;
  double counts   = 0.0;
  double mean     = 0.0;
  double variance = 0.0;
  int firstBin = a_histo->FindBin(lowBound);
  int lastBin  = a_histo->FindBin(highBound);
  if(firstBin < 1) firstBin = 1;
  if(lastBin > a_histo->GetNbinsX()) lastBin = a_histo->GetNbinsX();


  double binCount = 0; // not a sum of counts , just a temp storage variable
  for(int bin = firstBin; bin <= lastBin; bin++){
    if(a_isDifferential)  binCount = a_histo->GetBinContent(bin)*a_histo->GetBinWidth(bin);
    else                  binCount = a_histo->GetBinContent(bin);
    mean   += binCount*a_histo->GetBinCenter(bin);
    counts += binCount;
  }
  if(counts <= 0){
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " HistogramUtilities::getGausData:  The histogram is empty in the region specified..." << endl;
    #endif
    a_returnedData[0] = 0;
    a_returnedData[1] = 0;
    a_returnedData[2] = 1.0;
    return;
  }
  mean = mean/counts;


  for(int bin = firstBin; bin <= lastBin; bin++){
    if(a_isDifferential)  binCount = a_histo->GetBinContent(bin)*a_histo->GetBinWidth(bin);
    else                  binCount = a_histo->GetBinContent(bin);
    variance += binCount*pow(a_histo->GetBinCenter(bin) - mean, 2.0);
  }
  if(a_isntCounts) variance = variance/counts;
  else             variance = variance/(counts-1.0);

  if(!std::isfinite(variance)){
    cout << "WARNING: HistogramUtilities::getGausData gave a non finite number for the variance with counts: " << counts << " and mean " << mean << ", so setting to 1.0" << endl;
    variance = 1.0;
  }

  a_returnedData[0] = counts;
  a_returnedData[1] = mean;
  a_returnedData[2] = sqrt(variance);

  #ifdef _HISTO_UTIL_DEBUG_
    if(!std::isfinite(counts) || !std::isfinite(mean) || !std::isfinite(a_returnedData[2]) || a_returnedData[2] <= 0.0){
      cout << " WARNING in HistogramUtilities::getGausData: Some of the returned contents are not numbers or are infinite!" << endl;
      cout << "   N: " << a_returnedData[0] << " mean:" << a_returnedData[1] << " stdDev:" << a_returnedData[2] << endl;
    }
  #endif

  return;
}

//___________________________________________________________________________________________________________________________________________________________________________


#ifdef _MAC_OSX_
vector<double> HistogramUtilities::calcGausData(TH1D* a_histo, double a_guessMean, int a_numBinsToScan, int a_binScanWidth){
  vector< double > returnData;
  if(!a_histo) return returnData;
  int guessBin = a_histo->FindBin(a_guessMean);

  double sampleCount = 0;
  vector< double > N_samples;
  vector< double > mu_samples;
  vector< double > sig_samples;

  double N_sum   = 0.0;
  double mu_sum  = 0.0;
  double sig_sum = 0.0;

  int lowBin = guessBin - a_numBinsToScan/2;
  int highBin = guessBin + a_numBinsToScan/2;
  int binOfLastScan = highBin - 2 * a_binScanWidth;
  for(int binStart = lowBin; binStart <= binOfLastScan; binStart++ ){
    int bin1 = binStart;
    int bin2 = binStart + a_binScanWidth;
    int bin3 = binStart + 2*a_binScanWidth;

    double x1 = a_histo->GetBinCenter(bin1);
    double x2 = a_histo->GetBinCenter(bin2);
    double x3 = a_histo->GetBinCenter(bin3);

    double y1 = a_histo->GetBinContent(bin1);
    double y2 = a_histo->GetBinContent(bin2);
    double y3 = a_histo->GetBinContent(bin3);


    double trial_N = 0.0;
    double trial_mu = 0.0;
    double trial_sig = 0.0;

    PhysMath::get3PointGausParams(trial_N, trial_mu, trial_sig,x1,y1,x2,y2,x3,y3);

    N_samples.push_back( trial_N);
    mu_samples.push_back( trial_mu);
    sig_samples.push_back( trial_sig);

    N_sum   += trial_N;
    mu_sum  += trial_mu;
    sig_sum += trial_sig;
    sampleCount += 1.0;
  }

  N_sum   = N_sum/sampleCount;
  mu_sum  = mu_sum/sampleCount;
  sig_sum = sig_sum/sampleCount;

  double N_var = 0.0;
  double mu_var = 0.0;
  double sig_var = 0.0;

  for(int iii = 0; iii < N_samples.size(); iii++){
    N_var += (N_sum - N_samples[iii])*(N_sum - N_samples[iii]);
    mu_var += (mu_sum - mu_samples[iii])*(mu_sum - mu_samples[iii]);
    sig_var += (sig_sum - sig_samples[iii])*(sig_sum - sig_samples[iii]);
  }

  N_var   = sqrt(N_var/(sampleCount-1.0));
  mu_var  = sqrt(mu_var/(sampleCount-1.0));
  sig_var = sqrt(sig_var/(sampleCount-1.0));


  returnData.push_back(N_sum);
  returnData.push_back(mu_sum);
  returnData.push_back(sig_sum);
  returnData.push_back(N_var);
  returnData.push_back(mu_var);
  returnData.push_back(sig_var);

  return returnData;

}
#endif

//___________________________________________________________________________________________________________________________________________________________________________


void  HistogramUtilities::addToHistogram_binSplitting(TH1D* a_histo, TH1D* a_histoToAdd){
  if(!a_histo) return;
  for(int bin = 1; bin <= a_histoToAdd->GetNbinsX(); bin++){
    addToHistogram(a_histo, a_histoToAdd->GetBinContent(bin), a_histoToAdd->GetBinError(bin), a_histoToAdd->GetBinLowEdge(bin), a_histoToAdd->GetBinLowEdge(bin+1));
  }
}

//___________________________________________________________________________________________________________________________________________________________________________


void HistogramUtilities::addToHistogram(TH1D* a_histo, double a_value, double a_value_err, double a_lowRange, double a_highRange){
  if(!a_histo) return;
  int firstBin = a_histo->FindBin(a_lowRange);
  int lastBin  = a_histo->FindBin(a_highRange);
  if(firstBin == lastBin){ //if the entire range fits inside one bin
     a_histo->SetBinContent(firstBin,a_histo->GetBinContent(firstBin) + a_value);
     a_histo->SetBinError(firstBin,sqrt(a_histo->GetBinError(firstBin)*a_histo->GetBinError(firstBin) + a_value_err*a_value_err));
     return;
  }
  double width = a_highRange - a_lowRange;
  double partialWidthFrac = (a_histo->GetBinLowEdge(firstBin+1) - a_lowRange)/width;
  double newContent       = a_histo->GetBinContent(firstBin) + partialWidthFrac*a_value;
  double newVariance      = a_histo->GetBinError(firstBin)*a_histo->GetBinError(firstBin) + partialWidthFrac*a_value_err*a_value_err;
  a_histo->SetBinContent(firstBin,newContent);
  a_histo->SetBinError(  firstBin,sqrt(newVariance));
  for(int bin = firstBin+1; bin < lastBin; bin++){
    partialWidthFrac = a_histo->GetBinWidth(bin)/width;
    newContent       = a_histo->GetBinContent(bin) + partialWidthFrac*a_value;
    newVariance      = a_histo->GetBinError(bin)*a_histo->GetBinError(bin) + partialWidthFrac*a_value_err*a_value_err;
    a_histo->SetBinContent(bin,newContent);
    a_histo->SetBinError(bin,sqrt(newVariance));
  }
  partialWidthFrac = (a_highRange - a_histo->GetBinLowEdge(lastBin))/width;
  newContent       = a_histo->GetBinContent(lastBin) + partialWidthFrac*a_value;
  newVariance      = a_histo->GetBinError(lastBin)*a_histo->GetBinError(lastBin) + partialWidthFrac*a_value_err*a_value_err;
  a_histo->SetBinContent(lastBin,newContent);
  a_histo->SetBinError(lastBin,sqrt(newVariance));
}

//___________________________________________________________________________________________________________________________________________________________________________


TH1D* HistogramUtilities::truncateHistogram(TH1D* a_histo, double a_lowCut, double a_highCut){
  if(!a_histo) return nullptr;
  int lowBin = a_histo->FindBin(a_lowCut);
  int highBin = a_histo->FindBin(a_highCut);
  if(a_highCut == a_histo->GetBinLowEdge(highBin)){
    highBin--;
  }
  if(highBin-lowBin < 0){
    cout << "HistogramUtilities::truncateHistogram with lowCut " << a_lowCut << " and high cut " << a_highCut << " wouldn't have enough bins " << endl;
    TH1D* cloneHisto = (TH1D*) a_histo->Clone();
    cloneHisto->SetName(Form("%s_cloned",a_histo->GetName()));
    return cloneHisto;
  }
  vector<double> edges;
  for(int bin = lowBin; bin <= highBin+1; bin++){
    edges.push_back(a_histo->GetBinLowEdge(bin));
  }
  TH1D* truncHisto = new TH1D(
    Form("%s_trunc",a_histo->GetName()),
    a_histo->GetTitle(),
    (int)edges.size()-1,
    &edges[0]
    );
  truncHisto->GetXaxis()->SetTitle(a_histo->GetXaxis()->GetTitle());
  truncHisto->GetYaxis()->SetTitle(a_histo->GetYaxis()->GetTitle());
  for(int bin = lowBin; bin <= highBin; bin++){
    int newBin = truncHisto->FindBin(a_histo->GetBinCenter(bin));
    truncHisto->SetBinContent(newBin, a_histo->GetBinContent(bin));
    truncHisto->SetBinError(newBin, a_histo->GetBinError(bin));
  }
  return truncHisto;
}

//___________________________________________________________________________________________________________________________________________________________________________


TH2D* HistogramUtilities::truncateHistogram2D(TH2D* a_histo, double a_lowCut_X, double a_highCut_X, double a_lowCut_Y, double a_highCut_Y){
  if(!a_histo) return nullptr;
  int lowBin_X  = a_histo->GetXaxis()->FindBin(a_lowCut_X);
  int highBin_X = a_histo->GetXaxis()->FindBin(a_highCut_X);
  int lowBin_Y  = a_histo->GetYaxis()->FindBin(a_lowCut_Y);
  int highBin_Y = a_histo->GetYaxis()->FindBin(a_highCut_Y);
  if(a_highCut_X == a_histo->GetXaxis()->GetBinLowEdge(highBin_X)){
    highBin_X--;
  }
  if(a_highCut_Y == a_histo->GetYaxis()->GetBinLowEdge(highBin_Y)){
    highBin_Y--;
  }

  bool truncate_X = (highBin_X - lowBin_X >= 0);
  bool truncate_Y = (highBin_Y - lowBin_Y >= 0);

  if(!truncate_X && !truncate_Y){
    cout << "HistogramUtilities::truncateHistogram with lowCut " << a_lowCut_X << " and high cut " << a_highCut_X << " wouldn't have enough bins " << endl;
    TH2D* cloneHisto = (TH2D*) a_histo->Clone();
    cloneHisto->SetName(Form("%s_cloned",a_histo->GetName()));
    return cloneHisto;
  }

  vector<double> edges_X;
  if(truncate_X){
    for(int bin = lowBin_X; bin <= highBin_X+1; bin++){
      edges_X.push_back(a_histo->GetXaxis()->GetBinLowEdge(bin));
    }
  }else{
    for(int bin = 1; bin <= a_histo->GetNbinsX()+1; bin++){
      edges_X.push_back(a_histo->GetXaxis()->GetBinLowEdge(bin));
    }
  }

  vector<double> edges_Y;
  if(truncate_Y){
    for(int bin = lowBin_Y; bin <= highBin_Y+1; bin++){
      edges_Y.push_back(a_histo->GetYaxis()->GetBinLowEdge(bin));
    }
  }else{
    for(int bin = 1; bin <= a_histo->GetNbinsY()+1; bin++){
      edges_Y.push_back(a_histo->GetYaxis()->GetBinLowEdge(bin));
    }
  }

  TH2D* truncHisto = new TH2D(Form("%s_trunc",a_histo->GetName()),a_histo->GetTitle(),
    (int)edges_X.size()-1,&edges_X[0],(int)edges_Y.size()-1,&edges_Y[0]
    );
  truncHisto->GetXaxis()->SetTitle(a_histo->GetXaxis()->GetTitle());
  truncHisto->GetYaxis()->SetTitle(a_histo->GetYaxis()->GetTitle());

  for(int bin_X = lowBin_X; bin_X <= highBin_X; bin_X++){
    int newBin_X = truncHisto->GetXaxis()->FindBin(a_histo->GetXaxis()->GetBinCenter(bin_X));
    for(int bin_Y = lowBin_Y; bin_Y <= highBin_Y; bin_Y++){
      int newBin_Y = truncHisto->GetYaxis()->FindBin(a_histo->GetYaxis()->GetBinCenter(bin_Y));
      truncHisto->SetBinContent(newBin_X,newBin_Y, a_histo->GetBinContent(bin_X,bin_Y));
      truncHisto->SetBinError(newBin_X,newBin_Y, a_histo->GetBinError(bin_X,bin_Y));
    }
  }
  return truncHisto;
}

//___________________________________________________________________________________________________________________________________________________________________________



TGraphErrors* HistogramUtilities::makeGraph(TH1D* a_histo, bool a_removeZeros,double a_lowRange, double a_highRange, bool a_removeZeroErrors){
  if(!a_histo) return nullptr;
  vector<double> x_vect;
  vector<double> y_vect;
  vector<double> x_vect_err;
  vector<double> y_vect_err;

  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    if(a_lowRange != -999 && a_histo->GetBinCenter(bin) < a_lowRange) continue;
    if(a_highRange != -999 && a_histo->GetBinCenter(bin) > a_highRange) continue;
    if(a_removeZeros && a_histo->GetBinContent(bin) == 0) continue;
    if(a_removeZeroErrors && a_histo->GetBinError(bin) == 0) continue;
    x_vect.push_back(a_histo->GetBinCenter(bin));
    y_vect.push_back(a_histo->GetBinContent(bin));
    x_vect_err.push_back(a_histo->GetBinWidth(bin)/2.0);
    y_vect_err.push_back(a_histo->GetBinError(bin));
  }

  TGraphErrors* graph = new TGraphErrors(x_vect.size(), &x_vect[0], &y_vect[0], &x_vect_err[0], &y_vect_err[0]);
  graph->SetName(Form("%s_graph", a_histo->GetName()));
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(a_histo->GetLineColor());
  graph->SetLineColor(a_histo->GetLineColor());

  return graph;
}


//___________________________________________________________________________________________________________________________________________________________________________

TGraphErrors* HistogramUtilities::makeGraph(TH2D* a_histo, int a_binX){
  if(!a_histo) return nullptr;
  vector<double> x_vect;
  vector<double> y_vect;
  vector<double> x_vect_err;
  vector<double> y_vect_err;

  for(int bin = 1; bin <= a_histo->GetNbinsY(); bin++){
    x_vect.push_back(a_histo->GetYaxis()->GetBinCenter(bin));
    y_vect.push_back(a_histo->GetBinContent(a_binX,bin));
    x_vect_err.push_back(a_histo->GetYaxis()->GetBinWidth(bin)/2.0);
    y_vect_err.push_back(a_histo->GetBinError(a_binX,bin));
  }

  TGraphErrors* graph = new TGraphErrors(x_vect.size(), &x_vect[0], &y_vect[0], &x_vect_err[0], &y_vect_err[0]);
  graph->SetName(Form("%s_graph", a_histo->GetName()));
  return graph;

}


//___________________________________________________________________________________________________________________________________________________________________________


TGraphAsymmErrors* HistogramUtilities::makeGraphAsymm(TH1D* a_histo, bool a_removeZeros,double a_lowRange, double a_highRange, bool a_removeZeroErrors){
  if(!a_histo) return nullptr;
  vector<double> x_vect;
  vector<double> y_vect;
  vector<double> x_vect_low_err;
  vector<double> y_vect_low_err;
  vector<double> x_vect_high_err;
  vector<double> y_vect_high_err;

  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    if(a_lowRange != -999 && a_histo->GetBinCenter(bin) < a_lowRange) continue;
    if(a_highRange != -999 && a_histo->GetBinCenter(bin) > a_highRange) continue;
    if(a_removeZeros && a_histo->GetBinContent(bin) == 0.0) continue;
    if(a_removeZeroErrors && a_histo->GetBinError(bin) == 0.0) continue;
    x_vect.push_back(a_histo->GetBinCenter(bin));
    y_vect.push_back(a_histo->GetBinContent(bin));
    x_vect_low_err.push_back(a_histo->GetBinWidth(bin)/2.0);
    x_vect_high_err.push_back(a_histo->GetBinWidth(bin)/2.0);
    y_vect_low_err.push_back(a_histo->GetBinError(bin));
    y_vect_high_err.push_back(a_histo->GetBinError(bin));
  }

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(x_vect.size(), &x_vect[0], &y_vect[0], &x_vect_low_err[0], &x_vect_high_err[0], &y_vect_low_err[0],&y_vect_high_err[0]);
  graph->SetName(Form("%s_graph", a_histo->GetName()));
  graph->SetMarkerStyle(20);
  graph->SetMarkerColor(a_histo->GetLineColor());
  graph->SetLineColor(a_histo->GetLineColor());
  return graph;
}


//___________________________________________________________________________________________________________________________________________________________________________


TH1D* HistogramUtilities::makeFitRatioHistogram(TH1D* a_histo, TF1* a_funct, bool a_doIntegral){
  if(!a_histo) return nullptr;
  TH1D* ratioHisto = (TH1D*) a_histo->Clone();
  string newName = a_histo->GetName();
  newName = newName + "_ratioToFit";
  ratioHisto->SetName(newName.c_str());
  for(int bin = 1; bin <= ratioHisto->GetNbinsX(); bin++){
    double lowEdge = a_histo->GetBinLowEdge(bin);
    double center  = a_histo->GetBinCenter(bin);
    double highEdge = a_histo->GetBinLowEdge(bin + 1);
    double content  = a_histo->GetBinContent(bin);
    double error = a_histo->GetBinError(bin);
    double fit = a_funct->Eval(center);
    if(a_doIntegral){
      fit = a_funct->Integral(lowEdge,highEdge)/(highEdge - lowEdge);
    }
    double ratio = content/fit;
    double ratioErr = error/fit;
    ratioHisto->SetBinContent(bin,ratio);
    ratioHisto->SetBinError(bin,ratioErr);
  }
  return ratioHisto;
}


//___________________________________________________________________________________________________________________________________________________________________________


TH1D* HistogramUtilities::makeFitPullHistogram(TH1D* a_histo, TF1* a_funct, bool a_doIntegral){
  if(!a_histo) return nullptr;
  TH1D* pullHisto = HistogramUtilities::makeBinStructure(a_histo);
  string newName = a_histo->GetName();
  newName = newName + "_pullToFit";
  pullHisto->SetName(newName.c_str());
  pullHisto->GetYaxis()->SetTitle("Standardized Pulls");
  for(int bin = 1; bin <= pullHisto->GetNbinsX(); bin++){
    double lowEdge  = a_histo->GetBinLowEdge(bin);
    double center   = a_histo->GetBinCenter(bin);
    double highEdge = a_histo->GetBinLowEdge(bin + 1);
    double content  = a_histo->GetBinContent(bin);
    double error    = a_histo->GetBinError(bin);
    double fit      = a_funct->Eval(center);
    if(a_doIntegral){
      fit = a_funct->Integral(lowEdge,highEdge)/(highEdge - lowEdge);
    }
    double pull    = (content - fit)/error;
    double pullErr = 1.0;
    if(std::isnan(pull) || error == 0.0) pull = -999.999;
    pullHisto->SetBinContent(bin,pull);
    pullHisto->SetBinError(bin,pullErr);
  }
  return pullHisto;
}


//___________________________________________________________________________________________________________________________________________________________________________



void HistogramUtilities::truncateTGraph(TGraph* a_graph, double a_lowCut, double a_highCut){
  if(!a_graph) return;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
             double Y_ofPoint = 0.0;
             double X_ofPoint = 0.0;
             a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);

            if(X_ofPoint > a_highCut || X_ofPoint < a_lowCut){
               a_graph->RemovePoint(pointIndex);
               pointIndex--;
            }
          }//loop over points
  return;
}

//___________________________________________________________________________________________________________________________________________________________________________

void HistogramUtilities::scaleTGraph(TGraphErrors* a_graph, double a_scale){
  if(!a_graph) return;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    double errorY = a_graph->GetErrorY(pointIndex);
    a_graph->SetPoint(pointIndex,X_ofPoint,Y_ofPoint*a_scale);
    a_graph->SetPointError(pointIndex, a_graph->GetErrorX(pointIndex), errorY*a_scale);
   }//loop over points
}

//___________________________________________________________________________________________________________________________________________________________________________

void HistogramUtilities::scaleTGraph(TGraphAsymmErrors* a_graph, double a_scale){
  if(!a_graph) return;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    double errorYlow = a_graph->GetEYlow()[pointIndex];
    double errorYhigh = a_graph->GetEYhigh()[pointIndex];
    a_graph->SetPoint(pointIndex,X_ofPoint,Y_ofPoint*a_scale);
    a_graph->SetPointError(pointIndex, a_graph->GetEXlow()[pointIndex],a_graph->GetEXhigh()[pointIndex], errorYlow*a_scale,errorYhigh*a_scale);
   }//loop over points
}


TGraphErrors* HistogramUtilities::convertAsymToNoAsym(TGraphAsymmErrors* a_graph){
  TGraphErrors* errorGraph = new TGraphErrors();
  if(!a_graph) return errorGraph;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    double errorXlow = a_graph->GetEXlow()[pointIndex];
    double errorXhigh = a_graph->GetEXhigh()[pointIndex];
    double newErrorX = errorXlow + errorXhigh;
    double errorYlow = a_graph->GetEYlow()[pointIndex];
    double errorYhigh = a_graph->GetEYhigh()[pointIndex];
    double newErrorY = errorYlow + errorYhigh;
    int pointNumNew = errorGraph->GetN();

    errorGraph->SetPoint(pointIndex, X_ofPoint, Y_ofPoint);
    errorGraph->SetPointError(pointIndex,newErrorX,newErrorY);

  }//loop over points
  return errorGraph;
}

TGraphAsymmErrors* HistogramUtilities::convertNoAsymToAsym(TGraphErrors* a_graph){
  TGraphAsymmErrors* asymmGraph = new TGraphAsymmErrors();
  if(!a_graph) return asymmGraph;
  asymmGraph->SetName(Form("%s_noAsymm",a_graph->GetName()));
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    int pointNumNew = asymmGraph->GetN();
    asymmGraph->SetPoint(pointIndex, a_graph->GetX()[pointIndex], a_graph->GetY()[pointIndex]);
    asymmGraph->SetPointError(pointIndex,a_graph->GetEX()[pointIndex],a_graph->GetEX()[pointIndex],a_graph->GetEY()[pointIndex],a_graph->GetEY()[pointIndex]);
  }//loop over points
  return asymmGraph;
}

void HistogramUtilities::shiftTGraph(TGraphErrors* a_graph, double a_shift){
  if(!a_graph) return;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
             double Y_ofPoint = 0.0;
             double X_ofPoint = 0.0;
             a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
             a_graph->SetPoint(pointIndex, X_ofPoint + a_shift ,Y_ofPoint);
             //double errorX = a_graph->GetErrorX(pointIndex);
             //double errorY = a_graph->GetErrorY(pointIndex);
             //a_graph->SetPointError(pointIndex, errorX, errorY);
   }//loop over points
}

void HistogramUtilities::zeroTGraph(TGraphErrors* a_graph){
  if(!a_graph) return;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
             double Y_ofPoint = 0.0;
             double X_ofPoint = 0.0;
             a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
             a_graph->SetPoint(pointIndex, X_ofPoint,0);
             a_graph->SetPointError(pointIndex, a_graph->GetErrorX(pointIndex), 0);
   }//loop over points
}

void HistogramUtilities::removeBelowThreshold(TGraphErrors* a_graph,double a_threshold){
  if(!a_graph){cout << "NULL Ptr passed to HistogramUtilities::removeBelowThreshold" << endl; return;}
  int N = a_graph->GetN();
  for(int pointIndex = 0; pointIndex < N; pointIndex++){
    if(a_graph->GetY()[pointIndex] < a_threshold){
      a_graph->RemovePoint(pointIndex);
      pointIndex--;
    }
  }
}
void HistogramUtilities::removeBelowThreshold(TGraphAsymmErrors* a_graph,double a_threshold){
  if(!a_graph){cout << "NULL Ptr passed to HistogramUtilities::removeBelowThreshold" << endl; return;}
  int N = a_graph->GetN();
  for(int pointIndex = 0; pointIndex < N; pointIndex++){
    if(a_graph->GetY()[pointIndex] < a_threshold){
      a_graph->RemovePoint(pointIndex);
      pointIndex--;
    }
  }
}

void HistogramUtilities::discardNanData(TGraphErrors* a_graph, bool a_coupleDataAndError){
  if(!a_graph) return;
  if(a_coupleDataAndError){
    for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
      if( !std::isfinite(a_graph->GetY()[pointIndex]) || !std::isfinite(a_graph->GetEY()[pointIndex])){
        a_graph->RemovePoint(pointIndex);
        pointIndex--;
      }
    }
  }else{
    for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
      if( !std::isfinite(a_graph->GetY()[pointIndex])){
        a_graph->SetPoint(pointIndex, a_graph->GetX()[pointIndex], 0.0);
      }
      if( !std::isfinite(a_graph->GetEY()[pointIndex])){
        a_graph->SetPointError(pointIndex, a_graph->GetEX()[pointIndex], 0.0);
      }
    }
  }
}
void HistogramUtilities::discardNanData(TGraphAsymmErrors* a_graph, bool a_coupleDataAndError){
  if(!a_graph) return;
  if(a_coupleDataAndError){
    for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
      if( !std::isfinite(a_graph->GetY()[pointIndex]) || !std::isfinite(a_graph->GetEYlow()[pointIndex])|| !std::isfinite(a_graph->GetEYhigh()[pointIndex])){
        a_graph->RemovePoint(pointIndex);
        pointIndex--;
      }
    }
  }else{
    for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
      if( !std::isfinite(a_graph->GetY()[pointIndex])){
        a_graph->SetPoint(pointIndex, a_graph->GetX()[pointIndex], 0.0);
      }
      if( !std::isfinite(a_graph->GetEYlow()[pointIndex])){
        a_graph->SetPointError(pointIndex, a_graph->GetEXlow()[pointIndex],a_graph->GetEXlow()[pointIndex], 0.0,a_graph->GetEYhigh()[pointIndex]);
      }
      if( !std::isfinite(a_graph->GetEYhigh()[pointIndex])){
        a_graph->SetPointError(pointIndex, a_graph->GetEXlow()[pointIndex],a_graph->GetEXlow()[pointIndex],a_graph->GetEYlow()[pointIndex],0.0);
      }
    }
  }
}

TGraphErrors* HistogramUtilities::addTGraphErrors(TGraphErrors* a_graph1, TGraphErrors* a_graph2, TGraphErrors* a_graph3){
  TGraphErrors* totalGraph = new TGraphErrors(0);
  totalGraph->SetName("totalGraph");
  totalGraph->SetTitle("totalGraph");
  cout << "Created Total Graph" << endl;
  if(a_graph1){
    cout << "looping over graph 1 with " << a_graph1->GetN() << " points" << endl;
    for(int pointIndex = 0; pointIndex < a_graph1->GetN(); pointIndex++){
      double Y_ofPoint = 0.0;
      double X_ofPoint = 0.0;
      double errorX = a_graph1->GetErrorX(pointIndex);
      double errorY = a_graph1->GetErrorY(pointIndex);
      a_graph1->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);

      int totPointNum = totalGraph->GetN();
      totalGraph->SetPoint(totPointNum, X_ofPoint ,Y_ofPoint);
      totalGraph->SetPointError(totPointNum, errorX, errorY);
      cout << "     TotPt: " << totPointNum << " from Gr1Pt: " << pointIndex << " x: " << X_ofPoint << endl;
    }//loop over points
  }

  if(a_graph2){
    cout << "looping over graph 2 with " << a_graph2->GetN() << " points" << endl;

    for(int pointIndex = 0; pointIndex < a_graph2->GetN(); pointIndex++){
      double Y_ofPoint = 0.0;
      double X_ofPoint = 0.0;
      a_graph2->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
      double errorX = a_graph2->GetErrorX(pointIndex);
      double errorY = a_graph2->GetErrorY(pointIndex);
      int totPointNum = totalGraph->GetN();
      totalGraph->SetPoint(totPointNum, X_ofPoint ,Y_ofPoint);
      totalGraph->SetPointError(totPointNum, errorX, errorY);
      cout << "     TotPt: " << totPointNum << " from Gr2Pt: " << pointIndex << " x: " << X_ofPoint << endl;
    }//loop over points
  }

  if(a_graph3){
    cout << "looping over graph 3 with " << a_graph3->GetN() << " points" << endl;

    for(int pointIndex = 0; pointIndex < a_graph3->GetN(); pointIndex++){
      double Y_ofPoint = 0.0;
      double X_ofPoint = 0.0;
      a_graph3->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
      double errorX = a_graph3->GetErrorX(pointIndex);
      double errorY = a_graph3->GetErrorY(pointIndex);
      int totPointNum = totalGraph->GetN();
      totalGraph->SetPoint(totPointNum, X_ofPoint ,Y_ofPoint);
      totalGraph->SetPointError(totPointNum, errorX, errorY);
      cout << "     TotPt: " << totPointNum << " from Gr3Pt: " << pointIndex << " x: " << X_ofPoint << endl;
    }//loop over points
  }

  cout << "HistogramUtilities::addTGraphErrors Check Points Total: " << totalGraph->GetN() << endl;
  return totalGraph;
}

double HistogramUtilities::findSmallestYVal(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double smallestYVal = 1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) Y_ofPoint -= a_graph->GetEY()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint < smallestYVal) smallestYVal = Y_ofPoint;
    }
  }//loop over points
  return smallestYVal;
}

double HistogramUtilities::findLargestYVal(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double largestYVal = -1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) Y_ofPoint += a_graph->GetEY()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint > largestYVal) largestYVal = Y_ofPoint;
    }
  }//loop over points
  return largestYVal;
}

double HistogramUtilities::findSmallestXVal(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double smallestXVal = 1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) X_ofPoint -= a_graph->GetEX()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint < smallestXVal) smallestXVal = X_ofPoint;
    }
  }//loop over points
  return smallestXVal;
}

double HistogramUtilities::findLargestXVal(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double largestXVal = -1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) X_ofPoint += a_graph->GetEX()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint > largestXVal) largestXVal = X_ofPoint;
    }
  }//loop over points
  return largestXVal;
}


double HistogramUtilities::findSmallestYVal(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double smallestYVal = 1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) Y_ofPoint -= a_graph->GetEYlow()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint < smallestYVal) smallestYVal = Y_ofPoint;
    }
  }//loop over points
  return smallestYVal;
}

double HistogramUtilities::findLargestYVal(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double largestYVal = -1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) Y_ofPoint += a_graph->GetEYhigh()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint > largestYVal){
        largestYVal = Y_ofPoint;
      }
    }
  }//loop over points
  return largestYVal;
}

double HistogramUtilities::findSmallestXVal(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double smallestXVal = 1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) X_ofPoint -= a_graph->GetEXlow()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint < smallestXVal) smallestXVal = X_ofPoint;
    }
  }//loop over points
  return smallestXVal;
}

double HistogramUtilities::findLargestXVal(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange, bool a_includeErrorBar){
  if(!a_graph) return -999;
  double largestXVal = -1E200;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(a_includeErrorBar) X_ofPoint += a_graph->GetEXhigh()[pointIndex];
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint > largestXVal) largestXVal = X_ofPoint;
    }
  }//loop over points
  return largestXVal;
}



int HistogramUtilities::findSmallestYValPoint(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double smallestYVal = 1E200;
  int pointOfInterest = -1;

  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint < smallestYVal){
        smallestYVal = Y_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}

int HistogramUtilities::findLargestYValPoint(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double largestYVal = -1E200;
  int pointOfInterest = -1;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint > largestYVal){
        largestYVal = Y_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}

int HistogramUtilities::findSmallestXValPoint(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double smallestXVal = 1E200;
  int pointOfInterest = -1;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint < smallestXVal){
        smallestXVal = X_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}

int HistogramUtilities::findLargestXValPoint(TGraphErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double largestXVal = -1E200;
  int pointOfInterest = -1;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint > largestXVal){
        largestXVal = X_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}



//_________________________________________________________________________________________________________________________________________________


int HistogramUtilities::findSmallestYValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double smallestYVal = 1E200;
  int pointOfInterest = -1;

  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint < smallestYVal){
        smallestYVal = Y_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}

int HistogramUtilities::findLargestYValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double largestYVal = -1E200;
  int pointOfInterest = -1;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(Y_ofPoint > largestYVal){
        largestYVal = Y_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}

int HistogramUtilities::findSmallestXValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double smallestXVal = 1E200;
  int pointOfInterest = -1;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint < smallestXVal){
        smallestXVal = X_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}

int HistogramUtilities::findLargestXValPoint(TGraphAsymmErrors* a_graph, double a_lowXRange, double a_highXRange){
  if(!a_graph) return -999;
  double largestXVal = -1E200;
  int pointOfInterest = -1;
  bool useXRange = false;
  if(a_lowXRange != a_highXRange) useXRange = true;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double Y_ofPoint = 0.0;
    double X_ofPoint = 0.0;
    a_graph->GetPoint(pointIndex,X_ofPoint,Y_ofPoint);
    if(!useXRange||(useXRange && X_ofPoint > a_lowXRange && X_ofPoint < a_highXRange)){
      if(X_ofPoint > largestXVal){
        largestXVal = X_ofPoint;
        pointOfInterest = pointIndex;
      }
    }
  }//loop over points
  return pointOfInterest;
}




void HistogramUtilities::writeTgraphErrorsToFile(TGraphErrors *a_graph, TString a_outFile){
  if(!a_graph) return;
  int nP = a_graph->GetN();
  fstream out(a_outFile.Data(),ios::out);
  double x, y, ey;
  for(int i=0; i<nP; i++){
    a_graph->GetPoint(i,x,y);
    ey = a_graph->GetErrorY(i);
    out << setw(5) << x << "\t" << setw(7) << y << "\t" << setw(7) << ey << endl;
  }
  return;
}

void HistogramUtilities::writeTgraphErrorsWithSysToFile(TGraphErrors* a_graphNom, TGraphErrors* a_graphSys, TString a_outFile){
  if(!a_graphNom || !a_graphSys) return;
  int nP = a_graphNom->GetN();

  fstream out(a_outFile.Data(),ios::out);

  double x, y, ey, eySys;

  for(int i=0; i<nP; i++){
    a_graphNom->GetPoint(i,x,y);
    ey = a_graphNom->GetErrorY(i);
    eySys = a_graphSys->GetErrorY(i);
    out << setw(5) << x << "\t" << setw(7) << y << "\t" << setw(7) << ey << "\t" << setw(7) << eySys << endl;
  }
  return;
}



void HistogramUtilities::integrateGraph(TGraphErrors* a_graph, double a_lowRange, double a_highRange, double &a_integral, double &a_integralError, bool a_divideByRange){
  if(!a_graph || a_graph->GetN() == 0) return;
  a_graph->Sort();

  double range = a_highRange - a_lowRange;

  if(range <= 0.0){
    cout << " ERROR: The range you have supplied is not ordered..." << endl;
    return;
  }

  int closeLowPointIndex = getPointClosestToVal(a_graph, a_lowRange);
  int closeHighPointIndex = getPointClosestToVal(a_graph, a_lowRange);

  if(closeLowPointIndex == 0 && a_graph->GetX()[closeLowPointIndex] > a_lowRange){
    cout << " ERROR: low range is lower than lowest point on graph" << endl;
    return;
  }
  if(closeLowPointIndex == a_graph->GetN() - 1 && a_graph->GetX()[closeHighPointIndex] < a_highRange){
    cout << " ERROR: high range is higher than highest point on graph" << endl;
    return;
  }

  int startIndex = closeLowPointIndex;
  if(a_graph->GetX()[closeLowPointIndex] > a_lowRange) startIndex = closeLowPointIndex - 1;
  int endIndex = closeHighPointIndex;
  if(a_graph->GetX()[closeHighPointIndex] < a_highRange) endIndex = closeHighPointIndex + 1;


  double lowVal = 0;
  double lowValErr = 0;
  evalWithError(a_graph,a_lowRange,lowVal,lowValErr);


  double highVal = 0;
  double highValErr = 0;
  evalWithError(a_graph,a_highRange,highVal,highValErr);

  double firstVal = 0;
  double firstValErr = 0;
  double lastVal = 0;
  double lastValErr = 0;

  if(endIndex - startIndex == 1){ // the integral is between two points
    firstVal = range*0.5*(lowVal + highVal);
    firstValErr = range*0.5*sqrt(lowValErr*lowValErr + highValErr*highValErr);
  }else{
    double deltaX = a_graph->GetX()[startIndex+1] - a_lowRange;
    firstVal = deltaX*0.5*(lowVal + a_graph->GetY()[startIndex+1]);
    firstValErr = deltaX*0.5*sqrt(lowValErr*lowValErr + a_graph->GetEY()[startIndex+1]*a_graph->GetEY()[startIndex+1]);

    deltaX = a_highRange - a_graph->GetX()[endIndex-1];
    lastVal = deltaX*0.5*(highVal + a_graph->GetY()[endIndex-1]);
    lastValErr = deltaX*0.5*sqrt(highValErr*highValErr + a_graph->GetEY()[endIndex-1]*a_graph->GetEY()[endIndex-1]);
  }
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "   First and Last evaluations: " << firstVal << " +/- " << firstValErr << "  and " << lastVal << " +/- " << lastValErr << endl;
  #endif

  a_integral = 0.0;
  vector< double > integralErrors;
  for(int pointIndex = startIndex + 1; pointIndex < endIndex - 1; pointIndex++){
    a_integral +=             0.5*(a_graph->GetX()[pointIndex + 1] - a_graph->GetX()[pointIndex]) * (a_graph->GetY()[pointIndex + 1] + a_graph->GetY()[pointIndex]);
    integralErrors.push_back( 0.5*(a_graph->GetX()[pointIndex + 1] - a_graph->GetX()[pointIndex]) * sqrt(    a_graph->GetEY()[pointIndex + 1]*a_graph->GetEY()[pointIndex + 1] + a_graph->GetEY()[pointIndex]*a_graph->GetEY()[pointIndex]) );
    #ifdef _HISTO_UTIL_DEBUG_
      cout << "  Sub-Integral: " << a_integral << " +/- " << integralErrors[integralErrors.size()-1] << endl;
    #endif
  }
  a_integral += firstVal;
  a_integral += lastVal;


  a_integralError = firstValErr*firstValErr + lastValErr*lastValErr;
  for(unsigned int iii = 0; iii < integralErrors.size(); iii++){
    a_integralError += integralErrors[iii]*integralErrors[iii];
  }
  a_integralError = sqrt(a_integralError);

  if(a_divideByRange){
    a_integral = a_integral / range;
    a_integralError = a_integralError / range;
  }

}









//returns xlow,xhigh, ylow, yhigh
//it is used for making a smaller projection, cutting white space
void HistogramUtilities::binBounds2D(TH3D* a_histo, int a_zBin, int* a_returnVars, double a_minThreshold){
  if(!a_histo) return;
  int bounds[4] = {a_histo->GetNbinsX(),1,a_histo->GetNbinsY(),1};
  for(int xBin = 1; xBin < a_histo->GetNbinsX(); xBin++){
    for(int yBin = 1; yBin < a_histo->GetNbinsY(); yBin++){
      if(a_histo->GetBinContent(xBin,yBin, a_zBin) >= a_minThreshold){
        if(xBin < bounds[0]) bounds[0] = xBin;
        if(xBin > bounds[1]) bounds[1] = xBin;
        if(yBin < bounds[2]) bounds[2] = yBin;
        if(yBin > bounds[3]) bounds[3] = yBin;
      }
    }
  }
  //std::cout << "Bounds Found: " << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << std::endl;
  a_returnVars[0] = bounds[0];
  a_returnVars[1] = bounds[1];
  a_returnVars[2] = bounds[2];
  a_returnVars[3] = bounds[3];
  //cout << "copied to input parameters" << endl;
  return;
}

//returns xlow,xhigh, ylow, yhigh
//it is used for making a smaller projection, cutting white space
void HistogramUtilities::binBounds2D(TH3I* a_histo, int a_zBin, int* a_returnVars, double a_minThreshold){
  if(!a_histo) return;
  int bounds[4] = {a_histo->GetNbinsX(),1,a_histo->GetNbinsY(),1};
  for(int xBin = 1; xBin < a_histo->GetNbinsX(); xBin++){
    for(int yBin = 1; yBin < a_histo->GetNbinsY(); yBin++){
      if(a_histo->GetBinContent(xBin,yBin, a_zBin) >= a_minThreshold){
        if(xBin < bounds[0]) bounds[0] = xBin;
        if(xBin > bounds[1]) bounds[1] = xBin;
        if(yBin < bounds[2]) bounds[2] = yBin;
        if(yBin > bounds[3]) bounds[3] = yBin;
      }
    }
  }
  //std::cout << "Bounds Found: " << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << std::endl;
  a_returnVars[0] = bounds[0];
  a_returnVars[1] = bounds[1];
  a_returnVars[2] = bounds[2];
  a_returnVars[3] = bounds[3];
  //cout << "copied to input parameters" << endl;
  return;
}


double HistogramUtilities::sumOfContent(TH1D* a_histo, double a_threshold , bool a_isDifferential, double a_lowRange, double a_highRange){
  if(!a_histo) return -999;
  int numBins = a_histo->GetNbinsX();
  int firstBin = 1;
  int lastBin = numBins;
  if(a_lowRange != -999)  firstBin    = a_histo->FindBin(a_lowRange);
  if(a_highRange != -999) lastBin     = a_histo->FindBin(a_highRange);
  if(firstBin < 1) firstBin = 1;
  if(lastBin > numBins) lastBin = numBins;
  double sum = 0;
  for(int binIndex = firstBin; binIndex <= lastBin; binIndex++){
    double newContent = a_histo->GetBinContent(binIndex);
    if(a_isDifferential) newContent = newContent * a_histo->GetBinWidth(binIndex);
    if(newContent > a_threshold)  sum += newContent;
  }
  return sum;
}


double  HistogramUtilities::sumOfContent2D(TH2D* a_histo, double a_threshold, bool a_isDifferential){
  if(!a_histo) return -999;
  double sum = 0;
  for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      double newContent = a_histo->GetBinContent(binX,binY);
      if(a_isDifferential) newContent = newContent* a_histo->GetXaxis()->GetBinWidth(binX)* a_histo->GetYaxis()->GetBinWidth(binY);
      if(newContent > a_threshold){
        sum += newContent;
      }
    }
  }
  return sum;
}


double HistogramUtilities::integral(TH1D* a_histo, double a_lowRange, double a_highRange, double a_integral_error){
  if(!a_histo) return -999;
  int numBins = a_histo->GetNbinsX();
  a_integral_error = 0.0; // stored as variance until sqrt()
  int firstBin = 1; //a_histo->FindBin(a_lowRange);
  int lastBin = numBins; //a_histo->FindBin(a_highRange);

  if(a_lowRange != -999)  firstBin    = a_histo->FindBin(a_lowRange);
  else                    a_lowRange  = a_histo->GetBinLowEdge(firstBin);
  if(a_highRange != -999) lastBin     = a_histo->FindBin(a_highRange);
  else                    a_highRange = a_histo->GetBinLowEdge(lastBin+1);

  if(firstBin < 1){
    cout << "WARNING: HistogramUtilities::integral has limits beyond histogram size... [" << a_lowRange << "," << a_highRange << "]" << endl;
    firstBin     = 1;
    a_lowRange = a_histo->GetBinLowEdge(1);
  }
  if(lastBin > numBins){
    cout << "WARNING: HistogramUtilities::integral has limits beyond histogram size... [" << a_lowRange << "," << a_highRange << "]" << endl;
    lastBin     = numBins;
    a_highRange = a_histo->GetBinLowEdge(numBins+1);
  }

  double sum = 0;
  if(firstBin == lastBin){
    sum = a_histo->GetBinContent(lastBin)*(a_highRange - a_lowRange);
    a_integral_error += pow(a_histo->GetBinError(lastBin) * (a_highRange - a_lowRange),2);
  }else{
    sum += a_histo->GetBinContent(firstBin)*(a_histo->GetBinLowEdge(firstBin+1) - a_lowRange)
             + a_histo->GetBinContent(lastBin)*(a_highRange - a_histo->GetBinLowEdge(lastBin+1));
    a_integral_error += pow(a_histo->GetBinError(firstBin) * (a_histo->GetBinLowEdge(firstBin+1) - a_lowRange),2);
    a_integral_error += pow(a_histo->GetBinError(lastBin) * (a_highRange - a_histo->GetBinLowEdge(lastBin+1)),2);

  }
  for(int binIndex = lastBin - 1; binIndex > firstBin; binIndex--){
    sum += a_histo->GetBinContent(binIndex)*a_histo->GetBinWidth(binIndex);
    a_integral_error += pow(a_histo->GetBinError(binIndex)*a_histo->GetBinWidth(binIndex),2);
  }
  a_integral_error = sqrt(a_integral_error);
  return sum;
}


double HistogramUtilities::meanOfYVals(TH1D* a_histo,double a_lowRange , double a_highRange){
  if(!a_histo){
    cout << "HistogramUtilities::getGausData given a NULL histogram pointer..." << endl;
    return -999;
  }
  double lowBound = getLowEdge(a_histo);
  double highBound = getHighEdge(a_histo);
  if(a_lowRange != -999) lowBound = a_lowRange;
  if(a_highRange != -999) highBound = a_highRange;

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "  Call to HistogramUtilities::meanOfYVals with lowBound: " << lowBound << " and highBound: " << highBound << endl;
    //printHistoInfo(a_histo);
  #endif

  double weights = 0.0;
  double mean = 0.0;
  int firstBin = a_histo->FindBin(lowBound);
  int lastBin = a_histo->FindBin(highBound);


  double binCount = 0;
  double binWeight = 0;
  for(int bin = firstBin; bin <= lastBin; bin++){
    binCount = a_histo->GetBinContent(bin);
    binWeight = 1.0/(a_histo->GetBinError(bin)*a_histo->GetBinError(bin));
    if(std::isnan(binWeight) || std::isnan(binCount)) continue;
    mean += binWeight*binCount;
    weights += binWeight;
  }
  mean = mean/weights;

  return mean;
}


double HistogramUtilities::meanOfZValsFromProj(TH2D* a_histo, int a_xBin, double a_lowYRange, double a_highYRange){
  if(!a_histo){
    cout << "ERROR: HistogramUtilities::getGausData given a NULL histogram pointer..." << endl;
    return -999;
  }
  double lowBound  = a_histo->GetYaxis()->GetBinLowEdge(1);
  double highBound = a_histo->GetYaxis()->GetBinLowEdge(a_histo->GetNbinsY()+1);
  if(a_lowYRange  != -999) lowBound  = a_lowYRange;
  if(a_highYRange != -999) highBound = a_highYRange;

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "  Call to HistogramUtilities::meanOfZValsFromProj with lowBound: " << lowBound << " and highBound: " << highBound << endl;
  #endif

  double weights = 0.0;
  double mean    = 0.0;
  int firstBin = a_histo->GetYaxis()->FindBin(lowBound);
  int lastBin  = a_histo->GetYaxis()->FindBin(highBound);

  double binCount  = 0;
  double binWeight = 0;
  for(int bin = firstBin; bin <= lastBin; bin++){
    binCount = a_histo->GetBinContent(a_xBin,bin);
    binWeight = 1.0/(a_histo->GetBinError(a_xBin,bin)*a_histo->GetBinError(a_xBin,bin));
    #ifdef _HISTO_UTIL_DEBUG_
       cout << "     bin: " << bin << " y: " << binCount << "  weight: " << binWeight << endl;
    #endif
    if(!std::isfinite(binWeight) || !std::isfinite(binCount)) continue;
    mean    += binWeight*binCount;
    weights += binWeight;
  }
  mean = mean/weights;
  return mean;



}






void HistogramUtilities::makeSqrtNErrors(TH1D* a_histo){
  if(!a_histo){cout << "ERROR: NULL Ptr passed to HistogramUtilities::makeSqrtNErrors" << endl; return;}
  int numBins = a_histo->GetNbinsX();
  for(int binX = 1; binX <= numBins; binX++){
    a_histo->SetBinError(binX,sqrt(a_histo->GetBinContent(binX)));
  }
}

void HistogramUtilities::makeSqrtNErrors2D(TH2D* a_histo){
  if(!a_histo){cout << "ERROR: NULL Ptr passed to HistogramUtilities::makeSqrtNErrors" << endl; return;}
  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  for(int binX = 1; binX <= numBinsX; binX++){
    for(int binY = 1; binY <= numBinsY; binY++){
      a_histo->SetBinError(binX,binY,sqrt(a_histo->GetBinContent(binX,binY)));
    }
  }
}

void HistogramUtilities::makeSqrtNErrors3D(TH3D* a_histo){
  if(!a_histo){cout << "ERROR: NULL Ptr passed to HistogramUtilities::makeSqrtNErrors" << endl; return;}
  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  for(int binX = 1; binX <= numBinsX; binX++){
    for(int binY = 1; binY <= numBinsY; binY++){
      for(int binZ = 1; binZ <= numBinsY; binZ++){
        a_histo->SetBinError(binX,binY,binZ,sqrt(a_histo->GetBinContent(binX,binY,binZ)));
      }
    }
  }
}

void HistogramUtilities::zeroBinsBelowThreshold(TH1D* a_histo,double a_threshold){
  if(!a_histo){cout << "ERROR: NULL Ptr passed to HistogramUtilities::zeroBinsBelowThreshold" << endl; return;}
  int numBins = a_histo->GetNbinsX();
  for(int binIndex = numBins; binIndex > 0; binIndex--){
    if(a_histo->GetBinContent(binIndex) < a_threshold){
      a_histo->SetBinContent(binIndex,0);
    }
  }
}

TH1D* HistogramUtilities::projectionY(TH2D* a_histo, double a_lowVal, double a_highVal){
  if(!a_histo){cout << "ERROR: NULL Ptr passed to HistogramUtilities::projectionY" << endl; return NULL;}
  string name = Form("%s_%s",a_histo->GetName(),"_projY");
  int numBins = a_histo->GetNbinsY();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins;iii++)
  {
    edges[iii-1] = a_histo->GetYaxis()->GetBinLowEdge(iii);
    //cout << " Bin " << iii << " Edge " << edges[iii-1] << endl;
  }
  edges[numBins] = a_histo->GetYaxis()->GetBinLowEdge(numBins)+a_histo->GetYaxis()->GetBinWidth(numBins);
  //cout << " Bin " << numBins+1 << " Edge " << edges[numBins] << endl;
  TH1D* rtrPtr = new TH1D(name.c_str(),a_histo->GetTitle(),numBins,edges);

  int lowBin = 1;
  int highBin = numBins;
  if(a_lowVal != 999) lowBin = a_histo->GetXaxis()->FindBin(a_lowVal);
  if(a_highVal != 999) highBin = a_histo->GetXaxis()->FindBin(a_highVal);
  for(int binX = lowBin; binX <= highBin ; binX++ ){
    for(int binY = 1; binY <= a_histo->GetYaxis()->GetNbins() ; binY++ ){
      rtrPtr->SetBinContent(binY,
        a_histo->GetBinContent(binX,binY) + rtrPtr->GetBinContent(binY) );
    }
  }
  rtrPtr->Sumw2();
  //delete a_histo;
  return rtrPtr;
};

TH1I* HistogramUtilities::projectionY(TH2I* a_histo, double a_lowVal, double a_highVal){
  if(!a_histo) return nullptr;
  string name = Form("%s_%s",a_histo->GetName(),"_projY");
  int numBins = a_histo->GetNbinsY();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins;iii++){
    edges[iii-1] = a_histo->GetYaxis()->GetBinLowEdge(iii);
    //cout << " Bin " << iii << " Edge " << edges[iii-1] << endl;
  }
  edges[numBins] = a_histo->GetYaxis()->GetBinLowEdge(numBins) + a_histo->GetYaxis()->GetBinWidth(numBins);
  //cout << " Bin " << numBins+1 << " Edge " << edges[numBins] << endl;
  TH1I* rtrPtr = new TH1I(name.c_str(),a_histo->GetTitle(),numBins,edges);

  int lowBin = 1;
  int highBin = numBins;
  if(a_lowVal != 999) lowBin = a_histo->GetXaxis()->FindBin(a_lowVal);
  if(a_highVal != 999) highBin = a_histo->GetXaxis()->FindBin(a_highVal);
  for(int binX = lowBin; binX <= highBin ; binX++ ){
    for(int binY = 1; binY <= a_histo->GetYaxis()->GetNbins() ; binY++ ){
      rtrPtr->SetBinContent(binY,
        a_histo->GetBinContent(binX,binY) + rtrPtr->GetBinContent(binY) );
    }
  }
  //rtrPtr->Sumw2();
  //delete a_histo;
  return rtrPtr;
};



TH1D* HistogramUtilities::projectionZ(TH3D* a_histo, int a_binX, int a_binY, string a_uniqueIdentifier){
  if(!a_histo) return nullptr;
  string name = Form("%s_%s_projZ_xBin%04d_yBin%04d",a_uniqueIdentifier.c_str(),a_histo->GetName(),a_binX,a_binY);
  //if(a_binX > 99) name = name + "_hundX";
  //if(a_binY > 99) name = name + "_hundY";
  int numBins = a_histo->GetNbinsZ();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins + 1;iii++){
    edges[iii-1] = a_histo->GetZaxis()->GetBinLowEdge(iii);
  }
  TH1D* rtrPtr = new TH1D(name.c_str(),Form("%s_%s",a_uniqueIdentifier.c_str(),a_histo->GetTitle()),numBins,edges);
  for(int binZ = 1; binZ <= numBins; binZ++ ){
    rtrPtr->SetBinContent(binZ, a_histo->GetBinContent(a_binX,a_binY,binZ));
    rtrPtr->SetBinError(binZ, a_histo->GetBinError(a_binX,a_binY,binZ));
  }
  return rtrPtr;
};


vector< vector<TH1D*> > HistogramUtilities::projectionsZ(TH3D* a_histo, bool a_makeDifferential, bool a_normalize){
  vector< vector<TH1D*> > projections;
  if(!a_histo) return projections;

  for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
    vector<TH1D*> storage_vect;
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      TH1D* histo = HistogramUtilities::projectionZ(a_histo,binX,binY);
      histo->SetTitle(Form("%s x=[%f,%f] y=[%f,%f]", a_histo->GetTitle(),
        a_histo->GetXaxis()->GetBinLowEdge(binX),
        a_histo->GetXaxis()->GetBinLowEdge(binX+1),
        a_histo->GetYaxis()->GetBinLowEdge(binY),
        a_histo->GetYaxis()->GetBinLowEdge(binY+1)));
      if(a_normalize){
        double counts = HistogramUtilities::sumOfContent(histo);
        histo->Scale(1.0/counts);
      }
      if(a_makeDifferential){
        histo->Scale(1.0,"width");
      }
      storage_vect.push_back(histo);
    }
    projections.push_back(storage_vect);
  }
  return projections;

}



  //   binLow <= bin <= binHigh
  // only even binned histograms
  // NOTE: reduction only works for single bin at the moment
  TH2D* HistogramUtilities::projection2D(TH3D* a_histo, int a_binZLow, int a_binZHigh, bool a_reduceBounds, double a_boundThreshold){
    if(!a_histo) return nullptr;
    TH2D* returnHisto = NULL;
    if(a_reduceBounds){
      int bounds[4]={0,0,0,0};
      binBounds2D(a_histo, a_binZLow, &bounds[0], a_boundThreshold);
      //cout << "called binBounds2D " << endl;
      if(bounds[1]-bounds[0]+1 <2 || bounds[3]-bounds[2]+1 <2) return NULL;
      returnHisto = new TH2D("projectionXY", "projectionXY",
                        bounds[1]-bounds[0]+1,
                        a_histo->GetXaxis()->GetBinLowEdge(bounds[0]),
                        a_histo->GetXaxis()->GetBinLowEdge(bounds[1]+1),
                        bounds[3]-bounds[2]+1,
                        a_histo->GetYaxis()->GetBinLowEdge(bounds[2]),
                        a_histo->GetYaxis()->GetBinLowEdge(bounds[3]+1)
                        );
    }else{
      vector< double > binEdgesX;
      vector< double > binEdgesY;
      int counter = 1;
      for(int bin = 1; bin <= a_histo->GetNbinsX()+1; bin++){
        binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
      }
      for(int bin = 1; bin <= a_histo->GetNbinsY()+1; bin++){
        binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
      }
      returnHisto = new TH2D("projectionXY","projectionXY", (int)(binEdgesX.size()-1), &binEdgesX[0], (int)(binEdgesY.size()-1), &binEdgesY[0] );
    }

    int numBins = a_histo->GetNbinsZ();
    for(int binZ = a_binZLow; binZ <= a_binZHigh; binZ++){
      //cout << "On BinZ=" << binZ << endl;
      for(int binX = 1; binX <= a_histo->GetXaxis()->GetNbins(); binX++){
        for(int binY = 1; binY <= a_histo->GetYaxis()->GetNbins(); binY++){
          double content = returnHisto->GetBinContent(binX,binY);
          double newContent = a_histo->GetBinContent(binX,binY,binZ);
          returnHisto->SetBinContent(binX,binY, content + newContent);
        }
      }
    }
    return returnHisto;
  }

  TH2I* HistogramUtilities::projection2D(TH3I* a_histo, int a_binZLow, int a_binZHigh, bool a_reduceBounds, double a_boundThreshold){
    if(!a_histo) return nullptr;
    TH2I* returnHisto = NULL;
    if(a_reduceBounds){
      int bounds[4]={0,0,0,0};
      binBounds2D(a_histo, a_binZLow, &bounds[0], a_boundThreshold);
      //cout << "called binBounds2D " << endl;
      if(bounds[1]-bounds[0]+1 <2 || bounds[3]-bounds[2]+1 <2) return NULL;


      returnHisto = new TH2I("projectionXY", "projectionXY",
                        bounds[1]-bounds[0]+1,
                        a_histo->GetXaxis()->GetBinLowEdge(bounds[0]),
                        a_histo->GetXaxis()->GetBinLowEdge(bounds[1]+1),
                        bounds[3]-bounds[2]+1,
                        a_histo->GetYaxis()->GetBinLowEdge(bounds[2]),
                        a_histo->GetYaxis()->GetBinLowEdge(bounds[3]+1)
                        );

    }else{

      vector< double > binEdgesX;
      vector< double > binEdgesY;

      int counter = 1;
      for(int bin = 1; bin <= a_histo->GetNbinsX()+1; bin++){
        binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
      }
      for(int bin = 1; bin <= a_histo->GetNbinsY()+1; bin++){
        binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
      }

      returnHisto = new TH2I("projectionXY","projectionXY", (int)(binEdgesX.size()-1), &binEdgesX[0], (int)(binEdgesY.size()-1), &binEdgesY[0] );

    }


    //cout << "made the histogram " << returnHisto->GetNbinsX() << " " << returnHisto->GetNbinsY() << endl;

    int numBins = a_histo->GetNbinsZ();
    for(int binZ = a_binZLow; binZ <= a_binZHigh; binZ++){
      //cout << "On BinZ=" << binZ << endl;
      for(int binX = 1; binX <= a_histo->GetXaxis()->GetNbins(); binX++){
        for(int binY = 1; binY <= a_histo->GetYaxis()->GetNbins(); binY++){
          int content = returnHisto->GetBinContent(binX,binY);
          int newContent = a_histo->GetBinContent(binX,binY,binZ);
          returnHisto->SetBinContent(binX,binY, content + newContent);
        }
      }
    }

    return returnHisto;
  }


TH2D* HistogramUtilities::projection2D_XY(TH3D* a_histo, int a_binZ){
  if(!a_histo) return nullptr;
  vector< double > binEdgesX;
  vector< double > binEdgesY;

  for(int bin = 1; bin <= a_histo->GetNbinsX()+1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  for(int bin = 1; bin <= a_histo->GetNbinsY()+1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
  }


  string newName = a_histo->GetName();
  newName = newName + Form("_ZBin%02d",a_binZ);

  string newTitle = a_histo->GetTitle();
  newTitle = newTitle + Form(" z_var=[%2.3f,%2.3f]",
      a_histo->GetZaxis()->GetBinLowEdge(a_binZ),
      a_histo->GetZaxis()->GetBinLowEdge(a_binZ+1));
  TH2D* returnHisto = new TH2D(newName.c_str(),newTitle.c_str(), (int)(binEdgesX.size()-1), &binEdgesX[0], (int)(binEdgesY.size()-1), &binEdgesY[0] );

  for(int binX = 1; binX <= a_histo->GetXaxis()->GetNbins(); binX++){
    for(int binY = 1; binY <= a_histo->GetYaxis()->GetNbins(); binY++){
      returnHisto->SetBinContent(binX,binY, a_histo->GetBinContent(binX,binY,a_binZ));
      returnHisto->SetBinError(binX,binY, a_histo->GetBinError(binX,binY,a_binZ));
    }
  }
  return returnHisto;
}



TH2D* HistogramUtilities::projection2D_XZ(TH3D* a_histo, int a_binY){
  if(!a_histo) return nullptr;
  vector< double > binEdgesX;
  vector< double > binEdgesZ;

  for(int bin = 1; bin <= a_histo->GetNbinsX()+1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  for(int bin = 1; bin <= a_histo->GetNbinsZ()+1; bin++){
    binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
  }


  string newName = a_histo->GetName();
  newName = newName + Form("_YBin%02d",a_binY);

  string newTitle = a_histo->GetTitle();
  newTitle = newTitle + Form(" y_var=[%2.3f,%2.3f]",
      a_histo->GetYaxis()->GetBinLowEdge(a_binY),
      a_histo->GetYaxis()->GetBinLowEdge(a_binY+1));
  TH2D* returnHisto = new TH2D(newName.c_str(),newTitle.c_str(), (int)(binEdgesX.size()-1), &binEdgesX[0], (int)(binEdgesZ.size()-1), &binEdgesZ[0] );

  for(int binX = 1; binX <= a_histo->GetXaxis()->GetNbins(); binX++){
    for(int binZ = 1; binZ <= a_histo->GetZaxis()->GetNbins(); binZ++){
      returnHisto->SetBinContent(binX,binZ, a_histo->GetBinContent(binX,a_binY,binZ));
      returnHisto->SetBinError(binX,binZ, a_histo->GetBinError(binX,a_binY,binZ));
    }
  }


  return returnHisto;

}






TH2D* HistogramUtilities::projection2D_YZ(TH3D* a_histo, int a_binX){
  if(!a_histo) return nullptr;
  vector< double > binEdgesY;
  vector< double > binEdgesZ;

  for(int bin = 1; bin <= a_histo->GetNbinsY()+1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
  }
  for(int bin = 1; bin <= a_histo->GetNbinsZ()+1; bin++){
    binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
  }


  string newName = a_histo->GetName();
  newName = newName + Form("_XBin%02d",a_binX);

  string newTitle = a_histo->GetTitle();
  newTitle = newTitle + Form(" x_var=[%2.3f,%2.3f]",
      a_histo->GetXaxis()->GetBinLowEdge(a_binX),
      a_histo->GetXaxis()->GetBinLowEdge(a_binX+1));
  TH2D* returnHisto = new TH2D(newName.c_str(),newTitle.c_str(), (int)(binEdgesY.size()-1), &binEdgesY[0], (int)(binEdgesZ.size()-1), &binEdgesZ[0] );

  for(int binY = 1; binY <= a_histo->GetYaxis()->GetNbins(); binY++){
    for(int binZ = 1; binZ <= a_histo->GetZaxis()->GetNbins(); binZ++){
      returnHisto->SetBinContent(binY,binZ, a_histo->GetBinContent(a_binX,binY,binZ));
      returnHisto->SetBinError(binY,binZ, a_histo->GetBinError(a_binX,binY,binZ));
    }
  }


    return returnHisto;

}



vector< TH2I* > HistogramUtilities::projeciton2DSlices_X(TH3I* a_histo){ // steps through x producing all 2D slices
  vector<TH2I* > return_histo_vector;
  if(!a_histo) return return_histo_vector;
  string original_name = a_histo->GetName();
  string original_title = a_histo->GetTitle();

  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  int numBinsZ = a_histo->GetNbinsZ();

  vector< double > binEdgesY;
  vector< double > binEdgesZ;
  for(int bin = 1; bin <= numBinsY+1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
  }
  for(int bin = 1; bin <= numBinsZ+1; bin++){
    binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
  }

  TH2I* struct_histo = new TH2I("base_projectionYZ","base_projectionYZ", (int)(binEdgesY.size()-1), &binEdgesY[0], (int)(binEdgesZ.size()-1), &binEdgesZ[0] );
  zeroBins2D(struct_histo);


  for(int binX = 1; binX <= numBinsX; binX++){
    string newName = Form("%s_XSlice_Bin_%d",original_name.c_str(),binX);
    string newTitle = Form("%s x=[%1.2f,%1.2f]",original_name.c_str(),a_histo->GetXaxis()->GetBinLowEdge(binX), a_histo->GetXaxis()->GetBinLowEdge(binX+1) );
    TH2I* cloneHisto = (TH2I*) struct_histo->Clone();
    cloneHisto->SetName(newName.c_str());
    cloneHisto->SetTitle(newTitle.c_str());

    int numBins = a_histo->GetNbinsZ();

    for(int binY = 0; binY <= numBinsY + 1; binY++){
      for(int binZ = 0; binZ <= numBinsZ + 1; binZ++){
        cloneHisto->SetBinContent(binY, binZ, a_histo->GetBinContent(binX,binY,binZ));
        cloneHisto->SetBinError(binY, binZ, a_histo->GetBinError(binX,binY,binZ));
      }
    }
    return_histo_vector.push_back(cloneHisto);
  }
  #ifndef _SAVE_POINTERS_
    delete struct_histo;
  #endif
  return return_histo_vector;

}

vector< TH2D* > HistogramUtilities::projeciton2DSlices_X(TH3D* a_histo){ // steps through x producing all 2D slices
  vector<TH2D* > return_histo_vector;
  if(!a_histo) return return_histo_vector;
  string original_name = a_histo->GetName();
  string original_title = a_histo->GetTitle();

  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  int numBinsZ = a_histo->GetNbinsZ();

  vector< double > binEdgesY;
  vector< double > binEdgesZ;
  for(int bin = 1; bin <= numBinsY+1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
  }
  for(int bin = 1; bin <= numBinsZ+1; bin++){
    binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
  }

  TH2D* struct_histo = new TH2D("base_projectionYZ","base_projectionYZ", (int)(binEdgesY.size()-1), &binEdgesY[0], (int)(binEdgesZ.size()-1), &binEdgesZ[0] );
  zeroBins2D(struct_histo);

  for(int binX = 1; binX <= numBinsX; binX++){
    string newName = Form("%s_XSlice_Bin_%d",original_name.c_str(),binX);
    string newTitle = Form("%s x=[%1.2f,%1.2f]",original_name.c_str(),a_histo->GetXaxis()->GetBinLowEdge(binX), a_histo->GetXaxis()->GetBinLowEdge(binX+1) );
    TH2D* cloneHisto = (TH2D*) struct_histo->Clone();
    cloneHisto->SetName(newName.c_str());
    cloneHisto->SetTitle(newTitle.c_str());

    int numBins = a_histo->GetNbinsZ();

    for(int binY = 0; binY <= numBinsY + 1; binY++){
      for(int binZ = 0; binZ <= numBinsZ + 1; binZ++){
        cloneHisto->SetBinContent(binY, binZ, a_histo->GetBinContent(binX,binY,binZ));
        cloneHisto->SetBinError(binY, binZ,   a_histo->GetBinError(binX,binY,binZ));
      }
    }
    return_histo_vector.push_back(cloneHisto);
  }
  #ifndef _SAVE_POINTERS_
    delete struct_histo;
  #endif
  return return_histo_vector;

}

/*
vector< TH2I* > HistogramUtilities::projeciton2DSlices_Z(TH3I* a_histo, vector<double> a_zValueCuts){
  vector<TH2I* > return_histo_vector;
  if(!a_histo) return return_histo_vector;
  string original_name = a_histo->GetName();
  string original_title = a_histo->GetTitle();
  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  int numBinsZ = a_histo->GetNbinsZ();
  vector< double > binEdgesX;
  for(int bin = 1; bin <= numBinsX+1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  vector< double > binEdgesY;
  for(int bin = 1; bin <= numBinsY+1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
  }

  TH2I* struct_histo = new TH2I("base_projectionXY","base_projectionXY", (int)(binEdgesX.size()-1), &binEdgesX[0], (int)(binEdgesY.size()-1), &binEdgesY[0] );
  zeroBins2D(struct_histo);
  int numNewZProjections = a_zValueCuts.size();

  for(int zBinIndex = 0; zBinIndex < numNewZProjections; zBinIndex++){
    int    lowZBin  = a_histo->GetZaxis()->FindBin(a_zValueCuts[zBinIndex]);
    double lowZEdge = a_histo->GetZaxis()->GetBinLowEdge(lowZBin);
    int    highZBin  = 0;
    double highZEdge = 0.0;
    if(zBinIndex == 0){
      highZBin  = numBinsZ + 1;
    }else{
      highZBin = a_histo->GetZaxis()->FindBin(a_zValueCuts[zBinIndex - 1]);
    }
    highZEdge = a_histo->GetZaxis()->GetBinLowEdge(highZBin);

    string newName = Form("%s_ZSlice_Bin_%d",original_name.c_str(),zBinIndex);
    string newTitle = Form("%s z=[%1.2f,%1.2f]",original_name.c_str(),lowZEdge, highZEdge );
    TH2I* cloneHisto = (TH2I*) struct_histo->Clone();
    cloneHisto->SetName(newName.c_str());
    cloneHisto->SetTitle(newTitle.c_str());
    for(int binX = 0; binX <= numBinsX + 1; binX++){
      for(int binY = 0; binY <= numBinsY + 1; binY++){
        for(int binZ = lowZBin; binZ < highZBin; binZ++){
          double originalContent = cloneHisto->GetBinContent(binX, binY);
          cloneHisto->SetBinContent(binX, binY, originalContent + a_histo->GetBinContent(binX,binY,binZ));
        }
      }
    }
    return_histo_vector.push_back(cloneHisto);
  }
  #ifndef _SAVE_POINTERS_
    delete struct_histo;
  #endif
  return return_histo_vector;

}*/

vector< TH2I* > HistogramUtilities::projeciton2DSlices_Z(TH3I* a_histo, vector<double> a_zValueCuts){
  vector<TH2I* > return_histo_vector;
  if(!a_histo) return return_histo_vector;
  string original_name = a_histo->GetName();
  string original_title = a_histo->GetTitle();
  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  int numBinsZ = a_histo->GetNbinsZ();
  vector< double > binEdgesX;
  for(int bin = 1; bin <= numBinsX+1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  vector< double > binEdgesY;
  for(int bin = 1; bin <= numBinsY+1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
  }

  TH2I* struct_histo = new TH2I("base_projectionXY","base_projectionXY", (int)(binEdgesX.size()-1), &binEdgesX[0], (int)(binEdgesY.size()-1), &binEdgesY[0] );
  zeroBins2D(struct_histo);
  int numNewZProjections = a_zValueCuts.size() - 1; // only do the cuts specified, no under or overflow

  for(int zBinIndex = 0; zBinIndex < numNewZProjections; zBinIndex++){
    double value_low  = a_zValueCuts[zBinIndex];
    double value_high = a_zValueCuts[zBinIndex+1];
    if(value_low > value_high){
      value_low  = a_zValueCuts[zBinIndex+1];
      value_high = a_zValueCuts[zBinIndex];
    }
    int    lowZBin  = a_histo->GetZaxis()->FindBin(value_low);
    if(lowZBin <= 0) lowZBin = 1;
    double lowZEdge = a_histo->GetZaxis()->GetBinLowEdge(lowZBin);
    int    highZBin  = a_histo->GetZaxis()->FindBin(value_high);
    if(highZBin > numBinsZ) highZBin = numBinsZ;
    double highZEdge = a_histo->GetZaxis()->GetBinLowEdge(highZBin + 1);

    string newName = Form("%s_ZSlice_Bin_%d",original_name.c_str(),zBinIndex);
    string newTitle = Form("%s z=[%1.2f,%1.2f]",original_name.c_str(),lowZEdge, highZEdge );
    TH2I* cloneHisto = (TH2I*) struct_histo->Clone();
    cloneHisto->SetName(newName.c_str());
    cloneHisto->SetTitle(newTitle.c_str());
    for(int binX = 0; binX <= numBinsX + 1; binX++){
      for(int binY = 0; binY <= numBinsY + 1; binY++){
        for(int binZ = lowZBin; binZ <= highZBin; binZ++){
          double originalContent = cloneHisto->GetBinContent(binX, binY);
          cloneHisto->SetBinContent(binX, binY, originalContent + a_histo->GetBinContent(binX,binY,binZ));
        }
      }
    }
    return_histo_vector.push_back(cloneHisto);
  }
  #ifndef _SAVE_POINTERS_
    delete struct_histo;
  #endif
  return return_histo_vector;

}



vector< TH2D* > HistogramUtilities::projeciton2DSlices_Z(TH3D* a_histo, vector<double> a_zValueCuts){
  vector<TH2D* > return_histo_vector;
  if(!a_histo) return return_histo_vector;
  string original_name = a_histo->GetName();
  string original_title = a_histo->GetTitle();
  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  int numBinsZ = a_histo->GetNbinsZ();
  vector< double > binEdgesX;
  for(int bin = 1; bin <= numBinsX+1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  vector< double > binEdgesY;
  for(int bin = 1; bin <= numBinsY+1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
  }

  TH2D* struct_histo = new TH2D("base_projectionXY","base_projectionXY", (int)(binEdgesX.size()-1), &binEdgesX[0], (int)(binEdgesY.size()-1), &binEdgesY[0] );
  zeroBins2D(struct_histo);
  int numNewZProjections = a_zValueCuts.size() - 1; // only do the cuts specified, no under or overflow

  for(int zBinIndex = 0; zBinIndex < numNewZProjections; zBinIndex++){
    double value_low  = a_zValueCuts[zBinIndex];
    double value_high = a_zValueCuts[zBinIndex+1];
    if(value_low > value_high){
      value_low  = a_zValueCuts[zBinIndex+1];
      value_high = a_zValueCuts[zBinIndex];
    }
    int    lowZBin  = a_histo->GetZaxis()->FindBin(value_low);
    if(lowZBin <= 0) lowZBin = 1;
    double lowZEdge = a_histo->GetZaxis()->GetBinLowEdge(lowZBin);
    int    highZBin  = a_histo->GetZaxis()->FindBin(value_high);
    if(highZBin > numBinsZ) highZBin = numBinsZ;
    double highZEdge = a_histo->GetZaxis()->GetBinLowEdge(highZBin + 1);

    string newName = Form("%s_ZSlice_Bin_%d",original_name.c_str(),zBinIndex);
    string newTitle = Form("%s z=[%1.2f,%1.2f]",original_name.c_str(),lowZEdge, highZEdge );
    TH2D* cloneHisto = (TH2D*) struct_histo->Clone();
    cloneHisto->SetName(newName.c_str());
    cloneHisto->SetTitle(newTitle.c_str());
    for(int binX = 0; binX <= numBinsX + 1; binX++){
      for(int binY = 0; binY <= numBinsY + 1; binY++){
        for(int binZ = lowZBin; binZ <= highZBin; binZ++){
          double originalContent = cloneHisto->GetBinContent(binX, binY);
          cloneHisto->SetBinContent(binX, binY, originalContent + a_histo->GetBinContent(binX,binY,binZ));
        }
      }
    }
    return_histo_vector.push_back(cloneHisto);
  }
  #ifndef _SAVE_POINTERS_
    delete struct_histo;
  #endif
  return return_histo_vector;

}




  TH2D* HistogramUtilities::crop2DHisto( TH2D* a_histoToCrop,double a_lowX, double a_highX, double a_lowY, double a_highY){
    if(!a_histoToCrop) return nullptr;
    int binLowX = a_histoToCrop->GetXaxis()->FindBin(a_lowX);
    int binHighX = a_histoToCrop->GetXaxis()->FindBin(a_highX);
    int binLowY = a_histoToCrop->GetYaxis()->FindBin(a_lowY);
    int binHighY = a_histoToCrop->GetYaxis()->FindBin(a_highY);

    TH2D* croppedHisto = NULL;
    if(binHighX > binLowX && binHighY > binLowY){
      croppedHisto = new TH2D("cropped2D","cropped2D",binHighX - binLowX + 1,
          a_histoToCrop->GetXaxis()->GetBinLowEdge(binLowX),
          a_histoToCrop->GetXaxis()->GetBinLowEdge(binHighX+1),
          binHighY - binLowY + 1,
          a_histoToCrop->GetYaxis()->GetBinLowEdge(binLowY),
          a_histoToCrop->GetYaxis()->GetBinLowEdge(binHighY+1)
          );

      for(int binX = binLowX; binX <= binHighX; binX++){
        for(int binY = binLowY; binY <= binHighY; binY++){
           double content = a_histoToCrop->GetBinContent(binX,binY);
           croppedHisto->SetBinContent(binX-binLowX+1,binY-binLowY+1,content);
        }
      }

    }else{
      cerr <<"\n\n\nERROR: HistogramUtilities::crop2DHisto()  select different bounds to crop\n\n\n";
    }

    return croppedHisto;
  }


  ///calibrates the histogram passed to the function using a linear calibration of the form a*X +b
  void HistogramUtilities::applyCalibration(TH1*   a_histogram,///pointer to the histogram to apply the calibration to
                        double a_a,///<linear factor in the calibration
                        double a_b///<constant to change the x intercept
                        ){
    if(!a_histogram) return;
    TAxis* axis_1 = a_histogram->GetXaxis();
    axis_1->SetLimits(axis_1->GetXmin()*a_a+a_b,axis_1->GetXmax()*a_a+a_b);
  }

  TH1D* HistogramUtilities::shiftHisto(TH1D* a_histogram, double a_shift){
    if(!a_histogram) return nullptr;
    double numBins = a_histogram->GetNbinsX();
    double lowEnd = a_histogram->GetBinLowEdge(1);
    double highEnd = a_histogram->GetBinLowEdge(numBins+1);
    TH1D* histo = new TH1D("shifted","shifted",numBins,lowEnd+a_shift,highEnd+a_shift);
    for(int iii = 1; iii <= numBins; iii++){
      histo->SetBinContent(iii, a_histogram->GetBinContent(iii));
    }
    return histo;
  }


TH2D* HistogramUtilities::makeScaleAndStretchHisto(TH2D* a_histo, double a_scale, double a_stretch){
  if(!a_histo) return nullptr;
  string name = Form("%s_%s",a_histo->GetName(),"scaledAndStretched");
  int numBinsX = a_histo->GetNbinsX();
  int numBinsY = a_histo->GetNbinsY();
  double* edgesX = new double[numBinsX+1];
  double* edgesY = new double[numBinsY+1];
  for(int binX = 1;binX <= numBinsX + 1; binX++){
    edgesX[binX-1] = a_histo->GetXaxis()->GetBinLowEdge(binX);
  }
  for(int binY = 1;binY <= numBinsY + 1; binY++){
    edgesY[binY-1] = a_histo->GetYaxis()->GetBinLowEdge(binY)*a_stretch;
  }

  TH2D* rtrPtr = new TH2D(name.c_str(),Form("%s_%s",a_histo->GetTitle()," Scaled and Stretched"),numBinsX,edgesX,numBinsY,edgesY);
  for(int binX = 1; binX <= numBinsX ; binX++ ){
    for(int binY = 1; binY <= numBinsY ; binY++ ){
      rtrPtr->SetBinContent(binX,binY,a_histo->GetBinContent(binX,binY)*a_scale);
      rtrPtr->SetBinError(binX,binY,a_histo->GetBinError(binX,binY)*a_scale);
    }
  }
  //delete a_histo;
  return rtrPtr;

}



void HistogramUtilities::chi_squared(TH1D* a_histo_1, TH1D* a_histo_2,double &a_chi_sqr, double &a_ndf){

  int numBins = a_histo_1->GetNbinsX();
  if(numBins != a_histo_2->GetNbinsX()){
    cout << "ERROR: HistogramUtilities::chi_squared : you have two histograms with two different binnings..." << endl;
    return;
  }

  a_ndf = 0.0;
  a_chi_sqr = 0.0;
  double content_1 = 0.0;
  double content_2 = 0.0;
  double error_1   = 0.0;
  double error_2   = 0.0;
  for(int bin = 1; bin <= numBins; bin++){
    content_1 = a_histo_1->GetBinContent(bin);
    content_2 = a_histo_2->GetBinContent(bin);
    error_1 = a_histo_1->GetBinError(bin);
    error_2 = a_histo_2->GetBinError(bin);
    if(content_1 > 0.0 && content_2 > 0.0 && error_1 > 0.0 && error_2 > 0.0){
      a_ndf += 1.0;
      a_chi_sqr += pow(content_1 - content_2,2)/(pow(error_1,2) + pow(error_2,2));
    }




  }
}








TH1D* HistogramUtilities::logifyHisto(TH1D* a_histo, double a_minContent){
  if(!a_histo) return nullptr;
  string name = Form("%s_%s",a_histo->GetName(),"_logified");
  int numBins = a_histo->GetNbinsX();

  TH1D* returnHisto = (TH1D*) a_histo->Clone();
  returnHisto->SetName(name.c_str());
  returnHisto->GetYaxis()->SetTitle(Form("ln(%s)",returnHisto->GetYaxis()->GetTitle()));
  zeroBins(returnHisto);

  for(int bin = 1; bin <= numBins; bin++ ){
    double content = a_histo->GetBinContent(bin);
    double error = a_histo->GetBinError(bin);
    if(content < a_minContent){
      returnHisto->SetBinContent(bin,0.0);
      returnHisto->SetBinError(bin,0.0);
    }else{
      double newContent = log(content);
      double newError = 0.0;
      if(error <= 0.0){
        newError = 0.0;
      }else{
        newError = 1.0/error;
      }
      returnHisto->SetBinContent(bin,newContent);
      returnHisto->SetBinError(bin,newError);
    }

  }
  return returnHisto;
}

TH1D* HistogramUtilities::logifyHistoHorizontalAxis(TH1D* a_histo){
  if(!a_histo) return nullptr;
  string name = Form("%s_%s",a_histo->GetName(),"_logified");
  int numBins = a_histo->GetNbinsX();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins;iii++){
    edges[iii-1] = log(a_histo->GetBinLowEdge(iii));
  }
  edges[numBins] = log(a_histo->GetBinLowEdge(numBins)+a_histo->GetBinWidth(numBins));
  TH1D* rtrPtr = new TH1D(name.c_str(),a_histo->GetTitle(),numBins,edges);
  for(int iii = 1; iii <= numBins ; iii++ ){
    rtrPtr->SetBinContent(iii,a_histo->GetBinContent(iii));
    rtrPtr->SetBinError(iii,a_histo->GetBinError(iii));
  }
  //delete a_hist;
  string newAxisTitle = "ln(";
  newAxisTitle += a_histo->GetXaxis()->GetTitle();
  newAxisTitle += ")";
  rtrPtr->GetXaxis()->SetTitle(newAxisTitle.c_str());
  rtrPtr->GetYaxis()->SetTitle(a_histo->GetYaxis()->GetTitle());
  return rtrPtr;
}






TH1* HistogramUtilities::makeBinStructure(TH1* a_hist){
  if(!a_hist) return nullptr;
  string name = Form("%s_%s",a_hist->GetName(),"withBinStructure");
  int numBins = a_hist->GetNbinsX();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins;iii++)
  {
    edges[iii-1] = a_hist->GetBinLowEdge(iii);
    //cout << " Bin " << iii << " Edge " << edges[iii-1] << endl;
  }
  edges[numBins] = a_hist->GetBinLowEdge(numBins)+a_hist->GetBinWidth(numBins);
  //cout << " Bin " << numBins+1 << " Edge " << edges[numBins] << endl;
  TH1* rtrPtr = new TH1D(name.c_str(),a_hist->GetTitle(),numBins,edges);
  for(int iii = 1; iii <= numBins ; iii++ )
  {
    rtrPtr->SetBinContent(iii,a_hist->GetBinContent(iii));
    rtrPtr->SetBinError(iii,a_hist->GetBinError(iii));

  }
  //delete a_hist;
  return rtrPtr;
}

TH1D* HistogramUtilities::makeBinStructure(TH1D* a_hist){
  if(!a_hist) return nullptr;
  string name = Form("%s_%s",a_hist->GetName(),"withBinStructure");
  int numBins = a_hist->GetNbinsX();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins;iii++){
    edges[iii-1] = a_hist->GetBinLowEdge(iii);
    //cout << " Bin " << iii << " Edge " << edges[iii-1] << endl;
  }
  edges[numBins] = a_hist->GetBinLowEdge(numBins)+a_hist->GetBinWidth(numBins);
  //cout << " Bin " << numBins+1 << " Edge " << edges[numBins] << endl;
  TH1D* rtrPtr = new TH1D(name.c_str(),a_hist->GetTitle(),numBins,edges);
  for(int iii = 1; iii <= numBins ; iii++ ){
    rtrPtr->SetBinContent(iii,a_hist->GetBinContent(iii));
    rtrPtr->SetBinError(iii,a_hist->GetBinError(iii));
  }
  //delete a_hist;
  rtrPtr->GetXaxis()->SetTitle(a_hist->GetXaxis()->GetTitle());
  rtrPtr->GetYaxis()->SetTitle(a_hist->GetYaxis()->GetTitle());
  return rtrPtr;
}


vector< double > HistogramUtilities::getBinEdges(TH1* a_histo){
  vector<double> edges;
  if(!a_histo) return edges;
  for(int bin = 1; bin <= a_histo->GetNbinsX() + 1; bin++){
    edges.push_back(a_histo->GetBinLowEdge(bin));
  }
  return edges;
}


vector< double > HistogramUtilities::getBinEdges_X(TH2* a_histo){
  vector<double> edges;
  if(!a_histo) return edges;
  for(int bin = 1; bin <= a_histo->GetNbinsX() + 1; bin++){
    edges.push_back(a_histo->GetXaxis()->GetBinLowEdge(bin));
  }
  return edges;
}

vector< double > HistogramUtilities::getBinEdges_Y(TH2* a_histo){
  vector<double> edges;
  if(!a_histo) return edges;
  for(int bin = 1; bin <= a_histo->GetNbinsY() + 1; bin++){
    edges.push_back(a_histo->GetYaxis()->GetBinLowEdge(bin));
  }
  return edges;
}

vector< double > HistogramUtilities::getBinEdges_X(TH3* a_histo){
  vector<double> edges;
  if(!a_histo) return edges;
  for(int bin = 1; bin <= a_histo->GetNbinsX() + 1; bin++){
    edges.push_back(a_histo->GetXaxis()->GetBinLowEdge(bin));
  }
  return edges;
}
vector< double > HistogramUtilities::getBinEdges_Y(TH3* a_histo){
  vector<double> edges;
  if(!a_histo) return edges;
  for(int bin = 1; bin <= a_histo->GetNbinsY() + 1; bin++){
    edges.push_back(a_histo->GetYaxis()->GetBinLowEdge(bin));
  }
  return edges;
}
vector< double > HistogramUtilities::getBinEdges_Z(TH3* a_histo){
  vector<double> edges;
  if(!a_histo) return edges;
  for(int bin = 1; bin <= a_histo->GetNbinsZ() + 1; bin++){
    edges.push_back(a_histo->GetZaxis()->GetBinLowEdge(bin));
  }
  return edges;
}


vector< double > HistogramUtilities::getSmallerBinEdges(TH1* a_histo, int a_divNum){
  vector<double> edges;
  if(!a_histo) return edges;
  edges = getBinEdges(a_histo);
  for(int index = 0; index < edges.size() - 1; index++){
    double width = edges[index+1] - edges[index];
    double smallWidth = (edges[index+1] - edges[index])/((double)(a_divNum));
    for(int iii = 1; iii < a_divNum; iii++){
      edges.push_back(edges[index] + ((double) iii)*smallWidth);
    }
  }
  std::sort(edges.begin(), edges.end());
  return edges;
}


vector< double > HistogramUtilities::getBinCenters(TH1* a_histo){
  vector<double> centers;
  if(!a_histo) return centers;
  centers.push_back(a_histo->GetBinLowEdge(1) - a_histo->GetBinWidth(1));
  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    centers.push_back(a_histo->GetBinCenter(bin));
  }
  centers.push_back(a_histo->GetBinLowEdge(a_histo->GetNbinsX() + 1) + a_histo->GetBinWidth(a_histo->GetNbinsX()));
  return centers;
}

vector< double > HistogramUtilities::getBinContents(TH1D* a_histo){
  vector<double> contents;
  if(!a_histo) return contents;
  for(int bin = 0; bin <= a_histo->GetNbinsX()+1; bin++){
    contents.push_back(a_histo->GetBinContent(bin));
  }
  return contents;
}

vector< double > HistogramUtilities::getBinErrors(TH1D* a_histo){
  vector<double> errors;
  if(!a_histo) return errors;
  for(int bin = 0; bin <= a_histo->GetNbinsX()+1; bin++){
    errors.push_back(a_histo->GetBinError(bin));
  }
  return errors;
}


vector<double> HistogramUtilities::makeEvenBinEdges(int a_numBins, double a_lowBound, double a_highBound){
  double width = (a_highBound - a_lowBound)/((double) a_numBins);
  vector<double> edges;
  for(double step = 0.0; step < ((double) a_numBins) + 1; step++ ){
    edges.push_back(a_lowBound + step*width);
  }
  return edges;
}


vector<double> HistogramUtilities::makeLogBinEdges(int a_numBins, double a_lowBound, double a_highBound){
  vector<double> edges;
  for(int binIndex = 0; binIndex <= a_numBins; binIndex++){
     edges.push_back(pow(10, log10(a_lowBound) + log10(a_highBound/a_lowBound)*((double) binIndex)/((double) a_numBins)));
  }
  return edges;
}

double HistogramUtilities::normalize(TH1D* a_histo, bool a_isDifferential, double a_lowRange, double a_highRange){
  if(!a_histo) return -999;
  double normalizationFactor = 1.0;
  int lowBin = 1;
  int highBin = a_histo->GetNbinsX();
  if(a_lowRange != -999){
    lowBin = a_histo->FindBin(a_lowRange);
  }
  if(lowBin < 1) lowBin = 1;
  if(a_highRange == -999){
    highBin = a_histo->FindBin(a_highRange);
  }
  if(highBin > a_histo->GetNbinsX()) highBin = a_histo->GetNbinsX();
  if(a_isDifferential){
    normalizationFactor = 1.0/a_histo->Integral(lowBin,highBin,"width");
  }else{
    normalizationFactor = 1.0/a_histo->Integral(lowBin,highBin);
  }

  if(normalizationFactor == 0.0){
    cout << "WARNING: HistogramUtilities::normalize came up with a scale factor of 0. Not doing anything" << endl;
    return 0.0;
  }
  a_histo->Scale(normalizationFactor);
  return normalizationFactor;
}


void HistogramUtilities::normalizeEachVerticalSlice(TH2D* a_histo, bool a_isDifferential){
  if(!a_histo) return;
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
}






void HistogramUtilities::makeDifferential(TH1D* a_histo, double a_extraScale){
  if(!a_histo) return;
  for(int iii = 1; iii <= a_histo->GetNbinsX(); iii++){
    a_histo->SetBinContent(iii, a_histo->GetBinContent(iii)*a_extraScale/a_histo->GetBinWidth(iii));
    a_histo->SetBinError(iii, a_histo->GetBinError(iii)*a_extraScale/a_histo->GetBinWidth(iii));
  }
}


void HistogramUtilities::makeDifferential(TH2D* a_histo, double a_extraScale){
  if(!a_histo) return;
  for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      double factor = a_extraScale/(a_histo->GetXaxis()->GetBinWidth(binX) * a_histo->GetYaxis()->GetBinWidth(binY));
      a_histo->SetBinContent(binX,binY, a_histo->GetBinContent(binX,binY)*factor);
      a_histo->SetBinError(binX,binY, a_histo->GetBinError(binX,binY)*factor);
    }
  }
}

void HistogramUtilities::makeNotDifferential(TH1D* a_histo, double a_extraScale){
  if(!a_histo) return;
  for(int iii = 1; iii <= a_histo->GetNbinsX(); iii++){
    a_histo->SetBinContent(iii, a_histo->GetBinContent(iii)*a_extraScale*a_histo->GetBinWidth(iii));
    a_histo->SetBinError(iii, a_histo->GetBinError(iii)*a_extraScale*a_histo->GetBinWidth(iii));
  }
}


void HistogramUtilities::makeNotDifferential(TH2D* a_histo, double a_extraScale){
  if(!a_histo) return;
  for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
    for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
      double factor = a_extraScale* a_histo->GetXaxis()->GetBinWidth(binX) * a_histo->GetYaxis()->GetBinWidth(binY);
      a_histo->SetBinContent(binX,binY, a_histo->GetBinContent(binX,binY)*factor);
      a_histo->SetBinError(binX,binY, a_histo->GetBinError(binX,binY)*factor);
    }
  }
}


TH1D* HistogramUtilities::convertFromLogBins(TH1D* a_hist){
  if(!a_hist) return nullptr;
  string name = Form("%s_%s",a_hist->GetName(),"nowLog");
  int numBins = a_hist->GetNbinsX();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins;iii++){
    edges[iii-1] = log(a_hist->GetBinLowEdge(iii));
    //cout << " Bin " << iii << " Edge " << edges[iii-1] << endl;
  }
  edges[numBins] = log(a_hist->GetBinLowEdge(numBins) + a_hist->GetBinWidth(numBins));
  //cout << " Bin " << numBins+1 << " Edge " << edges[numBins] << endl;
  TH1D* rtrPtr = new TH1D(name.c_str(),a_hist->GetTitle(),numBins,edges);
  for(int iii = 1; iii <= numBins ; iii++ ){
    rtrPtr->SetBinContent(iii,a_hist->GetBinContent(iii));
  }
  //delete a_hist;
  return rtrPtr;
}

TH1D* HistogramUtilities::convertInt_to_DoubleHisto(TH1I* a_histo){
  if(!a_histo) return nullptr;
  string name = Form("%s_%s",a_histo->GetName(),"convToDouble");
  int numBins = a_histo->GetNbinsX();
  double* edges = new double[numBins+1];
  for(int iii = 1;iii <= numBins;iii++){
    edges[iii-1] = a_histo->GetBinLowEdge(iii);
  }
  edges[numBins] = a_histo->GetBinLowEdge(numBins) + a_histo->GetBinWidth(numBins);
  TH1D* rtrPtr = new TH1D(name.c_str(),a_histo->GetTitle(),numBins,edges);
  for(int iii = 1; iii <= numBins ; iii++ ){
    rtrPtr->SetBinContent(iii,(double) a_histo->GetBinContent(iii));
    rtrPtr->SetBinError(iii,sqrt( (double) a_histo->GetBinContent(iii) ));
  }
  return rtrPtr;
}

TH2D* HistogramUtilities::convertInt_to_DoubleHisto2D(TH2I* a_histo){
  if(!a_histo) return nullptr;
  string name = Form("%s_%s",a_histo->GetName(),"convToDouble");
  int numBinsX = a_histo->GetNbinsX();
  double* edgesX = new double[numBinsX+1];
  int numBinsY = a_histo->GetNbinsY();
  double* edgesY = new double[numBinsY+1];
  for(int binX = 1; binX <= numBinsX; binX++){
    edgesX[binX-1] = a_histo->GetXaxis()->GetBinLowEdge(binX);
  }
  edgesX[numBinsX] = a_histo->GetXaxis()->GetBinLowEdge(numBinsX) + a_histo->GetXaxis()->GetBinWidth(numBinsX);

  for(int binY = 1; binY <= numBinsY; binY++){
    edgesY[binY-1] = a_histo->GetYaxis()->GetBinLowEdge(binY);
  }
  edgesY[numBinsY] = a_histo->GetYaxis()->GetBinLowEdge(numBinsY) + a_histo->GetYaxis()->GetBinWidth(numBinsY);

  TH2D* rtrPtr = new TH2D(name.c_str(),a_histo->GetTitle(),numBinsX,edgesX,numBinsY,edgesY);
  for(int binX = 1; binX <= numBinsX ; binX++ ){
    for(int binY = 1; binY <= numBinsY; binY++){
      rtrPtr->SetBinContent(binX,binY,(double) a_histo->GetBinContent(binX,binY));
      rtrPtr->SetBinError(binX,binY,sqrt( (double) a_histo->GetBinContent(binX,binY) ));
    }
  }
  return rtrPtr;
}

TH3D* HistogramUtilities::convertInt_to_DoubleHisto3D(TH3I* a_histo){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " In HistogramUtilities::convertInt_to_DoubleHisto3D(" << a_histo << ")" << endl;
  #endif
  if(!a_histo){
    cout << "ERROR: You have handed HistogramUtilities::convertInt_to_DoubleHisto3D a NULL pointer! = " << a_histo << endl;
    return NULL;
  }

  string name = Form("%s_%s",a_histo->GetName(),"convToDouble");
  int numBinsX = a_histo->GetNbinsX();
  double* edgesX = new double[numBinsX+1];
  int numBinsY = a_histo->GetNbinsY();
  double* edgesY = new double[numBinsY+1];
  int numBinsZ = a_histo->GetNbinsZ();
  double* edgesZ = new double[numBinsZ+1];

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "   Number of Bins " << numBinsX << "  " << numBinsY << "  " << numBinsZ << endl;
  #endif

  for(int binX = 1; binX <= numBinsX; binX++){
    edgesX[binX-1] = a_histo->GetXaxis()->GetBinLowEdge(binX);
  }
  edgesX[numBinsX] = a_histo->GetXaxis()->GetBinLowEdge(numBinsX) + a_histo->GetXaxis()->GetBinWidth(numBinsX);

  for(int binY = 1; binY <= numBinsY; binY++){
    edgesY[binY-1] = a_histo->GetYaxis()->GetBinLowEdge(binY);
  }
  edgesY[numBinsY] = a_histo->GetYaxis()->GetBinLowEdge(numBinsY) + a_histo->GetYaxis()->GetBinWidth(numBinsY);

  for(int binZ = 1; binZ <= numBinsZ; binZ++){
    edgesZ[binZ-1] = a_histo->GetZaxis()->GetBinLowEdge(binZ);
  }
  edgesZ[numBinsZ] = a_histo->GetZaxis()->GetBinLowEdge(numBinsZ) + a_histo->GetZaxis()->GetBinWidth(numBinsZ);

  TH3D* rtrPtr = new TH3D(name.c_str(),a_histo->GetTitle(),numBinsX,edgesX,numBinsY,edgesY,numBinsZ,edgesZ);
  for(int binX = 1; binX <= numBinsX ; binX++){
    for(int binY = 1; binY <= numBinsY; binY++){
      for(int binZ = 1; binZ <= numBinsZ; binZ++){
        rtrPtr->SetBinContent(binX,binY,binZ,(double) a_histo->GetBinContent(binX,binY,binZ));
        rtrPtr->SetBinError(binX,binY,binZ, sqrt( (double) a_histo->GetBinContent(binX,binY,binZ) ));
      }
    }
  }
  return rtrPtr;
}


//must already be in order
TH1D* HistogramUtilities::stitchHistos(int a_numHistos, TH1D** a_histos){
  if(a_numHistos == 0 || !a_histos[0]) return nullptr;
  string ending = "_stitched";
  string name = a_histos[0]->GetName();
  name = name + ending;
  string title = a_histos[0]->GetTitle();

  //std::vector<int> ordering(a_numHistos,0);
  std::vector<double> binEdges;
  std::vector<double> binCounts;
  std::vector<double> binErrors;
  int binCounter = -1;
  int numBins = a_histos[0]->GetNbinsX();
  for(int binIndex = 1; binIndex <= numBins+1; binIndex++){
    binEdges.push_back(a_histos[0]->GetXaxis()->GetBinLowEdge(binIndex));
    if(binIndex != numBins+1){
      binCounts.push_back(a_histos[0]->GetBinContent(binIndex));
      binErrors.push_back(a_histos[0]->GetBinError(binIndex));
    }
    binCounter++;
  }
  for(int histoIndex = 1; histoIndex < a_numHistos; histoIndex++){
    //cout << "On Histo" << histoIndex<< endl;
    //see if overlaping
    double lastBinEdge = binEdges.at(binEdges.size()-1);
    //cout << "Last Bin Edge" << lastBinEdge << endl;
    if(a_histos[histoIndex]->GetXaxis()->GetBinLowEdge(1) < lastBinEdge) {
      //its overlapping, find bin to split
      //cout << "overlapping" << endl;
      int binToSplit = -1;
      numBins = a_histos[histoIndex]->GetNbinsX();
      for(int binIndex = 1; binIndex <= numBins+1; binIndex++){
        //cout << "comparing " << a_histos[histoIndex]->GetXaxis()->GetBinLowEdge(binIndex) << " to "
        //      << lastBinEdge << endl;
        if(a_histos[histoIndex]->GetXaxis()->GetBinLowEdge(binIndex) > lastBinEdge){
          binToSplit = binIndex -1;
          break;
        }
      }
      if(binToSplit < 0){
        cout << "ERROR IN HistogramUtilities::stitchHistos()  : The Histograms supplied are not in correct ordering..." << endl;
        cout << "last Bin edge: " << lastBinEdge << endl;
        return NULL;
      }

      //deal with the split bin
      double binWidth = a_histos[histoIndex]->GetBinWidth(binToSplit);
      double rightFraction = (a_histos[histoIndex]->GetBinLowEdge(binToSplit+1) - lastBinEdge)/binWidth;
      double leftFraction = (lastBinEdge - a_histos[histoIndex]->GetBinLowEdge(binToSplit))/binWidth;
      binCounts.push_back(rightFraction * a_histos[histoIndex]->GetBinContent(binToSplit));
      binErrors.push_back(TMath::Sqrt(rightFraction*a_histos[histoIndex]->GetBinContent(binToSplit)));


      //fill the rest of the bins
      for(int binIndex = binToSplit + 1; binIndex <= numBins+1; binIndex++){
        binEdges.push_back(a_histos[histoIndex]->GetXaxis()->GetBinLowEdge( binIndex ));
        if(binIndex != numBins+1){
          binCounts.push_back(a_histos[histoIndex]->GetBinContent(binIndex));
          binErrors.push_back(a_histos[histoIndex]->GetBinError(binIndex));
        }

      }

    }
    else{
      //cout << "Bins not overlapping" << endl;
      numBins = a_histos[histoIndex]->GetNbinsX();
      //deal with the empty bin
      binCounts.push_back(0);

      //fill the rest of the bins
      for(int binIndex = 1; binIndex <= numBins+1; binIndex++){
        binEdges.push_back(a_histos[histoIndex]->GetXaxis()->GetBinLowEdge(binIndex));
        if(binIndex != numBins+1){
          binCounts.push_back(a_histos[histoIndex]->GetBinContent(binIndex));
          binErrors.push_back(a_histos[histoIndex]->GetBinError(binIndex));
        }

      }
    }//end non overlap else
  }//end for loop over histograms

  //for(int iii = 0; iii < binEdges.size(); iii++) cout << " " << binEdges[iii];
  //cout << endl;
  //for(int iii = 0; iii < binCounts.size(); iii++) cout << " " << binCounts[iii];
  //cout << endl;
  //create the stitched histogram now knowing bin edges

  TH1D* rtrPtr = new TH1D(name.c_str(),title.c_str(),binEdges.size()-1,&binEdges[0]);
  for(int binIndex = 1; binIndex <= rtrPtr->GetNbinsX(); binIndex++ )
  {
    //cout << " " << binEdges[binIndex-1];
    rtrPtr->SetBinContent(binIndex, binCounts[binIndex-1]);
    rtrPtr->SetBinError(binIndex,binErrors[binIndex-1]);
  }
  //cout << "       Created Stitched Histo with " << rtrPtr->GetNbinsX() << " bins" << endl;
  return rtrPtr;



}//end stitch histos


vector<double> HistogramUtilities::makeAdaptiveBinEdgesEqualCounts(TH1D* a_histo, double a_minCounts, bool a_roundToMillionsPlace){

  std::vector<double> binEdges;
  if(!a_histo) return binEdges;
  int firstBin = firstBinWithContent(a_histo);
  int lastBin = lastBinWithContent(a_histo);
  int numBins = lastBin - firstBin + 1;
  if(numBins < 1){
    cout << "ERROR: HistogramUtilities::makeAdaptiveBinEdgesEqualCounts : NumBins is 0."
         << " Check to make sure your Histogram actually has counts " << endl;
  }

  double currentCount = 0.0;
  binEdges.push_back(a_histo->GetBinLowEdge(firstBin));
  int currentBin = firstBin;
  double currentLowEdge = a_histo->GetBinLowEdge(firstBin);
  while(currentBin <= lastBin){
    double remainingCountsInBin = (a_histo->GetBinLowEdge(currentBin+1) - currentLowEdge)*a_histo->GetBinContent(currentBin)/a_histo->GetBinWidth(currentBin);
    cout << " Bin: " << currentBin << " currentCount: " << currentCount << " currentLowEdge: " << currentLowEdge << "  remainingCountsInBin: " << remainingCountsInBin;
    if(currentCount + remainingCountsInBin < a_minCounts){
      currentCount += remainingCountsInBin;
      currentBin++;
      currentLowEdge = a_histo->GetBinLowEdge(currentBin);
      cout << " Not Enough " << endl;
    }else{
      double newWidth = (a_minCounts - currentCount)*a_histo->GetBinWidth(currentBin)/a_histo->GetBinContent(currentBin);
      currentLowEdge = currentLowEdge + newWidth;
      binEdges.push_back(currentLowEdge);
      currentCount = 0.0;
      cout << " Enough with new Width:  " << newWidth << endl;
    }
  }

  binEdges.push_back(a_histo->GetBinLowEdge(lastBin + 1));

  if(a_roundToMillionsPlace){
    for(int iii = 0; iii < binEdges.size(); iii++){
      binEdges[iii] = (double)(((int) (binEdges[iii] * 1000000)))/1000000.0;
    }
  }

  return binEdges;




}




TH1D* HistogramUtilities::rebinMinCount(TH1D* a_histo, double a_minCount){
  if(!a_histo) return nullptr;
  std::vector<double> binEdges;
  std::vector<double> binCounts;

  int firstBin = firstBinWithContent(a_histo);
  int lastBin = lastBinWithContent(a_histo);
  int numBins = lastBin - firstBin + 1;
  if(numBins < 1){
    cout << "ERROR: HistogramUtilities::rebinMinCount : NumBins is 0."
         << " Check to make sure your Histogram actually has counts " << endl;
  }

  //cout << "   num bins = " << numBins << endl;
  double currentCount = a_histo->GetBinContent(firstBin);
  binEdges.push_back(a_histo->GetBinLowEdge(firstBin));

  for(int binIndex = firstBin+1; binIndex <= lastBin; binIndex++){

    if(currentCount < a_minCount){
      currentCount += a_histo->GetBinContent(binIndex);

    }else{
      binCounts.push_back(currentCount);
      binEdges.push_back(a_histo->GetBinLowEdge(binIndex));
      currentCount = a_histo->GetBinContent(binIndex);
    }
  }
  binCounts.push_back(currentCount);
  binEdges.push_back(a_histo->GetBinLowEdge(lastBin+1));
  //for(int iii = 0; iii < binEdges.size(); iii++){
  // cout << " " << binEdges[iii];
  //}
  //binEdges.push_back(a_histo->GetBinLowEdge(numBins+1));
  string ending = "_reBinned";
  string name = a_histo->GetName();
  name = name + ending;
  string title = a_histo->GetTitle();
  TH1D* rtrPtr = new TH1D(name.c_str(),title.c_str(),binEdges.size()-1,&binEdges[0]);
  if(!rtrPtr){
    cout << "HistogramUtilities::RebinMinCount ERROR couldn't create histogram" << endl;
  }
  numBins = rtrPtr->GetNbinsX();
  for(int binIndex = 1; binIndex <= numBins; binIndex++){
    rtrPtr->SetBinContent(binIndex,binCounts[binIndex-1]);
  }
  rtrPtr->Sumw2();

  return rtrPtr;

}

TH1D* HistogramUtilities::rebinEqualCounts(TH1D* a_histo, int a_numBins){
  if(!a_histo) return nullptr;
  double counts = sumOfContent(a_histo);
  double equalBinCounts = counts/((double) a_numBins);
  return rebinMinCount(a_histo,equalBinCounts);


}



TH1D* HistogramUtilities::rebinEqualByVerticalChange(TH1D* a_histo, int a_numBins){
  if(!a_histo) return nullptr;
  TH1D* distHisto = (TH1D*) a_histo->Clone();
  distHisto->SetBinContent(0,1);
  double distanceCovered = 0;
  for(int bin = 2; bin <= a_histo->GetNbinsX(); bin++){
    distanceCovered += fabs(a_histo->GetBinContent(bin) - a_histo->GetBinContent(bin-1));
    distHisto->SetBinContent(bin,distanceCovered);
  }

  double distanceWidth = distanceCovered/((double) a_numBins);

  std::vector<double> binEdges;
  std::vector<double> binCounts;
  //int numBins = a_histo->GetNbinsX();

  int firstBin = firstBinWithContent(a_histo);
  int lastBin = lastBinWithContent(a_histo);
  int numBins = lastBin - firstBin + 1;
  if(numBins < 1){
    cout << "ERROR: HistogramUtilities::rebinMinCount : NumBins is 0."
         << " Check to make sure your Histogram actually has counts " << endl;
  }

  //cout << "   num bins = " << numBins << endl;
  double currentTravel = 0.0;
  double currentCount = a_histo->GetBinContent(firstBin);
  binEdges.push_back(a_histo->GetBinLowEdge(firstBin));

  for(int binIndex = firstBin+1; binIndex <= lastBin; binIndex++){
    currentTravel += fabs(a_histo->GetBinContent(binIndex) - a_histo->GetBinContent(binIndex - 1));
    currentCount += a_histo->GetBinContent(binIndex);
    if(currentTravel > distanceWidth){
      binCounts.push_back(currentCount);
      binEdges.push_back(a_histo->GetBinLowEdge(binIndex+1));
      currentTravel = 0.0;
      currentCount = 0.0;
    }
  }

  binEdges.push_back(a_histo->GetBinLowEdge(lastBin+1));
  binCounts.push_back(currentCount);

  string ending = "_reBinned";
  string name = a_histo->GetName();
  name = name + ending;
  string title = a_histo->GetTitle();
  TH1D* rtrPtr = new TH1D(name.c_str(),title.c_str(),binEdges.size()-1,&binEdges[0]);
  if(!rtrPtr){
    cout << "HistogramUtilities::RebinMinCount ERROR couldn't create histogram" << endl;
  }
  numBins = rtrPtr->GetNbinsX();
  for(int binIndex = 1; binIndex <= numBins; binIndex++){
    rtrPtr->SetBinContent(binIndex,binCounts[binIndex-1]);
  }
  rtrPtr->Sumw2();


  return rtrPtr;
}





vector< double > HistogramUtilities::getEdgesByVerticalChange(TH1D* a_histo, int a_numBins, bool a_weightByVariance, double a_lowRange,double a_highRange){
  std::vector<double> binEdges;
  if(!a_histo) return binEdges;
  double distanceCovered = 0;
  for(int bin = 2; bin <= a_histo->GetNbinsX(); bin++){
    if(a_weightByVariance) distanceCovered += fabs((a_histo->GetBinContent(bin) - a_histo->GetBinContent(bin-1))/(a_histo->GetBinError(bin)*a_histo->GetBinError(bin) + a_histo->GetBinContent(bin-1)*a_histo->GetBinContent(bin-1)));
    else distanceCovered += fabs(a_histo->GetBinContent(bin) - a_histo->GetBinContent(bin-1));
  }

  double distanceWidth = distanceCovered/((double) a_numBins);

  int firstBin = firstBinWithContent(a_histo);
  if(a_lowRange != -999) firstBin = a_histo->FindBin(a_lowRange);
  int lastBin = lastBinWithContent(a_histo);
  if(a_highRange != -999) lastBin = a_histo->FindBin(a_highRange);


  int numBins = lastBin - firstBin + 1;
  if(numBins < 1){
    cout << "ERROR: HistogramUtilities::rebinMinCount : NumBins is 0."
         << " Check to make sure your Histogram actually has counts " << endl;
  }

  //cout << "   num bins = " << numBins << endl;
  double currentTravel = 0.0;
  binEdges.push_back(a_histo->GetBinLowEdge(firstBin));

  for(int bin = firstBin+1; bin <= lastBin; bin++){
    if(a_weightByVariance) currentTravel += fabs((a_histo->GetBinContent(bin) - a_histo->GetBinContent(bin-1))/(a_histo->GetBinError(bin)*a_histo->GetBinError(bin) + a_histo->GetBinContent(bin-1)*a_histo->GetBinContent(bin-1)));
    else                   currentTravel += fabs(a_histo->GetBinContent(bin) - a_histo->GetBinContent(bin - 1));
    if(currentTravel > distanceWidth){
      binEdges.push_back(a_histo->GetBinLowEdge(bin+1));
      currentTravel = currentTravel - distanceWidth;
    }
  }

  if(binEdges[binEdges.size()-1] != a_histo->GetBinLowEdge(lastBin+1)) binEdges.push_back(a_histo->GetBinLowEdge(lastBin+1));
  else cout << "  WARNING:  HistogramUtilities::getEdgesByVerticalChange is changing the number of points by one, was " << a_numBins + 1 << " now " << binEdges.size() << " you might not have fine enough sampling for that many" << endl;

  return binEdges;
}

vector< double > HistogramUtilities::getEdgesByVerticalChange(TGraphAsymmErrors* a_graph, int a_numBins, bool a_weightByVariance, double a_lowRange,double a_highRange){
  double distanceCovered = 0;
  for(int pointIndex = 1; pointIndex < a_graph->GetN(); pointIndex++){
    double dist = a_graph->GetY()[pointIndex] - a_graph->GetY()[pointIndex-1];
    if(a_weightByVariance){
      double variance = 0.0;
      if(dist > 0.0){
        variance = a_graph->GetEYlow()[pointIndex]*a_graph->GetEYlow()[pointIndex] + a_graph->GetEYhigh()[pointIndex-1]*a_graph->GetEYhigh()[pointIndex-1];
      }else{
        variance = a_graph->GetEYhigh()[pointIndex]*a_graph->GetEYhigh()[pointIndex] + a_graph->GetEYlow()[pointIndex-1]*a_graph->GetEYlow()[pointIndex-1];
      }
      distanceCovered += fabs(dist/variance);
    }
    else distanceCovered += fabs(dist);
  }

  double distanceWidth = distanceCovered/((double) a_numBins);

  std::vector<double> binEdges;
  int lowPoint = 0;
  if(a_lowRange != -999) lowPoint = getPointClosestToVal(a_graph, a_lowRange);
  int highPoint = a_graph->GetN()-1;
  if(a_highRange != -999) highPoint = getPointClosestToVal(a_graph, a_highRange);

  double currentTravel = 0.0;
  binEdges.push_back(a_graph->GetX()[0]);

  for(int pointIndex = lowPoint + 1; pointIndex <= highPoint; pointIndex++){
    double dist = a_graph->GetY()[pointIndex] - a_graph->GetY()[pointIndex-1];
    if(a_weightByVariance){
      double variance = 0.0;
      if(dist > 0.0){
        variance = a_graph->GetEYlow()[pointIndex]*a_graph->GetEYlow()[pointIndex] + a_graph->GetEYhigh()[pointIndex-1]*a_graph->GetEYhigh()[pointIndex-1];
      }else{
        variance = a_graph->GetEYhigh()[pointIndex]*a_graph->GetEYhigh()[pointIndex] + a_graph->GetEYlow()[pointIndex-1]*a_graph->GetEYlow()[pointIndex-1];
      }
      currentTravel += fabs(dist/variance);
    }
    else currentTravel += fabs(dist);

    if(currentTravel > distanceWidth){
      binEdges.push_back(a_graph->GetX()[pointIndex]);
      currentTravel = currentTravel - distanceWidth;
    }
  }
  if(binEdges[binEdges.size()-1] != a_graph->GetX()[highPoint]) binEdges.push_back(a_graph->GetX()[highPoint]);
  else cout << "  WARNING:  HistogramUtilities::getEdgesByVerticalChange is changing the number of points by one, was " << a_numBins + 1 << " now " << binEdges.size() << " you might not have fine enough sampling for that many" << endl;
  return binEdges;
}


TH1D* HistogramUtilities::rebinNoZeroBins(TH1D* a_histo, bool a_isDifferential){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Call to HistogramUtilities::rebinNoZeroBins( " << a_histo << " ) " << endl;
  #endif
  if(!a_histo) return nullptr;

  int firstBin = firstBinWithContent(a_histo);
  int lastBin  = lastBinWithContent(a_histo);

  vector<double> binEdges;

  binEdges.push_back(a_histo->GetBinLowEdge(firstBin));
  for(int bin = firstBin; bin < lastBin; bin++){
    int nextBinWithContent = lastBin;
    for(int bin2 = bin+1; bin2 <= lastBin; bin2++){
      if(a_histo->GetBinContent(bin2) > 0.0){
        nextBinWithContent = bin2;
        break;
      }
    }
    if(nextBinWithContent == bin + 1){
      binEdges.push_back(a_histo->GetBinLowEdge(nextBinWithContent));
    }else{
      double newEdge = 0.5*(a_histo->GetBinLowEdge(bin+1) + a_histo->GetBinLowEdge(nextBinWithContent));
      binEdges.push_back(newEdge);
      bin = nextBinWithContent - 1;
    }
  }
  binEdges.push_back(a_histo->GetBinLowEdge(lastBin+1));

  string newName = a_histo->GetName();
  newName += "_noZeros";
  string newTitle = a_histo->GetTitle();
  newTitle += Form(";%s;%s",a_histo->GetXaxis()->GetTitle(),a_histo->GetYaxis()->GetTitle());
  TH1D* rebinHisto = new TH1D(newName.c_str(),newTitle.c_str(),((int) binEdges.size())-1,&binEdges[0]);
  for(int bin = firstBin; bin <= lastBin; bin++){
    double content = a_histo->GetBinContent(bin);
    if(content > 0.0){
      double error   = a_histo->GetBinError(bin);
      if(a_isDifferential){
        content = content*a_histo->GetBinWidth(bin);
        error   = error*a_histo->GetBinWidth(bin);
      }
      int newBin = rebinHisto->FindBin(a_histo->GetBinCenter(bin));
      rebinHisto->SetBinContent(newBin,content);
      rebinHisto->SetBinError(newBin,error);
    }
  }

  if(a_isDifferential){
    makeDifferential(rebinHisto);
  }

  return rebinHisto;
}




TH1D* HistogramUtilities::rebinEqualByCDF(TH1D* a_histo, int a_numBins, bool a_isDifferential, bool a_useLnWeight, bool a_makeNoZeroBins){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Call to HistogramUtilities::rebinEqualByCDF( " << a_histo << ", " << a_numBins << " , " << a_isDifferential << " , " << a_useLnWeight  << " ) " << endl;
  #endif
  if(!a_histo) return nullptr;
  discardNanHistoData(a_histo);
  TH1D* adaptedHisto = a_histo;
  if(a_makeNoZeroBins){
    adaptedHisto = rebinNoZeroBins(a_histo,a_isDifferential);
    if(a_numBins >= adaptedHisto->GetNbinsX()) return adaptedHisto;
  }
  double lnWeightCounts = 0.0;
  double noLnWeightCounts = 0.0;
  double totalCounts = 0.0;
  for(int binIndex = adaptedHisto->GetNbinsX(); binIndex > 0; binIndex--){
    if(adaptedHisto->GetBinContent(binIndex) > 0.0){
      if(a_isDifferential) lnWeightCounts += log(adaptedHisto->GetBinContent(binIndex) * adaptedHisto->GetBinWidth(binIndex) + 1.0);
      else                 lnWeightCounts += log(adaptedHisto->GetBinContent(binIndex) + 1.0);
    }
  }
  noLnWeightCounts = sumOfContent(adaptedHisto,-1,a_isDifferential);

  if(a_useLnWeight) totalCounts = lnWeightCounts;
  else              totalCounts = noLnWeightCounts;

  if(noLnWeightCounts <= 0.0){
    cout << "ERROR: HistogramUtilities::rebinEqualByCDF : Total Counts is " << noLnWeightCounts << " so the rebin is impossible, just returning a clone of the original histogram" << endl;
    cout << " Inspect the histogram for errors: " << endl;
    printHistoInfo(adaptedHisto);
    return (TH1D*) adaptedHisto->Clone();
  }

  std::vector<double> binEdges;

  int firstBin = firstBinWithContent(adaptedHisto);
  int lastBin  = lastBinWithContent(adaptedHisto);
  int numBins  = lastBin - firstBin + 1;
  double width = totalCounts/((double) a_numBins);
  if(numBins < 1){
    cout << "ERROR: HistogramUtilities::rebinEqualByCDF : NumBins is 0."
         << " Check to make sure your Histogram actually has counts : " << totalCounts << endl;
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "     CountsWithoutLN: " << noLnWeightCounts << "   CountsWithLN: " << lnWeightCounts << endl;
    cout << "     FirstBin: " << firstBin << "  lastBin: " << lastBin << "   counts: " << totalCounts << "   newBins: " << a_numBins <<  "   count Width: " << width << endl;
  #endif



  //#########   CALCULATE BIN EDGES ############
  double currentCDF = 0.0;
  double currentStep = 1.0;
  binEdges.push_back(adaptedHisto->GetBinLowEdge(firstBin));
  if(a_useLnWeight){
    if(adaptedHisto->GetBinContent(firstBin) > 0.0){
      if(a_isDifferential) currentCDF+= log(adaptedHisto->GetBinContent(firstBin) * adaptedHisto->GetBinWidth(firstBin) + 1.0);
      else                 currentCDF+= log(adaptedHisto->GetBinContent(firstBin) + 1.0);
    }
  }else{
    if(adaptedHisto->GetBinContent(firstBin) > 0.0){
      if(a_isDifferential) currentCDF+= adaptedHisto->GetBinContent(firstBin) * adaptedHisto->GetBinWidth(firstBin);
      else                 currentCDF+= adaptedHisto->GetBinContent(firstBin);
    }
  }



  for(int binIndex = firstBin+1; binIndex <= lastBin; binIndex++){
    if(currentCDF > currentStep*width){ // current bin is over the counts step, so set a bin edge
      currentStep += 1.0;
      binEdges.push_back(adaptedHisto->GetBinLowEdge(binIndex));
      #ifdef _HISTO_UTIL_DEBUG_
        cout << "        Bin Edge: " << adaptedHisto->GetBinLowEdge(binIndex) <<  " bc " << currentCDF/width << "  >  " << currentStep-1.0 <<  endl;
      #endif
    }
    if(adaptedHisto->GetBinContent(binIndex) > 0.0){ // add the content to the rolling sum
      if(a_useLnWeight){
        if(a_isDifferential) currentCDF += log(adaptedHisto->GetBinContent(binIndex) * adaptedHisto->GetBinWidth(binIndex) + 1.0);
        else                 currentCDF += log(adaptedHisto->GetBinContent(binIndex) + 1.0);
      }else{
        if(a_isDifferential) currentCDF += adaptedHisto->GetBinContent(binIndex) * adaptedHisto->GetBinWidth(binIndex);
        else                 currentCDF += adaptedHisto->GetBinContent(binIndex);
      }
    }
  }

  binEdges.push_back(adaptedHisto->GetBinLowEdge(lastBin+1));

  //#########   CREATE NEW HISTOGRAM ############

  string ending = "_reBinnedByCDF";
  string name = adaptedHisto->GetName();
  name = name + ending;
  string title = adaptedHisto->GetTitle();
  TH1D* rtrPtr = new TH1D(name.c_str(),title.c_str(),binEdges.size()-1,&binEdges[0]);
  rtrPtr->GetXaxis()->SetTitle( adaptedHisto->GetXaxis()->GetTitle());
  rtrPtr->GetYaxis()->SetTitle( adaptedHisto->GetYaxis()->GetTitle());
  if(!rtrPtr){
    cout << "ERROR: HistogramUtilities::rebinEqualByCDF couldn't create histogram" << endl;
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "    NumBins: " << numBins << "   newBins: " << rtrPtr->GetNbinsX() << endl;
  #endif


  //#########   FILL NEW HISTOGRAM ############
  int startNewBinSearchIndex = 1;
  int newNBins = rtrPtr->GetNbinsX();

  for(int binIndex = firstBin; binIndex <= lastBin; binIndex++){
    double binCenter = adaptedHisto->GetBinCenter(binIndex);
    for(int newBinIndex = startNewBinSearchIndex; newBinIndex <= newNBins; newBinIndex++){

      //#ifdef _HISTO_UTIL_DEBUG_
      //  cout << "      "   << binCenter << " > "<< rtrPtr->GetBinLowEdge(newBinIndex) << " &&" << binCenter<< "<" << rtrPtr->GetBinLowEdge(newBinIndex+1) << " = " << (binCenter > rtrPtr->GetBinLowEdge(newBinIndex) && binCenter < rtrPtr->GetBinLowEdge(newBinIndex+1)) << endl;
      //#endif
      if(binCenter > rtrPtr->GetBinLowEdge(newBinIndex) && binCenter < rtrPtr->GetBinLowEdge(newBinIndex+1)){
        if(adaptedHisto->GetBinContent(binIndex) > 0.0){
          if(a_isDifferential){
            rtrPtr->SetBinContent(newBinIndex, rtrPtr->GetBinContent(newBinIndex) + (adaptedHisto->GetBinContent(binIndex) * adaptedHisto->GetBinWidth(binIndex)));
          }else{
            rtrPtr->SetBinContent(newBinIndex, rtrPtr->GetBinContent(newBinIndex) + adaptedHisto->GetBinContent(binIndex) );
          }
        }
        startNewBinSearchIndex = newBinIndex;
      }
    }
  }


  makeSqrtNErrors(rtrPtr);
  if(a_isDifferential){
    rtrPtr->Scale(1.0,"width");
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "           Count double check: " << totalCounts << "  =? "  << sumOfContent(rtrPtr,-1,a_isDifferential) << endl;;
  #endif

  if(a_makeNoZeroBins && adaptedHisto){
    delete adaptedHisto;
  }

  return rtrPtr;
}



TH2D* HistogramUtilities::addTH2D(TH2D* a_histo1, TH2D* a_histo2){
  if(!a_histo1 || !a_histo2) return nullptr;
  TH2D* retHist = (TH2D*) a_histo1->Clone();
  for(int binx = 1; binx <= a_histo1->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo1->GetNbinsY(); biny++){
      retHist->SetBinContent(binx,biny,
        a_histo1->GetBinContent(binx,biny) + a_histo2->GetBinContent(binx,biny)
      );
      retHist->SetBinError(binx,biny,
       sqrt( pow( a_histo1->GetBinError(binx,biny)/a_histo1->GetBinContent(binx,biny),2.0)
          + pow( a_histo2->GetBinError(binx,biny)/a_histo2->GetBinContent(binx,biny),2.0)
        )
      );
    }
  }

  return retHist;
}


void HistogramUtilities::zeroBins(TH1D* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    a_histo->SetBinContent(binx,0.0);
    a_histo->SetBinError(binx,0.0);
  }

}

void HistogramUtilities::zeroBins2D(TH2D* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      a_histo->SetBinContent(binx,biny,0.0);
      a_histo->SetBinError(binx,biny,0.0);
    }
  }
}

void HistogramUtilities::zeroBins3D(TH3D* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      for(int binz = 1; binz <= a_histo->GetNbinsZ(); binz++){
        a_histo->SetBinContent(binx,biny,binz,0.0);
        a_histo->SetBinError(binx,biny,binz,0.0);
      }
    }
  }
}

void HistogramUtilities::zeroErrorBars2D(TH2D* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      a_histo->SetBinError(binx,biny,0.0);
    }
  }
}

void HistogramUtilities::zeroBins(TH1I* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    a_histo->SetBinContent(binx,0);
    a_histo->SetBinError(binx,0);
  }

}

void HistogramUtilities::zeroBins2D(TH2I* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      a_histo->SetBinContent(binx,biny,0);
      a_histo->SetBinError(binx,biny,0);
    }
  }
}

void HistogramUtilities::zeroBins3D(TH3I* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      for(int binz = 1; binz <= a_histo->GetNbinsZ(); binz++){
        a_histo->SetBinContent(binx,biny,binz,0);
        a_histo->SetBinError(binx,biny,binz,0);
      }
    }
  }
}

void HistogramUtilities::zeroInfBins(TH1D* a_histo){
  if(!a_histo) return;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    if(a_histo->GetBinWidth(binx) <=0.0){
      cout << "Bin " << binx << " has bin width of 0 or less" << endl;
    }
    if( std::isnan(a_histo->GetBinContent(binx))){
      a_histo->SetBinContent(binx, 0.0);
    }
  }

}

void HistogramUtilities::discardNanHistoData(TH1D* a_histo, bool a_coupleDataAndError){
  if(!a_histo) return;
  if(a_coupleDataAndError){
    for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
      if( !std::isfinite(a_histo->GetBinContent(binx)) || !std::isfinite(a_histo->GetBinError(binx))){
        a_histo->SetBinContent(binx, 0.0);
        a_histo->SetBinError(binx,0.0);
      }
    }
  }else{
    for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
      if( !std::isfinite(a_histo->GetBinContent(binx))){
        a_histo->SetBinContent(binx, 0.0);
      }
      if( !std::isfinite(a_histo->GetBinError(binx))){
        a_histo->SetBinError(binx,0.0);
      }
    }
  }


}

void HistogramUtilities::discardNanHistoData(TH2D* a_histo, bool a_coupleDataAndError){
  if(!a_histo) return;
  if(a_coupleDataAndError){
    for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
      for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
        if( !std::isfinite(a_histo->GetBinContent(binx,biny)) || !std::isfinite(a_histo->GetBinError(binx,biny))){
          a_histo->SetBinContent(binx,biny, 0.0);
          a_histo->SetBinError(binx,biny,0.0);
        }
      }
    }
  }else{
    for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
      for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
        if( !std::isfinite(a_histo->GetBinContent(binx,biny))){
          a_histo->SetBinContent(binx,biny, 0.0);
        }
        if( !std::isfinite(a_histo->GetBinError(binx,biny))){
          a_histo->SetBinError(binx,biny,0.0);
        }
      }
    }
  }
}

void HistogramUtilities::removeHighCountCorruptData(TH1D* a_histo, double a_thresholdFactor){
  discardNanHistoData(a_histo,true);
  double highestBin = getMaxBinInRange(a_histo);
  double highestBinCounts = a_histo->GetBinContent(highestBin);
  double secondHighestBin = -1;
  double secondHighestBinCounts = -1E100;

  double averageBinContent = 0.0;
  int endBin = a_histo->GetNbinsX();
  for(int bin = 1; bin <= endBin; bin++){
    averageBinContent += a_histo->GetBinContent(bin);
    if(bin == highestBin) continue;
    if(a_histo->GetBinContent(bin) > secondHighestBinCounts){
      secondHighestBinCounts = a_histo->GetBinContent(bin);
      secondHighestBin = bin;
    }
  }
  averageBinContent = averageBinContent / ((double) endBin);

  if(highestBinCounts > a_thresholdFactor*secondHighestBinCounts){
    a_histo->SetBinContent(highestBin,0.0);
    a_histo->SetBinError(highestBin,0.0);
  }


  //remove the ones that are greater than right or left bin above average greater than threshold
  for(int bin = 1; bin <= endBin; bin++){
    if(bin != 1){
      if(a_histo->GetBinContent(bin) > averageBinContent && a_histo->GetBinContent(bin-1) != 0.0
          && a_histo->GetBinContent(bin)/a_histo->GetBinContent(bin-1) > a_thresholdFactor){
        a_histo->SetBinContent(bin,0.0);
        a_histo->SetBinError(bin,0.0);
      }
    }
    if(bin != endBin){
      if(a_histo->GetBinContent(bin) > averageBinContent && a_histo->GetBinContent(bin+1) != 0.0
          && a_histo->GetBinContent(bin)/a_histo->GetBinContent(bin+1) > a_thresholdFactor){
        a_histo->SetBinContent(bin,0.0);
        a_histo->SetBinError(bin,0.0);
      }
    }
  }


}




void HistogramUtilities::symmeterizeHistoX(TH1D* a_histo){
  if(!a_histo) return;
  int zeroBinX = a_histo->GetXaxis()->FindBin(0.0);
  for(int binx = 1; binx < zeroBinX; binx++){
    if(binx > a_histo->GetNbinsX()) continue;
    double x_center  = a_histo->GetXaxis()->GetBinCenter(binx);
    int    mirrorBin = a_histo->GetXaxis()->FindBin(-x_center);
    if(mirrorBin > a_histo->GetNbinsX()) continue;
    double newContent = 0.5*(a_histo->GetBinContent(binx) + a_histo->GetBinContent(mirrorBin));
    double newError   = 0.5*sqrt(pow(a_histo->GetBinError(binx),2) + pow(a_histo->GetBinError(mirrorBin),2));
    a_histo->SetBinContent(binx, newContent);
    a_histo->SetBinContent(mirrorBin, newContent);
    a_histo->SetBinError(binx, newError);
    a_histo->SetBinError(mirrorBin, newError);
  }
}


void HistogramUtilities::symmeterizeHistoX(TH2D* a_histo){
  if(!a_histo) return;
  int zeroBinX = a_histo->GetXaxis()->FindBin(0.0);
  for(int binx = 1; binx < zeroBinX; binx++){
    if(binx > a_histo->GetNbinsX()) continue;
    double x_center  = a_histo->GetXaxis()->GetBinCenter(binx);
    int    mirrorBin = a_histo->GetXaxis()->FindBin(-x_center);
    if(mirrorBin > a_histo->GetNbinsX()) continue;
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      double newContent = 0.5*(a_histo->GetBinContent(binx,biny) + a_histo->GetBinContent(mirrorBin,biny));
      double newError   = 0.5*sqrt(pow(a_histo->GetBinError(binx,biny),2) + pow(a_histo->GetBinError(mirrorBin,biny),2));
      a_histo->SetBinContent(binx,      biny, newContent);
      a_histo->SetBinContent(mirrorBin, biny, newContent);
      a_histo->SetBinError(binx,      biny, newError);
      a_histo->SetBinError(mirrorBin, biny, newError);
    }
  }
}


void HistogramUtilities::setYAxisRanges_collectionOfHistos(vector<TH1D*> a_histo_vect, bool a_logScale){
  if(a_histo_vect.size() == 0){
    cout << "ERROR: HistogramUtilities::setYAxisRanges_collectionOfHistos was handed an empty vector..." << endl;
    return;
  }
  double smallestVal_overall = a_histo_vect[0]->GetBinContent(1);
  double biggestVal_overall  = a_histo_vect[0]->GetBinContent(1);
  double smallestPositiveVal_overall = (a_histo_vect[0]->GetBinContent(1) > 0 ? a_histo_vect[0]->GetBinContent(1) : 1E200);

  for(unsigned int index = 0; index < a_histo_vect.size(); index++){
    int    startBin = 1; //a_histo_vect[index]->FindBin(a_lowRangeX);
    double highestVal = a_histo_vect[index]->GetBinContent(startBin);
    double lowestVal  = highestVal;
    double lowestPositiveValue = highestVal;
    if(lowestPositiveValue <= 0.0) lowestPositiveValue = 1.0E99;
    int endBin = a_histo_vect[index]->GetNbinsX(); //a_histo_vect[index]->FindBin(a_highRangeX);
    for(int bin = startBin+1; bin <= endBin; bin++){
      if(a_histo_vect[index]->GetBinContent(bin) > highestVal){
        highestVal = a_histo_vect[index]->GetBinContent(bin);
      }
      if(a_histo_vect[index]->GetBinContent(bin) < lowestVal){
        lowestVal = a_histo_vect[index]->GetBinContent(bin);
      }
      if(a_histo_vect[index]->GetBinContent(bin) > 0.0 && a_histo_vect[index]->GetBinContent(bin) < lowestPositiveValue){
          lowestPositiveValue = a_histo_vect[index]->GetBinContent(bin);
      }
    }

    if(lowestVal < smallestVal_overall) smallestVal_overall = lowestVal;
    if(highestVal > biggestVal_overall) biggestVal_overall = highestVal;
    if(lowestPositiveValue < smallestPositiveVal_overall) smallestPositiveVal_overall = lowestPositiveValue;

  }


  double distance = biggestVal_overall - smallestVal_overall;
  //a_histo->GetXaxis()->SetRangeUser(a_lowRangeX,a_highRangeX);
  double lowYRange = smallestVal_overall - 0.1*distance;
  if(a_logScale){
    if(smallestPositiveVal_overall < 1.0E99){
      lowYRange = smallestPositiveVal_overall*0.5;
    }else{
      lowYRange = 1E-10;
    }
  }
  if(distance <= 0.0){
    cout << "WARNING: setYAxisRanges_collectionOfHistos got empty histograms or you have negative values in log scale..." << endl;
    return;
  }

  for(unsigned int index = 0; index < a_histo_vect.size(); index++){
    a_histo_vect[index]->GetYaxis()->SetRangeUser(lowYRange,biggestVal_overall + 0.1*distance);
  }
}



void HistogramUtilities::setAxisRanges(TH1D* a_histo, double a_lowRangeX, double a_highRangeX, bool a_logScale, bool a_includeErrorBars){
  int startBin      = a_histo->FindBin(a_lowRangeX);
  double highestVal = a_histo->GetBinContent(startBin);
  double lowestVal  = highestVal;
  double lowestPositiveValue = highestVal;
  if(a_includeErrorBars){
    highestVal += a_histo->GetBinError(startBin);
    lowestVal -= a_histo->GetBinError(startBin);
    lowestPositiveValue -= a_histo->GetBinError(startBin);
  }
  if(lowestPositiveValue <= 0.0) lowestPositiveValue = 1.0E99;
  int endBin = a_histo->FindBin(a_highRangeX);
  for(int bin = startBin+1; bin <= endBin; bin++){
    double valueLow  = a_histo->GetBinContent(bin);
    double valueHigh = a_histo->GetBinContent(bin);
    if(a_includeErrorBars){
      valueLow  -= a_histo->GetBinError(bin);
      valueHigh += a_histo->GetBinError(bin);
    }
    if( valueHigh > highestVal)   highestVal = valueHigh;
    if( valueLow  < lowestVal)    lowestVal  = valueLow;
    if( valueLow > 0.0 && valueLow < lowestPositiveValue){
        lowestPositiveValue = valueLow;
    }
  }

  double distance = highestVal - lowestVal;
  a_histo->GetXaxis()->SetRangeUser(a_lowRangeX,a_highRangeX);
  double lowYRange = lowestVal - 0.1*distance;
  if(a_logScale){
    if(lowestPositiveValue < 1.0E99){
      lowYRange = lowestPositiveValue*0.5;
    }else{
      lowYRange = 1E-10;
    }
  }
  if(distance > 0){
    a_histo->GetYaxis()->SetRangeUser(lowYRange,highestVal + 0.1*distance);
  }

  return;
}



// used for interpolating errors
TGraph2D* HistogramUtilities::makeErrorGraph2D(TH2D* a_histo){
  TGraph2D* graph = new TGraph2D(a_histo->GetNbinsX() * a_histo->GetNbinsY());
  int count = 0;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      graph->SetPoint(count,
            a_histo->GetXaxis()->GetBinCenter(binx),
            a_histo->GetYaxis()->GetBinCenter(biny),
            a_histo->GetBinError(binx,biny)
            );
      count++;
    }
  }
  return graph;
}



// used for making
TGraph2DErrors* HistogramUtilities::makeGraph2D(TH2D* a_histo, TGraph2DErrors* a_graph){

  TGraph2DErrors* graph = NULL;
  if(!a_graph){
     graph = new TGraph2DErrors(a_histo->GetNbinsX() * a_histo->GetNbinsY());
     graph->SetName(Form("%s_%s",a_histo->GetName(),"_graph2D"));
     graph->SetTitle(a_histo->GetTitle());
  }else{
    graph = a_graph;
  }
  int count = 0;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      graph->SetPoint(count,
            a_histo->GetXaxis()->GetBinCenter(binx),
            a_histo->GetYaxis()->GetBinCenter(biny),
            a_histo->GetBinContent(binx,biny)
            );
      graph->SetPointError(count,
            a_histo->GetXaxis()->GetBinWidth(binx),
            a_histo->GetYaxis()->GetBinWidth(biny),
            a_histo->GetBinError(binx,biny)
         );
      count++;
    }
  }
  return graph;
}


  TH1D* HistogramUtilities::make1DHistLogBins(string a_name, string a_title, int a_numBins, double a_start, double a_stop){
    if(a_start <= 0 || a_stop <=0){
      cerr << "In HistogramUtilities::make1DHistLogBins: start and stop both need to be greater than zero" << endl;
    }
    vector<double> binEdges;
    for(int binIndex = 0; binIndex <= a_numBins; binIndex++){
      binEdges.push_back(pow(10, log10(a_start) + log10(a_stop/a_start)*((double) binIndex)/((double) a_numBins)));
    }
    TH1D* logHisto = new TH1D(a_name.c_str(),a_title.c_str(),a_numBins,&binEdges[0]);
    return logHisto;
  }

  TH1I* HistogramUtilities::make1DHistLogBinsInt(string a_name, string a_title, int a_numBins, double a_start, double a_stop){
    if(a_start <= 0 || a_stop <=0){
      cerr << "In HistogramUtilities::make1DHistLogBins: start and stop both need to be greater than zero" << endl;
    }
    vector<double> binEdges;
    for(int binIndex = 0; binIndex <= a_numBins; binIndex++){
      binEdges.push_back(pow(10, log10(a_start) + log10(a_stop/a_start)*((double) binIndex)/((double) a_numBins)));
    }
    TH1I* logHisto = new TH1I(a_name.c_str(),a_title.c_str(),a_numBins,&binEdges[0]);
    return logHisto;

  }

  TH2D* HistogramUtilities::make2DHistLogBins(string a_name, string a_title,
                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY){
    if((a_makeLogX && (a_startX <= 0 || a_stopX <=0)) ||(a_makeLogY && (a_startY <= 0 || a_stopY <=0))){
      cerr << "In HistogramUtilities::make2DHistLogBins: start and stop both need to be greater than zero" << endl;
    }

    vector<double> binEdgesX;
    vector<double> binEdgesY;

    if(!a_makeLogX){
      double constXWidth = (a_stopX-a_startX)/((double) a_numBinsX);
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(constXWidth*((double) binIndex) + a_startX);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(pow(10, log10(a_startX) + log10(a_stopX/a_startX)*((double) binIndex)/((double) a_numBinsX)));
      }
    }

    if(!a_makeLogY){
      double constYWidth = (a_stopY-a_startY)/((double) a_numBinsY);
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(constYWidth*((double) binIndex) + a_startY);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(pow(10, log10(a_startY) + log10(a_stopY/a_startY)*((double) binIndex)/((double) a_numBinsY)));
      }
    }

    TH2D* logHisto = new TH2D(a_name.c_str(),a_title.c_str(),a_numBinsX,&binEdgesX[0],a_numBinsY,&binEdgesY[0]);
    return logHisto;
  }

  TH2I* HistogramUtilities::make2DHistLogBinsInt(string a_name, string a_title,
                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY){
    if((a_makeLogX && (a_startX <= 0 || a_stopX <=0)) ||(a_makeLogY && (a_startY <= 0 || a_stopY <=0))){
      cerr << "In HistogramUtilities::make2DHistLogBins: start and stop both need to be greater than zero" << endl;
    }

    vector<double> binEdgesX;
    vector<double> binEdgesY;

    if(!a_makeLogX){
      double constXWidth = (a_stopX-a_startX)/((double) a_numBinsX);
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(constXWidth*((double) binIndex) + a_startX);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(pow(10, log10(a_startX) + log10(a_stopX/a_startX)*((double) binIndex)/((double) a_numBinsX)));
      }
    }

    if(!a_makeLogY){
      double constYWidth = (a_stopY-a_startY)/((double) a_numBinsY);
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(constYWidth*((double) binIndex) + a_startY);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(pow(10, log10(a_startY) + log10(a_stopY/a_startY)*((double) binIndex)/((double) a_numBinsY)));
      }
    }


    TH2I* logHisto = new TH2I(a_name.c_str(),a_title.c_str(),a_numBinsX,&binEdgesX[0],a_numBinsY,&binEdgesY[0]);
    return logHisto;
  }

  TH3D* HistogramUtilities::make3DHistLogBins(string a_name, string a_title,
                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY,
                          bool a_makeLogZ, int a_numBinsZ, double a_startZ, double a_stopZ){
    if((a_makeLogX && (a_startX <= 0 || a_stopX <=0)) ||
       (a_makeLogY && (a_startY <= 0 || a_stopY <=0)) ||
       (a_makeLogZ && (a_startZ <= 0 || a_stopZ <=0))){
      cerr << "In HistogramUtilities::make2DHistLogBins: start and stop both need to be greater than zero" << endl;
    }
    vector<double> binEdgesX;
    vector<double> binEdgesY;
    vector<double> binEdgesZ;

    if(!a_makeLogX){
      double constXWidth = (a_stopX-a_startX)/((double) a_numBinsX);
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(constXWidth*((double) binIndex) + a_startX);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(pow(10, log10(a_startX) + log10(a_stopX/a_startX)*((double) binIndex)/((double) a_numBinsX)));
      }
    }

    if(!a_makeLogY){
      double constYWidth = (a_stopY-a_startY)/((double) a_numBinsY);
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(constYWidth*((double) binIndex) + a_startY);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(pow(10, log10(a_startY) + log10(a_stopY/a_startY)*((double) binIndex)/((double) a_numBinsY)));
      }
    }

    if(!a_makeLogZ){
      double constZWidth = (a_stopZ-a_startZ)/((double) a_numBinsZ);
      for(int binIndex = 0; binIndex <= a_numBinsZ; binIndex++){
        binEdgesZ.push_back(constZWidth*((double) binIndex) + a_startZ);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsZ; binIndex++){
        binEdgesZ.push_back(pow(10, log10(a_startZ) + log10(a_stopZ/a_startZ)*((double) binIndex)/((double) a_numBinsZ)));
      }
    }

    TH3D* logHisto = new TH3D(a_name.c_str(),a_title.c_str(),a_numBinsX,&binEdgesX[0],a_numBinsY,&binEdgesY[0],a_numBinsZ,&binEdgesZ[0]);
  return logHisto;

  }



 TH3I* HistogramUtilities::make3DHistLogBinsInt(string a_name, string a_title,

                          bool a_makeLogX, int a_numBinsX, double a_startX, double a_stopX,
                          bool a_makeLogY, int a_numBinsY, double a_startY, double a_stopY,
                          bool a_makeLogZ, int a_numBinsZ, double a_startZ, double a_stopZ){
    if((a_makeLogX && (a_startX <= 0 || a_stopX <=0)) ||
       (a_makeLogY && (a_startY <= 0 || a_stopY <=0)) ||
       (a_makeLogZ && (a_startZ <= 0 || a_stopZ <=0))){
      cerr << "In HistogramUtilities::make2DHistLogBins: start and stop both need to be greater than zero" << endl;
    }
    vector<double> binEdgesX;
    vector<double> binEdgesY;
    vector<double> binEdgesZ;

    if(!a_makeLogX){
      double constXWidth = (a_stopX-a_startX)/((double) a_numBinsX);
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(constXWidth*((double) binIndex) + a_startX);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsX; binIndex++){
        binEdgesX.push_back(pow(10, log10(a_startX) + log10(a_stopX/a_startX)*((double) binIndex)/((double) a_numBinsX)));
      }
    }

    if(!a_makeLogY){
      double constYWidth = (a_stopY-a_startY)/((double) a_numBinsY);
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(constYWidth*((double) binIndex) + a_startY);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsY; binIndex++){
        binEdgesY.push_back(pow(10, log10(a_startY) + log10(a_stopY/a_startY)*((double) binIndex)/((double) a_numBinsY)));
      }
    }

    if(!a_makeLogZ){
      double constZWidth = (a_stopZ-a_startZ)/((double) a_numBinsZ);
      for(int binIndex = 0; binIndex <= a_numBinsZ; binIndex++){
        binEdgesZ.push_back(constZWidth*((double) binIndex) + a_startZ);
      }
    }else{
      for(int binIndex = 0; binIndex <= a_numBinsZ; binIndex++){
        binEdgesZ.push_back(pow(10, log10(a_startZ) + log10(a_stopZ/a_startZ)*((double) binIndex)/((double) a_numBinsZ)));
      }
    }

    TH3I* logHisto = new TH3I(a_name.c_str(),a_title.c_str(),a_numBinsX,&binEdgesX[0],a_numBinsY,&binEdgesY[0],a_numBinsZ,&binEdgesZ[0]);
  return logHisto;

  }

vector<TH1D* > HistogramUtilities::HorizontalProjectionHistos(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol){
  vector< TH1D* > projections;
  if(!a_histo) return projections;
  for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
    TH1D* graph = a_histo->ProjectionX(Form("%s_%sIndex_%02d", a_histo->GetName(), a_axisNamePlain.c_str(), biny-1),biny,biny);
    graph->SetTitle( Form("%s with %s=[%f,%f]", a_histo->GetTitle(), a_axisSymbol.c_str(),
                  a_histo->GetYaxis()->GetBinLowEdge(biny),
                  a_histo->GetYaxis()->GetBinLowEdge(biny+1) ) );
    projections.push_back(graph);
  }
  return projections;
}



vector<TH1D* > HistogramUtilities::VerticalProjectionHistos(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol){
  vector< TH1D* > projections;
  if(!a_histo) return projections;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    TH1D* graph = a_histo->ProjectionY(Form("%s_%sIndex_%02d", a_histo->GetName(), a_axisNamePlain.c_str(), binx-1),binx,binx);
    graph->SetTitle( Form("%s with %s=[%f,%f]", a_histo->GetTitle(), a_axisSymbol.c_str(),
                  a_histo->GetXaxis()->GetBinLowEdge(binx),
                  a_histo->GetXaxis()->GetBinLowEdge(binx+1) ) );
    projections.push_back(graph);
  }
  return projections;
}


vector< TGraphErrors* > HistogramUtilities::VerticalProjectionGraphs(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol, bool a_removeZeros){
  vector< TGraphErrors* > projections;
  if(!a_histo) return projections;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    TGraphErrors* graph = new TGraphErrors();
    graph->SetName(  Form("%s_%sIndex_%02d", a_histo->GetName(), a_axisNamePlain.c_str(), binx-1));
    graph->SetTitle( Form("%s with %s=[%f,%f]", a_histo->GetTitle(), a_axisSymbol.c_str(),
                  a_histo->GetXaxis()->GetBinLowEdge(binx),
                  a_histo->GetXaxis()->GetBinLowEdge(binx+1) ) );

    // add the points into the TGraph
    int pointCounter = 0;
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      if(a_removeZeros && a_histo->GetBinContent(binx,biny) == 0.0) continue;
      graph->SetPoint(pointCounter,
            a_histo->GetYaxis()->GetBinCenter(biny),
            a_histo->GetBinContent(binx,biny));
      graph->SetPointError(pointCounter,
            a_histo->GetYaxis()->GetBinWidth(biny)*0.5,
            a_histo->GetBinError(binx,biny));
      pointCounter++;
    }
    projections.push_back(graph);
  }
  return projections;
}

vector<TGraph* > HistogramUtilities::VerticalProjectionGraphsNoErrors(TH2D* a_histo, string a_axisNamePlain, string a_axisSymbol, bool a_removeZeros){
  vector<TGraph*> projections;
  for(int binx = 1; binx <= a_histo->GetNbinsX(); binx++){
    TGraph* graph = new TGraph();
    string name = Form("%s_%sIndex_%02d", a_histo->GetName(), a_axisNamePlain.c_str(), binx-1);
    string title = Form("%s with %s=[%f,%f]", a_histo->GetTitle(), a_axisSymbol.c_str(),
                  a_histo->GetXaxis()->GetBinLowEdge(binx),
                  a_histo->GetXaxis()->GetBinLowEdge(binx+1)
                  );

    graph->SetName(name.c_str());
    graph->SetTitle(title.c_str());
    int count = 0;
    cout << "  HistogramUtilities::VerticalProjectionGraphsNoErrors making : " << name << "  " << title << endl;
    for(int biny = 1; biny <= a_histo->GetNbinsY(); biny++){
      if(a_removeZeros && a_histo->GetBinContent(binx,biny) == 0.0) continue;
      graph->SetPoint(count,
            a_histo->GetYaxis()->GetBinCenter(biny),
            a_histo->GetBinContent(binx,biny)
            );
      count++;
    }
    projections.push_back(graph);
  }
  return projections;

}

TGraphAsymmErrors* HistogramUtilities::convertSpectra_mTm0_to_pT(TGraphErrors* a_graph, double a_mass){
  if(a_graph == NULL || a_mass <=0.0) cout << "ERROR: HistogramUtilities::convertSpectra_mTm0_to_pT Handed ptr: " << a_graph << " ptr and mass=" << a_mass << endl;
  vector<double> exl;
  vector<double> exh;
  vector<double> eyl;
  vector<double> eyh;
  vector<double> xc;
  vector<double> yc;

  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double mTm0Cent;
    double yieldCent;
    a_graph->GetPoint(pointIndex,mTm0Cent,yieldCent);
    double yieldError = a_graph->GetErrorY(pointIndex);
    double xError = a_graph->GetErrorX(pointIndex);
    double mTm0Low = mTm0Cent - xError;
    double mTm0High = mTm0Cent + xError;


    double pTLow = TMath::Sqrt((mTm0Low+a_mass)*(mTm0Low+a_mass)-a_mass*a_mass);
    double pTHigh = TMath::Sqrt((mTm0High+a_mass)*(mTm0High+a_mass)-a_mass*a_mass);
    double pTCent = TMath::Sqrt((mTm0Cent+a_mass)*(mTm0Cent+a_mass)-a_mass*a_mass);
    exl.push_back( pTCent - pTLow);
    exh.push_back( pTHigh - pTCent);
    xc.push_back(pTCent);
    yc.push_back(yieldCent);
    eyl.push_back(yieldError);
    eyh.push_back(yieldError);
  }



  TGraphAsymmErrors* mappedGraph = new TGraphAsymmErrors((int)xc.size(),&xc[0],&yc[0],&exl[0],&exh[0],&eyl[0],&eyh[0]);
  string name = a_graph->GetName();
  string ending = "_pTSpace";
  string newName = name + ending;
  mappedGraph->SetName(newName.c_str());
  string title = a_graph->GetTitle();
  string newTitle = title + ending;
  mappedGraph->SetTitle(newTitle.c_str());
  return mappedGraph;
}

TGraphAsymmErrors* HistogramUtilities::convertSpectra_mTm0_to_pT(TGraphAsymmErrors* a_graph, double a_mass){
  if(a_graph == NULL || a_mass <=0.0) cout << "Handed: " << a_graph << " ptr and mass=" << a_mass << endl;
  vector<double> exl;
  vector<double> exh;
  vector<double> eyl;
  vector<double> eyh;
  vector<double> xc;
  vector<double> yc;

  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double mTm0Cent = a_graph->GetX()[pointIndex];
    double mTm0Low = mTm0Cent - a_graph->GetEXlow()[pointIndex];
    double mTm0High = mTm0Cent + a_graph->GetEXhigh()[pointIndex];
    double pTLow = TMath::Sqrt((mTm0Low+a_mass)*(mTm0Low+a_mass)-a_mass*a_mass);
    double pTHigh = TMath::Sqrt((mTm0High+a_mass)*(mTm0High+a_mass)-a_mass*a_mass);
    double pTCent = TMath::Sqrt((mTm0Cent+a_mass)*(mTm0Cent+a_mass)-a_mass*a_mass);
    exl.push_back( pTCent - pTLow);
    exh.push_back( pTHigh - pTCent);
    xc.push_back(pTCent);
    yc.push_back(a_graph->GetY()[pointIndex]);
    eyl.push_back(a_graph->GetEYlow()[pointIndex]);
    eyh.push_back(a_graph->GetEYhigh()[pointIndex]);
  }



  TGraphAsymmErrors* mappedGraph = new TGraphAsymmErrors((int)xc.size(),&xc[0],&yc[0],&exl[0],&exh[0],&eyl[0],&eyh[0]);
  string name = a_graph->GetName();
  string ending = "_pTSpace";
  string newName = name + ending;
  mappedGraph->SetName(newName.c_str());
  string title = a_graph->GetTitle();
  string newTitle = title + ending;
  mappedGraph->SetTitle(newTitle.c_str());
  return mappedGraph;
}

//____________________________________________________________________________________________________________________________________________________________________________________
TGraphAsymmErrors* HistogramUtilities::convertSpectra_pT_to_mTm0(TGraphErrors* a_graph, double a_mass){
  cout << "Handed: " << a_graph << " ptr and mass=" << a_mass << endl;
  vector<double> exl;
  vector<double> exh;
  vector<double> eyl;
  vector<double> eyh;
  vector<double> xc;
  vector<double> yc;

  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double pTCent;
    double yieldCent;
    a_graph->GetPoint(pointIndex,pTCent,yieldCent);
    double yieldError = a_graph->GetEY()[pointIndex];
    double xErrorLow  = a_graph->GetEX()[pointIndex];
    double xErrorHigh = a_graph->GetEX()[pointIndex];

    if(xErrorLow  == 0.0 && a_graph->GetN() > 1){
      if(pointIndex == 0){
        xErrorLow = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
        xErrorHigh = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
      }else if(pointIndex == a_graph->GetN() - 1){
        xErrorLow = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]);
        xErrorHigh = xErrorLow;
      }else{
        xErrorLow = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]);
        xErrorHigh = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
      }
    }
    double pTLow  = pTCent - xErrorLow;
    double pTHigh = pTCent + xErrorHigh;

    double mTm0Low  = PhysMath::mTm0(a_mass, pTLow);
    double mTm0High = PhysMath::mTm0(a_mass, pTHigh);
    double mTm0Cent = PhysMath::mTm0(a_mass, pTCent);

    exl.push_back( fabs(mTm0Cent - mTm0Low));
    exh.push_back( fabs(mTm0High - mTm0Cent));
    xc.push_back(mTm0Cent);
    yc.push_back(yieldCent);
    eyl.push_back(yieldError);
    eyh.push_back(yieldError);
  }



  TGraphAsymmErrors* mappedGraph = new TGraphAsymmErrors((int)xc.size(),&xc[0],&yc[0],&exl[0],&exh[0],&eyl[0],&eyh[0]);
  string name = a_graph->GetName();
  string ending = "_mTm0Space";
  string newName = name + ending;
  mappedGraph->SetName(newName.c_str());
  string title = a_graph->GetTitle();
  string newTitle = title + ending;
  mappedGraph->SetTitle(newTitle.c_str());
  return mappedGraph;
}

//____________________________________________________________________________________________________________________________________________________________________________________
TGraphAsymmErrors* HistogramUtilities::convertSpectra_pT_to_mTm0(TGraphAsymmErrors* a_graph, double a_mass){
  if(!a_graph || a_mass < 0.0){
    cout << "ERROR: HistogramUtilities::convertSpectra_pT_to_mTm0( TGraphAsymmErrors* " << a_graph << ", mass = " << a_mass << " )" << endl;
    return nullptr;
  }
  vector<double> exl;
  vector<double> exh;
  vector<double> eyl;
  vector<double> eyh;
  vector<double> xc;
  vector<double> yc;

  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double pTCent = a_graph->GetX()[pointIndex];
    double xErrorLow  = a_graph->GetEXlow()[pointIndex];
    double xErrorHigh = a_graph->GetEXhigh()[pointIndex];

    if(xErrorLow  == 0.0 && a_graph->GetN() > 1){
      if(pointIndex == 0){
        xErrorLow = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
        xErrorHigh = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
      }else if(pointIndex == a_graph->GetN() - 1){
        xErrorLow = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]);
        xErrorHigh = xErrorLow;
      }else{
        xErrorLow = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]);
        xErrorHigh = 0.5*fabs(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
      }
    }
    double pTLow  = pTCent - xErrorLow;
    double pTHigh = pTCent + xErrorHigh;

    double mTm0Low  = PhysMath::mTm0(a_mass, pTLow);
    double mTm0High = PhysMath::mTm0(a_mass, pTHigh);
    double mTm0Cent = PhysMath::mTm0(a_mass, pTCent);
    exl.push_back( mTm0Cent - mTm0Low);
    exh.push_back( mTm0High - mTm0Cent);
    xc.push_back(mTm0Cent);
    yc.push_back(a_graph->GetY()[pointIndex]);
    eyl.push_back(a_graph->GetEYlow()[pointIndex]);
    eyh.push_back(a_graph->GetEYhigh()[pointIndex]);
  }



  TGraphAsymmErrors* mappedGraph = new TGraphAsymmErrors((int)xc.size(),&xc[0],&yc[0],&exl[0],&exh[0],&eyl[0],&eyh[0]);
  string name = a_graph->GetName();
  string ending = "_mTm0Space";
  string newName = name + ending;
  mappedGraph->SetName(newName.c_str());
  string title = a_graph->GetTitle();
  string newTitle = title + ending;
  mappedGraph->SetTitle(newTitle.c_str());
  return mappedGraph;
}



// dmT/dpT * dN/dmT = dN/dpT
TGraphAsymmErrors* HistogramUtilities::convertPlot_mT_to_pT(TGraphErrors* a_graph, double a_mass){
  cout << "Handed: " << a_graph << " ptr and mass=" << a_mass << endl;
  vector<double> exl;
  vector<double> exh;
  vector<double> eyl;
  vector<double> eyh;
  vector<double> xc;
  vector<double> yc;

  //mT = sqrt(pT^2+m^2)
  //dmT/dpT = pT/mT


  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double mT;
    double yieldCent;
    a_graph->GetPoint(pointIndex,mT,yieldCent);
    double yieldError = a_graph->GetErrorY(pointIndex);
    double xError = a_graph->GetErrorX(pointIndex);
    double mTLow = mT - xError;
    double mTHigh = mT + xError;



    double pTLow = TMath::Sqrt(mTLow*mTLow-a_mass*a_mass);
    double pTHigh = TMath::Sqrt(mTHigh*mTLow-a_mass*a_mass);
    double pTCent = TMath::Sqrt(mT*mT-a_mass*a_mass);
    double scale = pTCent/mT;


    exl.push_back( pTCent - pTLow);
    exh.push_back( pTHigh - pTCent);
    xc.push_back(pTCent);
    yc.push_back(yieldCent*scale);
    eyl.push_back(yieldError*scale);
    eyh.push_back(yieldError);
  }



  TGraphAsymmErrors* mappedGraph = new TGraphAsymmErrors((int)xc.size(),&xc[0],&yc[0],&exl[0],&exh[0],&eyl[0],&eyh[0]);
  string name = a_graph->GetName();
  string ending = "_pTSpace";
  string newName = name + ending;
  mappedGraph->SetName(newName.c_str());
  string title = a_graph->GetTitle();
  string newTitle = title + ending;
  mappedGraph->SetTitle(newTitle.c_str());
  return mappedGraph;
}

//____________________________________________________________________________________________________________________________________________________________________________________
TGraphErrors* HistogramUtilities::convertFromYieldToSpectra(TGraphErrors* a_mTm0YieldGraph,double a_centEvents, double a_rapWidth, double a_mass, double a_efficiency, double a_efficiencyErr){
  TGraphErrors* invariantGraph = (TGraphErrors*) a_mTm0YieldGraph->Clone();
  string newName = a_mTm0YieldGraph->GetName();
  newName += "_convToSpectra";
  invariantGraph->SetName(newName.c_str());
  for(int pointIndex = 0; pointIndex < a_mTm0YieldGraph->GetN(); pointIndex++){
    double mTm0         = a_mTm0YieldGraph->GetX()[pointIndex];
    double mTm0Width    = 2.0*a_mTm0YieldGraph->GetEX()[pointIndex];
    double mT           = mTm0 + a_mass;
    double scale        = 1.0 /(a_rapWidth * mTm0Width * 2.0 * 3.14159265359 *
                          mT * a_centEvents * a_efficiency);
    double content      = a_mTm0YieldGraph->GetY()[pointIndex];
    double contentError = a_mTm0YieldGraph->GetEY()[pointIndex];
    double percentError = TMath::Sqrt(
                          1.0/a_centEvents
                          +(contentError*contentError/(content*content))
                          +(a_efficiencyErr*a_efficiencyErr/(a_efficiency*a_efficiency)));
    double newContent      = content * scale;
    double newContentError = percentError * newContent;
    invariantGraph->SetPoint(pointIndex,mTm0,newContent);
    invariantGraph->SetPointError(pointIndex,a_mTm0YieldGraph->GetEX()[pointIndex],newContentError);

  } //point Index loop
  return invariantGraph;
}

//____________________________________________________________________________________________________________________________________________________________________________________
TGraphErrors* HistogramUtilities::convertFrom_Invariant_To_NotInvariant_Spectra(TGraphErrors* a_mTm0SpectraGraph, double a_mass){
  TGraphErrors* yieldGraph = (TGraphErrors*) a_mTm0SpectraGraph->Clone();
  string newName = a_mTm0SpectraGraph->GetName();
  newName += "_convToNotInvariant";
  yieldGraph->SetName(newName.c_str());
  for(int pointIndex = 0; pointIndex < a_mTm0SpectraGraph->GetN(); pointIndex++){
    double mTm0            = a_mTm0SpectraGraph->GetX()[pointIndex];
    double scale           = 2.0 * 3.14159265359 * (mTm0 + a_mass); // 2 Pi Mt
    double newContent      = scale*a_mTm0SpectraGraph->GetY()[pointIndex];
    double newContentError = scale*a_mTm0SpectraGraph->GetEY()[pointIndex];
    yieldGraph->SetPoint(pointIndex,mTm0,newContent);
    yieldGraph->SetPointError(pointIndex,a_mTm0SpectraGraph->GetEX()[pointIndex],newContentError);
  } //mTm0 point index
  return yieldGraph;
}
//____________________________________________________________________________________________________________________________________________________________________________________
TGraphAsymmErrors* HistogramUtilities::convertFrom_Invariant_To_NotInvariant_Spectra(TGraphAsymmErrors* a_mTm0SpectraGraph, double a_mass){
  TGraphAsymmErrors* yieldGraph = (TGraphAsymmErrors*) a_mTm0SpectraGraph->Clone();
  string newName = a_mTm0SpectraGraph->GetName();
  newName += "_convToNotInvariant";
  yieldGraph->SetName(newName.c_str());
  for(int pointIndex = 0; pointIndex < a_mTm0SpectraGraph->GetN(); pointIndex++){
    double mTm0            = a_mTm0SpectraGraph->GetX()[pointIndex];
    double scale           = 2.0 * 3.14159265359 * (mTm0 + a_mass); // 2 Pi Mt
    double newContent      = scale*a_mTm0SpectraGraph->GetY()[pointIndex];
    yieldGraph->SetPoint(pointIndex,mTm0,newContent);
    yieldGraph->SetPointError(pointIndex,a_mTm0SpectraGraph->GetEXlow()[pointIndex],a_mTm0SpectraGraph->GetEXhigh()[pointIndex],scale*a_mTm0SpectraGraph->GetEYlow()[pointIndex],scale*a_mTm0SpectraGraph->GetEYhigh()[pointIndex]);
  } //mTm0 point index
  return yieldGraph;
}

//____________________________________________________________________________________________________________________________________________________________________________________

TGraphAsymmErrors* HistogramUtilities::convertFrom_pT_To_betaGamma(TGraphAsymmErrors* a_graph, double a_mass){
  TGraphAsymmErrors* graph = (TGraphAsymmErrors*) a_graph->Clone();
  string newName = a_graph->GetName();
  newName += "_convToBetaGamma";
  graph->SetName(newName.c_str());
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double pT        = a_graph->GetX()[pointIndex];
    double pT_low    = pT - a_graph->GetEXlow()[pointIndex];
    double pT_high   = pT + a_graph->GetEXhigh()[pointIndex];
    double betaGamma        = PhysMath::betagamma(pT,a_mass);
    double betaGamma_low    = PhysMath::betagamma(pT_low,a_mass);
    double betaGamma_high   = PhysMath::betagamma(pT_high,a_mass);
    graph->SetPoint(pointIndex, betaGamma, a_graph->GetY()[pointIndex]);
    graph->SetPointError(pointIndex, betaGamma-betaGamma_low, betaGamma_high-betaGamma, a_graph->GetEYlow()[pointIndex], a_graph->GetEYhigh()[pointIndex]);
  } //point index
  return graph;
}


//____________________________________________________________________________________________________________________________________________________________________________________
void HistogramUtilities::calculateFiducialDnDy(TGraphErrors* a_NotInvariantSpectra, double &a_dNdy, double &a_dNdyErr){
  a_dNdy = 0.0;
  double sum_squared_error = 0.0;
  for(int pointIndex = 0; pointIndex < a_NotInvariantSpectra->GetN(); pointIndex++){
    double mTm0Width    = 2.0*a_NotInvariantSpectra->GetEX()[pointIndex];
    a_dNdy += mTm0Width * a_NotInvariantSpectra->GetY()[pointIndex];
    sum_squared_error += mTm0Width * mTm0Width * a_NotInvariantSpectra->GetEY()[pointIndex]* a_NotInvariantSpectra->GetEY()[pointIndex];
  }
  a_dNdyErr = sqrt(sum_squared_error);
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::calculateFiducialDnDy,  dNdy = " << a_dNdy << " +/- " << a_dNdyErr << endl;
  #endif
  return;
}

//____________________________________________________________________________________________________________________________________________________________________________________
void HistogramUtilities::calculateFiducialDnDy(TGraphAsymmErrors* a_NotInvariantSpectra, double &a_dNdy, double &a_dNdyErr){
  a_dNdy = 0.0;
  double sum_squared_error = 0.0;
  for(int pointIndex = 0; pointIndex < a_NotInvariantSpectra->GetN(); pointIndex++){
    double mTm0Width    = a_NotInvariantSpectra->GetEXlow()[pointIndex] + a_NotInvariantSpectra->GetEXhigh()[pointIndex];
    a_dNdy += mTm0Width * a_NotInvariantSpectra->GetY()[pointIndex];
    sum_squared_error += mTm0Width * mTm0Width * a_NotInvariantSpectra->GetEYlow()[pointIndex]* a_NotInvariantSpectra->GetEYlow()[pointIndex];
  }
  a_dNdyErr = sqrt(sum_squared_error);
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::calculateFiducialDnDy,  dNdy = " << a_dNdy << " +/- " << a_dNdyErr << endl;
  #endif
  return;
}
//____________________________________________________________________________________________________________________________________________________________________________________
vector< vector<double> > HistogramUtilities::findSpectraHoles(TGraphErrors* a_spectra, double a_tolerance){
  a_spectra->Sort();
  vector< vector<double> > mTm0_hole_vect;
  for(int pointIndex = 0; pointIndex < a_spectra->GetN() - 1; pointIndex++){
    double edgeHighFromLowPoint = a_spectra->GetX()[pointIndex]   + a_spectra->GetEX()[pointIndex];
    double edgeLowFromHighPoint = a_spectra->GetX()[pointIndex+1] - a_spectra->GetEX()[pointIndex+1];
    if(edgeLowFromHighPoint < edgeHighFromLowPoint){
      cout << "WARNING: HistogramUtilities::findSpectraHoles has found an overlap region... " << endl;
      continue;
    }
    if(a_spectra->GetX()[pointIndex+1] < a_spectra->GetX()[pointIndex] ){
      cout << "ERROR: HistogramUtilities::findSpectraHoles is only valid with TGrpah Points In Order! " << endl;
      return mTm0_hole_vect;
    }
    if(fabs(edgeLowFromHighPoint - edgeHighFromLowPoint) > a_tolerance){
      vector< double > hole;
      hole.push_back(edgeHighFromLowPoint);
      hole.push_back(edgeLowFromHighPoint);
      mTm0_hole_vect.push_back(hole);
      #ifdef _HISTO_UTIL_DEBUG_
        cout << "Hole in spectra found for " << a_spectra->GetName() << " from " << edgeHighFromLowPoint << " to " << edgeLowFromHighPoint << endl;
      #endif
    }
  }
  return mTm0_hole_vect;
}

vector< vector<double> > HistogramUtilities::findSpectraHoles(TGraphAsymmErrors* a_spectra, double a_tolerance){
  a_spectra->Sort();
  vector< vector<double> > mTm0_hole_vect;
  for(int pointIndex = 0; pointIndex < a_spectra->GetN() - 1; pointIndex++){
    double edgeHighFromLowPoint = a_spectra->GetX()[pointIndex] + a_spectra->GetEXhigh()[pointIndex];
    double edgeLowFromHighPoint = a_spectra->GetX()[pointIndex+1] - a_spectra->GetEXlow()[pointIndex+1];
    if(edgeLowFromHighPoint < edgeHighFromLowPoint){
      cout << "WARNING: HistogramUtilities::findSpectraHoles has found an overlap region... " << endl;
      continue;
    }
    if(a_spectra->GetX()[pointIndex+1] < a_spectra->GetX()[pointIndex] ){
      cout << "ERROR: HistogramUtilities::findSpectraHoles is only valid with TGrpah Points In Order! " << endl;
      return mTm0_hole_vect;
    }
    if(fabs(edgeLowFromHighPoint - edgeHighFromLowPoint) > a_tolerance){
      vector< double > hole;
      hole.push_back(edgeHighFromLowPoint);
      hole.push_back(edgeLowFromHighPoint);
      mTm0_hole_vect.push_back(hole);
      #ifdef _HISTO_UTIL_DEBUG_
        cout << "Hole in spectra found for " << a_spectra->GetName() << " from " << edgeHighFromLowPoint << " to " << edgeLowFromHighPoint << endl;
      #endif
    }
  }
  return mTm0_hole_vect;
}


//____________________________________________________________________________________________________________________________________________________________________________________
double HistogramUtilities::scaledEntropy(TH1D* a_histo){
  TH1D* histoClone = (TH1D*) a_histo->Clone();
  double counts = sumOfContent(a_histo);
  histoClone->Scale(1.0/counts);

  double entropy = 0;
  double entropyMax = 0;
  for(int bin = 1; bin <= (int) a_histo->GetNbinsX(); bin++){
    double content = a_histo->GetBinContent(bin);
    double width = a_histo->GetBinWidth(bin);
    entropy += content*log(content)*width;
  }
  entropy = -1.0*entropy/log(2);
  double numBins = (double) a_histo->GetNbinsX();
  entropy = entropy/log(numBins);

  return entropy;
}

//____________________________________________________________________________________________________________________________________________________________________________________
TH1D* HistogramUtilities::convolution(TH1D* a_histo_1, TH1D* a_histo_2){
  //need histo with twice the limits
  double low_edge_1  = a_histo_1->GetBinLowEdge(1);
  double high_edge_1 =  a_histo_1->GetBinLowEdge(a_histo_1->GetNbinsX()+1);
  double low_edge_2  = a_histo_2->GetBinLowEdge(1);
  double high_edge_2 =  a_histo_2->GetBinLowEdge(a_histo_2->GetNbinsX()+1);

  double lowest_edge = low_edge_1;
  if(low_edge_1 > low_edge_2) lowest_edge = low_edge_2;
  double highest_edge = high_edge_1;
  if(high_edge_1 < high_edge_2) highest_edge = high_edge_2;
  double width = highest_edge - lowest_edge;

  double binWidth = a_histo_1->GetBinWidth(1);
  int numBins = (int)(width/binWidth);
  double low_edge = lowest_edge - width;
  double high_edge = highest_edge + width;

  vector<double> edges;

  for(double edge = low_edge_1; low_edge < edge; edge -= binWidth){
    edges.push_back(edge);
  }
  for(double edge = low_edge_1 + binWidth; edge < high_edge; edge += binWidth){
    edges.push_back(edge);
  }
  std::sort(edges.begin(),edges.end());

  TH1D* conv_histo = new TH1D(Form("%s_conv", a_histo_1->GetName()),
              Form("%s Convlution", a_histo_1->GetTitle()),
              ((int)edges.size()) - 1, &edges[0]);

  for(int bin_conv = 1; bin_conv <= conv_histo->GetNbinsX(); bin_conv++){
    //  h(t) = sum[ f(tau)*g(t-tau) ]
    double t = conv_histo->GetBinCenter(bin_conv);
    double newContent = 0.0;
    for(int bin2 = 1; bin2 <= a_histo_1->GetNbinsX(); bin2++){
      double tau             = a_histo_1->GetBinCenter(bin2);
      double f_tau           = a_histo_1->GetBinContent(bin2);
      int    t_minus_tau_bin = a_histo_2->FindBin(t - tau);
      if(t_minus_tau_bin < 1 || t_minus_tau_bin > a_histo_2->GetNbinsX()) continue;
      double f_t_m_tau = a_histo_2->GetBinContent(t_minus_tau_bin);
      newContent += f_tau * f_t_m_tau;
    }
    conv_histo->SetBinContent(bin_conv, newContent);
    conv_histo->SetBinError(bin_conv, 0);
  }
  return conv_histo;
}


//____________________________________________________________________________________________________________________________________________________________________________________
TH1D* HistogramUtilities::selfConvolution(TH1D* a_histo){

  return convolution(a_histo,a_histo);

  /*
  //need histo with twice the limits
  double low_edge  = a_histo->GetBinLowEdge(1);
  double high_edge =  a_histo->GetBinLowEdge(a_histo->GetNbinsX()+1);
  double width = high_edge - low_edge;
  TH1D* conv_histo = new TH1D(Form("%s_conv", a_histo->GetName()),
              Form("%s Self-Convlution", a_histo->GetTitle()),
              3*a_histo->GetNbinsX(), low_edge - width, high_edge + width);

  for(int bin_conv = 1; bin_conv <= conv_histo->GetNbinsX(); bin_conv++){
    //  h(t) = sum[ f(tau)*f(t-tau) ]
    double t = conv_histo->GetBinCenter(bin_conv);
    double newContent = 0.0;
    for(int bin2 = 1; bin2 <= a_histo->GetNbinsX(); bin2++){
      double tau             = a_histo->GetBinCenter(bin2);
      double f_tau           = a_histo->GetBinContent(bin2);
      int    t_minus_tau_bin = a_histo->FindBin(t - tau);
      if(t_minus_tau_bin < 1 || t_minus_tau_bin > a_histo->GetNbinsX()) continue;
      double f_t_m_tau = a_histo->GetBinContent(t_minus_tau_bin);
      newContent += f_tau * f_t_m_tau;
    }
    conv_histo->SetBinContent(bin_conv, newContent);
    conv_histo->SetBinError(bin_conv, 0);
  }
  return conv_histo;
  */
}


TH3D* HistogramUtilities::reBinVar3D(TH3D* a_histo, int a_rebinX, int a_rebinY, int a_rebinZ){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " In HistogramUtilities::reBinVar3D(" << a_histo << ", " << a_rebinX << ", " << a_rebinY << ", " << a_rebinZ << ") " << endl;
  #endif

  if(!a_histo) cout << "HistogramUtilities::reBinVar3D received NULL histogram Pointer" << endl;
  if(a_rebinX == 1 && a_rebinY == 1 && a_rebinZ == 1) return (TH3D*) a_histo->Clone();

  vector< double > binEdgesX;
  vector< double > binEdgesY;
  vector< double > binEdgesZ;

  int counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    if( counter % a_rebinX == 0){
      binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(a_histo->GetNbinsX()+1) );

  counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsY(); bin++){
    if( counter % a_rebinY == 0){
      binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(a_histo->GetNbinsY()+1) );

  counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsZ(); bin++){
    if( counter % a_rebinZ == 0){
      binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(a_histo->GetNbinsZ()+1) );

  TH3D* rebinHisto = new TH3D(Form("%s_rebin",a_histo->GetName()),
        Form("%s Rebin",a_histo->GetTitle()),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0],
        (int) (binEdgesZ.size()-1), &binEdgesZ[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      for(int binZ = 1; binZ <= a_histo->GetNbinsZ(); binZ++){
        double binCentZ = a_histo->GetZaxis()->GetBinCenter(binZ);
        double newbinContent = a_histo->GetBinContent(binX,binY,binZ);
        double newbinError = a_histo->GetBinError(binX,binY,binZ);

        int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
        int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
        int binFoundZ = rebinHisto->GetZaxis()->FindBin(binCentZ);
        double currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY,binFoundZ);
        double currentError = rebinHisto->GetBinError(binFoundX,binFoundY,binFoundZ);

        double sum = currentContent + newbinContent;
        double error = sqrt(newbinError*newbinError + currentError*currentError);

        rebinHisto->SetBinContent(binFoundX,binFoundY,binFoundZ,sum);
        rebinHisto->SetBinError(binFoundX,binFoundY,binFoundZ,error);
      }
    }
  }

  //rebinHisto->Sumw2();
  return rebinHisto;


}
//____________________________________________________________________________________________________________________________________________________________________________________
TH2D* HistogramUtilities::reBinVar2D(TH2D* a_histo, int a_rebinX, int a_rebinY){
  if(!a_histo){
    cout << "Null Histo passed to HistogramUtilities::reBinVar2D" << endl;
    return NULL;
  }

  if(a_rebinX == 1 && a_rebinY == 1) return (TH2D*) a_histo->Clone();
  vector< double > binEdgesX;
  vector< double > binEdgesY;

  int counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    if( counter % a_rebinX == 0){
      binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(a_histo->GetNbinsX()+1) );

  counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsY(); bin++){
    if( counter % a_rebinY == 0){
      binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(a_histo->GetNbinsY()+1) );

  TH2D* rebinHisto = new TH2D(Form("%s_rebin",a_histo->GetName()),
        a_histo->GetTitle(),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      double newbinContent = a_histo->GetBinContent(binX,binY);
      double newbinError = a_histo->GetBinError(binX,binY);

      int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
      int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
      double currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY);
      double currentError = rebinHisto->GetBinError(binFoundX,binFoundY);

      double sum = currentContent + newbinContent;
      double error = sqrt(newbinError*newbinError + currentError*currentError);

      rebinHisto->SetBinContent(binFoundX,binFoundY,sum);
      rebinHisto->SetBinError(binFoundX,binFoundY,error);


    }
  }

  //rebinHisto->Sumw2();
  return rebinHisto;


}
//____________________________________________________________________________________________________________________________________________________________________________________
TH3I* HistogramUtilities::reBinVar3D(TH3I* a_histo, int a_rebinX, int a_rebinY, int a_rebinZ){
  if(!a_histo) cout << "HistogramUtilities::reBinVar3D received NULL histogram Pointer" << endl;
  if(a_rebinX == 1 && a_rebinY == 1 && a_rebinZ == 1) return (TH3I*) a_histo->Clone();

  vector< double > binEdgesX;
  vector< double > binEdgesY;
  vector< double > binEdgesZ;

  int counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    if( counter % a_rebinX == 0){
      binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(a_histo->GetNbinsX()+1) );

  counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsY(); bin++){
    if( counter % a_rebinY == 0){
      binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(a_histo->GetNbinsY()+1) );

  counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsZ(); bin++){
    if( counter % a_rebinZ == 0){
      binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(a_histo->GetNbinsZ()+1) );

  TH3I* rebinHisto = new TH3I(Form("%s_rebin",a_histo->GetName()),
        Form("%s Rebin",a_histo->GetTitle()),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0],
        (int) (binEdgesZ.size()-1), &binEdgesZ[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      for(int binZ = 1; binZ <= a_histo->GetNbinsZ(); binZ++){
        double binCentZ = a_histo->GetZaxis()->GetBinCenter(binZ);

        int newbinContent = a_histo->GetBinContent(binX,binY,binZ);
        int newbinError = a_histo->GetBinError(binX,binY,binZ);

        int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
        int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
        int binFoundZ = rebinHisto->GetZaxis()->FindBin(binCentZ);
        int currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY,binFoundZ);
        int currentError = rebinHisto->GetBinError(binFoundX,binFoundY,binFoundZ);

        int sum = currentContent + newbinContent;
        int error = sqrt(newbinError*newbinError + currentError*currentError);

        rebinHisto->SetBinContent(binFoundX,binFoundY,binFoundZ,sum);
        rebinHisto->SetBinError(binFoundX,binFoundY,binFoundZ,error);

      }
    }
  }

  //rebinHisto->Sumw2();
  return rebinHisto;


}

 TH2I* HistogramUtilities::reBinVar2D(TH2I* a_histo, int a_rebinX, int a_rebinY){
  if(!a_histo){
    cout << "Null Histo passed to HistogramUtilities::reBinVar2D" << endl;
    return NULL;
  }

  if(a_rebinX == 1 && a_rebinY == 1) return (TH2I*) a_histo->Clone();
  vector< double > binEdgesX;
  vector< double > binEdgesY;

  int counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsX(); bin++){
    if( counter % a_rebinX == 0){
      binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(a_histo->GetNbinsX()+1) );

  counter = 0;
  for(int bin = 1; bin <= a_histo->GetNbinsY(); bin++){
    if( counter % a_rebinY == 0){
      binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
      counter = 1;
    }else{
      counter++;
    }
  }
  binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(a_histo->GetNbinsY()+1) );

  TH2I* rebinHisto = new TH2I(Form("%s_rebin",a_histo->GetName()),
        a_histo->GetTitle(),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      int newbinContent = a_histo->GetBinContent(binX,binY);
      int newbinError = a_histo->GetBinError(binX,binY);

      int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
      int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
      int currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY);
      int currentError = rebinHisto->GetBinError(binFoundX,binFoundY);

      int sum = currentContent + newbinContent;
      int error = sqrt(newbinError*newbinError + currentError*currentError);

      rebinHisto->SetBinContent(binFoundX,binFoundY,sum);
      rebinHisto->SetBinError(binFoundX,binFoundY,error);

    }
  }

  //rebinHisto->Sumw2();
  return rebinHisto;


}




TH2I* HistogramUtilities::reBinBar2D_specificY(TH2I* a_histo, vector< vector<double> > a_rebinInfo){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::reBinBar3D_specificY(" << a_histo << ")" << endl;
    cout << "   with rebin info: " << endl;
    for(unsigned int iii = 0; iii < a_rebinInfo.size(); iii++){
      cout << "        " << iii << " : ";
      for(unsigned int jjj = 0; jjj < a_rebinInfo[iii].size(); jjj++){
        cout << " " << a_rebinInfo[iii][jjj];
      }
      cout << endl;
    }
    cout << endl;
  #endif
  if(!a_histo) cout << "HistogramUtilities::reBinBar2D_specificY received NULL histogram Pointer" << endl;
  double nan_double = 1.0/0.0;

  vector< double > binEdgesX;
  vector< double > binEdgesY;

  for(int bin = 1; bin <= a_histo->GetNbinsX() + 1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Initial Y bins:  ";
  #endif
  for(int bin = 1; bin <= a_histo->GetNbinsY() + 1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
    #ifdef _HISTO_UTIL_DEBUG_
      cout << "\t" << a_histo->GetYaxis()->GetBinLowEdge(bin);
    #endif
  }
  #ifdef _HISTO_UTIL_DEBUG_
    cout << endl;
  #endif

  for(unsigned int rangeIndex = 0; rangeIndex < a_rebinInfo.size(); rangeIndex++){
    int    rebinNum  = (int) (a_rebinInfo[rangeIndex][0]);
    double minRange  = a_rebinInfo[rangeIndex][1];
    double maxRange  = a_rebinInfo[rangeIndex][2];
    int lowBinIndex  = a_histo->GetYaxis()->FindBin(minRange)-1;
    int highBinIndex = a_histo->GetYaxis()->FindBin(maxRange)-1;
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Attempting a " << rebinNum << " rebin from [" << minRange << "->" << maxRange << "] which corresponds to binIndex [" << lowBinIndex << "," << highBinIndex << "]" << endl;
    #endif
    for(int binIndex = lowBinIndex+1; binIndex < highBinIndex; binIndex++){
      if((binIndex - lowBinIndex)%rebinNum != 0) binEdgesY[binIndex] = nan_double;
    }
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Modified Edges in range:  ";
      for(int index = lowBinIndex; index < highBinIndex; index++){
        cout << "\t" << binEdgesY[index];
      }
      cout << endl;
    #endif
  }
  for(int binIndex = binEdgesY.size() -1; binIndex >= 0; binIndex--){
    if(!std::isfinite(binEdgesY[binIndex])){
      binEdgesY.erase(binEdgesY.begin() + binIndex);
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Modified Edges in range after removal:  ";
    for(int index = 0; index < binEdgesY.size(); index++){
      cout << "\t" << binEdgesY[index];
    }
    cout << endl;
  #endif

  TH2I* rebinHisto = new TH2I(Form("%s_rebin",a_histo->GetName()),
        Form("%s Rebin",a_histo->GetTitle()),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      int newbinContent = a_histo->GetBinContent(binX,binY);
      int newbinError   = a_histo->GetBinError(binX,binY);
      int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
      int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
      int currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY);
      int currentError = rebinHisto->GetBinError(binFoundX,binFoundY);
      int sum   = currentContent + newbinContent;
      int error = sqrt(newbinError*newbinError + currentError*currentError);
      rebinHisto->SetBinContent(binFoundX,binFoundY,sum);
      rebinHisto->SetBinError(binFoundX,binFoundY,error);
    }
  }
  //rebinHisto->Sumw2();
  return rebinHisto;
}



TH2D* HistogramUtilities::reBinBar2D_specificY(TH2D* a_histo, vector< vector<double> > a_rebinInfo){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::reBinBar3D_specificY(" << a_histo << ")" << endl;
    cout << "   with rebin info: " << endl;
    for(unsigned int iii = 0; iii < a_rebinInfo.size(); iii++){
      cout << "        " << iii << " : ";
      for(unsigned int jjj = 0; jjj < a_rebinInfo[iii].size(); jjj++){
        cout << " " << a_rebinInfo[iii][jjj];
      }
      cout << endl;
    }
    cout << endl;
  #endif
  if(!a_histo) cout << "HistogramUtilities::reBinBar2D_specificY received NULL histogram Pointer" << endl;
  double nan_double = 1.0/0.0;

  vector< double > binEdgesX;
  vector< double > binEdgesY;

  for(int bin = 1; bin <= a_histo->GetNbinsX() + 1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Initial Y bins:  ";
  #endif
  for(int bin = 1; bin <= a_histo->GetNbinsY() + 1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
    #ifdef _HISTO_UTIL_DEBUG_
      cout << "\t" << a_histo->GetYaxis()->GetBinLowEdge(bin);
    #endif
  }
  #ifdef _HISTO_UTIL_DEBUG_
    cout << endl;
  #endif

  for(unsigned int rangeIndex = 0; rangeIndex < a_rebinInfo.size(); rangeIndex++){
    int    rebinNum  = (int) (a_rebinInfo[rangeIndex][0]);
    double minRange  = a_rebinInfo[rangeIndex][1];
    double maxRange  = a_rebinInfo[rangeIndex][2];
    int lowBinIndex  = a_histo->GetYaxis()->FindBin(minRange)-1;
    int highBinIndex = a_histo->GetYaxis()->FindBin(maxRange)-1;
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Attempting a " << rebinNum << " rebin from [" << minRange << "->" << maxRange << "] which corresponds to binIndex [" << lowBinIndex << "," << highBinIndex << "]" << endl;
    #endif
    for(int binIndex = lowBinIndex+1; binIndex < highBinIndex; binIndex++){
      if((binIndex - lowBinIndex)%rebinNum != 0) binEdgesY[binIndex] = nan_double;
    }
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Modified Edges in range:  ";
      for(int index = lowBinIndex; index < highBinIndex; index++){
        cout << "\t" << binEdgesY[index];
      }
      cout << endl;
    #endif
  }
  for(int binIndex = binEdgesY.size() -1; binIndex >= 0; binIndex--){
    if(!std::isfinite(binEdgesY[binIndex])){
      binEdgesY.erase(binEdgesY.begin() + binIndex);
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Modified Edges in range after removal:  ";
    for(int index = 0; index < binEdgesY.size(); index++){
      cout << "\t" << binEdgesY[index];
    }
    cout << endl;
  #endif

  TH2D* rebinHisto = new TH2D(Form("%s_rebin",a_histo->GetName()),
        Form("%s Rebin",a_histo->GetTitle()),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      double newbinContent = a_histo->GetBinContent(binX,binY);
      double newbinError   = a_histo->GetBinError(binX,binY);
      int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
      int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
      double currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY);
      double currentError = rebinHisto->GetBinError(binFoundX,binFoundY);
      double sum   = currentContent + newbinContent;
      double error = sqrt(newbinError*newbinError + currentError*currentError);
      rebinHisto->SetBinContent(binFoundX,binFoundY,sum);
      rebinHisto->SetBinError(binFoundX,binFoundY,error);
    }
  }
  //rebinHisto->Sumw2();
  return rebinHisto;
}







TH3I* HistogramUtilities::reBinBar3D_specificY(TH3I* a_histo, vector< vector<double> > a_rebinInfo){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::reBinBar3D_specificY(" << a_histo << ")" << endl;
    cout << "   with rebin info: " << endl;
    for(unsigned int iii = 0; iii < a_rebinInfo.size(); iii++){
      cout << "        " << iii << " : ";
      for(unsigned int jjj = 0; jjj < a_rebinInfo[iii].size(); jjj++){
        cout << " " << a_rebinInfo[iii][jjj];
      }
      cout << endl;
    }
    cout << endl;
  #endif
  if(!a_histo) cout << "HistogramUtilities::reBinBar3D_specificY received NULL histogram Pointer" << endl;


  double nan_double = 1.0/0.0;

  vector< double > binEdgesX;
  vector< double > binEdgesY;
  vector< double > binEdgesZ;

  for(int bin = 1; bin <= a_histo->GetNbinsX() + 1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Initial Y bins:  ";
  #endif
  for(int bin = 1; bin <= a_histo->GetNbinsY() + 1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
    #ifdef _HISTO_UTIL_DEBUG_
      cout << "\t" << a_histo->GetYaxis()->GetBinLowEdge(bin);
    #endif
  }
  #ifdef _HISTO_UTIL_DEBUG_
      cout << endl;
  #endif
  for(int bin = 1; bin <= a_histo->GetNbinsZ()+1; bin++){
    binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
  }

  for(unsigned int rangeIndex = 0; rangeIndex < a_rebinInfo.size(); rangeIndex++){
    int rebinNum = (int) (a_rebinInfo[rangeIndex][0]);
    double minRange = a_rebinInfo[rangeIndex][1];
    double maxRange = a_rebinInfo[rangeIndex][2];
    int lowBinIndex = a_histo->GetYaxis()->FindBin(minRange)-1;
    int highBinIndex = a_histo->GetYaxis()->FindBin(maxRange)-1;
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Attempting a " << rebinNum << " rebin from [" << minRange << "->" << maxRange << "] which corresponds to binIndex [" << lowBinIndex << "," << highBinIndex << "]" << endl;
    #endif
    for(int binIndex = lowBinIndex+1; binIndex < highBinIndex; binIndex++){
      if((binIndex - lowBinIndex)%rebinNum != 0) binEdgesY[binIndex] = nan_double;
    }
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Modified Edges in range:  ";
      for(int index = lowBinIndex; index < highBinIndex; index++){
        cout << "\t" << binEdgesY[index];
      }
      cout << endl;
    #endif
  }

  for(int binIndex = binEdgesY.size() -1; binIndex >= 0; binIndex--){
    if(!std::isfinite(binEdgesY[binIndex])){
      binEdgesY.erase(binEdgesY.begin() + binIndex);
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Modified Edges in range after removal:  ";
    for(int index = 0; index < binEdgesY.size(); index++){
      cout << "\t" << binEdgesY[index];
    }
    cout << endl;
  #endif

  TH3I* rebinHisto = new TH3I(Form("%s_rebin",a_histo->GetName()),
        Form("%s Rebin",a_histo->GetTitle()),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0],
        (int) (binEdgesZ.size()-1), &binEdgesZ[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      for(int binZ = 1; binZ <= a_histo->GetNbinsZ(); binZ++){
        double binCentZ = a_histo->GetZaxis()->GetBinCenter(binZ);

        int newbinContent = a_histo->GetBinContent(binX,binY,binZ);
        int newbinError = a_histo->GetBinError(binX,binY,binZ);

        int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
        int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
        int binFoundZ = rebinHisto->GetZaxis()->FindBin(binCentZ);
        int currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY,binFoundZ);
        int currentError = rebinHisto->GetBinError(binFoundX,binFoundY,binFoundZ);

        int sum = currentContent + newbinContent;
        int error = sqrt(newbinError*newbinError + currentError*currentError);

        rebinHisto->SetBinContent(binFoundX,binFoundY,binFoundZ,sum);
        rebinHisto->SetBinError(binFoundX,binFoundY,binFoundZ,error);

      }
    }
  }

  //rebinHisto->Sumw2();
  return rebinHisto;


 }



TH3D* HistogramUtilities::reBinBar3D_specificY(TH3D* a_histo, vector< vector<double> > a_rebinInfo){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::reBinBar3D_specificY(" << a_histo << ")" << endl;
    cout << "   with rebin info: " << endl;
    for(unsigned int iii = 0; iii < a_rebinInfo.size(); iii++){
      cout << "        " << iii << " : ";
      for(unsigned int jjj = 0; jjj < a_rebinInfo[iii].size(); jjj++){
        cout << " " << a_rebinInfo[iii][jjj];
      }
      cout << endl;
    }
    cout << endl;
  #endif
  if(!a_histo) cout << "HistogramUtilities::reBinBar3D_specificY received NULL histogram Pointer" << endl;


  double nan_double = 1.0/0.0;

  vector< double > binEdgesX;
  vector< double > binEdgesY;
  vector< double > binEdgesZ;

  for(int bin = 1; bin <= a_histo->GetNbinsX() + 1; bin++){
    binEdgesX.push_back( a_histo->GetXaxis()->GetBinLowEdge(bin) );
  }
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Initial Y bins:  ";
  #endif
  for(int bin = 1; bin <= a_histo->GetNbinsY() + 1; bin++){
    binEdgesY.push_back( a_histo->GetYaxis()->GetBinLowEdge(bin) );
    #ifdef _HISTO_UTIL_DEBUG_
      cout << "\t" << a_histo->GetYaxis()->GetBinLowEdge(bin);
    #endif
  }
  #ifdef _HISTO_UTIL_DEBUG_
      cout << endl;
  #endif
  for(int bin = 1; bin <= a_histo->GetNbinsZ()+1; bin++){
    binEdgesZ.push_back( a_histo->GetZaxis()->GetBinLowEdge(bin) );
  }

  for(unsigned int rangeIndex = 0; rangeIndex < a_rebinInfo.size(); rangeIndex++){
    int rebinNum = (int) (a_rebinInfo[rangeIndex][0]);
    double minRange = a_rebinInfo[rangeIndex][1];
    double maxRange = a_rebinInfo[rangeIndex][2];
    int lowBinIndex = a_histo->GetYaxis()->FindBin(minRange)-1;
    int highBinIndex = a_histo->GetYaxis()->FindBin(maxRange)-1;
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Attempting a " << rebinNum << " rebin from [" << minRange << "->" << maxRange << "] which corresponds to binIndex [" << lowBinIndex << "," << highBinIndex << "]" << endl;
    #endif
    for(int binIndex = lowBinIndex+1; binIndex < highBinIndex; binIndex++){
      if((binIndex - lowBinIndex)%rebinNum != 0) binEdgesY[binIndex] = nan_double;
    }
    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Modified Edges in range:  ";
      for(int index = lowBinIndex; index < highBinIndex; index++){
        cout << "\t" << binEdgesY[index];
      }
      cout << endl;
    #endif
  }

  for(int binIndex = binEdgesY.size() -1; binIndex >= 0; binIndex--){
    if(!std::isfinite(binEdgesY[binIndex])){
      binEdgesY.erase(binEdgesY.begin() + binIndex);
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Modified Edges in range after removal:  ";
    for(int index = 0; index < binEdgesY.size(); index++){
      cout << "\t" << binEdgesY[index];
    }
    cout << endl;
  #endif

  TH3D* rebinHisto = new TH3D(Form("%s_rebin",a_histo->GetName()),
        Form("%s Rebin",a_histo->GetTitle()),
        (int) (binEdgesX.size()-1), &binEdgesX[0],
        (int) (binEdgesY.size()-1), &binEdgesY[0],
        (int) (binEdgesZ.size()-1), &binEdgesZ[0]);

  for(int binY = 1; binY <= a_histo->GetNbinsY(); binY++){
    double binCentY = a_histo->GetYaxis()->GetBinCenter(binY);
    for(int binX = 1; binX <= a_histo->GetNbinsX(); binX++){
      double binCentX = a_histo->GetXaxis()->GetBinCenter(binX);
      for(int binZ = 1; binZ <= a_histo->GetNbinsZ(); binZ++){
        double binCentZ = a_histo->GetZaxis()->GetBinCenter(binZ);

        int newbinContent = a_histo->GetBinContent(binX,binY,binZ);
        int newbinError = a_histo->GetBinError(binX,binY,binZ);

        int binFoundX = rebinHisto->GetXaxis()->FindBin(binCentX);
        int binFoundY = rebinHisto->GetYaxis()->FindBin(binCentY);
        int binFoundZ = rebinHisto->GetZaxis()->FindBin(binCentZ);
        int currentContent = rebinHisto->GetBinContent(binFoundX,binFoundY,binFoundZ);
        int currentError = rebinHisto->GetBinError(binFoundX,binFoundY,binFoundZ);

        int sum = currentContent + newbinContent;
        int error = sqrt(newbinError*newbinError + currentError*currentError);

        rebinHisto->SetBinContent(binFoundX,binFoundY,binFoundZ,sum);
        rebinHisto->SetBinError(binFoundX,binFoundY,binFoundZ,error);

      }
    }
  }

  //rebinHisto->Sumw2();
  return rebinHisto;


 }

TH1D* HistogramUtilities::rebinHistoToChebyshevSpacing(TH1D* a_histo,double a_lowRange, double a_highRange, int a_reBinFactor, int a_numPoints){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::rebinHistoToChebyshevSpacing(TH1D* " << a_histo << ",double " << a_lowRange << ", double " << a_highRange << ", int " << a_reBinFactor << ")" << endl;
  #endif
  int lowBin  = a_histo->GetXaxis()->FindBin(a_lowRange);
  int highBin = a_histo->GetXaxis()->FindBin(a_highRange);
  if(lowBin < 1) cout << "ERROR: in HistogramUtilities::rebinHistoToChebyshevSpacing... Low Range is too low for histogram: " << a_lowRange << " < " << a_histo->GetBinLowEdge(1) << endl;
  if(highBin > a_histo->GetNbinsX()) cout << "ERROR: in HistogramUtilities::rebinHistoToChebyshevSpacing... High Range is too high for histogram: " << a_highRange << " < " << a_histo->GetBinLowEdge(a_histo->GetNbinsX()+1) << endl;

  int numBinsInRange = highBin - lowBin + 1; //inclusive ex. 3 to 7 is 5 bins
  int numBinsWithRebin = numBinsInRange/a_reBinFactor;
  if(a_numPoints > 0){
    numBinsWithRebin = a_numPoints;
  }
  double lowBinCenter  = a_histo->GetXaxis()->GetBinCenter(lowBin);
  double highBinCenter = a_histo->GetXaxis()->GetBinCenter(highBin);
  double radius = 0.5*(highBinCenter - lowBinCenter);
  double center = 0.5*(highBinCenter + lowBinCenter);
  double angularStep = TMath::Pi()/numBinsWithRebin;
  vector<double> chebyshevPoints;
  chebyshevPoints.push_back(lowBinCenter);
  for(unsigned int iii = 1; iii < numBinsWithRebin; iii++){
    double index = (double) iii;
    chebyshevPoints.push_back(center - radius*cos(angularStep*index));
  }
  chebyshevPoints.push_back(highBinCenter);


  vector<double> binEdges;
  double firstWidth = chebyshevPoints[1]-chebyshevPoints[0];
  binEdges.push_back(chebyshevPoints[0] - 0.4*firstWidth);
  for(unsigned int iii = 0; iii <= numBinsWithRebin; iii++){
    double halfWidth = chebyshevPoints[iii] - binEdges[iii];
    binEdges.push_back(chebyshevPoints[iii] + halfWidth);
  }

  TH1D* rebinHisto = rebinHistoToAlternateBinStructure(a_histo, binEdges);

  #ifdef _HISTO_UTIL_DEBUG_
    cout <<" Are the bin centers at the Chebyshev Points: " << endl;;
    for(unsigned int iii = 0; iii < chebyshevPoints.size(); iii++){
      cout << chebyshevPoints[iii] << "  =? " << rebinHisto->GetBinCenter(iii+1) << endl;
    }
    cout << endl;
  #endif
  return rebinHisto;
}


TH2D* HistogramUtilities::rebin2DHistoToVertChebyshevSpacing(TH2D* a_histo,double a_lowRange, double a_highRange, int a_reBinFactor, int a_numPoints){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "Call to HistogramUtilities::rebin2DHistoToVertChebyshevSpacing(TH1D* " << a_histo << ",double " << a_lowRange << ", double " << a_highRange << ", int " << a_reBinFactor << ")" << endl;
  #endif
  int lowBin  = a_histo->GetYaxis()->FindBin(a_lowRange);
  int highBin = a_histo->GetYaxis()->FindBin(a_highRange);
  if(lowBin < 1) cout << "ERROR: in HistogramUtilities::rebin2DHistoToVertChebyshevSpacing... Low Range is too low for histogram: " << a_lowRange << " < " << a_histo->GetBinLowEdge(1) << endl;
  if(highBin > a_histo->GetNbinsY()) cout << "ERROR: in HistogramUtilities::rebin2DHistoToVertChebyshevSpacing... High Range is too high for histogram: " << a_highRange << " < " << a_histo->GetBinLowEdge(a_histo->GetNbinsX()+1) << endl;

  int numBinsInRange = highBin - lowBin + 1; //inclusive ex. 3 to 7 is 5 bins
  int numBinsWithRebin = numBinsInRange/a_reBinFactor;
  if(a_numPoints > 0){
    numBinsWithRebin = a_numPoints;
  }
  double lowBinCenter  = a_histo->GetYaxis()->GetBinCenter(lowBin);
  double highBinCenter = a_histo->GetYaxis()->GetBinCenter(highBin);
  double radius = 0.5*(highBinCenter - lowBinCenter);
  double center = 0.5*(highBinCenter + lowBinCenter);
  double angularStep = TMath::Pi()/numBinsWithRebin;
  vector<double> chebyshevPoints;
  chebyshevPoints.push_back(lowBinCenter);
  for(unsigned int iii = 1; iii < numBinsWithRebin; iii++){
    double index = (double) iii;
    chebyshevPoints.push_back(center - radius*cos(angularStep*index));
  }
  chebyshevPoints.push_back(highBinCenter);


  vector<double> binEdges;
  double firstWidth = chebyshevPoints[1]-chebyshevPoints[0];
  binEdges.push_back(chebyshevPoints[0] - 0.4*firstWidth);
  for(unsigned int iii = 0; iii <= numBinsWithRebin; iii++){
    double halfWidth = chebyshevPoints[iii] - binEdges[iii];
    binEdges.push_back(chebyshevPoints[iii] + halfWidth);
  }

  TH2D* rebinHisto = rebin2DHistoToAlternateVertBinStructure(a_histo, binEdges);

  #ifdef _HISTO_UTIL_DEBUG_
    cout <<" Are the bin centers at the Chebyshev Points: " << endl;;
    for(unsigned int iii = 0; iii < chebyshevPoints.size(); iii++){
      cout << chebyshevPoints[iii] << "  =? " << rebinHisto->GetYaxis()->GetBinCenter(iii+1) << endl;
    }
    cout << endl;
  #endif
  return rebinHisto;
}





TH1D* HistogramUtilities::rebinHistoToAlternateBinStructure(TH1D* a_histo, vector<double> a_binEdges){
  TH1D* rebinHisto = new TH1D(Form("%s_rebin",a_histo->GetName()),
        a_histo->GetTitle(), (int) (a_binEdges.size()-1), &a_binEdges[0]);
  HistogramUtilities::rebinHistoToAlternateBinStructure(a_histo,rebinHisto);
  return rebinHisto;
}

void HistogramUtilities::rebinHistoToAlternateBinStructure(TH1D* a_histo, TH1D* a_histoToPopulate){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Call to rebinHistoToAlternateBinStructure::rebinHistoToAlternateBinStructure( " << a_histo << ", " << a_histoToPopulate << " ) " << endl;
  #endif

  double new_lowEdge  = HistogramUtilities::getLowEdge(a_histoToPopulate);
  double new_highEdge = HistogramUtilities::getHighEdge(a_histoToPopulate);
  int new_nBins = a_histo->GetNbinsX();
  int old_nBins = a_histo->GetNbinsX();

  vector<double> edges; // combined edges of both histograms
  //Old Histogram
  for(int bin = 1; bin <= old_nBins+1; bin++){ // loop over old histo bins
    double lowBinEdge = a_histo->GetBinLowEdge(bin);
    if(lowBinEdge > new_lowEdge && lowBinEdge < new_highEdge){
      edges.push_back(lowBinEdge);
    }
  }
  for(int bin = 1; bin <= new_nBins+1; bin++){ // loop over old histo bins
    double lowBinEdge = a_histoToPopulate->GetBinLowEdge(bin);
    edges.push_back(lowBinEdge);
  }
  std::sort(edges.begin(), edges.end());

  //remove duplicates in the edges
  for(int iii = 1; iii < edges.size(); iii++){
    if(edges[iii-1] == edges[iii]){
      edges.erase(edges.begin()+iii);
      iii--;
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout<< " Sorted Edge Vector (duplicates removed): " << endl;
    for(int iii = 0; iii < edges.size(); iii++){
      cout << edges[iii] << "  ";
    }
    cout << endl;
  #endif


  //Determine Contents of Each Section And Add to new histogram
  for(int index = 0; index < edges.size(); index++){
    double section_center = 0.5*(edges[index]+edges[index+1]);
    int old_bin = a_histo->FindBin(section_center);
    double new_content = a_histo->GetBinContent(old_bin)*(edges[index+1] - edges[index])/a_histo->GetBinWidth(old_bin);
    int new_bin = a_histoToPopulate->FindBin(section_center);
    a_histoToPopulate->SetBinContent(new_bin, a_histoToPopulate->GetBinContent(new_bin) + new_content);
  }
  makeSqrtNErrors(a_histoToPopulate);
  a_histoToPopulate->SetTitle(a_histo->GetTitle());
  a_histoToPopulate->GetXaxis()->SetTitle(a_histo->GetXaxis()->GetTitle());
  a_histoToPopulate->GetYaxis()->SetTitle(a_histo->GetYaxis()->GetTitle());
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " End Call to rebinHistoToAlternateBinStructure::rebinHistoToAlternateBinStructure( ) " << endl;
  #endif
}


TH2D* HistogramUtilities::rebin2DHistoToAlternateVertBinStructure(TH2D* a_histo, vector<double> a_binEdges){
  vector<double> x_edges = HistogramUtilities::getBinEdges_X(a_histo);
  TH2D* rebinHisto = new TH2D(Form("%s_rebin",a_histo->GetName()),
        a_histo->GetTitle(), (int) (x_edges.size()-1), &x_edges[0], (int) (a_binEdges.size()-1), &a_binEdges[0]);
  HistogramUtilities::rebin2DHistoToAlternateVertBinStructure(a_histo,rebinHisto);
  return rebinHisto;
}

void HistogramUtilities::rebin2DHistoToAlternateVertBinStructure(TH2D* a_histo, TH2D* a_histoToPopulate){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Call to rebinHistoToAlternateBinStructure::rebin2DHistoToAlternateVertBinStructure( " << a_histo << ", " << a_histoToPopulate << " ) " << endl;
  #endif
  int new_nBins = a_histo->GetNbinsY();
  int old_nBins = a_histo->GetNbinsY();
  double new_lowEdge  = a_histoToPopulate->GetYaxis()->GetBinLowEdge(1);
  double new_highEdge = a_histoToPopulate->GetYaxis()->GetBinLowEdge(new_nBins);

  vector<double> edges; // combined edges of both histograms
  //Old Histogram
  for(int bin = 1; bin <= old_nBins+1; bin++){ // loop over old histo bins
    double lowBinEdge = a_histo->GetYaxis()->GetBinLowEdge(bin);
    if(lowBinEdge > new_lowEdge && lowBinEdge < new_highEdge){
      edges.push_back(lowBinEdge);
    }
  }
  for(int bin = 1; bin <= new_nBins+1; bin++){ // loop over old histo bins
    double lowBinEdge = a_histoToPopulate->GetYaxis()->GetBinLowEdge(bin);
    edges.push_back(lowBinEdge);
  }
  std::sort(edges.begin(), edges.end());

  //remove duplicates in the edges
  for(int iii = 1; iii < edges.size(); iii++){
    if(edges[iii-1] == edges[iii]){
      edges.erase(edges.begin()+iii);
      iii--;
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout<< " Sorted Edge Vector (duplicates removed): " << endl;
    for(int iii = 0; iii < edges.size(); iii++){
      cout << edges[iii] << "  ";
    }
    cout << endl;
  #endif


  //Determine Contents of Each Section And Add to new histogram
  for(int binX = 1; binX <= a_histo->GetNbinsX();binX++){
    for(int index = 0; index < edges.size(); index++){
      double section_center = 0.5*(edges[index]+edges[index+1]);
      int old_bin = a_histo->GetYaxis()->FindBin(section_center);
      double new_content = a_histo->GetBinContent(binX,old_bin)*(edges[index+1] - edges[index])/a_histo->GetYaxis()->GetBinWidth(old_bin);
      int new_bin = a_histoToPopulate->GetYaxis()->FindBin(section_center);
      a_histoToPopulate->SetBinContent(binX,new_bin, a_histoToPopulate->GetBinContent(binX,new_bin) + new_content);
    }
  }
  makeSqrtNErrors2D(a_histoToPopulate);
}







TH1D* HistogramUtilities::convertFromYieldToSpectra(TH1D* a_mTm0YieldHisto,double a_centEvents, double a_rapWidth, double a_mass, double a_efficiency, double a_efficiencyErr){
  TH1D* invariantHisto = (TH1D*) a_mTm0YieldHisto->Clone();
  string newName = a_mTm0YieldHisto->GetName();
  newName += "_convToSpectra";
  invariantHisto->SetName(newName.c_str());
  for(int mTm0Bin = 1; mTm0Bin <= a_mTm0YieldHisto->GetNbinsX(); mTm0Bin++){
    double mTm0         = a_mTm0YieldHisto->GetXaxis()->GetBinCenter(mTm0Bin);
    double mTm0Width    = a_mTm0YieldHisto->GetXaxis()->GetBinWidth(mTm0Bin);
    double mT           = mTm0 + a_mass;
    double scale        = 1.0 /(a_rapWidth * mTm0Width * 2.0 * 3.14159265359 *
                          mT * a_centEvents * a_efficiency);
    double content      = a_mTm0YieldHisto->GetBinContent( mTm0Bin);
    double contentError = a_mTm0YieldHisto->GetBinError( mTm0Bin);
    double percentError = TMath::Sqrt(
                          1.0/a_centEvents
                          +(contentError*contentError/(content*content))
                          +(a_efficiencyErr*a_efficiencyErr/(a_efficiency*a_efficiency)));
    double newContent      = content * scale;
    double newContentError = percentError * newContent;
    invariantHisto->SetBinContent(mTm0Bin, newContent);
    invariantHisto->SetBinError(mTm0Bin, newContentError);

  } //mTm0 Bin
  return invariantHisto;
}


TH1D* HistogramUtilities::convertFromSpectraToYield(TH1D* a_mTm0SpectraHisto, double a_centEvents, double a_rapWidth, double a_mass){
  TH1D* yieldHisto = (TH1D*) a_mTm0SpectraHisto->Clone();
  string newName = a_mTm0SpectraHisto->GetName();
  newName += "_convToYield";
  yieldHisto->SetName(newName.c_str());
  for(int mTm0Bin = 1; mTm0Bin <= a_mTm0SpectraHisto->GetNbinsX(); mTm0Bin++){
    double mTm0         = a_mTm0SpectraHisto->GetXaxis()->GetBinCenter(mTm0Bin);
    double mTm0Width    = a_mTm0SpectraHisto->GetXaxis()->GetBinWidth(mTm0Bin);
    double mT           = mTm0 + a_mass;
    double scale        = a_rapWidth * mTm0Width * 2.0 * 3.14159265359 * mT * a_centEvents;
    double content      = a_mTm0SpectraHisto->GetBinContent( mTm0Bin);
    double contentError = a_mTm0SpectraHisto->GetBinError( mTm0Bin);
    double percentError = TMath::Sqrt(
                            (contentError*contentError/(content*content))
                            -1.0/a_centEvents);
    double newContent      = content * scale;
    double newContentError = percentError * newContent;
    yieldHisto->SetBinContent(mTm0Bin, newContent);
    yieldHisto->SetBinError(mTm0Bin, newContentError);

  } //mTm0 Bin
  return yieldHisto;
}



TH2D* HistogramUtilities::convertFromYieldToSpectra2D(TH2D* a_mTm0YieldHisto,double a_centEvents, double a_mass, double a_efficiency, double a_efficiencyErr){
  TH2D* invariantHisto = (TH2D*) a_mTm0YieldHisto->Clone();
  string newName = a_mTm0YieldHisto->GetName();
  newName += "_convToSpectra";
  invariantHisto->SetName(newName.c_str());

  for(int rapBin = 1; rapBin <= a_mTm0YieldHisto->GetNbinsX(); rapBin++){
    double rapWidth = a_mTm0YieldHisto->GetXaxis()->GetBinWidth(rapBin);
    for(int mTm0Bin = 1; mTm0Bin <= a_mTm0YieldHisto->GetNbinsY(); mTm0Bin++){
      double mTm0         = a_mTm0YieldHisto->GetYaxis()->GetBinCenter(mTm0Bin);
      double mTm0Width    = a_mTm0YieldHisto->GetYaxis()->GetBinWidth(mTm0Bin);
      double mT           = mTm0 + a_mass;
      double scale        = 1.0 /(rapWidth * mTm0Width * 2.0 * 3.14159265359 *
                            mT * a_centEvents * a_efficiency);
      double content      = a_mTm0YieldHisto->GetBinContent(rapBin, mTm0Bin);
      double contentError = a_mTm0YieldHisto->GetBinError(rapBin, mTm0Bin);
      double percentError = TMath::Sqrt(
                            1.0/a_centEvents
                            +(contentError*contentError/(content*content))
                            +(a_efficiencyErr*a_efficiencyErr/(a_efficiency*a_efficiency)));
      double newContent      = content * scale;
      double newContentError = percentError * newContent;
      invariantHisto->SetBinContent(rapBin,mTm0Bin, newContent);
      invariantHisto->SetBinError(rapBin,mTm0Bin, newContentError);

    } //mTm0 Bin
  }// rap Bin
  return invariantHisto;
}


TH1D* HistogramUtilities::convertZTOFToMassSqr(TH1D* a_invBetaHisto,double a_p, double a_m_calib, bool a_ZTOF_isDifferential, bool a_massSqr_isDifferential){
  vector<double> orig_edges = HistogramUtilities::getBinEdges(a_invBetaHisto);
  vector<double> new_edges;
  #ifdef _HISTO_UTIL_DEBUG_
    cout << "HistogramUtilities::convertInvBetaToMassSqr(" << a_invBetaHisto << ", " << a_p << " )" << endl;
    cout << "  Number of InvBeta Edges: " << orig_edges.size() << endl;
    cout << "\n\n\n\n\n  Inverse Beta Histogram before conversion: " << endl;
    printHistoInfo(a_invBetaHisto);
  #endif
  for(int index = 0; index < orig_edges.size(); index++){
    double invBeta = orig_edges[index] + PhysMath::invBeta(a_m_calib,a_p);
    new_edges.push_back(PhysMath::massSquared_signed_fromInvBeta(a_p,invBeta));
  }
  string newName = a_invBetaHisto->GetName();
  newName += "_convToM2";
  string newTitle  = a_invBetaHisto->GetTitle();
  string yAxisName = a_invBetaHisto->GetYaxis()->GetTitle();
  TH1D* newHisto   = new TH1D(newName.c_str(),newTitle.c_str(),new_edges.size()-1,&new_edges[0]);
  newHisto->GetYaxis()->SetTitle(yAxisName.c_str());
  newHisto->GetXaxis()->SetTitle("m^{2} (GeV/c)^{2}");
  if(a_ZTOF_isDifferential){
    if(a_massSqr_isDifferential){
      for(int bin = 1; bin <= a_invBetaHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_invBetaHisto->GetBinContent(bin)*a_invBetaHisto->GetBinWidth(bin)/newHisto->GetBinWidth(bin));
        newHisto->SetBinError(bin,  a_invBetaHisto->GetBinError(bin)  *a_invBetaHisto->GetBinWidth(bin)/newHisto->GetBinWidth(bin));
      } //mTm0 Bin
    }else{
      for(int bin = 1; bin <= a_invBetaHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_invBetaHisto->GetBinContent(bin)*a_invBetaHisto->GetBinWidth(bin));
        newHisto->SetBinError(bin,  a_invBetaHisto->GetBinError(bin)  *a_invBetaHisto->GetBinWidth(bin));
      } //mTm0 Bin
    }
  }else{
    if(a_massSqr_isDifferential){
      for(int bin = 1; bin <= a_invBetaHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_invBetaHisto->GetBinContent(bin)/newHisto->GetBinWidth(bin));
        newHisto->SetBinError(bin,  a_invBetaHisto->GetBinError(bin)  /newHisto->GetBinWidth(bin));
      } //mTm0 Bin
    }else{
      for(int bin = 1; bin <= a_invBetaHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_invBetaHisto->GetBinContent(bin));
        newHisto->SetBinError(bin,  a_invBetaHisto->GetBinError(bin));
      } //mTm0 Bin
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "\n\n\n\n\n  Number of m^2 Edges: " << new_edges.size() << endl;
    cout << "  Mass Squared Histogram after conversion: " << endl;
    printHistoInfo(newHisto);
  #endif

  return newHisto;
}


TH1D* HistogramUtilities::convertMassSqrToInvBeta(TH1D* a_massSqrHisto,double a_p, bool a_massSqr_isDifferential, bool a_ZTOF_isDifferential, double a_m_calib){
  vector<double> orig_edges = HistogramUtilities::getBinEdges(a_massSqrHisto);
  vector<double> new_edges;
  double invBeta_calib = PhysMath::invBeta_from_MassSquaredSigned(a_m_calib*a_m_calib,a_p);

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "HistogramUtilities::convertMassSqrToInvBeta(" << a_massSqrHisto << ", " << a_p << " )" << endl;
    cout << "  Number of m2 Edges: " << orig_edges.size() << endl;
    cout << "\n\n\n\n\n  m2 Beta Histogram before conversion: " << endl;
    printHistoInfo(a_massSqrHisto);
    cout << "Inv Beta Calib: " << invBeta_calib << endl;
  #endif

  for(int index = 0; index < orig_edges.size(); index++){
    double invBetaEdge = PhysMath::invBeta_from_MassSquaredSigned(orig_edges[index],a_p);
    if(a_m_calib >= 0.0) invBetaEdge = invBetaEdge - invBeta_calib;
    new_edges.push_back(invBetaEdge);
  }
  string newName = a_massSqrHisto->GetName();
  newName += "_convToInvBeta";
  string newTitle  = a_massSqrHisto->GetTitle();
  string yAxisName = a_massSqrHisto->GetYaxis()->GetTitle();
  TH1D* newHisto   = new TH1D(newName.c_str(),newTitle.c_str(),new_edges.size()-1,&new_edges[0]);
  newHisto->GetYaxis()->SetTitle(yAxisName.c_str());
  string xLabel = "1/#beta";
  if(a_m_calib >= 0 ) xLabel += " - 1/#beta_{prediction}";
  newHisto->GetXaxis()->SetTitle(xLabel.c_str());
  if(a_massSqr_isDifferential){
    if(a_ZTOF_isDifferential){
      for(int bin = 1; bin <= a_massSqrHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_massSqrHisto->GetBinContent(bin)*a_massSqrHisto->GetBinWidth(bin)/newHisto->GetBinWidth(bin));
        newHisto->SetBinError(bin,  a_massSqrHisto->GetBinError(bin)  *a_massSqrHisto->GetBinWidth(bin)/newHisto->GetBinWidth(bin));
      } //mTm0 Bin
    }else{
      for(int bin = 1; bin <= a_massSqrHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_massSqrHisto->GetBinContent(bin)*a_massSqrHisto->GetBinWidth(bin));
        newHisto->SetBinError(bin,  a_massSqrHisto->GetBinError(bin)  *a_massSqrHisto->GetBinWidth(bin));
      } //mTm0 Bin
    }
  }else{
    if(a_ZTOF_isDifferential){
      for(int bin = 1; bin <= a_massSqrHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_massSqrHisto->GetBinContent(bin)/newHisto->GetBinWidth(bin));
        newHisto->SetBinError(bin,  a_massSqrHisto->GetBinError(bin)  /newHisto->GetBinWidth(bin));
      } //mTm0 Bin
    }else{
      for(int bin = 1; bin <= a_massSqrHisto->GetNbinsX(); bin++){
        newHisto->SetBinContent(bin,a_massSqrHisto->GetBinContent(bin));
        newHisto->SetBinError(bin,  a_massSqrHisto->GetBinError(bin));
      } //mTm0 Bin
    }
  }

  #ifdef _HISTO_UTIL_DEBUG_
    cout << "\n\n\n\n\n  Number of invBeta Edges: " << new_edges.size() << endl;
    cout << "  InvBeta Histogram after conversion: " << endl;
    printHistoInfo(newHisto);
  #endif

  return newHisto;
}


/*TH1D* HistogramUtilities::convertMassSqrToInvBeta(TH1D* a_massSqrHisto,double a_p){
  vector<double> orig_edges = HistogramUtilities::getBinEdges(a_massSqrHisto);
  vector<double> new_edges;
  for(int index = 0; index < orig_edges.size(); index++){
    new_edges.push_back(PhysMath::invBeta(orig_edges[index],a_p));
  }
  string newName = a_massSqrHisto->GetName();
  newName += "_convToInvBeta";
  string newTitle  = a_massSqrHisto->GetTitle();
  string yAxisName = a_massSqrHisto->GetYaxis()->GetTitle();
  TH1D* newHisto   = new TH1D(newName.c_str(),newTitle.c_str(),new_edges.size()-1,&new_edges[0]);
  newHisto->GetYaxis()->SetTitle(yAxisName.c_str());
  for(int bin = 1; bin <= a_massSqrHisto->GetNbinsX(); bin++){
    newHisto->SetBinContent(bin,a_massSqrHisto->GetBinContent(bin));
    newHisto->SetBinError(bin,a_massSqrHisto->GetBinError(bin));
  } //mTm0 Bin
  return newHisto;
}*/


TH1D* HistogramUtilities::convertSpectraHisto_mTm0_to_pT(TH1D* a_histo, double a_mass){
  if(a_histo == NULL || a_mass <=0.0 || a_histo->GetNbinsX() <= 0){
    cout << "ERROR: HistogramUtilities... Handed: " << a_histo << " ptr and mass=" << a_mass << " numbins: " << a_histo->GetNbinsX() <<  endl;
    return (TH1D*) NULL;
  }
  vector<double> edges = HistogramUtilities::getBinEdges(a_histo);
  for(unsigned int index = 0; index < edges.size(); index++){
    edges[index] = PhysMath::pT_from_mTm0(edges[index],a_mass);
  }
  string name = a_histo->GetName();
  string newName = name + "_pTSpace";
  string title = a_histo->GetTitle();
  TH1D* pTHisto = new TH1D(newName.c_str(),title.c_str(),edges.size()-1,&edges[0]);
  pTHisto->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  pTHisto->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{2#pi p_{T}} #frac{d^{2}N}{dp_{T}dy}");
  for(int bin = 1; bin <= pTHisto->GetNbinsX(); bin++){
    pTHisto->SetBinContent(bin,a_histo->GetBinContent(bin));
    pTHisto->SetBinError(bin,a_histo->GetBinError(bin));
  }
  return pTHisto;
}

TH1D* HistogramUtilities::convertSpectraHisto_pT_to_mTm0(TH1D* a_histo, double a_mass){
  if(a_histo == NULL || a_mass <=0.0 || a_histo->GetNbinsX() <= 0){
    cout << "ERROR: HistogramUtilities... Handed: " << a_histo << " ptr and mass=" << a_mass << " numbins: " << a_histo->GetNbinsX() <<  endl;
    return (TH1D*) NULL;
  }
  vector<double> edges = HistogramUtilities::getBinEdges(a_histo);
  for(unsigned int index = 0; index < edges.size(); index++){
    edges[index] = PhysMath::mTm0(a_mass,edges[index]);
  }
  string name = a_histo->GetName();
  string newName = name + "_mTm0Space";
  string title = a_histo->GetTitle();
  TH1D* mTHisto = new TH1D(newName.c_str(),title.c_str(),edges.size()-1,&edges[0]);
  mTHisto->GetXaxis()->SetTitle("m_{T}-m_{0} [GeV/c^{2}]");
  mTHisto->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{1}{2#pi m_{T}} #frac{d^{2}N}{dm_{T}dy}");
  for(int bin = 1; bin <= mTHisto->GetNbinsX(); bin++){
    mTHisto->SetBinContent(bin,a_histo->GetBinContent(bin));
    mTHisto->SetBinError(bin,a_histo->GetBinError(bin));
  }
  return mTHisto;
}














void HistogramUtilities::printGraphInfo(TGraphAsymmErrors* a_graph){
  cout << "TGraphAsymmErrors Ptr: " << a_graph << endl;
  if(!a_graph) return;
  cout << "Name:  " << a_graph->GetName() << endl;
  cout << "Title: " << a_graph->GetTitle() << endl;
  cout << "Point by Point Values: " << endl;
  a_graph->Print();
}
void HistogramUtilities::printGraphInfo(TGraphErrors* a_graph){
  cout << "TGraphErrors Ptr: " << a_graph << endl;
  if(!a_graph) return;
  cout << "Name:  " << a_graph->GetName() << endl;
  cout << "Title: " << a_graph->GetTitle() << endl;
  cout << "Point by Point Values: " << endl;
  a_graph->Print();
}

void HistogramUtilities::printReconstrucitonCode_TGraphErrors(TGraphErrors* a_graph){
  if(!a_graph) return;
  int numPoints = a_graph->GetN();
  string name = a_graph->GetName();
  if(name == "") name = "temp";

  cout << "double "<< name << "_x[" << numPoints << "] = {";
  for(int iii = 0; iii < numPoints - 1; iii++){
    cout << a_graph->GetX()[iii] << ",";
  }
  cout << a_graph->GetX()[numPoints - 1] << "};" << endl;

  cout << "double "<< name << "_y[" << numPoints << "] = {";
  for(int iii = 0; iii < numPoints - 1; iii++){
    cout << a_graph->GetY()[iii] << ",";
  }
  cout << a_graph->GetY()[numPoints - 1] << "};" << endl;

  cout << "double "<< name << "_ex[" << numPoints << "] = {";
  for(int iii = 0; iii < numPoints - 1; iii++){
    cout << a_graph->GetEX()[iii] << ",";
  }
  cout << a_graph->GetEX()[numPoints - 1] << "};" << endl;

  cout << "double "<< name << "_ey[" << numPoints << "] = {";
  for(int iii = 0; iii < numPoints - 1; iii++){
    cout << a_graph->GetEY()[iii] << ",";
  }
  cout << a_graph->GetEY()[numPoints - 1] << "};" << endl;
  cout << "TGraphErrors* " << name << " = new TGraphErrors(" << numPoints << "," << name << "_x,"<< name << "_y,"<< name << "_ex,"<< name << "_ey);" << endl;
  cout << name << "->SetName(\"" << a_graph->GetName() << "\");" << endl;
  cout << name << "->SetTitle(\"" << a_graph->GetTitle() << "\");" << endl;
  cout << name << "->SetMarkerStyle(" << a_graph->GetMarkerStyle() << ");" << endl;
  cout << name << "->SetMarkerColor(" << a_graph->GetMarkerColor() << ");" << endl;
  cout << name << "->SetLineColor(" << a_graph->GetLineColor() << ");" << endl;
  cout << name << "->SetMarkerSize(" << a_graph->GetMarkerSize() << ");" << endl;

}


TGraphErrors* HistogramUtilities::makeGraphFromDataTheif(string a_name, int a_numPoints, double* a_xVal, double* a_yVal, double* a_yValPlusErr, double* a_xValPlusErr){
  vector<double> xErr;
  vector<double> yErr;
  for(int iii = 0; iii < a_numPoints; iii++){
    if(a_xValPlusErr) xErr.push_back(fabs(a_xValPlusErr[iii] - a_xVal[iii]));
    yErr.push_back(fabs(a_yValPlusErr[iii] - a_yVal[iii]));
  }
  TGraphErrors* graph = NULL;
  if(!a_xValPlusErr){
    vector<double> zeros(a_numPoints);
    graph = new TGraphErrors(a_numPoints,a_xVal,a_yVal,&zeros[0],&yErr[0]);
  }else{
    graph = new TGraphErrors(a_numPoints,a_xVal,a_yVal,&xErr[0],&yErr[0]);
  }
  graph->SetName(a_name.c_str());
  return graph;
}





TH1D* HistogramUtilities::makeHistogram(TGraphErrors* a_graph){
  vector<double> binEdges;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    binEdges.push_back(a_graph->GetX()[pointIndex] - a_graph->GetEX()[pointIndex]);
  }
  binEdges.push_back(a_graph->GetX()[a_graph->GetN() - 1] + a_graph->GetEX()[a_graph->GetN() - 1]);
  string name = a_graph->GetName();
  name =+ "_graph";
  TH1D* histo = new TH1D(name.c_str(),a_graph->GetTitle(),((int) binEdges.size()) - 1, &binEdges[0]);
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    histo->SetBinContent(pointIndex + 1,a_graph->GetY()[pointIndex]);
    histo->SetBinError(pointIndex + 1,a_graph->GetEY()[pointIndex]);
  }
  return histo;
}

TH1D* HistogramUtilities::makeHistogram_Asymm(TGraphAsymmErrors* a_graph){
  vector<double> binEdges;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double edge = a_graph->GetX()[pointIndex] - a_graph->GetEXlow()[pointIndex];
    binEdges.push_back(edge);
  }
  double edge = a_graph->GetX()[a_graph->GetN() - 1] + a_graph->GetEXhigh()[a_graph->GetN() - 1];
  binEdges.push_back(edge);
  string name = a_graph->GetName();
  name =+ "_graph";
  TH1D* histo = new TH1D(name.c_str(),a_graph->GetTitle(),((int)binEdges.size()) - 1, &binEdges[0]);
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    int bin = pointIndex + 1;
    histo->SetBinContent(bin,a_graph->GetY()[pointIndex]);
    histo->SetBinError(bin,a_graph->GetEYhigh()[pointIndex]);
  }
  return histo;
}





int HistogramUtilities::getPointClosestToVal(TGraph* a_graph, double a_x){

  int pointClosest = -1;
  double smallestDiff = 1e100;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double difference = fabs(a_x - a_graph->GetX()[pointIndex]);
    if(difference < smallestDiff){
      smallestDiff = difference;
      pointClosest = pointIndex;
    }
  }//loop over points

  if(pointClosest == -1){
    cout << "ERROR IN HistogramUtilities::getPointClosestToVal, couldnt find closest x, probably because no points in TGraph" << endl;
  }
  return pointClosest;
}


int HistogramUtilities::getPointClosestToVal(TGraphErrors* a_graph, double a_x){
  int pointClosest = -1;
  double smallestDiff = 1e100;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double difference = fabs(a_x - a_graph->GetX()[pointIndex]);
    if(difference < smallestDiff){
      smallestDiff = difference;
      pointClosest = pointIndex;
    }
  }//loop over points

  if(pointClosest == -1){
    cout << "ERROR IN HistogramUtilities::getPointClosestToVal, couldnt find closest x, probably because no points in TGraph" << endl;
  }
  return pointClosest;
}


int HistogramUtilities::getPointClosestToVal(TGraphAsymmErrors* a_graph, double a_x){
  int pointClosest = -1;
  double smallestDiff = 1e100;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double difference = fabs(a_x - a_graph->GetX()[pointIndex]);
    if(difference < smallestDiff){
      smallestDiff = difference;
      pointClosest = pointIndex;
    }
  }//loop over points

  if(pointClosest == -1){
    cout << "ERROR IN HistogramUtilities::getPointClosestToVal, couldnt find closest x, probably because no points in TGraph" << endl;
  }
  return pointClosest;
}




void HistogramUtilities::evalWithError(TGraphErrors* a_graph, double a_evalX, double &a_yVal, double &a_yErr){
  if( !a_graph || a_graph->GetN() < 2 ) return;
  a_graph->Sort();
  int point = TMath::BinarySearch(a_graph->GetN(),a_graph->GetX(),a_evalX);
  int firstPoint  = 0;
  int secondPoint = 0;
  if(point < 1 ){
    firstPoint  = 0;
    secondPoint = 1;
  }else{
    if( point >= a_graph->GetN() - 1){
      secondPoint = a_graph->GetN() - 1;
      firstPoint  = secondPoint - 1;
    }else{
      if(a_evalX < a_graph->GetX()[point]){
        firstPoint  = point - 1;
        secondPoint = point;
      }else{
        firstPoint  = point;
        secondPoint = point + 1;
      }
    }
  }

  double x1    = a_graph->GetX()[firstPoint];
  double y1    = a_graph->GetY()[firstPoint];
  double yErr1 = a_graph->GetEY()[firstPoint];
  double x2    = a_graph->GetX()[secondPoint];
  double y2    = a_graph->GetY()[secondPoint];
  double yErr2 = a_graph->GetEY()[secondPoint];
  PhysMath::interpolate(a_evalX, x1, y1, yErr1, x2, y2, yErr2, a_yVal, a_yErr);

  /*
  int low = -1;
  int up = -1;
  int nPoints = a_graph->GetN();
  double xLow = -999;
  double xUp = -999;
  double yLow = -999;
  double yUp = -999;

  a_graph->Sort();
  low = TMath::BinarySearch(nPoints,a_graph->GetX(),a_evalX);
  if(low == -1) low = 0;
  if(low == nPoints-1) low--;
  if(a_graph->GetX()[low] == a_evalX){
    yVal = a_graph->GetY()[low];
    yErr = a_graph->GetEY()[low];
    return;
  }
  up = low + 1;
  if(a_graph->GetX()[up] == a_evalX){
    yVal = a_graph->GetY()[up];
    yErr = a_graph->GetEY()[up];
    return;
  }
  if(a_graph->GetX()[low] == a_graph->GetX()[up]){
    yVal = a_graph->GetY()[low];
    yErr = a_graph->GetEY()[low];
    return;
  }

  xLow = a_graph->GetX()[low];
  xUp = a_graph->GetX()[up];
  yLow = a_graph->GetY()[low];
  yUp = a_graph->GetY()[up];

  yVal = yUp + (a_evalX - xUp)*(yLow - yUp)/(xLow - xUp);
  yErr = sqrt(pow(1-(a_evalX - xUp)/(xLow - xUp),2)*pow(a_graph->GetEY()[up],2) + pow((a_evalX - xUp)/(xLow - xUp),2)*pow(a_graph->GetEY()[low],2));
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " HistogramUtilities::evalWithError  lowPt: " << low << " xLow: " << xLow << " yLow: " << yLow << " yErr: " <<  a_graph->GetEY()[low] << yLow
                                            << " highPt: " << up << " xUp: " << xUp << " yUp: " << yUp << " yErr: " << a_graph->GetEY()[up] << endl;
  #endif

  return;*/
}

TGraphErrors* HistogramUtilities::makeDerivativeGraph(TGraphAsymmErrors* a_graph){
  vector< double > x_vect;
  vector< double > xErr_vect;
  vector< double > yp_vect;
  vector< double > ypErr_vect;

  for(int pointIndex = 1; pointIndex < a_graph->GetN()-1; pointIndex++){
     double yM1Err = 0.5*(a_graph->GetEYlow()[pointIndex-1] + a_graph->GetEYlow()[pointIndex-1]);
     double yErr = 0.5*(a_graph->GetEYlow()[pointIndex] + a_graph->GetEYlow()[pointIndex]);
     double yP1Err = 0.5*(a_graph->GetEYlow()[pointIndex+1] + a_graph->GetEYlow()[pointIndex+1]);
     double derivative_one = (a_graph->GetY()[pointIndex] - a_graph->GetY()[pointIndex-1])/(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]);
     double derivative_one_sqrErr = (yM1Err*yM1Err + yErr*yErr)/((a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1])*(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]));
     double derivative_two = (a_graph->GetY()[pointIndex] - a_graph->GetY()[pointIndex+1])/(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
     double derivative_two_sqrErr = (yP1Err*yP1Err + yErr*yErr)/((a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1])*(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]));

     double derivative = 0.5*(derivative_one + derivative_two);
     double derivativeError = 0.5*sqrt(derivative_one_sqrErr + derivative_two_sqrErr);

     x_vect.push_back(a_graph->GetX()[pointIndex]);
     xErr_vect.push_back(0.0);
     yp_vect.push_back(derivative);
     ypErr_vect.push_back(derivativeError);

  }
  return new TGraphErrors(x_vect.size(),&x_vect[0],&yp_vect[0],&xErr_vect[0],&ypErr_vect[0]);

}






TGraphErrors* HistogramUtilities::makeDerivativeGraph(TGraphErrors* a_graph){

  vector< double > x_vect;
  vector< double > xErr_vect;
  vector< double > yp_vect;
  vector< double > ypErr_vect;

  for(int pointIndex = 1; pointIndex < a_graph->GetN()-1; pointIndex++){
     double yM1Err = a_graph->GetEY()[pointIndex-1];
     double yErr = a_graph->GetEY()[pointIndex];
     double yP1Err = a_graph->GetEY()[pointIndex+1];
     double derivative_one = (a_graph->GetY()[pointIndex] - a_graph->GetY()[pointIndex-1])/(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]);
     double derivative_one_sqrErr = (yM1Err*yM1Err + yErr*yErr)/((a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1])*(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex-1]));
     double derivative_two = (a_graph->GetY()[pointIndex] - a_graph->GetY()[pointIndex+1])/(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]);
     double derivative_two_sqrErr = (yP1Err*yP1Err + yErr*yErr)/((a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1])*(a_graph->GetX()[pointIndex] - a_graph->GetX()[pointIndex+1]));

     double derivative = 0.5*(derivative_one + derivative_two);
     double derivativeError = 0.5*sqrt(derivative_one_sqrErr + derivative_two_sqrErr);

     x_vect.push_back(a_graph->GetX()[pointIndex]);
     xErr_vect.push_back(0.0);
     yp_vect.push_back(derivative);
     ypErr_vect.push_back(derivativeError);

  }


  return new TGraphErrors(x_vect.size(),&x_vect[0],&yp_vect[0],&xErr_vect[0],&ypErr_vect[0]);

}


TGraphAsymmErrors* HistogramUtilities::makeFitPullGraph(TGraphAsymmErrors* a_graph, TF1* a_funct, bool a_doIntegral){
  TGraphAsymmErrors* pullGraph = new TGraphAsymmErrors();
  string newName = a_graph->GetName();
  newName = newName + "_pullToFit";
  pullGraph->SetName(newName.c_str());
  double lowFunctRange = 0.0;
  double highFunctRange = 0.0;
  a_funct->GetRange(lowFunctRange,highFunctRange);
  //pullGraph->GetYaxis()->SetTitle("Pulls");
  int numPointsCurrent = 0;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double center  = a_graph->GetX()[pointIndex];
    if(center < lowFunctRange || center > highFunctRange) continue;
    double lowEdge = center - a_graph->GetEXlow()[pointIndex];
    double highEdge = center - a_graph->GetEXhigh()[pointIndex];
    double content  = a_graph->GetY()[pointIndex];
    double errorYhigh = a_graph->GetEYhigh()[pointIndex];
    double errorYlow = a_graph->GetEYlow()[pointIndex];
    double fit = a_funct->Eval(center);
    if(a_doIntegral){
      fit = a_funct->Integral(lowEdge,highEdge)/(highEdge - lowEdge);
    }
    if(content > fit){
      double pull = (content - fit)/errorYlow;
      if(!std::isfinite(pull)) pull = -999.999;
      double pullErrLow = 1.0;
      double pullErrHigh = errorYhigh/errorYlow;
      if(!std::isfinite(pullErrHigh)) pullErrHigh = 0.0;
      pullGraph->SetPoint(numPointsCurrent,center,pull);
      pullGraph->SetPointError(numPointsCurrent,a_graph->GetEXlow()[pointIndex],a_graph->GetEXhigh()[pointIndex],pullErrLow,pullErrHigh);
    }else{
      double pull = (content - fit)/errorYhigh;
      if(!std::isfinite(pull)) pull = -999.999;
      double pullErrLow = errorYlow/errorYhigh;
      double pullErrHigh = 1.0;
      if(!std::isfinite(pullErrLow)) pullErrLow = 0.0;
      pullGraph->SetPoint(numPointsCurrent,center,pull);
      pullGraph->SetPointError(numPointsCurrent,a_graph->GetEXlow()[pointIndex],a_graph->GetEXhigh()[pointIndex],pullErrLow,pullErrHigh);
    }
    numPointsCurrent += 1;
  }
  return pullGraph;
}


TGraphAsymmErrors* HistogramUtilities::makeFitPullGraph(TGraphErrors* a_graph, TF1* a_funct, bool a_doIntegral){
  TGraphAsymmErrors* convGraph = convertNoAsymToAsym(a_graph);
  TGraphAsymmErrors* returnGraph = HistogramUtilities::makeFitPullGraph(convGraph,a_funct,a_doIntegral);
  #ifndef _SAVE_POINTERS_
    delete convGraph;
  #endif
  return returnGraph;
}

TGraphAsymmErrors* HistogramUtilities::makeFitPullGraph(TGraphAsymmErrors* a_graph, TGraphAsymmErrors* a_graph2){
  TGraphAsymmErrors* pullGraph = new TGraphAsymmErrors();
  string newName = a_graph->GetName();
  newName = newName + "_pullToFit";
  pullGraph->SetName(newName.c_str());
  //pullGraph->GetYaxis()->SetTitle("Pulls");
  int numPointsCurrent = 0;
  for(int pointIndex = 0; pointIndex < a_graph->GetN(); pointIndex++){
    double center  = a_graph->GetX()[pointIndex];
    double lowEdge = center - a_graph->GetEXlow()[pointIndex];
    double highEdge = center - a_graph->GetEXhigh()[pointIndex];
    double content  = a_graph->GetY()[pointIndex];
    double errorYhigh = a_graph->GetEYhigh()[pointIndex];
    double errorYlow = a_graph->GetEYlow()[pointIndex];
    double content2  = a_graph2->GetY()[pointIndex];
    double errorYhigh2 = a_graph2->GetEYhigh()[pointIndex];
    double errorYlow2 = a_graph2->GetEYlow()[pointIndex];
    if(content > content2){
      double pull = (content - content2)/sqrt(errorYlow*errorYlow + errorYhigh2*errorYhigh2);
      if(!std::isfinite(pull)) pull = -999.999;
      pullGraph->SetPoint(numPointsCurrent,center,pull);
      pullGraph->SetPointError(numPointsCurrent,a_graph->GetEXlow()[pointIndex],a_graph->GetEXhigh()[pointIndex],1.0,1.0);
    }else{
      double pull = (content - content2)/sqrt(errorYhigh*errorYhigh + errorYlow2*errorYlow2);
      if(!std::isfinite(pull)) pull = -999.999;
      pullGraph->SetPoint(numPointsCurrent,center,pull);
      pullGraph->SetPointError(numPointsCurrent,a_graph->GetEXlow()[pointIndex],a_graph->GetEXhigh()[pointIndex],1.0,1.0);
    }
    numPointsCurrent += 1;
  }
  return pullGraph;
}



  TGraphErrors* HistogramUtilities::makeErrorBandGraph(TF1* a_funct, TMatrixDSym* a_covMatrix, int a_num_sample_points){
    if(!a_funct || !a_covMatrix || a_num_sample_points < 1){
      cout << "ERROR: HistogramUtilities::makeErrorBandGraph  check passed arguments: " << a_funct << "  " << a_covMatrix << "  " << a_num_sample_points << endl;
    }
    TGraphErrors* error_graph = new TGraphErrors();
    string name = a_funct->GetName();
    name = name + "_ErrorBandGraph";
    error_graph->SetName(name.c_str());
    error_graph->SetFillColor(kViolet);
    error_graph->SetFillStyle(3002);
    double lowRange, highRange;
    a_funct->GetRange(lowRange,highRange);
    double width = highRange - lowRange;
    double small_width = width/((double) a_num_sample_points);
    double half_small_width = 0.5*small_width;

    a_covMatrix->Print();

    for(int iii = 0; iii < a_num_sample_points; iii++){
      double x_low  = ((double) iii) * small_width;
      double x_high = x_low + small_width;
      double value = a_funct->Integral(x_low,x_high)/small_width;
      double error = a_funct->IntegralError(x_low, x_high,a_funct->GetParameters(), a_covMatrix->GetMatrixArray())/small_width;
      cout << " point: " << iii << "   " <<  x_low+half_small_width << "  = " << value << " +/- " << error << endl;
      error_graph->SetPoint(iii, x_low+half_small_width,value);
      error_graph->SetPointError(iii, half_small_width,error);
    }

    return error_graph;





  }




void HistogramUtilities::evalWithError(TGraphAsymmErrors* a_graph, double a_evalX, double &a_yVal, double &a_yErr){
  if( !a_graph || a_graph->GetN() < 2 ) return;
  a_graph->Sort();
  int point = TMath::BinarySearch(a_graph->GetN(),a_graph->GetX(),a_evalX);
  int firstPoint  = 0;
  int secondPoint = 0;
  if(point < 1 ){
    firstPoint  = 0;
    secondPoint = 1;
  }else{
    if( point >= a_graph->GetN() - 1){
      secondPoint = a_graph->GetN() - 1;
      firstPoint  = secondPoint - 1;
    }else{
      if(a_evalX < a_graph->GetX()[point]){
        firstPoint  = point - 1;
        secondPoint = point;
      }else{
        firstPoint  = point;
        secondPoint = point + 1;
      }
    }
  }

  double x1    = a_graph->GetX()[firstPoint];
  double y1    = a_graph->GetY()[firstPoint];
  double yErr1 = 0.5*(a_graph->GetEYlow()[firstPoint] + a_graph->GetEYhigh()[firstPoint]);
  double x2    = a_graph->GetX()[secondPoint];
  double y2    = a_graph->GetY()[secondPoint];
  double yErr2 = 0.5*(a_graph->GetEYlow()[secondPoint] + a_graph->GetEYhigh()[secondPoint]);
  PhysMath::interpolate(a_evalX, x1, y1, yErr1, x2, y2, yErr2, a_yVal, a_yErr);




  /*
  int low = -1;
  int up = -1;
  int nPoints = a_graph->GetN();
  double xLow = -999;
  double xUp = -999;
  double yLow = -999;
  double yUp = -999;

  double yErrLow = -999;
  double yErrUp = -999;

  a_graph->Sort();
  low = TMath::BinarySearch(nPoints,a_graph->GetX(),a_evalX);
  if(low == -1) low = 0;
  if(low == nPoints-1) low--;
  if(a_graph->GetX()[low] == a_evalX){
    yVal = a_graph->GetY()[low];
    yErr = (a_graph->GetEYhigh()[low] + a_graph->GetEYlow()[low])/2.0;
    return;
  }
  up = low + 1;
  if(a_graph->GetX()[up] == a_evalX){
    yVal = a_graph->GetY()[up];
    yErr = (a_graph->GetEYhigh()[up] + a_graph->GetEYlow()[up])/2.0;
    return;
  }
  if(a_graph->GetX()[low] == a_graph->GetX()[up]){
    yVal = a_graph->GetY()[low];
    yErr = (a_graph->GetEYhigh()[low] + a_graph->GetEYlow()[low])/2.0;
    return;
  }

  xLow = a_graph->GetX()[low];
  xUp = a_graph->GetX()[up];
  yLow = a_graph->GetY()[low];
  yUp = a_graph->GetY()[up];

  yErrLow = (a_graph->GetEYhigh()[low] + a_graph->GetEYlow()[low])/2.0;
  yErrUp = (a_graph->GetEYhigh()[up] + a_graph->GetEYlow()[up])/2.0;

  yVal = yUp + (a_evalX - xUp)*(yLow - yUp)/(xLow - xUp);
  yErr = sqrt(pow(1-(a_evalX - xUp)/(xLow - xUp),2)*pow(yErrUp,2) + pow((a_evalX - xUp)/(xLow - xUp),2)*pow(yErrLow,2));
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " HistogramUtilities::evalWithError  lowPt: " << low << " xLow: " << xLow << " yLow: " << yLow << " yErr: " <<  yErrLow << yLow
                                            << " highPt: " << up << " xUp: " << xUp << " yUp: " << yUp << " yErr: " << yErrUp << endl;
  #endif
  return;*/
}



void  HistogramUtilities::evalWithAsymmError(TGraphAsymmErrors* a_graph, double a_evalX, double &a_yVal, double &a_yErrLow, double &a_yErrHigh){
  if( !a_graph || a_graph->GetN() < 2 ) return;
  a_graph->Sort();
  int point = TMath::BinarySearch(a_graph->GetN(),a_graph->GetX(),a_evalX);
  int firstPoint  = 0;
  int secondPoint = 0;
  if(point < 1 ){
    firstPoint  = 0;
    secondPoint = 1;
  }else{
    if( point >= a_graph->GetN() - 1){
      secondPoint = a_graph->GetN() - 1;
      firstPoint  = secondPoint - 1;
    }else{
      if(a_evalX < a_graph->GetX()[point]){
        firstPoint  = point - 1;
        secondPoint = point;
      }else{
        firstPoint  = point;
        secondPoint = point + 1;
      }
    }
  }

  double x1    = a_graph->GetX()[firstPoint];
  double y1    = a_graph->GetY()[firstPoint];
  double yErrLow1  = a_graph->GetEYlow()[firstPoint];
  double yErrHigh1 = a_graph->GetEYhigh()[firstPoint];

  double x2    = a_graph->GetX()[secondPoint];
  double y2    = a_graph->GetY()[secondPoint];
  double yErrLow2  = a_graph->GetEYlow()[secondPoint];
  double yErrHigh2 = a_graph->GetEYlow()[secondPoint];

  PhysMath::interpolateAsymmError(a_evalX, x1, y1, yErrLow1, yErrHigh1, x2, y2, yErrLow2, yErrHigh2, a_yVal, a_yErrLow, a_yErrHigh);
}




void HistogramUtilities::functionEvalWithError(TF1* a_funct, double a_x, double &a_eval, double &a_evalErr, TMatrixDSym* a_covMatrix){
  #ifdef _HISTO_UTIL_DEBUG_
    cout << " Starting HistogramUtilities::functionEvalWithError" << endl;
  #endif
  double fracOfParErr = 1.0/20.0; // how far to sample from par value for d/dx calc
  double deltaFunctionSqrd = 0.0;
  int numPar = a_funct->GetNpar();

  vector<double> slopes;  //  df/da for


  for(int parIndex = 0; parIndex < numPar; parIndex++){
    slopes.push_back(0.0);
    double par = a_funct->GetParameter(parIndex);
    double parErr = a_funct->GetParError(parIndex);
    if(parErr <= 0.0) continue;
    double par_m = par - fracOfParErr*parErr;
    double par_p = par + fracOfParErr*parErr;

    a_eval = a_funct->Eval(a_x);
    a_funct->SetParameter(parIndex,par_m);
    double eval_m = a_funct->Eval(a_x);
    a_funct->SetParameter(parIndex,par_p);
    double eval_p = a_funct->Eval(a_x);
    a_funct->SetParameter(parIndex,par);   ///  SET IT BACK TO ORIGINAL PARAMETER

    double slope_m = (a_eval - eval_m)/(par - par_m);
    double slope_p = (eval_p - a_eval)/(par_p - par);
    double slope = 0.5*(slope_m + slope_p);
    slopes[parIndex] = slope;



    double contribution = pow(slope*parErr,2.0);

    #ifdef _HISTO_UTIL_DEBUG_
      cout << " Par " << parIndex << " val " << par << " parStep: " << fracOfParErr*parErr << endl
      << "     slope1: " << slope_m << " slope2: " << slope_p << " slope avg:" << slope << endl
      << "     contribution: " << contribution << endl;
    #endif

    deltaFunctionSqrd += contribution;
  }


  if(a_covMatrix){
    for(int parIndex = 0; parIndex < numPar; parIndex++){
      for(int parIndexTwo = parIndex + 1; parIndexTwo < numPar; parIndexTwo++){
        double contribution = 2.0*slopes[parIndex]*slopes[parIndexTwo]* a_covMatrix->operator()(parIndex,parIndexTwo);

        #ifdef _HISTO_UTIL_DEBUG_
          cout << "  " << parIndex << "  " << parIndexTwo << "     slope1: " << slopes[parIndex] << " slope2: " << slopes[parIndexTwo] << " contribution:" << contribution << endl;
        #endif

        deltaFunctionSqrd += contribution;
      }
    }

  }else{
    #ifdef _HISTO_UTIL_DEBUG_
      cout << "    No covariance matrix, not necessarily a problem" << endl;
    #endif
  }





  a_evalErr = sqrt(deltaFunctionSqrd);
  return;
}


void HistogramUtilities::printTF1Info(TF1* a_funct){
  cout <<"TF1 Info:" << endl;
  cout << "  Ptr:    \t" << a_funct << endl;
  if(a_funct){
    cout << "  Name:   \t" << a_funct->GetName() << endl;
    cout << "  Title:  \t" << a_funct->GetTitle() << endl;
    cout << "  Formula:\t" << a_funct->GetExpFormula() << endl;
    double x_low = 0; double x_high = 0;
    a_funct->GetRange(x_low,x_high);
    cout << "  Range:  \t[" << x_low << ", " << x_high << "]" << endl;
    cout << "  " << a_funct->GetNpar() << "Parameters:" << endl;
    for(int iii = 0; iii < a_funct->GetNpar(); iii++){
      double v_low = 0; double v_high = 0;
      a_funct->GetParLimits(iii,v_low,v_high);
      cout << iii << ": " << a_funct->GetParName(iii) << " = " << a_funct->GetParameter(iii)
          << "  limits: [" << v_low << ", " << v_high << "]   +/- " << a_funct->GetParError(iii) << endl;
    }
  }
}

bool HistogramUtilities::parametersAreFinite(TF1* a_funct){
  bool finiteParams = true;
  for(int paramIndex = 0; paramIndex < a_funct->GetNpar(); paramIndex++){
    if(!std::isfinite(a_funct->GetParameter(paramIndex))) finiteParams = false;
    if(!std::isfinite(a_funct->GetParError(paramIndex))) finiteParams = false;
  }
  return finiteParams;
}


TObject* HistogramUtilities::FileGet(TFile* a_file, string a_name){
  TObject* ptr = a_file->Get(a_name.c_str());
  if(!ptr){
    cout << "ERROR: " << a_name << " not read from file! Address: " << ptr << endl;
  }
  #ifdef _HISTO_UTIL_DEBUG_
  else{
    cout << "Successful reading " << a_name << " from file (" << a_file << ") yielded a pointer (" << ptr << ")" << endl;
  }
  #endif

  return ptr;
};


TPaveText* HistogramUtilities::makeTF1TextBox(TF1* a_funct, double a_x1, double a_y1, double a_x2, double a_y2){
  TPaveText* box = new TPaveText(a_x1, a_y1, a_x2, a_y2,"NDC");
  int numPar = a_funct->GetNpar();
  for(int paramIndex = 0; paramIndex < numPar; paramIndex++){
    double lowLimit;
    double highLimit;
    a_funct->GetParLimits(paramIndex,lowLimit,highLimit);
    double width = highLimit - lowLimit;
    bool atLimit = (fabs(a_funct->GetParameter(paramIndex) - lowLimit) < 0.01*width) || (fabs(highLimit - a_funct->GetParameter(paramIndex)) < 0.01*width);
    string qualifier = Form("[ %f, %f ]%s",lowLimit,highLimit,(atLimit ? " AT LIMIT! " : ""));
    if(highLimit == lowLimit){
      qualifier = "(fixed)";
    }
    box->AddText(Form("%s \t %f  %s",a_funct->GetParName(paramIndex),a_funct->GetParameter(paramIndex),qualifier.c_str()));
  }
  box->SetFillStyle(0);
  box->SetBorderSize(5);
  return box;
}



void HistogramUtilities::testRebinEqualByCDF(){
  srand (time(NULL));

  TH1D* histo = new TH1D("histoThing","HistogamUtilities:rebinEqualByCDF Test;x; dN/dx",200,0,200);
  TRandom3* randObj = new TRandom3();
  randObj->SetSeed();

  double counts_1 = 4000;
  double mean_1 = 75;
  double sigma_1 = 25;
  double counts_2 = 400;
  double mean_2 = 125;
  double sigma_2 = 15;

  for(int trialIndex = 0; trialIndex < (int) counts_1; trialIndex++){
    histo->Fill(randObj->Gaus(mean_1,sigma_1));
  }
  for(int trialIndex = 0; trialIndex < (int) counts_2; trialIndex++){
    histo->Fill(randObj->Gaus(mean_2,sigma_2));
  }

  makeSqrtNErrors(histo);
  histo->Scale(1.0,"width");
  histo->SetLineColor(kBlack);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(kBlack);
  histo->SetMarkerSize(0.8);

  TH1D* rebinHisto = rebinEqualByCDF(histo,25,true,false);
  rebinHisto->SetLineColor(kRed);
  rebinHisto->SetMarkerStyle(20);
  rebinHisto->SetMarkerColor(kRed);
  rebinHisto->SetMarkerSize(0.8);

  TH1D* rebinHistoLN = rebinEqualByCDF(histo,25,true,true);
  rebinHistoLN->SetLineColor(kBlue);
  rebinHistoLN->SetMarkerStyle(20);
  rebinHistoLN->SetMarkerColor(kBlue);
  rebinHistoLN->SetMarkerSize(0.8);

  TH1D* rebinHistoNoZeros = rebinNoZeroBins(histo, true);
  rebinHistoNoZeros->SetLineColor(kGreen+2);
  rebinHistoNoZeros->SetMarkerStyle(20);
  rebinHistoNoZeros->SetMarkerColor(kGreen+2);
  rebinHistoNoZeros->SetMarkerSize(0.8);


  TCanvas* canvas = new TCanvas("c1");
  histo->Draw();
  rebinHistoNoZeros->Draw("SAME");
  rebinHisto->Draw("same");
  rebinHistoLN->Draw("same");
  canvas->Update();


  double counts_A = sumOfContent(histo,-1,true);
  double counts_B = sumOfContent(rebinHisto,-1,true);
  double counts_C = sumOfContent(rebinHistoLN,-1,true);
  double counts_D = sumOfContent(rebinHistoNoZeros,-1,true);

  cout << "Gaus     Counts: " << counts_1 + counts_2 << endl;
  cout << "Original Counts: " << counts_A << endl;
  cout << "Rebin    Counts: " << counts_B << endl;
  cout << "RebinLN  Counts: " << counts_C << endl;
  cout << "RebinNoZeros Counts: " << counts_D << endl;



}



void HistogramUtilities::testRebinNoZeros(){
  srand (time(NULL));

  TH1D* histo = new TH1D("histoThing","HistogamUtilities Test;x; dN/dx",200,0,200);
  TRandom3* randObj = new TRandom3();
  randObj->SetSeed();

  double counts_1 = 30;
  double mean_1 = 100;
  double sigma_1 = 15;

  for(int trialIndex = 0; trialIndex < (int) counts_1; trialIndex++){
    histo->Fill(randObj->Gaus(mean_1,sigma_1));
  }

  setBinValues(histo,1,1.0,1.0);
  setBinValues(histo,200,1.0,1.0);

  makeSqrtNErrors(histo);
  histo->Scale(1.0,"width");
  histo->SetLineColor(kBlack);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(kBlack);
  histo->SetMarkerSize(0.8);

  TH1D* rebinHisto = rebinNoZeroBins(histo, true);
  rebinHisto->SetLineColor(kRed);
  rebinHisto->SetMarkerStyle(20);
  rebinHisto->SetMarkerColor(kRed);
  rebinHisto->SetMarkerSize(0.8);


  TCanvas* canvas = new TCanvas("c1");
  histo->Draw();
  rebinHisto->Draw("same");
  canvas->Update();


  double counts_A = sumOfContent(histo,-1,true);
  double counts_B = sumOfContent(rebinHisto,-1,true);

  cout << "Gaus     Counts: " << counts_1 + 2 << endl;
  cout << "Original Counts: " << counts_A << endl;
  cout << "Rebin    Counts: " << counts_B << endl;



}



void HistogramUtilities::test_functionEvalWithError(){

  TF1* funct = new TF1("sinFunct","[0]*sin([1]*x+[2])",-15,15);
  funct->SetParameter(0,1.0);
  funct->SetParameter(1,1.0);
  funct->SetParameter(2,0.0);

  funct->SetParError(0,0.01);
  funct->SetParError(1,0.01);
  funct->SetParError(2,0.01);

  TGraphErrors* graph = new TGraphErrors();

  double numPoints = 200;
  double step = 20.0/numPoints;

  for(int index = 0; index < numPoints; index++){
    double x = -10.0 + step*((double) index);
    double y = 0.0;
    double yErr = 0.0;
    functionEvalWithError(funct,x, y, yErr);
    graph->SetPoint(index,x,y);
    graph->SetPointError(index, 0.0,yErr);
  }



  TCanvas* canv = new TCanvas();
  funct->Draw();
  graph->Draw("sameP");



}

void HistogramUtilities::test_rebinHistoToAlternateBinStructure(){

  TH1D* histo_even = new TH1D("histo_even","Counts; x; N",20,0,20);
  for(int bin = 1; bin <= histo_even->GetNbinsX(); bin++){
    double content = 100.0 - 10.0*((double) abs(10 - bin));
    histo_even->SetBinContent(bin,content);
    histo_even->SetBinError(bin,sqrt(content));
  }

  vector<double> newEdges;
  newEdges.push_back(-0.5);
  newEdges.push_back(1.5);
  newEdges.push_back(4.8);
  newEdges.push_back(6.5);
  newEdges.push_back(7.0);
  newEdges.push_back(13.0);
  newEdges.push_back(15.2);
  newEdges.push_back(16.7);
  newEdges.push_back(20.0);


  TH1D* histo_unEven = HistogramUtilities::rebinHistoToAlternateBinStructure(histo_even, newEdges);
  histo_even->SetLineColor(kRed);
  histo_unEven->SetLineColor(kBlue);

  TH1D* histo_even_diff   = (TH1D*) histo_even->Clone();
  TH1D* histo_unEven_diff = (TH1D*) histo_unEven->Clone();
  histo_even_diff->Scale(1.0,"width");
  histo_unEven_diff->Scale(1.0,"width");
  histo_even_diff->GetYaxis()->SetTitle("dN/dx");
  histo_unEven_diff->GetYaxis()->SetTitle("dN/dx");

  TCanvas* canv = new TCanvas();
  histo_even->Draw("E");
  histo_unEven->Draw("E same");

  TCanvas* canv_diff = new TCanvas();
  histo_even_diff->Draw("E");
  histo_unEven_diff->Draw("E same");

  cout << "Sum of inital histo: " << HistogramUtilities::sumOfContent(histo_even) << endl;
  cout << "Sum of new histo: " << HistogramUtilities::sumOfContent(histo_unEven) << endl;
  return;
}


void HistogramUtilities::test_rebinHistoToChebyshevSpacing(){
  TH1D* histo_even = new TH1D("histo_even","Counts; x; N",20,0,20);
  for(int bin = 1; bin <= histo_even->GetNbinsX(); bin++){
    double content = 100.0 - 10.0*((double) abs(10 - bin));
    histo_even->SetBinContent(bin,content);
    histo_even->SetBinError(bin,sqrt(content));
  }

  TH1D* histo_cheb = HistogramUtilities::rebinHistoToChebyshevSpacing(histo_even,3,18,1);
  histo_even->SetLineColor(kRed);
  histo_cheb->SetLineColor(kBlue);
  histo_even->SetLineWidth(2);
  histo_cheb->SetLineWidth(2);
  histo_even->Scale(1.0,"width");
  histo_cheb->Scale(1.0,"width");
  histo_even->GetYaxis()->SetTitle("dN/dx");
  histo_cheb->GetYaxis()->SetTitle("dN/dx");

  TCanvas* canv = new TCanvas();
  histo_even->Draw("E");
  histo_cheb->Draw("E same");
  return;
}

void HistogramUtilities::test_rebin2DHistoToVertChebyshevSpacing(){
  TH2D* histo_even = new TH2D("histo_even","Counts; x;y;N",20,0,20,20,0,20);
  for(int bin = 1; bin <= histo_even->GetNbinsX(); bin++){
    for(int binY = 1; binY <= histo_even->GetNbinsX(); binY++){
      double content = 100.0 - 10.0*((double) abs(10 - bin) - 10.0*((double) abs(10 - binY)));
      histo_even->SetBinContent(bin,binY,content);
      histo_even->SetBinError(bin,binY,sqrt(content));
    }
  }

  TH2D* histo_cheb = HistogramUtilities::rebin2DHistoToVertChebyshevSpacing(histo_even,3,18,1);
  histo_even->SetLineColor(kRed);
  histo_cheb->SetLineColor(kBlue);
  histo_even->SetLineWidth(2);
  histo_cheb->SetLineWidth(2);
  makeDifferential(histo_even);
  makeDifferential(histo_cheb);
  histo_even->GetZaxis()->SetTitle("d^{2}N/dxdy");
  histo_cheb->GetZaxis()->SetTitle("d^{2}N/dxdy");

  TCanvas* canv = new TCanvas();
  histo_even->Draw("colz");
  TCanvas* canv2 = new TCanvas();
  histo_cheb->Draw("colz");
  return;
}



void HistogramUtilities::test_histoAddWithBinSplitting(){
  TH1D* histo = new TH1D("histo","Counts; x; dN/dx",20,0,20);
  for(int bin = 1; bin <= histo->GetNbinsX(); bin++){
    histo->SetBinContent(bin,bin);
    histo->SetBinError(bin,sqrt(bin));
  }
  histo->SetLineColor(kBlue);

  double binEdges[10] = {0.0,1.1,3.2,5.5,5.9,6,8.8,15,17.5,20};
  TH1D* histo_offset = new TH1D("histo_offset","Counts; x; dN/dx",9,binEdges);
  addToHistogram_binSplitting(histo_offset,histo);
  histo->Scale(1.0,"width");
  histo_offset->Scale(1.0,"width");
  histo_offset->SetLineColor(kRed);


  TCanvas* canv = new TCanvas();
  histo->Draw("E");
  histo_offset->Draw("E SAME");

}

void HistogramUtilities::test_m2_invBeta_conversion(){
  srand (time(NULL));
  TRandom3* randObj = new TRandom3();
  randObj->SetSeed();
  double counts_1 = 4000;
  double mean_1 = 0.0;//1.18;
  double sigma_1 = 0.0133;

  TH1D* histo = new TH1D("invBetaHisto","Original ZTOF; 1/#beta - 1/#beta_{predict}; dN/d(1/#beta)",200,-0.05,0.05);
  TH1D* histo2 = new TH1D("invBetaHisto","Original 1/#beta; 1/#beta ; dN/d(1/#beta)",200,1.05,1.25);

  for(int trialIndex = 0; trialIndex < (int) counts_1; trialIndex++){
    double point = randObj->Gaus(mean_1,sigma_1);
    histo->Fill(point);
    histo2->Fill(point + 1.18);
  }

  makeSqrtNErrors(histo);
  histo->Scale(1.0,"width");
  histo->SetLineColor(kBlack);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(kBlack);
  histo->SetMarkerSize(0.8);

  makeSqrtNErrors(histo2);
  histo2->Scale(1.0,"width");
  histo2->SetLineColor(kViolet+1);
  histo2->SetMarkerStyle(20);
  histo2->SetMarkerColor(kViolet+1);
  histo2->SetMarkerSize(0.8);

  TH1D* histo_m2 = convertZTOFToMassSqr(histo, 1.5, 0.938, true, true);
  histo_m2->SetLineColor(kGreen+2);
  TH1D* histo_invBeta = convertMassSqrToInvBeta(histo_m2, 1.5,true, true, 0.938);
  histo_invBeta->SetLineColor(kRed);

  TH1D* histo_ZTOF = convertMassSqrToInvBeta(histo_m2, 1.5,true, true);
  histo_ZTOF->SetLineColor(kRed);


  histo->Draw("E");
  histo_invBeta->Draw("SAME");

  TCanvas* canv2 = new TCanvas();
  histo_m2->Draw("E");

  TCanvas* canv3 = new TCanvas();
  histo2->Draw("E");
  histo_ZTOF->Draw("SAME");

}



void HistogramUtilities::test_self_convolution(){

  srand (time(NULL));
  TRandom3* randObj = new TRandom3();
  randObj->SetSeed();
  int counts_1 = 100000;
  double mean_1 = 1.0;
  double sigma_1 = 1.0;

  TH1D* histo = new TH1D("base_histo","Gaussian;x;1/N dN/dx",200,-5,5);
  TH1D* histo_2 = new TH1D("base_histo_2","Gaussian;x;1/N dN/dx",130,-3,8);

  for(int trialIndex = 0; trialIndex < counts_1; trialIndex++){
    double point = randObj->Gaus(mean_1,sigma_1);
    histo->Fill(point);
    histo_2->Fill(point);
  }
  makeSqrtNErrors(histo);
  histo->Scale(1.0/histo->Integral(),"width");
  histo->SetLineColor(kBlack);
  histo->SetMarkerStyle(20);
  histo->SetMarkerColor(kBlack);
  histo->SetMarkerSize(0.8);

  makeSqrtNErrors(histo_2);
  histo_2->Scale(1.0/histo_2->Integral(),"width");
  histo_2->SetLineColor(kGreen+1);
  histo_2->SetMarkerStyle(20);
  histo_2->SetMarkerColor(kGreen+1);
  histo_2->SetMarkerSize(0.8);


  TH1D* histo_self_conv =  selfConvolution(histo); // must have equal bin widths
  histo_self_conv->SetLineColor(kRed+1);
  histo_self_conv->SetMarkerStyle(20);
  histo_self_conv->SetMarkerColor(kRed+1);
  histo_self_conv->SetMarkerSize(0.8);

  //TH1D* histo_conv =  convolution(histo,histo_2); // must have equal bin widths
  //histo_conv->SetLineColor(kViolet+1);
  //histo_conv->SetMarkerStyle(20);
  //histo_conv->SetMarkerColor(kViolet+1);
  //histo_conv->SetMarkerSize(0.8);

  //histo_conv->Scale(1.0/histo_conv->Integral(),"width");
  histo_self_conv->Scale(1.0/histo_self_conv->Integral(),"width");
  histo_self_conv->GetYaxis()->SetRangeUser(0,0.45);
  histo_self_conv->Draw("E");
  histo->Draw("SAME");
  //histo_2->Draw("SAME");
  //histo_conv->Draw("SAME");


}


void HistogramUtilities::test_convolution(){

  srand (time(NULL));
  TRandom3* randObj = new TRandom3();
  randObj->SetSeed();
  int counts_1 = 100000;
  double mean_1 = 1.0;
  double sigma_1 = 1.0;

  TH1D* histo_gaus = new TH1D("base_histo","Gaussian;x;1/N dN/dx",200,-5,5);
  TH1D* histo_uniform = new TH1D("base_histo_2","Gaussian;x;1/N dN/dx",130,-3,8);

  for(int trialIndex = 0; trialIndex < counts_1; trialIndex++){
    double point = randObj->Gaus(mean_1,sigma_1);
    histo_gaus->Fill(point);
  }

  for(int trialIndex = 0; trialIndex < counts_1; trialIndex++){
    double point = randObj->Uniform(-4,6);
    histo_uniform->Fill(point);
  }

  makeSqrtNErrors(histo_gaus);
  histo_gaus->Scale(1.0/histo_gaus->Integral(),"width");
  histo_gaus->SetLineColor(kBlack);
  histo_gaus->SetMarkerStyle(20);
  histo_gaus->SetMarkerColor(kBlack);
  histo_gaus->SetMarkerSize(0.8);

  makeSqrtNErrors(histo_uniform);
  histo_uniform->Scale(1.0/histo_uniform->Integral(),"width");
  histo_uniform->SetLineColor(kGreen+1);
  histo_uniform->SetMarkerStyle(20);
  histo_uniform->SetMarkerColor(kGreen+1);
  histo_uniform->SetMarkerSize(0.8);


  TH1D* histo_conv =  convolution(histo_gaus,histo_uniform); // must have equal bin widths
  histo_conv->SetLineColor(kRed+1);
  histo_conv->SetMarkerStyle(20);
  histo_conv->SetMarkerColor(kRed+1);
  histo_conv->SetMarkerSize(0.8);

  TH1D* histo_conv_2 =  convolution(histo_uniform,histo_gaus); // must have equal bin widths
  histo_conv_2->SetLineColor(kViolet+1);
  histo_conv_2->SetMarkerStyle(20);
  histo_conv_2->SetMarkerColor(kViolet+1);
  histo_conv_2->SetMarkerSize(0.8);


  histo_conv->Scale(1.0/histo_conv->Integral(),"width");
  histo_conv->GetYaxis()->SetRangeUser(0,0.45);
  histo_conv->Draw("E");
  histo_conv->Draw("SAME");
  histo_gaus->Draw("SAME");
  histo_uniform->Draw("SAME");


}


#endif
