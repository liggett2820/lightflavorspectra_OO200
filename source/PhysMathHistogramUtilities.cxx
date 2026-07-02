#ifndef _PHYSMATHHISTOGRAMUTILITIES_CXX_
#define _PHYSMATHHISTOGRAMUTILITIES_CXX_

#include "../headers/PhysMathHistogramUtilities.h"


void PhysMathHistogramUtilities::initializeLegendreExpansion(TF1* a_funct, TGraphAsymmErrors* a_graph, double a_lowRange, double a_highRange, bool a_fixRange, int a_numIntegrationSteps){

  int order = a_funct->GetNpar()-2;
  if(order < 0){
    cout << "ERROR: PhysMathHistogramUtilities::initializeLegendreExpansion: you requested an order < 0. I don't know what to do."<< endl;
    return;
  }

  int indexForLowRange = order;
  int indexForHighRange = order+1;

  int numSteps = a_graph->GetN();
  if(a_numIntegrationSteps > 0) numSteps = a_numIntegrationSteps;

  //Determine Function Range
  double lowRange = a_funct->GetParameter(indexForLowRange);
  double highRange = a_funct->GetParameter(indexForHighRange);
  if(a_lowRange != -999.999){
    lowRange = a_lowRange;
    a_funct->SetParameter(indexForLowRange,lowRange);
  }
  if(a_highRange != -999.999){
    highRange = a_highRange;
    a_funct->SetParameter(indexForHighRange,highRange);
  }

  // Allow the range to vary. This lowers edge effects from taylor expansion truncation.
  double smallestX = HistogramUtilities::findSmallestXVal(a_graph,a_lowRange,a_highRange);
  double largestX = HistogramUtilities::findLargestXVal(a_graph,a_lowRange,a_highRange);
  double deltaX = largestX - smallestX;
  if(a_fixRange){
    a_funct->FixParameter(indexForLowRange,lowRange);
    a_funct->FixParameter(indexForHighRange,highRange);
  }else{
    a_funct->SetParLimits(indexForLowRange, lowRange - 1.5*deltaX,lowRange);
    a_funct->SetParLimits(indexForHighRange,highRange, highRange + 1.5*deltaX);
  }

  a_funct->SetRange(lowRange,highRange);
  double lowPlusHighRange = lowRange + highRange;
  double highMinusLowRange = highRange - lowRange;
  double smallestY = HistogramUtilities::findSmallestYVal(a_graph,a_lowRange,a_highRange);
  double largestY = HistogramUtilities::findLargestYVal(a_graph,a_lowRange,a_highRange);
  double deltaY = largestY - smallestY;

  //#### Determine coefficients
  double integrationStep = 2.0/((double) numSteps);
  for(int orderIndex = 0; orderIndex < order; orderIndex++){

    //First way...
    /*double A_i = 0.0;
    double previousEval = PhysMath::legendre(orderIndex,lowRange)*a_graph->Eval(lowRange);
    for(int integrationIndex = 0; integrationIndex < numSteps - 1; integrationIndex++){
      double highEdge = -1.0 + ((double) (integrationIndex + 1))*integrationStep;
      double highEdgeEval = PhysMath::legendre(orderIndex,highEdge)*a_graph->Eval(0.5*(highEdge*highMinusLowRange + lowPlusHighRange));
      A_i += previousEval + highEdgeEval;
      previousEval = highEdgeEval;
    }
    A_i = (((double)orderIndex)+0.5)*0.5*integrationStep*A_i; //0.5 from trapazoidal approx and prefactor of delta function
    a_funct->SetParameter(orderIndex,A_i);
    */

    //Second Way...
    double A_i = 0.5*PhysMath::legendre(orderIndex,-1.0)*a_graph->Eval(lowRange);
    for(int integrationIndex = 1; integrationIndex < numSteps; integrationIndex++){
      double sample_x = -1.0 + ((double) (integrationIndex))*integrationStep;
      A_i += PhysMath::legendre(orderIndex,sample_x) * a_graph->Eval(0.5*(sample_x*highMinusLowRange + lowPlusHighRange));
    }
    A_i += 0.5*PhysMath::legendre(orderIndex,1.0)*a_graph->Eval(highRange);
    A_i = (((double)orderIndex)+0.5)*integrationStep*A_i; //prefactor of delta function  Int(PmPm) = 1/(m+0.5)
    a_funct->SetParameter(orderIndex,A_i);




    if(orderIndex == 0){
      a_funct->SetParLimits(0,smallestY - 2.0*deltaY,largestY + 2.0*deltaY);
    }else{
      if(A_i > 0.0){
        if(A_i < deltaY ){
           a_funct->SetParLimits(orderIndex,-1000.0*deltaY,1000.0*deltaY);
        }else{
           a_funct->SetParLimits(orderIndex,-1000.0*A_i,1000.0*A_i);
        }
      }else if(A_i < 0.0){
        if(fabs(A_i) < deltaY ){
           a_funct->SetParLimits(orderIndex,-1000.0*deltaY,1000.0*deltaY);
        }else{
           a_funct->SetParLimits(orderIndex,1000.0*A_i,-1000.0*A_i);
        }
      }else{
        a_funct->SetParLimits(orderIndex,-1000.0*deltaY,1000.0*deltaY);
        cout << " ERROR: Weird initialization of legendre polinomial coefficient for order " << orderIndex << ". The initialization in zero, so the parameter limits are difficult to set." << endl;
      }
    }
  }

}


TF1* PhysMathHistogramUtilities::determineGoodLegendreOrder(TGraphAsymmErrors* a_graph, double a_lowRange, double a_highRange , int a_lowOrderRange, int a_highOrderRange,  int a_leaveOutEvery){
  TVirtualFitter::SetMaxIterations(1000000);

  double lowRange = HistogramUtilities::findSmallestXVal(a_graph);
  double highRange =  HistogramUtilities::findLargestXVal(a_graph);
  if(a_lowRange != -999.999){
    lowRange = a_lowRange;
  }
  if(a_highRange != -999.999){
    highRange = a_highRange;
  }

  TF1* legendreFunct[30];
  for(int orderIndex = 0; orderIndex < 30; ++orderIndex){
    if( orderIndex < a_lowOrderRange || orderIndex > a_highOrderRange){
      legendreFunct[orderIndex] = NULL;
    }else{
      legendreFunct[orderIndex] = PhysMath::getLegendreExpansion(Form("legendreFunct_%d",orderIndex), orderIndex,false,lowRange, highRange);
    }
  }

  int numTrialsForEachOrder = 3;
  double chiSqrdPartial[30][3]; // same as above number
  double chiSqrdFull[30][3];
  double chiSqrd_Ndf_Partial[30][3];
  double chiSqrd_Ndf_Full[30][3];
  int    ndf_Partial[30][3];
  int    ndf_Full[30][3];

  for(int trialIndex = 0; trialIndex < numTrialsForEachOrder; trialIndex++){

    // Randomly remove points to estimate amount of overfitting when they are put back into chi^2
    TGraphAsymmErrors* graph_with_holes = (TGraphAsymmErrors*) a_graph->Clone();
    TRandom3* randObj = new TRandom3();
    randObj->SetSeed();
    int numPoints = a_graph->GetN();
    int numPoints_deleted = 0;
    double fracExclude = 1.0/((double) a_leaveOutEvery);
    for(int pointIndex = numPoints-1; pointIndex > 0; pointIndex--){
      if(randObj->Uniform() < fracExclude){
        numPoints_deleted += 1;
        graph_with_holes->RemovePoint(pointIndex);
      }
    }
    int numPoints_with_holes = numPoints - numPoints_deleted;

    #ifdef _PHYS_MATH_DEBUG_
      cout << " #### NumPoints: " << numPoints << "   NewNumPoints: " << numPoints_with_holes << "   so " << numPoints_deleted << "   deleted" << endl;
    #endif

    numPoints = HistogramUtilities::getNPointsInRange(a_graph,lowRange,highRange);
    numPoints_with_holes = HistogramUtilities::getNPointsInRange(graph_with_holes,lowRange,highRange);

    for(int orderIndex = 0; orderIndex < 30; orderIndex++){
      if(!legendreFunct[orderIndex]){
        chiSqrdPartial[orderIndex][trialIndex] = -999;
        chiSqrdFull[orderIndex][trialIndex] = -999;
        chiSqrd_Ndf_Partial[orderIndex][trialIndex] = -999;
        chiSqrd_Ndf_Full[orderIndex][trialIndex] = -999;
        ndf_Partial[orderIndex][trialIndex] = -999;
        ndf_Full[orderIndex][trialIndex] = -999;

        continue;
      }
      PhysMathHistogramUtilities::initializeLegendreExpansion(legendreFunct[orderIndex], a_graph, lowRange, highRange,false);
      int status = graph_with_holes->Fit(legendreFunct[orderIndex],"","EX0",lowRange,highRange);
      TF1* fitResut = graph_with_holes->GetFunction(Form("legendreFunct_%d",orderIndex));
      int ndf_partial = numPoints_with_holes - orderIndex - 2;
      int ndf_full = numPoints - orderIndex - 2;

      ndf_Partial[orderIndex][trialIndex] = ndf_partial;
      ndf_Full[orderIndex][trialIndex] = ndf_full;
      chiSqrdFull[orderIndex][trialIndex] = a_graph->Chisquare(fitResut,"R EX0");
      chiSqrdPartial[orderIndex][trialIndex] = graph_with_holes->Chisquare(fitResut,"R EX0");
      chiSqrd_Ndf_Partial[orderIndex][trialIndex] = chiSqrdPartial[orderIndex][trialIndex]/((double) ndf_partial);
      chiSqrd_Ndf_Full[orderIndex][trialIndex] = chiSqrdFull[orderIndex][trialIndex] /((double) ndf_full);

      #ifdef _PHYS_MATH_DEBUG_
        cout << "   Order: " << orderIndex << "  Fit Status: " << status << "   chi^2 partial: " <<  chiSqrdPartial[orderIndex][trialIndex] << "   chi^2 full: " << chiSqrdFull[orderIndex][trialIndex] << endl;
      #endif
    }
  }

  double avg_chiSqrdPartial[30];
  double avg_chiSqrdFull[30];
  double avg_chiSqrd_Ndf_Partial[30];
  double avg_chiSqrd_Ndf_Full[30];

  for(int orderIndex = 0; orderIndex < 30; orderIndex++){
    if(legendreFunct[orderIndex]){
      for(int trialIndex = 0; trialIndex < numTrialsForEachOrder; trialIndex++){
        avg_chiSqrdPartial[orderIndex] += chiSqrdPartial[orderIndex][trialIndex];
        avg_chiSqrdFull[orderIndex] += chiSqrdFull[orderIndex][trialIndex];
        avg_chiSqrd_Ndf_Partial[orderIndex] += chiSqrd_Ndf_Partial[orderIndex][trialIndex];
        avg_chiSqrd_Ndf_Full[orderIndex] += chiSqrd_Ndf_Full[orderIndex][trialIndex];
      }
        avg_chiSqrdPartial[orderIndex] = avg_chiSqrdPartial[orderIndex] / ((double) numTrialsForEachOrder);
        avg_chiSqrdFull[orderIndex] = avg_chiSqrdFull[orderIndex] / ((double) numTrialsForEachOrder);
        avg_chiSqrd_Ndf_Partial[orderIndex] = avg_chiSqrd_Ndf_Partial[orderIndex] / ((double) numTrialsForEachOrder);
        avg_chiSqrd_Ndf_Full[orderIndex] = avg_chiSqrd_Ndf_Full[orderIndex] / ((double) numTrialsForEachOrder);

    }
  }

  #ifdef _PHYS_MATH_DEBUG_
    TGraph* chi2Graph_full = new TGraph();
    TGraph* chi2Graph_partial = new TGraph();
    TGraph* chi2Graph_ndf_full = new TGraph();
    TGraph* chi2Graph_ndf_partial = new TGraph();
    chi2Graph_full->SetMarkerStyle(20);
    chi2Graph_partial->SetMarkerStyle(20);
    chi2Graph_full->SetLineColor(kRed);
    chi2Graph_partial->SetLineColor(kBlue);
    chi2Graph_full->SetMarkerColor(kRed);
    chi2Graph_partial->SetMarkerColor(kBlue);

    chi2Graph_ndf_full->SetMarkerStyle(20);
    chi2Graph_ndf_partial->SetMarkerStyle(20);
    chi2Graph_ndf_full->SetLineColor(kRed);
    chi2Graph_ndf_partial->SetLineColor(kBlue);
    chi2Graph_ndf_full->SetMarkerColor(kRed);
    chi2Graph_ndf_partial->SetMarkerColor(kBlue);



    for(int orderIndex = 0; orderIndex < 30; orderIndex++){
      if(!legendreFunct[orderIndex]) continue;
      chi2Graph_full->SetPoint(chi2Graph_full->GetN(),orderIndex,avg_chiSqrdFull[orderIndex]);
      chi2Graph_partial->SetPoint(chi2Graph_partial->GetN(),orderIndex,avg_chiSqrdPartial[orderIndex]);
      chi2Graph_ndf_full->SetPoint(chi2Graph_ndf_full->GetN(),orderIndex,avg_chiSqrd_Ndf_Full[orderIndex]);
      chi2Graph_ndf_partial->SetPoint(chi2Graph_ndf_partial->GetN(),orderIndex,avg_chiSqrd_Ndf_Partial[orderIndex]);
    }

    TCanvas* canv = new TCanvas();
    chi2Graph_full->Draw("APL");
    chi2Graph_partial->Draw("PL same");

    gPad->Update();

    TCanvas* canv2 = new TCanvas();
    chi2Graph_ndf_full->Draw("APL");
    chi2Graph_ndf_partial->Draw("PL same");


    gSystem->ProcessEvents();

    //delete chi2Graph_full;
    //delete chi2Graph_partial;

  #endif

  int bestOrder = a_lowOrderRange;
  for(int orderIndex = a_lowOrderRange; orderIndex <= a_highOrderRange; orderIndex++){
    if(avg_chiSqrd_Ndf_Full[orderIndex] < avg_chiSqrd_Ndf_Full[bestOrder]) bestOrder = orderIndex;
  }

  for(int orderIndex = 0; orderIndex < 30; orderIndex++){
    if(!legendreFunct[orderIndex]) continue;
    if(orderIndex != bestOrder) delete legendreFunct[orderIndex];
  }


  //the previous fit used the graph with holes
  int status = a_graph->Fit(legendreFunct[bestOrder],"","EX0",lowRange,highRange);
  TF1* bestFitResult = a_graph->GetFunction(Form("legendreFunct_%d",bestOrder));

  return bestFitResult;




}




#ifdef _MAC_OSX_

void PhysMathHistogramUtilities::testLambdaDecay(){

  double y1 = 0;
  double y2 = 0;
  double m1 = 0;
  double m2 = 0;
  double y3 = 0;
  double m3 = 0;

  TH2D* proton_y_mTm0_histo = new TH2D("proton_y_mTm0_histo","Proton Phase Space Probability; y ;m_{T}-m_{p}", 300, 0.35,0.65,300,0,0.2);
  TH2D* pion_y_mTm0_histo = new TH2D("pion_y_mTm0_histo","Pion Phase Space Probability; y ;m_{T}-m_{#pi}", 300, 0,1.4,300,0,0.2);
  double lambda_rapidity = 0.5;
  double lambda_mTm0 = 0.1;
  TRandom3 randomThing;
  randomThing.SetSeed(time(NULL));
  int numTrials = 100000;
  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    PhysMath::lambdaDecay(&randomThing, lambda_rapidity , lambda_mTm0, y1, m1, y2, m2);
    proton_y_mTm0_histo->Fill(y1,m1);
    pion_y_mTm0_histo->Fill(y2,m2);
  }

  HistogramUtilities::makeDifferential(proton_y_mTm0_histo, 1.0/((double) numTrials));
  HistogramUtilities::makeDifferential(pion_y_mTm0_histo, 1.0/((double) numTrials));

  TCanvas* canv = new TCanvas();
  proton_y_mTm0_histo->Draw("colz");
  TCanvas* canv2 = new TCanvas();
  pion_y_mTm0_histo->Draw("colz");

  return;
}


void PhysMathHistogramUtilities::testNewLambdaDecay(){

  TH2D* proton_y_mTm0_histo = new TH2D("proton_y_mTm0_histo","Proton Phase Space Probability; y ;m_{T}-m_{p} [GeV/c^{2}]", 300, 0.35,0.65,300,0,0.3);
  TH2D* pion_y_mTm0_histo   = new TH2D("pion_y_mTm0_histo","Pion Phase Space Probability; y ;m_{T}-m_{#pi} [GeV/c^{2}]", 300, -0.2,1.4,300,0,0.15);
  TH1D* proton_dca          = new TH1D("proton_dca","DCA of Proton from Lambda Decay;DCA [cm];Counts", 300, 0, 30);
  proton_dca->SetLineColor(kGreen+1);
  TH1D* pion_dca = new TH1D("pion_dca","DCA of Pion from Lambda Decay;DCA [cm];Counts", 300, 0, 30);
  pion_dca->SetLineColor(kRed+1);
  TH1D* lambda_decay_length = new TH1D("lambda_decay_length","Lambda Decay Length; Decay Length [cm]; Counts", 300, 0, 21);
  TH1D* lambda_pX = new TH1D("lambda_pX","Lambda p_x;p_x", 200, -1, 1);
  TH1D* lambda_pY = new TH1D("lambda_pY","Lambda p_y;p_y", 200, -1, 1);
  TH1D* lambda_pZ = new TH1D("lambda_pZ","Lambda p_z;p_z", 200, -1, 1);
  TH1D* proton_pX = new TH1D("proton_pX","proton p_x;p_x", 200, -1, 1);
  TH1D* proton_pY = new TH1D("proton_pY","proton p_y;p_y", 200, -1, 1);
  TH1D* proton_pZ = new TH1D("proton_pZ","proton p_z;p_z", 200, -1, 1);
  TH1D* secondVertex_X = new TH1D("secondVertex_pX","vertex_x;vertex_x", 200, -10, 10);
  TH1D* secondVertex_Y = new TH1D("secondVertex_pY","vertex_y;vertex_y", 200, -10, 10);
  TH1D* secondVertex_Z = new TH1D("secondVertex_pZ","vertex_z;vertex_z", 200, -10, 10);
  TH1D* proton_dcaX = new TH1D("proton_dcaX","DCA of Proton from Lambda Decay;DCA X [cm];Counts", 300, 0, 30);
  TH1D* proton_dcaY = new TH1D("proton_dcaY","DCA of Proton from Lambda Decay;DCA Y [cm];Counts", 300, 0, 30);
  TH1D* proton_dcaZ = new TH1D("proton_dcaZ","DCA of Proton from Lambda Decay;DCA Z [cm];Counts", 300, 0, 30);
  lambda_pX->SetLineColor(kBlue);
  lambda_pY->SetLineColor(kRed);
  lambda_pZ->SetLineColor(kViolet);
  proton_pX->SetLineColor(kBlue);
  proton_pY->SetLineColor(kRed);
  proton_pZ->SetLineColor(kViolet);
  proton_dcaX->SetLineColor(kBlue);
  proton_dcaY->SetLineColor(kRed);
  proton_dcaZ->SetLineColor(kViolet);
  secondVertex_X->SetLineColor(kBlue);
  secondVertex_Y->SetLineColor(kRed);
  secondVertex_Z->SetLineColor(kViolet);

  double lambda_rapidity = 0.5;//0.5;
  double lambda_mTm0 = 0.1;//0.1;
  TRandom3 randomThing;
  randomThing.SetSeed(time(NULL));
  ROOT::Math::Boost boostObj;
  int numTrials = 1000000;
  //int numTrials = 10000;
  ROOT::Math::XYZVector primaryVertex(0,0,0);
  double lambda_mass = 1.115683;
  double pT   = PhysMath::pT_from_mTm0(lambda_mTm0,lambda_mass);
  double pX = 0;
  double pY = 0;
  randomThing.Circle(pX,pY,pT);
  double pZ     = PhysMath::pZ_from_y_mTm0(lambda_rapidity,   lambda_mTm0,lambda_mass);
  double pTot   = PhysMath::pTot_from_y_mTm0(lambda_rapidity, lambda_mTm0,lambda_mass);
  double energy = sqrt(lambda_mass*lambda_mass + pTot*pTot);
  ROOT::Math::PxPyPzEVector lambda_vect(pX,pY,pZ,energy);
  ROOT::Math::PxPyPzEVector* proton_vect = NULL;
  ROOT::Math::PxPyPzEVector* pion_vect   = NULL;
  ROOT::Math::XYZVector* secondaryVertex = NULL;

  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    double pX = 0;
    double pY = 0;
    randomThing.Circle(pX,pY,pT);
    lambda_vect.SetPx(pX);
    lambda_vect.SetPy(pY);
    lambda_vect.SetPz(pZ);
    lambda_pX->Fill(lambda_vect.Px());
    lambda_pY->Fill(lambda_vect.Py());
    lambda_pZ->Fill(lambda_vect.Pz());
    if(trialIndex%100 == 0){
      cout << "####Trial " << trialIndex << "  Lambda Check: m = " << lambda_mass << " =? " << sqrt(lambda_vect.E()*lambda_vect.E() - lambda_vect.P()*lambda_vect.P()) <<  endl;
      cout << "    pT = " << sqrt(pX*pX + pY*pY) << "=?" << pT << endl;
    }



    PhysMath::lambdaDecay(&boostObj,&randomThing, &primaryVertex, &lambda_vect, proton_vect, pion_vect, secondaryVertex);
    double displacement[3] = {primaryVertex.X() - secondaryVertex->X(), primaryVertex.Y() - secondaryVertex->Y(), primaryVertex.Z() - secondaryVertex->Z()};
    double distance = sqrt(displacement[0]*displacement[0] + displacement[1]*displacement[1] + displacement[2]*displacement[2]);


    double proton_y    = 0.0;
    double proton_mTm0 = 0.0;
    PhysMath::rapidity_mTm0(proton_vect,proton_y,proton_mTm0);
    double pion_y    = 0.0;
    double pion_mTm0 = 0.0;
    PhysMath::rapidity_mTm0(pion_vect,pion_y,pion_mTm0);
    double protonDCA = PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, proton_vect);
    double pionDCA   = PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, pion_vect);
    double protonDCA_X, protonDCA_Y,protonDCA_Z;
    //cout << "Proton Vector: " << proton_vect[0] << "  " << proton_vect[1] << "  " << proton_vect[2] << endl;
    PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, proton_vect, protonDCA_X, protonDCA_Y,protonDCA_Z);
    //cout << "dca check: " << sqrt(protonDCA_X*protonDCA_X+ protonDCA_Y*protonDCA_Y+ protonDCA_Z*protonDCA_Z) << "   =?   " << protonDCA << endl;
    proton_y_mTm0_histo->Fill(proton_y,proton_mTm0);
    pion_y_mTm0_histo->Fill(pion_y,pion_mTm0);
    proton_dca->Fill(protonDCA);
    pion_dca->Fill(pionDCA);
    lambda_decay_length->Fill(distance);
    proton_pX->Fill(proton_vect->Px());
    proton_pY->Fill(proton_vect->Py());
    proton_pZ->Fill(proton_vect->Pz());

    proton_dcaX->Fill(fabs(protonDCA_X));
    proton_dcaY->Fill(fabs(protonDCA_Y));
    proton_dcaZ->Fill(fabs(protonDCA_Z));
    secondVertex_X->Fill(secondaryVertex->X());
    secondVertex_Y->Fill(secondaryVertex->Y());
    secondVertex_Z->Fill(secondaryVertex->Z());
    #ifdef _PHYS_MATH_DEBUG_
      cout << " Proton    y: " << proton_y << " mTm0: " << proton_mTm0 << "  DCA: " << protonDCA
               <<  "    (" << protonDCA_X << ", " << protonDCA_Y << ", " << protonDCA_Z << ")" <<   endl;
      cout << " Pion      y: " << pion_y << " mTm0: " << pion_mTm0 << "  DCA: " << pionDCA << endl;
      cout << " Lambda Decay Length: " << distance << endl;
    #endif

  }

  HistogramUtilities::makeDifferential(proton_y_mTm0_histo, 1.0/((double) numTrials));
  HistogramUtilities::makeDifferential(pion_y_mTm0_histo, 1.0/((double) numTrials));

  TCanvas* canv = new TCanvas();
  proton_y_mTm0_histo->Draw("colz");
  TCanvas* canv2 = new TCanvas();
  pion_y_mTm0_histo->Draw("colz");
  TCanvas* canv3 = new TCanvas();
  gPad->SetLogy();
  lambda_decay_length->Draw("E");
  TCanvas* canv4 = new TCanvas();
  gPad->SetLogy();
  pion_dca->Draw("E");
  proton_dca->Draw("E Same");

  TCanvas* canv5 = new TCanvas();
  gPad->SetLogy();
  proton_dcaZ->Draw("E");
  proton_dcaY->Draw("E SAME");
  proton_dcaX->Draw("E SAME");

  TCanvas* canv6 = new TCanvas();
  gPad->SetLogy(false);
  proton_pZ->Draw("E");
  proton_pY->Draw("E SAME");
  proton_pX->Draw("E SAME");

  TCanvas* canv7 = new TCanvas();
  gPad->SetLogy(false);
  lambda_pX->Draw("E");
  lambda_pY->Draw("E SAME");
  lambda_pZ->Draw("E SAME");

  TCanvas* canv8 = new TCanvas();
  gPad->SetLogy(false);
  secondVertex_X->Draw("E");
  secondVertex_Y->Draw("E SAME");
  secondVertex_Z->Draw("E SAME");


  return;
}





void PhysMathHistogramUtilities::testXiMinusDecay(){

  double y1 = 0;
  double y2 = 0;
  double m1 = 0;
  double m2 = 0;
  double y3 = 0;
  double m3 = 0;

  TH2D* proton_y_mTm0_histo = new TH2D("proton_y_mTm0_histo","Proton Phase Space Probability from #Xi^{-} Decay; y ;m_{T}-m_{p}", 300, 0.1,0.9,300,0,0.23);
  TH2D* pion_y_mTm0_histo = new TH2D("pion_y_mTm0_histo","#pi^{-} Space Probability from #Xi^{-} Decay; y ;m_{T}-m_{#pi}", 300, 0,1.5,300,0,0.14);
  double lambda_rapidity = 0.5;
  double lambda_mTm0 = 0.1;
  TRandom3 randomThing;
  randomThing.SetSeed(time(NULL));
  int numTrials = 10000000;
  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    PhysMath::xiMinusDecay(&randomThing, lambda_rapidity , lambda_mTm0, y1, m1, y2, m2, y3,m3);
    proton_y_mTm0_histo->Fill(y1,m1);
    pion_y_mTm0_histo->Fill(y2,m2);
    pion_y_mTm0_histo->Fill(y3,m3);
  }

  HistogramUtilities::makeDifferential(proton_y_mTm0_histo, 1.0/((double) numTrials));
  HistogramUtilities::makeDifferential(pion_y_mTm0_histo, 0.5/((double) numTrials));

  TCanvas* canv = new TCanvas();
  proton_y_mTm0_histo->Draw("colz");
  TCanvas* canv2 = new TCanvas();
  pion_y_mTm0_histo->Draw("colz");

  return;
}


void PhysMathHistogramUtilities::testNewXiMinusDecay(){

  TH2D* proton_y_mTm0_histo = new TH2D("proton_y_mTm0_histo","Proton Phase Space Probability; y ;m_{T}-m_{p} [GeV/c^{2}]", 300, -0.2,1.2,300,0,0.25);
  TH2D* pion_y_mTm0_histo = new TH2D("pion_y_mTm0_histo","Pion Phase Space Probability; y ;m_{T}-m_{#pi} [GeV/c^{2}]", 300, -0.5,1.4,300,0,0.15);

  TH1D* proton_dca = new TH1D("proton_dca","DCA of Proton from Lambda Decay;DCA [cm];Counts", 300, 0, 30);
  proton_dca->SetLineColor(kGreen+1);
  TH1D* pion_dca = new TH1D("pion_dca","DCA of Pion from Lambda Decay;DCA [cm];Counts", 300, 0, 30);
  pion_dca->SetLineColor(kRed+1);
  TH1D* lambda_decay_length = new TH1D("lambda_decay_length","Lambda Decay Length; Decay Length [cm]; Counts", 300, 0, 35);

  double xi_rapidity = 0.5;
  double xi_mTm0 = 0.1;
  TRandom3 randomThing;
  randomThing.SetSeed(time(NULL));
  ROOT::Math::Boost boostObj;
  int numTrials = 10000000;
  ROOT::Math::XYZVector primaryVertex(0,0,0);
  double xi_mass = 1.32171;
  double pT     = PhysMath::pT_from_mTm0(xi_mTm0,xi_mass);
  double pZ     = PhysMath::pZ_from_y_mTm0(xi_rapidity, xi_mTm0,xi_mass);
  double pTot   = PhysMath::pTot_from_y_mTm0(xi_rapidity, xi_mTm0,xi_mass);
  double energy = sqrt(xi_mass*xi_mass + pTot*pTot);
  ROOT::Math::PxPyPzEVector xi_vect(0,pT,pZ,energy);
  ROOT::Math::PxPyPzEVector* proton_vect = NULL;
  ROOT::Math::PxPyPzEVector* pion_vect   = NULL;
  ROOT::Math::PxPyPzEVector* second_pion_vect   = NULL;
  ROOT::Math::XYZVector*     secondaryVertex = NULL;
  ROOT::Math::XYZVector*     thirdVertex = NULL;

  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    double pX = 0;
    double pY = 0;
    randomThing.Circle(pX,pY,pT);
    xi_vect.SetPx(pX);
    xi_vect.SetPy(pY);
    xi_vect.SetPz(pZ);
    #ifdef _PHYS_MATH_DEBUG_
      if( trialIndex%100 == 0 ){
        cout << "####Trial " << trialIndex << "  Xi Check: m = " << xi_mass << " =? " << sqrt(xi_vect.E()*xi_vect.E() - xi_vect.P()*xi_vect.P()) <<  endl;
        cout << "    pT = " << sqrt(pX*pX + pY*pY) << "=?" << pT << endl;
      }
    #endif

    PhysMath::xiMinusDecay( &boostObj,&randomThing, &primaryVertex, &xi_vect, proton_vect, pion_vect,second_pion_vect, secondaryVertex, thirdVertex); // XiMinus -> p + piMinus + piMinus

    double displacement[3]     = {primaryVertex.X() - secondaryVertex->X(), primaryVertex.Y() - secondaryVertex->Y(), primaryVertex.Z() - secondaryVertex->Z()};
    double displacement_two[3] = {primaryVertex.X() - thirdVertex->X(), primaryVertex.Y() - thirdVertex->Y(), primaryVertex.Z() - thirdVertex->Z()};
    double distance            = sqrt(displacement[0]*displacement[0] + displacement[1]*displacement[1] + displacement[2]*displacement[2]);
    double distance_two        = sqrt(displacement_two[0]*displacement_two[0] + displacement_two[1]*displacement_two[1] + displacement_two[2]*displacement_two[2]);

    double proton_y    = 0.0;
    double proton_mTm0 = 0.0;
    PhysMath::rapidity_mTm0(proton_vect,proton_y,proton_mTm0);
    double pion_y      = 0.0;
    double pion_mTm0   = 0.0;
    PhysMath::rapidity_mTm0(pion_vect,pion_y,pion_mTm0);
    double second_pion_y      = 0.0;
    double second_pion_mTm0   = 0.0;
    PhysMath::rapidity_mTm0(second_pion_vect,second_pion_y,second_pion_mTm0);

    double pionDCA           = PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, pion_vect);
    double protonDCA         = PhysMath::dca_from_decay(&primaryVertex, thirdVertex, proton_vect);
    double secondPionDCA     = PhysMath::dca_from_decay(&primaryVertex, thirdVertex, second_pion_vect);

    proton_y_mTm0_histo->Fill(proton_y,proton_mTm0);
    pion_y_mTm0_histo->Fill(pion_y,pion_mTm0);
    pion_y_mTm0_histo->Fill(second_pion_y,second_pion_mTm0);
    proton_dca->Fill(protonDCA);
    pion_dca->Fill(pionDCA);
    pion_dca->Fill(secondPionDCA);
    lambda_decay_length->Fill(distance_two);

    #ifdef _PHYS_MATH_DEBUG_
      cout << " Proton    y: " << proton_y << " mTm0: " << proton_mTm0 << "  DCA: " << protonDCA << endl;
      cout << " Pion      y: " << pion_y << " mTm0: " << pion_mTm0 << "  DCA: " << pionDCA << endl;
      cout << " Pion      y: " << second_pion_y << " mTm0: " << second_pion_mTm0 << "  DCA: " << secondPionDCA << endl;
      cout << " Lambda Decay Length: " << distance_two << "   from initial Xi " << distance <<  endl;
    #endif

  }

  HistogramUtilities::makeDifferential(proton_y_mTm0_histo, 1.0/((double) numTrials));
  HistogramUtilities::makeDifferential(pion_y_mTm0_histo, 1.0/((double) numTrials));

  TCanvas* canv = new TCanvas();
  proton_y_mTm0_histo->Draw("colz");
  TCanvas* canv2 = new TCanvas();
  pion_y_mTm0_histo->Draw("colz");
  TCanvas* canv3 = new TCanvas();
  gPad->SetLogy();
  lambda_decay_length->Draw("E");
  TCanvas* canv4 = new TCanvas();
  gPad->SetLogy();
  pion_dca->Draw("E");
  proton_dca->Draw("E Same");

  return;
}







void PhysMathHistogramUtilities::testOmegaMinusDecay(){
  double y_p = 0;
  double y_K = 0;
  double m_p = 0;
  double m_K = 0;
  double y_pi = 0;
  double m_pi = 0;
  double y_pi2 = 0;
  double m_pi2 = 0;

  TH2D* proton_y_mTm0_histo = new TH2D("proton_y_mTm0_histo","Proton Phase Space Probability from #Omega^{-} Decay; y ;m_{T}-m_{p}", 300, -0.1,1.2,300,0,0.33);
  TH2D* kaon_y_mTm0_histo = new TH2D("kaon_y_mTm0_histo","Kaon Phase Space Probability from #Omega^{-} Decay; y ;m_{T}-m_{K}", 300, -0.1,1.1,300,0,0.18);
  TH2D* pion_y_mTm0_histo = new TH2D("pion_y_mTm0_histo","#pi^{-} Phase Space Probability from #Omega^{-} Decay; y ;m_{T}-m_{#pi}", 300, 0,2.3,300,0,0.4);
  double omega_rapidity = 0.5;
  double omega_mTm0 = 0.1;
  TRandom3 randomThing;
  randomThing.SetSeed(time(NULL));
  int numTrials = 10000000;
  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    double y_p = -1;
    double y_K = -1;
    double m_p = -1;
    double m_K = -1;
    double y_pi = -1;
    double m_pi = -1;
    double y_pi2 = -1;
    double m_pi2 = -1;
    PhysMath::omegaMinusDecay(&randomThing, omega_rapidity , omega_mTm0, y_p, m_p, y_K, m_K, y_pi,m_pi, y_pi2,m_pi2);
    if(m_p > 0.0) proton_y_mTm0_histo->Fill(y_p,m_p);
    if(m_pi > 0.0) pion_y_mTm0_histo->Fill(y_pi,m_pi);
    if(m_pi2 > 0.0) pion_y_mTm0_histo->Fill(y_pi2,m_pi2);
    if(m_K > 0.0) kaon_y_mTm0_histo->Fill(y_K,m_K);;
  }

  HistogramUtilities::makeDifferential(proton_y_mTm0_histo, 1.0/((double) numTrials));
  HistogramUtilities::makeDifferential(kaon_y_mTm0_histo, 0.5/((double) numTrials));
  HistogramUtilities::makeDifferential(pion_y_mTm0_histo, 0.5/((double) numTrials));

  TCanvas* canv = new TCanvas();
  proton_y_mTm0_histo->Draw("colz");
  TCanvas* canv2 = new TCanvas();
  kaon_y_mTm0_histo->Draw("colz");
  TCanvas* canv3 = new TCanvas();
  pion_y_mTm0_histo->Draw("colz");

  return;
}


void PhysMathHistogramUtilities::testNewOmegaMinusDecay(){

  TH2D* proton_y_mTm0_histo = new TH2D("proton_y_mTm0_histo","Proton Phase Space Probability; y ;m_{T}-m_{p} [GeV/c^{2}]", 300, -0.2,1.2,300,0,0.35);
  TH2D* kaon_y_mTm0_histo = new TH2D("kaon_y_mTm0_histo","Kaon Phase Space Probability from #Omega^{-} Decay; y ;m_{T}-m_{K}", 300, -0.1,1.1,300,0,0.18);
  TH2D* pion_y_mTm0_histo = new TH2D("pion_y_mTm0_histo","Pion Phase Space Probability; y ;m_{T}-m_{#pi} [GeV/c^{2}]", 300, -1.5,2.5,300,0,0.4);

  TH1D* proton_dca = new TH1D("proton_dca","DCA of Proton from Lambda Decay;DCA [cm];Counts", 300, 0, 30);
  proton_dca->SetLineColor(kGreen+1);
  TH1D* pion_dca = new TH1D("pion_dca","DCA of Pion from Lambda Decay;DCA [cm];Counts", 300, 0, 30);
  pion_dca->SetLineColor(kRed+1);
  TH1D* kaon_dca = new TH1D("kaon_dca","DCA of kaon from Lambda Decay;DCA [cm];Counts", 300, 0, 30);
  kaon_dca->SetLineColor(kBlue+1);

  TH1D* decay_length = new TH1D("decay_length","Decay Length; Decay Length [cm]; Counts", 300, 0, 35);

  double omega_rapidity = 0.5;
  double omega_mTm0 = 0.1;
  TRandom3 randomThing;
  randomThing.SetSeed(time(NULL));
  ROOT::Math::Boost boostObj;
  int numTrials = 10000000;
  ROOT::Math::XYZVector primaryVertex(0,0,0);
  double omega_mass = 1.67245;
  double pT     = PhysMath::pT_from_mTm0(omega_mTm0,omega_mass);
  double pZ     = PhysMath::pZ_from_y_mTm0(omega_rapidity, omega_mTm0,omega_mass);
  double pTot   = PhysMath::pTot_from_y_mTm0(omega_rapidity, omega_mTm0,omega_mass);
  double energy = sqrt(omega_mass*omega_mass + pTot*pTot);
  ROOT::Math::PxPyPzEVector  omega_vect(0,pT,pZ,energy);
  ROOT::Math::PxPyPzEVector* proton_vect = NULL;
  ROOT::Math::PxPyPzEVector* pion_vect   = NULL;
  ROOT::Math::PxPyPzEVector* kaon_vect   = NULL;
  ROOT::Math::PxPyPzEVector* second_pion_vect   = NULL;
  ROOT::Math::PxPyPzEVector* pion_zero_vect   = NULL;
  ROOT::Math::XYZVector*     secondaryVertex = NULL;
  ROOT::Math::XYZVector*     thirdVertex = NULL;
  ROOT::Math::XYZVector*     fourthVertex = NULL;

  for(int trialIndex = 0; trialIndex < numTrials; trialIndex++){
    double pX = 0;
    double pY = 0;
    randomThing.Circle(pX,pY,pT);
    omega_vect.SetPx(pX);
    omega_vect.SetPy(pY);
    omega_vect.SetPz(pZ);
    #ifdef _PHYS_MATH_DEBUG_
      if( trialIndex%100 == 0 ){
        cout << "####Trial " << trialIndex << "  Xi Check: m = " << omega_mass << " =? " << sqrt(omega_vect.E()*omega_vect.E() - omega_vect.P()*omega_vect.P()) <<  endl;
        cout << "    pT = " << sqrt(pX*pX + pY*pY) << "=?" << pT << endl;
      }
    #endif

    int proton_vertex_index;
    int kaon_vertex_index;
    int pion_vertex_index;
    int second_pion_vertex_index;
    PhysMath::omegaMinusDecay( &boostObj,&randomThing, &primaryVertex, &omega_vect, proton_vect, pion_zero_vect, kaon_vect,pion_vect,second_pion_vect, secondaryVertex, thirdVertex, fourthVertex,proton_vertex_index, kaon_vertex_index, pion_vertex_index, second_pion_vertex_index);


    double proton_y    = 0.0;
    double proton_mTm0 = -999;
    if(proton_vect) PhysMath::rapidity_mTm0(proton_vect,proton_y,proton_mTm0);
    double pion_y      = 0.0;
    double pion_mTm0   = -999;
    if(pion_vect) PhysMath::rapidity_mTm0(pion_vect,pion_y,pion_mTm0);
    double kaon_y      = 0.0;
    double kaon_mTm0   = -999;
    if(kaon_vect) PhysMath::rapidity_mTm0(kaon_vect,kaon_y,kaon_mTm0);
    double second_pion_y      = 0.0;
    double second_pion_mTm0   = -999;
    if(kaon_vect) PhysMath::rapidity_mTm0(kaon_vect,kaon_y,kaon_mTm0);

    double pionDCA       = -999;
    double protonDCA     = -999;
    double secondPionDCA = -999;
    double kaonDCA       = -999;

    if(pion_vertex_index == 2)          pionDCA = PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, pion_vect);
    else if(pion_vertex_index == 3)     pionDCA = PhysMath::dca_from_decay(&primaryVertex, thirdVertex, pion_vect);
    else if(pion_vertex_index == 4)     pionDCA = PhysMath::dca_from_decay(&primaryVertex, fourthVertex, pion_vect);

    if(kaon_vertex_index == 2)          kaonDCA = PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, kaon_vect);

    if(proton_vertex_index == 2)        protonDCA = PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, proton_vect);
    else if(proton_vertex_index == 3)   protonDCA = PhysMath::dca_from_decay(&primaryVertex, thirdVertex, proton_vect);
    else if(proton_vertex_index == 4)   protonDCA = PhysMath::dca_from_decay(&primaryVertex, fourthVertex, proton_vect);

    if(second_pion_vertex_index == 2)        secondPionDCA = PhysMath::dca_from_decay(&primaryVertex, secondaryVertex, second_pion_vect);
    else if(second_pion_vertex_index == 3)   secondPionDCA = PhysMath::dca_from_decay(&primaryVertex, thirdVertex, second_pion_vect);
    else if(second_pion_vertex_index == 4)   secondPionDCA = PhysMath::dca_from_decay(&primaryVertex, fourthVertex, second_pion_vect);


    proton_y_mTm0_histo->Fill(proton_y,proton_mTm0);
    pion_y_mTm0_histo->Fill(pion_y,pion_mTm0);
    pion_y_mTm0_histo->Fill(second_pion_y,second_pion_mTm0);
    kaon_y_mTm0_histo->Fill(kaon_y,kaon_mTm0);

    proton_dca->Fill(protonDCA);
    pion_dca->Fill(pionDCA);
    pion_dca->Fill(secondPionDCA);
    kaon_dca->Fill(kaonDCA);

  }

  HistogramUtilities::makeDifferential(proton_y_mTm0_histo, 1.0/((double) numTrials));
  HistogramUtilities::makeDifferential(pion_y_mTm0_histo, 1.0/((double) numTrials));
  HistogramUtilities::makeDifferential(kaon_y_mTm0_histo, 1.0/((double) numTrials));


  TCanvas* canv = new TCanvas();
  proton_y_mTm0_histo->Draw("colz");
  TCanvas* canv2 = new TCanvas();
  pion_y_mTm0_histo->Draw("colz");
  TCanvas* canv5 = new TCanvas();
  kaon_y_mTm0_histo->Draw("colz");
  //TCanvas* canv3 = new TCanvas();
  //gPad->SetLogy();
  //lambda_decay_length->Draw("E");
  TCanvas* canv4 = new TCanvas();
  gPad->SetLogy();
  pion_dca->Draw("E");
  kaon_dca->Draw("E Same");
  proton_dca->Draw("E Same");

  return;
}













#endif



void PhysMathHistogramUtilities::testLegendreExpansionInitialization(int a_order, bool a_freeLegendreBoundsInFit){

  TVirtualFitter::SetMaxIterations(100000);
  double x[190];
  double y[190];
  double exl[190];
  double exh[190];
  double eyl[190];
  double eyh[190];

  x[0]=0.005; y[0]=0; exl[0]=0.005; exh[0]=0.005; eyl[0]=0; eyh[0]=0.000910597;
  x[1]=0.015; y[1]=0.00675676; exl[1]=0.005; exh[1]=0.005; eyl[1]=0.00287919; eyh[1]=0.00402919;
  x[2]=0.025; y[2]=0.228797; exl[2]=0.005; exh[2]=0.005; eyl[2]=0.01827; eyh[2]=0.0189811;
  x[3]=0.035; y[3]=0.331034; exl[3]=0.005; exh[3]=0.005; eyl[3]=0.0222562; eyh[3]=0.022772;
  x[4]=0.045; y[4]=0.303109; exl[4]=0.005; exh[4]=0.005; eyl[4]=0.0229991; eyh[4]=0.0236762;
  x[5]=0.055; y[5]=0.398204; exl[5]=0.005; exh[5]=0.005; eyl[5]=0.0265053; eyh[5]=0.0269096;
  x[6]=0.065; y[6]=0.637363; exl[6]=0.005; exh[6]=0.005; eyl[6]=0.0253827; eyh[6]=0.0248819;
  x[7]=0.075; y[7]=0.884013; exl[7]=0.005; exh[7]=0.005; eyl[7]=0.0187084; eyh[7]=0.0171079;
  x[8]=0.085; y[8]=0.866013; exl[8]=0.005; exh[8]=0.005; eyl[8]=0.0202368; eyh[8]=0.0186476;
  x[9]=0.095; y[9]=0.891374; exl[9]=0.005; exh[9]=0.005; eyl[9]=0.0184032; eyh[9]=0.0167404;
  x[10]=0.105; y[10]=0.911032; exl[10]=0.005; exh[10]=0.005; eyl[10]=0.0179528; eyh[10]=0.0160059;
  x[11]=0.115; y[11]=0.90301; exl[11]=0.005; exh[11]=0.005; eyl[11]=0.0180016; eyh[11]=0.0162083;
  x[12]=0.125; y[12]=0.881967; exl[12]=0.005; exh[12]=0.005; eyl[12]=0.0192844; eyh[12]=0.0176196;
  x[13]=0.135; y[13]=0.89726; exl[13]=0.005; exh[13]=0.005; eyl[13]=0.0186583; eyh[13]=0.0168488;
  x[14]=0.145; y[14]=0.916149; exl[14]=0.005; exh[14]=0.005; eyl[14]=0.0163065; eyh[14]=0.0145854;
  x[15]=0.155; y[15]=0.870036; exl[15]=0.005; exh[15]=0.005; eyl[15]=0.0210579; eyh[15]=0.0192834;
  x[16]=0.165; y[16]=0.897764; exl[16]=0.005; exh[16]=0.005; eyl[16]=0.0179568; eyh[16]=0.0162663;
  x[17]=0.175; y[17]=0.926316; exl[17]=0.005; exh[17]=0.005; eyl[17]=0.0164823; eyh[17]=0.0144888;
  x[18]=0.185; y[18]=0.919861; exl[18]=0.005; exh[18]=0.005; eyl[18]=0.0170045; eyh[18]=0.0150561;
  x[19]=0.195; y[19]=0.909091; exl[19]=0.005; exh[19]=0.005; eyl[19]=0.0191424; eyh[19]=0.0169911;
  x[20]=0.205; y[20]=0.918919; exl[20]=0.005; exh[20]=0.005; eyl[20]=0.0168105; eyh[20]=0.0149256;
  x[21]=0.215; y[21]=0.926641; exl[21]=0.005; exh[21]=0.005; eyl[21]=0.0173103; eyh[21]=0.0151151;
  x[22]=0.225; y[22]=0.942953; exl[22]=0.005; exh[22]=0.005; eyl[22]=0.0144561; eyh[22]=0.0124711;
  x[23]=0.235; y[23]=0.929889; exl[23]=0.005; exh[23]=0.005; eyl[23]=0.0165828; eyh[23]=0.014468;
  x[24]=0.245; y[24]=0.92963; exl[24]=0.005; exh[24]=0.005; eyl[24]=0.0166411; eyh[24]=0.0145198;
  x[25]=0.255; y[25]=0.931174; exl[25]=0.005; exh[25]=0.005; eyl[25]=0.0172909; eyh[25]=0.0149633;
  x[26]=0.265; y[26]=0.936759; exl[26]=0.005; exh[26]=0.005; eyl[26]=0.0164803; eyh[26]=0.0141766;
  x[27]=0.275; y[27]=0.915058; exl[27]=0.005; exh[27]=0.005; eyl[27]=0.0183893; eyh[27]=0.016256;
  x[28]=0.285; y[28]=0.906137; exl[28]=0.005; exh[28]=0.005; eyl[28]=0.0184896; eyh[28]=0.0165388;
  x[29]=0.295; y[29]=0.921502; exl[29]=0.005; exh[29]=0.005; eyl[29]=0.0166757; eyh[29]=0.0147593;
  x[30]=0.305; y[30]=0.916968; exl[30]=0.005; exh[30]=0.005; eyl[30]=0.0175823; eyh[30]=0.015578;
  x[31]=0.315; y[31]=0.948207; exl[31]=0.005; exh[31]=0.005; eyl[31]=0.0152276; eyh[31]=0.0128396;
  x[32]=0.325; y[32]=0.927835; exl[32]=0.005; exh[32]=0.005; eyl[32]=0.01616; eyh[32]=0.0142003;
  x[33]=0.335; y[33]=0.936759; exl[33]=0.005; exh[33]=0.005; eyl[33]=0.0164803; eyh[33]=0.0141766;
  x[34]=0.345; y[34]=0.949416; exl[34]=0.005; exh[34]=0.005; eyl[34]=0.0148859; eyh[34]=0.0125469;
  x[35]=0.355; y[35]=0.952569; exl[35]=0.005; exh[35]=0.005; eyl[35]=0.0146147; eyh[35]=0.0122199;
  x[36]=0.365; y[36]=0.925311; exl[36]=0.005; exh[36]=0.005; eyl[36]=0.0181215; eyh[36]=0.0157701;
  x[37]=0.375; y[37]=0.957692; exl[37]=0.005; exh[37]=0.005; eyl[37]=0.0137261; eyh[37]=0.0113661;
  x[38]=0.385; y[38]=0.960784; exl[38]=0.005; exh[38]=0.005; eyl[38]=0.0134403; eyh[38]=0.0110149;
  x[39]=0.395; y[39]=0.959259; exl[39]=0.005; exh[39]=0.005; eyl[39]=0.0132341; eyh[39]=0.0109531;
  x[40]=0.405; y[40]=0.911017; exl[40]=0.005; exh[40]=0.005; eyl[40]=0.0196867; eyh[40]=0.0173695;
  x[41]=0.415; y[41]=0.928; exl[41]=0.005; exh[41]=0.005; eyl[41]=0.0175036; eyh[41]=0.0152218;
  x[42]=0.425; y[42]=0.934866; exl[42]=0.005; exh[42]=0.005; eyl[42]=0.016408; eyh[42]=0.0141853;
  x[43]=0.435; y[43]=0.946619; exl[43]=0.005; exh[43]=0.005; eyl[43]=0.0145082; eyh[43]=0.012384;
  x[44]=0.445; y[44]=0.938406; exl[44]=0.005; exh[44]=0.005; eyl[44]=0.0155566; eyh[44]=0.0134366;
  x[45]=0.455; y[45]=0.934545; exl[45]=0.005; exh[45]=0.005; eyl[45]=0.0159886; eyh[45]=0.0138807;
  x[46]=0.465; y[46]=0.935622; exl[46]=0.005; exh[46]=0.005; eyl[46]=0.0173534; eyh[46]=0.0148588;
  x[47]=0.475; y[47]=0.899614; exl[47]=0.005; exh[47]=0.005; eyl[47]=0.0196837; eyh[47]=0.0176319;
  x[48]=0.485; y[48]=0.956989; exl[48]=0.005; exh[48]=0.005; eyl[48]=0.0132984; eyh[48]=0.0111038;
  x[49]=0.495; y[49]=0.950382; exl[49]=0.005; exh[49]=0.005; eyl[49]=0.0146127; eyh[49]=0.0123129;
  x[50]=0.505; y[50]=0.964427; exl[50]=0.005; exh[50]=0.005; eyl[50]=0.0129619; eyh[50]=0.0104939;
  x[51]=0.515; y[51]=0.946154; exl[51]=0.005; exh[51]=0.005; eyl[51]=0.0151848; eyh[51]=0.0128913;
  x[52]=0.525; y[52]=0.951852; exl[52]=0.005; exh[52]=0.005; eyl[52]=0.0141957; eyh[52]=0.0119562;
  x[53]=0.535; y[53]=0.928571; exl[53]=0.005; exh[53]=0.005; eyl[53]=0.0164239; eyh[53]=0.0143836;
  x[54]=0.545; y[54]=0.936567; exl[54]=0.005; exh[54]=0.005; eyl[54]=0.0159994; eyh[54]=0.0138258;
  x[55]=0.555; y[55]=0.932271; exl[55]=0.005; exh[55]=0.005; eyl[55]=0.0170291; eyh[55]=0.0147324;
  x[56]=0.565; y[56]=0.96124; exl[56]=0.005; exh[56]=0.005; eyl[56]=0.0132889; eyh[56]=0.0108891;
  x[57]=0.575; y[57]=0.933594; exl[57]=0.005; exh[57]=0.005; eyl[57]=0.0167128; eyh[57]=0.0144537;
  x[58]=0.585; y[58]=0.944649; exl[58]=0.005; exh[58]=0.005; eyl[58]=0.0150215; eyh[58]=0.0128293;
  x[59]=0.595; y[59]=0.942652; exl[59]=0.005; exh[59]=0.005; eyl[59]=0.0150098; eyh[59]=0.0128909;
  x[60]=0.605; y[60]=0.942966; exl[60]=0.005; exh[60]=0.005; eyl[60]=0.0154589; eyh[60]=0.0132092;
  x[61]=0.615; y[61]=0.943548; exl[61]=0.005; exh[61]=0.005; eyl[61]=0.0158882; eyh[61]=0.0134987;
  x[62]=0.625; y[62]=0.935743; exl[62]=0.005; exh[62]=0.005; eyl[62]=0.0167324; eyh[62]=0.0143975;
  x[63]=0.635; y[63]=0.934363; exl[63]=0.005; exh[63]=0.005; eyl[63]=0.0165286; eyh[63]=0.0142914;
  x[64]=0.645; y[64]=0.941667; exl[64]=0.005; exh[64]=0.005; eyl[64]=0.0163943; eyh[64]=0.0139364;
  x[65]=0.655; y[65]=0.885827; exl[65]=0.005; exh[65]=0.005; eyl[65]=0.0209367; eyh[65]=0.0189183;
  x[66]=0.665; y[66]=0.945312; exl[66]=0.005; exh[66]=0.005; eyl[66]=0.0154122; eyh[66]=0.0130876;
  x[67]=0.675; y[67]=0.941176; exl[67]=0.005; exh[67]=0.005; eyl[67]=0.0165259; eyh[67]=0.0140503;
  x[68]=0.685; y[68]=0.928839; exl[68]=0.005; exh[68]=0.005; eyl[68]=0.0168185; eyh[68]=0.0146775;
  x[69]=0.695; y[69]=0.911647; exl[69]=0.005; exh[69]=0.005; eyl[69]=0.0190811; eyh[69]=0.0168811;
  x[70]=0.705; y[70]=0.931915; exl[70]=0.005; exh[70]=0.005; eyl[70]=0.0176787; eyh[70]=0.0152278;
  x[71]=0.715; y[71]=0.935484; exl[71]=0.005; exh[71]=0.005; eyl[71]=0.0180468; eyh[71]=0.015369;
  x[72]=0.725; y[72]=0.945312; exl[72]=0.005; exh[72]=0.005; eyl[72]=0.0154122; eyh[72]=0.0130876;
  x[73]=0.735; y[73]=0.955947; exl[73]=0.005; exh[73]=0.005; eyl[73]=0.015039; eyh[73]=0.0123459;
  x[74]=0.745; y[74]=0.908088; exl[74]=0.005; exh[74]=0.005; eyl[74]=0.0185083; eyh[74]=0.0165119;
  x[75]=0.755; y[75]=0.925795; exl[75]=0.005; exh[75]=0.005; eyl[75]=0.0165926; eyh[75]=0.0145876;
  x[76]=0.765; y[76]=0.908; exl[76]=0.005; exh[76]=0.005; eyl[76]=0.019357; eyh[76]=0.0171859;
  x[77]=0.775; y[77]=0.918728; exl[77]=0.005; exh[77]=0.005; eyl[77]=0.0172311; eyh[77]=0.0152606;
  x[78]=0.785; y[78]=0.886447; exl[78]=0.005; exh[78]=0.005; eyl[78]=0.0201189; eyh[78]=0.0182375;
  x[79]=0.795; y[79]=0.912698; exl[79]=0.005; exh[79]=0.005; eyl[79]=0.0188682; eyh[79]=0.0166886;
  x[80]=0.805; y[80]=0.909434; exl[80]=0.005; exh[80]=0.005; eyl[80]=0.0186515; eyh[80]=0.0165956;
  x[81]=0.815; y[81]=0.924779; exl[81]=0.005; exh[81]=0.005; eyl[81]=0.0188084; eyh[81]=0.0163043;
  x[82]=0.825; y[82]=0.916031; exl[82]=0.005; exh[82]=0.005; eyl[82]=0.0181914; eyh[82]=0.0160774;
  x[83]=0.835; y[83]=0.896825; exl[83]=0.005; exh[83]=0.005; eyl[83]=0.0201904; eyh[83]=0.0180969;
  x[84]=0.845; y[84]=0.910714; exl[84]=0.005; exh[84]=0.005; eyl[84]=0.020266; eyh[84]=0.0178268;
  x[85]=0.855; y[85]=0.894958; exl[85]=0.005; exh[85]=0.005; eyl[85]=0.0209555; eyh[85]=0.01875;
  x[86]=0.865; y[86]=0.910569; exl[86]=0.005; exh[86]=0.005; eyl[86]=0.0192988; eyh[86]=0.0170781;
  x[87]=0.875; y[87]=0.884921; exl[87]=0.005; exh[87]=0.005; eyl[87]=0.0210893; eyh[87]=0.0190599;
  x[88]=0.885; y[88]=0.946058; exl[88]=0.005; exh[88]=0.005; eyl[88]=0.0158333; eyh[88]=0.0133591;
  x[89]=0.895; y[89]=0.912548; exl[89]=0.005; exh[89]=0.005; eyl[89]=0.0184601; eyh[89]=0.0163723;
  x[90]=0.905; y[90]=0.888476; exl[90]=0.005; exh[90]=0.005; eyl[90]=0.0201294; eyh[90]=0.0182099;
  x[91]=0.915; y[91]=0.901887; exl[91]=0.005; exh[91]=0.005; eyl[91]=0.019269; eyh[91]=0.0172519;
  x[92]=0.925; y[92]=0.921569; exl[92]=0.005; exh[92]=0.005; eyl[92]=0.0179433; eyh[92]=0.0157413;
  x[93]=0.935; y[93]=0.906137; exl[93]=0.005; exh[93]=0.005; eyl[93]=0.0184896; eyh[93]=0.0165388;
  x[94]=0.945; y[94]=0.913223; exl[94]=0.005; exh[94]=0.005; eyl[94]=0.0192293; eyh[94]=0.0169568;
  x[95]=0.955; y[95]=0.925197; exl[95]=0.005; exh[95]=0.005; eyl[95]=0.0176326; eyh[95]=0.0154021;
  x[96]=0.965; y[96]=0.903084; exl[96]=0.005; exh[96]=0.005; eyl[96]=0.0208014; eyh[96]=0.0184405;
  x[97]=0.975; y[97]=0.888889; exl[97]=0.005; exh[97]=0.005; eyl[97]=0.0207976; eyh[97]=0.0187469;
  x[98]=0.985; y[98]=0.887446; exl[98]=0.005; exh[98]=0.005; eyl[98]=0.0218788; eyh[98]=0.0196507;
  x[99]=0.995; y[99]=0.903915; exl[99]=0.005; exh[99]=0.005; eyl[99]=0.0185268; eyh[99]=0.0166145;
  x[100]=1.005; y[100]=0.894531; exl[100]=0.005; exh[100]=0.005; eyl[100]=0.0202018; eyh[100]=0.0181531;
  x[101]=1.015; y[101]=0.921811; exl[101]=0.005; exh[101]=0.005; eyl[101]=0.0183854; eyh[101]=0.0160734;
  x[102]=1.025; y[102]=0.912351; exl[102]=0.005; exh[102]=0.005; eyl[102]=0.0189386; eyh[102]=0.0167523;
  x[103]=1.035; y[103]=0.928; exl[103]=0.005; exh[103]=0.005; eyl[103]=0.0175036; eyh[103]=0.0152218;
  x[104]=1.045; y[104]=0.928839; exl[104]=0.005; exh[104]=0.005; eyl[104]=0.0168185; eyh[104]=0.0146775;
  x[105]=1.055; y[105]=0.903509; exl[105]=0.005; exh[105]=0.005; eyl[105]=0.0207166; eyh[105]=0.0183634;
  x[106]=1.065; y[106]=0.899642; exl[106]=0.005; exh[106]=0.005; eyl[106]=0.0189281; eyh[106]=0.0170229;
  x[107]=1.075; y[107]=0.92; exl[107]=0.005; exh[107]=0.005; eyl[107]=0.0182814; eyh[107]=0.0160441;
  x[108]=1.085; y[108]=0.91791; exl[108]=0.005; exh[108]=0.005; eyl[108]=0.0178081; eyh[108]=0.0157317;
  x[109]=1.095; y[109]=0.910781; exl[109]=0.005; exh[109]=0.005; eyl[109]=0.0183917; eyh[109]=0.0163595;
  x[110]=1.105; y[110]=0.912281; exl[110]=0.005; exh[110]=0.005; eyl[110]=0.0199332; eyh[110]=0.0175272;
  x[111]=1.115; y[111]=0.88755; exl[111]=0.005; exh[111]=0.005; eyl[111]=0.0210282; eyh[111]=0.0189601;
  x[112]=1.125; y[112]=0.894977; exl[112]=0.005; exh[112]=0.005; eyl[112]=0.0218906; eyh[112]=0.0194943;
  x[113]=1.135; y[113]=0.907336; exl[113]=0.005; exh[113]=0.005; eyl[113]=0.0190551; eyh[113]=0.0169628;
  x[114]=1.145; y[114]=0.896414; exl[114]=0.005; exh[114]=0.005; eyl[114]=0.0202649; eyh[114]=0.0181653;
  x[115]=1.155; y[115]=0.917603; exl[115]=0.005; exh[115]=0.005; eyl[115]=0.0178708; eyh[115]=0.0157882;
  x[116]=1.165; y[116]=0.926641; exl[116]=0.005; exh[116]=0.005; eyl[116]=0.0173103; eyh[116]=0.0151151;
  x[117]=1.175; y[117]=0.891892; exl[117]=0.005; exh[117]=0.005; eyl[117]=0.0202786; eyh[117]=0.0182674;
  x[118]=1.185; y[118]=0.906383; exl[118]=0.005; exh[118]=0.005; eyl[118]=0.0201413; eyh[118]=0.0178413;
  x[119]=1.195; y[119]=0.896947; exl[119]=0.005; exh[119]=0.005; eyl[119]=0.019773; eyh[119]=0.0177586;
  x[120]=1.205; y[120]=0.904215; exl[120]=0.005; exh[120]=0.005; eyl[120]=0.0192354; eyh[120]=0.0171754;
  x[121]=1.215; y[121]=0.923664; exl[121]=0.005; exh[121]=0.005; eyl[121]=0.0174903; eyh[121]=0.015336;
  x[122]=1.225; y[122]=0.915058; exl[122]=0.005; exh[122]=0.005; eyl[122]=0.0183893; eyh[122]=0.016256;
  x[123]=1.235; y[123]=0.909091; exl[123]=0.005; exh[123]=0.005; eyl[123]=0.019597; eyh[123]=0.017348;
  x[124]=1.245; y[124]=0.887719; exl[124]=0.005; exh[124]=0.005; eyl[124]=0.0195839; eyh[124]=0.0177755;
  x[125]=1.255; y[125]=0.907631; exl[125]=0.005; exh[125]=0.005; eyl[125]=0.0194296; eyh[125]=0.0172518;
  x[126]=1.265; y[126]=0.85; exl[126]=0.005; exh[126]=0.005; eyl[126]=0.0229894; eyh[126]=0.0212028;
  x[127]=1.275; y[127]=0.903084; exl[127]=0.005; exh[127]=0.005; eyl[127]=0.0208014; eyh[127]=0.0184405;
  x[128]=1.285; y[128]=0.904; exl[128]=0.005; exh[128]=0.005; eyl[128]=0.0196942; eyh[128]=0.0175449;
  x[129]=1.295; y[129]=0.928302; exl[129]=0.005; exh[129]=0.005; eyl[129]=0.0169388; eyh[129]=0.0147845;
  x[130]=1.305; y[130]=0.909483; exl[130]=0.005; exh[130]=0.005; eyl[130]=0.0200039; eyh[130]=0.0176558;
  x[131]=1.315; y[131]=0.896552; exl[131]=0.005; exh[131]=0.005; eyl[131]=0.0198432; eyh[131]=0.0178231;
  x[132]=1.325; y[132]=0.928; exl[132]=0.005; exh[132]=0.005; eyl[132]=0.0175036; eyh[132]=0.0152218;
  x[133]=1.335; y[133]=0.875; exl[133]=0.005; exh[133]=0.005; eyl[133]=0.0219687; eyh[133]=0.0199606;
  x[134]=1.345; y[134]=0.883721; exl[134]=0.005; exh[134]=0.005; eyl[134]=0.0209174; eyh[134]=0.0189412;
  x[135]=1.355; y[135]=0.914729; exl[135]=0.005; exh[135]=0.005; eyl[135]=0.0184563; eyh[135]=0.0163165;
  x[136]=1.365; y[136]=0.890977; exl[136]=0.005; exh[136]=0.005; eyl[136]=0.0200651; eyh[136]=0.0181113;
  x[137]=1.375; y[137]=0.899614; exl[137]=0.005; exh[137]=0.005; eyl[137]=0.0196837; eyh[137]=0.0176319;
  x[138]=1.385; y[138]=0.903614; exl[138]=0.005; exh[138]=0.005; eyl[138]=0.0197679; eyh[138]=0.0176121;
  x[139]=1.395; y[139]=0.893617; exl[139]=0.005; exh[139]=0.005; eyl[139]=0.0212027; eyh[139]=0.0189768;
  x[140]=1.405; y[140]=0.876596; exl[140]=0.005; exh[140]=0.005; eyl[140]=0.0224812; eyh[140]=0.0203534;
  x[141]=1.415; y[141]=0.883408; exl[141]=0.005; exh[141]=0.005; eyl[141]=0.022598; eyh[141]=0.0203148;
  x[142]=1.425; y[142]=0.901235; exl[142]=0.005; exh[142]=0.005; eyl[142]=0.0202215; eyh[142]=0.018026;
  x[143]=1.435; y[143]=0.895197; exl[143]=0.005; exh[143]=0.005; eyl[143]=0.0213649; eyh[143]=0.0190715;
  x[144]=1.445; y[144]=0.886792; exl[144]=0.005; exh[144]=0.005; eyl[144]=0.020409; eyh[144]=0.0184692;
  x[145]=1.455; y[145]=0.904167; exl[145]=0.005; exh[145]=0.005; eyl[145]=0.0201082; eyh[145]=0.0178686;
  x[146]=1.465; y[146]=0.900383; exl[146]=0.005; exh[146]=0.005; eyl[146]=0.0195435; eyh[146]=0.0175034;
  x[147]=1.475; y[147]=0.91635; exl[147]=0.005; exh[147]=0.005; eyl[147]=0.0181264; eyh[147]=0.0160187;
  x[148]=1.485; y[148]=0.882353; exl[148]=0.005; exh[148]=0.005; eyl[148]=0.021143; eyh[148]=0.0191509;
  x[149]=1.495; y[149]=0.897436; exl[149]=0.005; exh[149]=0.005; eyl[149]=0.0209421; eyh[149]=0.0186845;
  x[150]=1.505; y[150]=0.910156; exl[150]=0.005; exh[150]=0.005; eyl[150]=0.0189325; eyh[150]=0.0168006;
  x[151]=1.515; y[151]=0.90681; exl[151]=0.005; exh[151]=0.005; eyl[151]=0.0183658; eyh[151]=0.0164257;
  x[152]=1.525; y[152]=0.859438; exl[152]=0.005; exh[152]=0.005; eyl[152]=0.0229382; eyh[152]=0.0210222;
  x[153]=1.535; y[153]=0.873077; exl[153]=0.005; exh[153]=0.005; eyl[153]=0.0215628; eyh[153]=0.019657;
  x[154]=1.545; y[154]=0.925197; exl[154]=0.005; exh[154]=0.005; eyl[154]=0.0176326; eyh[154]=0.0154021;
  x[155]=1.555; y[155]=0.917603; exl[155]=0.005; exh[155]=0.005; eyl[155]=0.0178708; eyh[155]=0.0157882;
  x[156]=1.565; y[156]=0.870036; exl[156]=0.005; exh[156]=0.005; eyl[156]=0.0210579; eyh[156]=0.0192834;
  x[157]=1.575; y[157]=0.898113; exl[157]=0.005; exh[157]=0.005; eyl[157]=0.0195653; eyh[157]=0.0175676;
  x[158]=1.585; y[158]=0.880309; exl[158]=0.005; exh[158]=0.005; eyl[158]=0.0211149; eyh[158]=0.0191642;
  x[159]=1.595; y[159]=0.897338; exl[159]=0.005; exh[159]=0.005; eyl[159]=0.0197033; eyh[159]=0.0176944;
  x[160]=1.605; y[160]=0.85159; exl[160]=0.005; exh[160]=0.005; eyl[160]=0.0219154; eyh[160]=0.0202659;
  x[161]=1.615; y[161]=0.855019; exl[161]=0.005; exh[161]=0.005; eyl[161]=0.0222994; eyh[161]=0.0205473;
  x[162]=1.625; y[162]=0.872428; exl[162]=0.005; exh[162]=0.005; eyl[162]=0.0223797; eyh[162]=0.0203448;
  x[163]=1.635; y[163]=0.872881; exl[163]=0.005; exh[163]=0.005; eyl[163]=0.0226918; eyh[163]=0.0205941;
  x[164]=1.645; y[164]=0.897638; exl[164]=0.005; exh[164]=0.005; eyl[164]=0.020043; eyh[164]=0.0179616;
  x[165]=1.655; y[165]=0.9; exl[165]=0.005; exh[165]=0.005; eyl[165]=0.0204562; eyh[165]=0.0182403;
  x[166]=1.665; y[166]=0.884462; exl[166]=0.005; exh[166]=0.005; eyl[166]=0.0211664; eyh[166]=0.0191314;
  x[167]=1.675; y[167]=0.872428; exl[167]=0.005; exh[167]=0.005; eyl[167]=0.0223797; eyh[167]=0.0203448;
  x[168]=1.685; y[168]=0.89916; exl[168]=0.005; exh[168]=0.005; eyl[168]=0.0206156; eyh[168]=0.018386;
  x[169]=1.695; y[169]=0.88716; exl[169]=0.005; exh[169]=0.005; eyl[169]=0.0207118; eyh[169]=0.0187099;
  x[170]=1.705; y[170]=0.889734; exl[170]=0.005; exh[170]=0.005; eyl[170]=0.0202762; eyh[170]=0.0183066;
  x[171]=1.715; y[171]=0.900794; exl[171]=0.005; exh[171]=0.005; eyl[171]=0.019874; eyh[171]=0.0177591;
  x[172]=1.725; y[172]=0.889344; exl[172]=0.005; exh[172]=0.005; eyl[172]=0.0211174; eyh[172]=0.0189971;
  x[173]=1.735; y[173]=0.896825; exl[173]=0.005; exh[173]=0.005; eyl[173]=0.0201904; eyh[173]=0.0180969;
  x[174]=1.745; y[174]=0.844961; exl[174]=0.005; exh[174]=0.005; eyl[174]=0.0233689; eyh[174]=0.0215946;
  x[175]=1.755; y[175]=0.880658; exl[175]=0.005; exh[175]=0.005; eyl[175]=0.0218042; eyh[175]=0.0197236;
  x[176]=1.765; y[176]=0.886792; exl[176]=0.005; exh[176]=0.005; eyl[176]=0.020409; eyh[176]=0.0184692;
  x[177]=1.775; y[177]=0.875502; exl[177]=0.005; exh[177]=0.005; eyl[177]=0.0218883; eyh[177]=0.0198855;
  x[178]=1.785; y[178]=0.897321; exl[178]=0.005; exh[178]=0.005; eyl[178]=0.0214384; eyh[178]=0.0190811;
  x[179]=1.795; y[179]=0.926471; exl[179]=0.005; exh[179]=0.005; eyl[179]=0.0168814; eyh[179]=0.0147919;
  x[180]=1.805; y[180]=0.842105; exl[180]=0.005; exh[180]=0.005; eyl[180]=0.0240635; eyh[180]=0.0222261;
  x[181]=1.815; y[181]=0.912351; exl[181]=0.005; exh[181]=0.005; eyl[181]=0.0189386; eyh[181]=0.0167523;
  x[182]=1.825; y[182]=0.881679; exl[182]=0.005; exh[182]=0.005; eyl[182]=0.0208933; eyh[182]=0.0189578;
  x[183]=1.835; y[183]=0.904762; exl[183]=0.005; exh[183]=0.005; eyl[183]=0.0187424; eyh[183]=0.0167699;
  x[184]=1.845; y[184]=0.869919; exl[184]=0.005; exh[184]=0.005; eyl[184]=0.0224058; eyh[184]=0.0204093;
  x[185]=1.855; y[185]=0.875; exl[185]=0.005; exh[185]=0.005; eyl[185]=0.0212649; eyh[185]=0.0193781;
  x[186]=1.865; y[186]=0.857708; exl[186]=0.005; exh[186]=0.005; eyl[186]=0.0228557; eyh[186]=0.020979;
  x[187]=1.875; y[187]=0.866667; exl[187]=0.005; exh[187]=0.005; eyl[187]=0.0215528; eyh[187]=0.0197492;
  x[188]=1.885; y[188]=0.87395; exl[188]=0.005; exh[188]=0.005; eyl[188]=0.0225183; eyh[188]=0.0204321;
  x[189]=1.895; y[189]=0.892704; exl[189]=0.005; exh[189]=0.005; eyl[189]=0.0213707; eyh[189]=0.019131;


  TGraphAsymmErrors* graph = new TGraphAsymmErrors(190,x,y,exl,exh,eyl,eyh);
  graph->SetTitle(Form("Fit with Legendre Polynomial (order = %d)",a_order));
  TF1* function = PhysMath::getLegendreExpansion("legendreFunct", a_order,false,0.0, 1.5);
  PhysMathHistogramUtilities::initializeLegendreExpansion(function, graph, 0.0, 1.5, false);
  if(a_freeLegendreBoundsInFit){
    function->SetParLimits(function->GetNpar()-1,-50,50);
    function->SetParLimits(function->GetNpar()-2,-50,50);
  }
  TF1* initialFunct = (TF1*) function->Clone();
  initialFunct->SetName("initialFunct");
  initialFunct->SetLineColor(kYellow);


  int status = graph->Fit(function,"","EX0",0,1.5);
  function->SetLineColor(kRed);
  cout << " Fit Status: " << status << endl;




  TCanvas* canv = new TCanvas();
  graph->Draw("AP");
  initialFunct->Draw("same");
  function->Draw("same");


}

void PhysMathHistogramUtilities::testDetermineGoodLegendreOrder(){

  double x[190];
  double y[190];
  double exl[190];
  double exh[190];
  double eyl[190];
  double eyh[190];

  x[0]=0.005; y[0]=0; exl[0]=0.005; exh[0]=0.005; eyl[0]=0; eyh[0]=0.000910597;
  x[1]=0.015; y[1]=0.00675676; exl[1]=0.005; exh[1]=0.005; eyl[1]=0.00287919; eyh[1]=0.00402919;
  x[2]=0.025; y[2]=0.228797; exl[2]=0.005; exh[2]=0.005; eyl[2]=0.01827; eyh[2]=0.0189811;
  x[3]=0.035; y[3]=0.331034; exl[3]=0.005; exh[3]=0.005; eyl[3]=0.0222562; eyh[3]=0.022772;
  x[4]=0.045; y[4]=0.303109; exl[4]=0.005; exh[4]=0.005; eyl[4]=0.0229991; eyh[4]=0.0236762;
  x[5]=0.055; y[5]=0.398204; exl[5]=0.005; exh[5]=0.005; eyl[5]=0.0265053; eyh[5]=0.0269096;
  x[6]=0.065; y[6]=0.637363; exl[6]=0.005; exh[6]=0.005; eyl[6]=0.0253827; eyh[6]=0.0248819;
  x[7]=0.075; y[7]=0.884013; exl[7]=0.005; exh[7]=0.005; eyl[7]=0.0187084; eyh[7]=0.0171079;
  x[8]=0.085; y[8]=0.866013; exl[8]=0.005; exh[8]=0.005; eyl[8]=0.0202368; eyh[8]=0.0186476;
  x[9]=0.095; y[9]=0.891374; exl[9]=0.005; exh[9]=0.005; eyl[9]=0.0184032; eyh[9]=0.0167404;
  x[10]=0.105; y[10]=0.911032; exl[10]=0.005; exh[10]=0.005; eyl[10]=0.0179528; eyh[10]=0.0160059;
  x[11]=0.115; y[11]=0.90301; exl[11]=0.005; exh[11]=0.005; eyl[11]=0.0180016; eyh[11]=0.0162083;
  x[12]=0.125; y[12]=0.881967; exl[12]=0.005; exh[12]=0.005; eyl[12]=0.0192844; eyh[12]=0.0176196;
  x[13]=0.135; y[13]=0.89726; exl[13]=0.005; exh[13]=0.005; eyl[13]=0.0186583; eyh[13]=0.0168488;
  x[14]=0.145; y[14]=0.916149; exl[14]=0.005; exh[14]=0.005; eyl[14]=0.0163065; eyh[14]=0.0145854;
  x[15]=0.155; y[15]=0.870036; exl[15]=0.005; exh[15]=0.005; eyl[15]=0.0210579; eyh[15]=0.0192834;
  x[16]=0.165; y[16]=0.897764; exl[16]=0.005; exh[16]=0.005; eyl[16]=0.0179568; eyh[16]=0.0162663;
  x[17]=0.175; y[17]=0.926316; exl[17]=0.005; exh[17]=0.005; eyl[17]=0.0164823; eyh[17]=0.0144888;
  x[18]=0.185; y[18]=0.919861; exl[18]=0.005; exh[18]=0.005; eyl[18]=0.0170045; eyh[18]=0.0150561;
  x[19]=0.195; y[19]=0.909091; exl[19]=0.005; exh[19]=0.005; eyl[19]=0.0191424; eyh[19]=0.0169911;
  x[20]=0.205; y[20]=0.918919; exl[20]=0.005; exh[20]=0.005; eyl[20]=0.0168105; eyh[20]=0.0149256;
  x[21]=0.215; y[21]=0.926641; exl[21]=0.005; exh[21]=0.005; eyl[21]=0.0173103; eyh[21]=0.0151151;
  x[22]=0.225; y[22]=0.942953; exl[22]=0.005; exh[22]=0.005; eyl[22]=0.0144561; eyh[22]=0.0124711;
  x[23]=0.235; y[23]=0.929889; exl[23]=0.005; exh[23]=0.005; eyl[23]=0.0165828; eyh[23]=0.014468;
  x[24]=0.245; y[24]=0.92963; exl[24]=0.005; exh[24]=0.005; eyl[24]=0.0166411; eyh[24]=0.0145198;
  x[25]=0.255; y[25]=0.931174; exl[25]=0.005; exh[25]=0.005; eyl[25]=0.0172909; eyh[25]=0.0149633;
  x[26]=0.265; y[26]=0.936759; exl[26]=0.005; exh[26]=0.005; eyl[26]=0.0164803; eyh[26]=0.0141766;
  x[27]=0.275; y[27]=0.915058; exl[27]=0.005; exh[27]=0.005; eyl[27]=0.0183893; eyh[27]=0.016256;
  x[28]=0.285; y[28]=0.906137; exl[28]=0.005; exh[28]=0.005; eyl[28]=0.0184896; eyh[28]=0.0165388;
  x[29]=0.295; y[29]=0.921502; exl[29]=0.005; exh[29]=0.005; eyl[29]=0.0166757; eyh[29]=0.0147593;
  x[30]=0.305; y[30]=0.916968; exl[30]=0.005; exh[30]=0.005; eyl[30]=0.0175823; eyh[30]=0.015578;
  x[31]=0.315; y[31]=0.948207; exl[31]=0.005; exh[31]=0.005; eyl[31]=0.0152276; eyh[31]=0.0128396;
  x[32]=0.325; y[32]=0.927835; exl[32]=0.005; exh[32]=0.005; eyl[32]=0.01616; eyh[32]=0.0142003;
  x[33]=0.335; y[33]=0.936759; exl[33]=0.005; exh[33]=0.005; eyl[33]=0.0164803; eyh[33]=0.0141766;
  x[34]=0.345; y[34]=0.949416; exl[34]=0.005; exh[34]=0.005; eyl[34]=0.0148859; eyh[34]=0.0125469;
  x[35]=0.355; y[35]=0.952569; exl[35]=0.005; exh[35]=0.005; eyl[35]=0.0146147; eyh[35]=0.0122199;
  x[36]=0.365; y[36]=0.925311; exl[36]=0.005; exh[36]=0.005; eyl[36]=0.0181215; eyh[36]=0.0157701;
  x[37]=0.375; y[37]=0.957692; exl[37]=0.005; exh[37]=0.005; eyl[37]=0.0137261; eyh[37]=0.0113661;
  x[38]=0.385; y[38]=0.960784; exl[38]=0.005; exh[38]=0.005; eyl[38]=0.0134403; eyh[38]=0.0110149;
  x[39]=0.395; y[39]=0.959259; exl[39]=0.005; exh[39]=0.005; eyl[39]=0.0132341; eyh[39]=0.0109531;
  x[40]=0.405; y[40]=0.911017; exl[40]=0.005; exh[40]=0.005; eyl[40]=0.0196867; eyh[40]=0.0173695;
  x[41]=0.415; y[41]=0.928; exl[41]=0.005; exh[41]=0.005; eyl[41]=0.0175036; eyh[41]=0.0152218;
  x[42]=0.425; y[42]=0.934866; exl[42]=0.005; exh[42]=0.005; eyl[42]=0.016408; eyh[42]=0.0141853;
  x[43]=0.435; y[43]=0.946619; exl[43]=0.005; exh[43]=0.005; eyl[43]=0.0145082; eyh[43]=0.012384;
  x[44]=0.445; y[44]=0.938406; exl[44]=0.005; exh[44]=0.005; eyl[44]=0.0155566; eyh[44]=0.0134366;
  x[45]=0.455; y[45]=0.934545; exl[45]=0.005; exh[45]=0.005; eyl[45]=0.0159886; eyh[45]=0.0138807;
  x[46]=0.465; y[46]=0.935622; exl[46]=0.005; exh[46]=0.005; eyl[46]=0.0173534; eyh[46]=0.0148588;
  x[47]=0.475; y[47]=0.899614; exl[47]=0.005; exh[47]=0.005; eyl[47]=0.0196837; eyh[47]=0.0176319;
  x[48]=0.485; y[48]=0.956989; exl[48]=0.005; exh[48]=0.005; eyl[48]=0.0132984; eyh[48]=0.0111038;
  x[49]=0.495; y[49]=0.950382; exl[49]=0.005; exh[49]=0.005; eyl[49]=0.0146127; eyh[49]=0.0123129;
  x[50]=0.505; y[50]=0.964427; exl[50]=0.005; exh[50]=0.005; eyl[50]=0.0129619; eyh[50]=0.0104939;
  x[51]=0.515; y[51]=0.946154; exl[51]=0.005; exh[51]=0.005; eyl[51]=0.0151848; eyh[51]=0.0128913;
  x[52]=0.525; y[52]=0.951852; exl[52]=0.005; exh[52]=0.005; eyl[52]=0.0141957; eyh[52]=0.0119562;
  x[53]=0.535; y[53]=0.928571; exl[53]=0.005; exh[53]=0.005; eyl[53]=0.0164239; eyh[53]=0.0143836;
  x[54]=0.545; y[54]=0.936567; exl[54]=0.005; exh[54]=0.005; eyl[54]=0.0159994; eyh[54]=0.0138258;
  x[55]=0.555; y[55]=0.932271; exl[55]=0.005; exh[55]=0.005; eyl[55]=0.0170291; eyh[55]=0.0147324;
  x[56]=0.565; y[56]=0.96124; exl[56]=0.005; exh[56]=0.005; eyl[56]=0.0132889; eyh[56]=0.0108891;
  x[57]=0.575; y[57]=0.933594; exl[57]=0.005; exh[57]=0.005; eyl[57]=0.0167128; eyh[57]=0.0144537;
  x[58]=0.585; y[58]=0.944649; exl[58]=0.005; exh[58]=0.005; eyl[58]=0.0150215; eyh[58]=0.0128293;
  x[59]=0.595; y[59]=0.942652; exl[59]=0.005; exh[59]=0.005; eyl[59]=0.0150098; eyh[59]=0.0128909;
  x[60]=0.605; y[60]=0.942966; exl[60]=0.005; exh[60]=0.005; eyl[60]=0.0154589; eyh[60]=0.0132092;
  x[61]=0.615; y[61]=0.943548; exl[61]=0.005; exh[61]=0.005; eyl[61]=0.0158882; eyh[61]=0.0134987;
  x[62]=0.625; y[62]=0.935743; exl[62]=0.005; exh[62]=0.005; eyl[62]=0.0167324; eyh[62]=0.0143975;
  x[63]=0.635; y[63]=0.934363; exl[63]=0.005; exh[63]=0.005; eyl[63]=0.0165286; eyh[63]=0.0142914;
  x[64]=0.645; y[64]=0.941667; exl[64]=0.005; exh[64]=0.005; eyl[64]=0.0163943; eyh[64]=0.0139364;
  x[65]=0.655; y[65]=0.885827; exl[65]=0.005; exh[65]=0.005; eyl[65]=0.0209367; eyh[65]=0.0189183;
  x[66]=0.665; y[66]=0.945312; exl[66]=0.005; exh[66]=0.005; eyl[66]=0.0154122; eyh[66]=0.0130876;
  x[67]=0.675; y[67]=0.941176; exl[67]=0.005; exh[67]=0.005; eyl[67]=0.0165259; eyh[67]=0.0140503;
  x[68]=0.685; y[68]=0.928839; exl[68]=0.005; exh[68]=0.005; eyl[68]=0.0168185; eyh[68]=0.0146775;
  x[69]=0.695; y[69]=0.911647; exl[69]=0.005; exh[69]=0.005; eyl[69]=0.0190811; eyh[69]=0.0168811;
  x[70]=0.705; y[70]=0.931915; exl[70]=0.005; exh[70]=0.005; eyl[70]=0.0176787; eyh[70]=0.0152278;
  x[71]=0.715; y[71]=0.935484; exl[71]=0.005; exh[71]=0.005; eyl[71]=0.0180468; eyh[71]=0.015369;
  x[72]=0.725; y[72]=0.945312; exl[72]=0.005; exh[72]=0.005; eyl[72]=0.0154122; eyh[72]=0.0130876;
  x[73]=0.735; y[73]=0.955947; exl[73]=0.005; exh[73]=0.005; eyl[73]=0.015039; eyh[73]=0.0123459;
  x[74]=0.745; y[74]=0.908088; exl[74]=0.005; exh[74]=0.005; eyl[74]=0.0185083; eyh[74]=0.0165119;
  x[75]=0.755; y[75]=0.925795; exl[75]=0.005; exh[75]=0.005; eyl[75]=0.0165926; eyh[75]=0.0145876;
  x[76]=0.765; y[76]=0.908; exl[76]=0.005; exh[76]=0.005; eyl[76]=0.019357; eyh[76]=0.0171859;
  x[77]=0.775; y[77]=0.918728; exl[77]=0.005; exh[77]=0.005; eyl[77]=0.0172311; eyh[77]=0.0152606;
  x[78]=0.785; y[78]=0.886447; exl[78]=0.005; exh[78]=0.005; eyl[78]=0.0201189; eyh[78]=0.0182375;
  x[79]=0.795; y[79]=0.912698; exl[79]=0.005; exh[79]=0.005; eyl[79]=0.0188682; eyh[79]=0.0166886;
  x[80]=0.805; y[80]=0.909434; exl[80]=0.005; exh[80]=0.005; eyl[80]=0.0186515; eyh[80]=0.0165956;
  x[81]=0.815; y[81]=0.924779; exl[81]=0.005; exh[81]=0.005; eyl[81]=0.0188084; eyh[81]=0.0163043;
  x[82]=0.825; y[82]=0.916031; exl[82]=0.005; exh[82]=0.005; eyl[82]=0.0181914; eyh[82]=0.0160774;
  x[83]=0.835; y[83]=0.896825; exl[83]=0.005; exh[83]=0.005; eyl[83]=0.0201904; eyh[83]=0.0180969;
  x[84]=0.845; y[84]=0.910714; exl[84]=0.005; exh[84]=0.005; eyl[84]=0.020266; eyh[84]=0.0178268;
  x[85]=0.855; y[85]=0.894958; exl[85]=0.005; exh[85]=0.005; eyl[85]=0.0209555; eyh[85]=0.01875;
  x[86]=0.865; y[86]=0.910569; exl[86]=0.005; exh[86]=0.005; eyl[86]=0.0192988; eyh[86]=0.0170781;
  x[87]=0.875; y[87]=0.884921; exl[87]=0.005; exh[87]=0.005; eyl[87]=0.0210893; eyh[87]=0.0190599;
  x[88]=0.885; y[88]=0.946058; exl[88]=0.005; exh[88]=0.005; eyl[88]=0.0158333; eyh[88]=0.0133591;
  x[89]=0.895; y[89]=0.912548; exl[89]=0.005; exh[89]=0.005; eyl[89]=0.0184601; eyh[89]=0.0163723;
  x[90]=0.905; y[90]=0.888476; exl[90]=0.005; exh[90]=0.005; eyl[90]=0.0201294; eyh[90]=0.0182099;
  x[91]=0.915; y[91]=0.901887; exl[91]=0.005; exh[91]=0.005; eyl[91]=0.019269; eyh[91]=0.0172519;
  x[92]=0.925; y[92]=0.921569; exl[92]=0.005; exh[92]=0.005; eyl[92]=0.0179433; eyh[92]=0.0157413;
  x[93]=0.935; y[93]=0.906137; exl[93]=0.005; exh[93]=0.005; eyl[93]=0.0184896; eyh[93]=0.0165388;
  x[94]=0.945; y[94]=0.913223; exl[94]=0.005; exh[94]=0.005; eyl[94]=0.0192293; eyh[94]=0.0169568;
  x[95]=0.955; y[95]=0.925197; exl[95]=0.005; exh[95]=0.005; eyl[95]=0.0176326; eyh[95]=0.0154021;
  x[96]=0.965; y[96]=0.903084; exl[96]=0.005; exh[96]=0.005; eyl[96]=0.0208014; eyh[96]=0.0184405;
  x[97]=0.975; y[97]=0.888889; exl[97]=0.005; exh[97]=0.005; eyl[97]=0.0207976; eyh[97]=0.0187469;
  x[98]=0.985; y[98]=0.887446; exl[98]=0.005; exh[98]=0.005; eyl[98]=0.0218788; eyh[98]=0.0196507;
  x[99]=0.995; y[99]=0.903915; exl[99]=0.005; exh[99]=0.005; eyl[99]=0.0185268; eyh[99]=0.0166145;
  x[100]=1.005; y[100]=0.894531; exl[100]=0.005; exh[100]=0.005; eyl[100]=0.0202018; eyh[100]=0.0181531;
  x[101]=1.015; y[101]=0.921811; exl[101]=0.005; exh[101]=0.005; eyl[101]=0.0183854; eyh[101]=0.0160734;
  x[102]=1.025; y[102]=0.912351; exl[102]=0.005; exh[102]=0.005; eyl[102]=0.0189386; eyh[102]=0.0167523;
  x[103]=1.035; y[103]=0.928; exl[103]=0.005; exh[103]=0.005; eyl[103]=0.0175036; eyh[103]=0.0152218;
  x[104]=1.045; y[104]=0.928839; exl[104]=0.005; exh[104]=0.005; eyl[104]=0.0168185; eyh[104]=0.0146775;
  x[105]=1.055; y[105]=0.903509; exl[105]=0.005; exh[105]=0.005; eyl[105]=0.0207166; eyh[105]=0.0183634;
  x[106]=1.065; y[106]=0.899642; exl[106]=0.005; exh[106]=0.005; eyl[106]=0.0189281; eyh[106]=0.0170229;
  x[107]=1.075; y[107]=0.92; exl[107]=0.005; exh[107]=0.005; eyl[107]=0.0182814; eyh[107]=0.0160441;
  x[108]=1.085; y[108]=0.91791; exl[108]=0.005; exh[108]=0.005; eyl[108]=0.0178081; eyh[108]=0.0157317;
  x[109]=1.095; y[109]=0.910781; exl[109]=0.005; exh[109]=0.005; eyl[109]=0.0183917; eyh[109]=0.0163595;
  x[110]=1.105; y[110]=0.912281; exl[110]=0.005; exh[110]=0.005; eyl[110]=0.0199332; eyh[110]=0.0175272;
  x[111]=1.115; y[111]=0.88755; exl[111]=0.005; exh[111]=0.005; eyl[111]=0.0210282; eyh[111]=0.0189601;
  x[112]=1.125; y[112]=0.894977; exl[112]=0.005; exh[112]=0.005; eyl[112]=0.0218906; eyh[112]=0.0194943;
  x[113]=1.135; y[113]=0.907336; exl[113]=0.005; exh[113]=0.005; eyl[113]=0.0190551; eyh[113]=0.0169628;
  x[114]=1.145; y[114]=0.896414; exl[114]=0.005; exh[114]=0.005; eyl[114]=0.0202649; eyh[114]=0.0181653;
  x[115]=1.155; y[115]=0.917603; exl[115]=0.005; exh[115]=0.005; eyl[115]=0.0178708; eyh[115]=0.0157882;
  x[116]=1.165; y[116]=0.926641; exl[116]=0.005; exh[116]=0.005; eyl[116]=0.0173103; eyh[116]=0.0151151;
  x[117]=1.175; y[117]=0.891892; exl[117]=0.005; exh[117]=0.005; eyl[117]=0.0202786; eyh[117]=0.0182674;
  x[118]=1.185; y[118]=0.906383; exl[118]=0.005; exh[118]=0.005; eyl[118]=0.0201413; eyh[118]=0.0178413;
  x[119]=1.195; y[119]=0.896947; exl[119]=0.005; exh[119]=0.005; eyl[119]=0.019773; eyh[119]=0.0177586;
  x[120]=1.205; y[120]=0.904215; exl[120]=0.005; exh[120]=0.005; eyl[120]=0.0192354; eyh[120]=0.0171754;
  x[121]=1.215; y[121]=0.923664; exl[121]=0.005; exh[121]=0.005; eyl[121]=0.0174903; eyh[121]=0.015336;
  x[122]=1.225; y[122]=0.915058; exl[122]=0.005; exh[122]=0.005; eyl[122]=0.0183893; eyh[122]=0.016256;
  x[123]=1.235; y[123]=0.909091; exl[123]=0.005; exh[123]=0.005; eyl[123]=0.019597; eyh[123]=0.017348;
  x[124]=1.245; y[124]=0.887719; exl[124]=0.005; exh[124]=0.005; eyl[124]=0.0195839; eyh[124]=0.0177755;
  x[125]=1.255; y[125]=0.907631; exl[125]=0.005; exh[125]=0.005; eyl[125]=0.0194296; eyh[125]=0.0172518;
  x[126]=1.265; y[126]=0.85; exl[126]=0.005; exh[126]=0.005; eyl[126]=0.0229894; eyh[126]=0.0212028;
  x[127]=1.275; y[127]=0.903084; exl[127]=0.005; exh[127]=0.005; eyl[127]=0.0208014; eyh[127]=0.0184405;
  x[128]=1.285; y[128]=0.904; exl[128]=0.005; exh[128]=0.005; eyl[128]=0.0196942; eyh[128]=0.0175449;
  x[129]=1.295; y[129]=0.928302; exl[129]=0.005; exh[129]=0.005; eyl[129]=0.0169388; eyh[129]=0.0147845;
  x[130]=1.305; y[130]=0.909483; exl[130]=0.005; exh[130]=0.005; eyl[130]=0.0200039; eyh[130]=0.0176558;
  x[131]=1.315; y[131]=0.896552; exl[131]=0.005; exh[131]=0.005; eyl[131]=0.0198432; eyh[131]=0.0178231;
  x[132]=1.325; y[132]=0.928; exl[132]=0.005; exh[132]=0.005; eyl[132]=0.0175036; eyh[132]=0.0152218;
  x[133]=1.335; y[133]=0.875; exl[133]=0.005; exh[133]=0.005; eyl[133]=0.0219687; eyh[133]=0.0199606;
  x[134]=1.345; y[134]=0.883721; exl[134]=0.005; exh[134]=0.005; eyl[134]=0.0209174; eyh[134]=0.0189412;
  x[135]=1.355; y[135]=0.914729; exl[135]=0.005; exh[135]=0.005; eyl[135]=0.0184563; eyh[135]=0.0163165;
  x[136]=1.365; y[136]=0.890977; exl[136]=0.005; exh[136]=0.005; eyl[136]=0.0200651; eyh[136]=0.0181113;
  x[137]=1.375; y[137]=0.899614; exl[137]=0.005; exh[137]=0.005; eyl[137]=0.0196837; eyh[137]=0.0176319;
  x[138]=1.385; y[138]=0.903614; exl[138]=0.005; exh[138]=0.005; eyl[138]=0.0197679; eyh[138]=0.0176121;
  x[139]=1.395; y[139]=0.893617; exl[139]=0.005; exh[139]=0.005; eyl[139]=0.0212027; eyh[139]=0.0189768;
  x[140]=1.405; y[140]=0.876596; exl[140]=0.005; exh[140]=0.005; eyl[140]=0.0224812; eyh[140]=0.0203534;
  x[141]=1.415; y[141]=0.883408; exl[141]=0.005; exh[141]=0.005; eyl[141]=0.022598; eyh[141]=0.0203148;
  x[142]=1.425; y[142]=0.901235; exl[142]=0.005; exh[142]=0.005; eyl[142]=0.0202215; eyh[142]=0.018026;
  x[143]=1.435; y[143]=0.895197; exl[143]=0.005; exh[143]=0.005; eyl[143]=0.0213649; eyh[143]=0.0190715;
  x[144]=1.445; y[144]=0.886792; exl[144]=0.005; exh[144]=0.005; eyl[144]=0.020409; eyh[144]=0.0184692;
  x[145]=1.455; y[145]=0.904167; exl[145]=0.005; exh[145]=0.005; eyl[145]=0.0201082; eyh[145]=0.0178686;
  x[146]=1.465; y[146]=0.900383; exl[146]=0.005; exh[146]=0.005; eyl[146]=0.0195435; eyh[146]=0.0175034;
  x[147]=1.475; y[147]=0.91635; exl[147]=0.005; exh[147]=0.005; eyl[147]=0.0181264; eyh[147]=0.0160187;
  x[148]=1.485; y[148]=0.882353; exl[148]=0.005; exh[148]=0.005; eyl[148]=0.021143; eyh[148]=0.0191509;
  x[149]=1.495; y[149]=0.897436; exl[149]=0.005; exh[149]=0.005; eyl[149]=0.0209421; eyh[149]=0.0186845;
  x[150]=1.505; y[150]=0.910156; exl[150]=0.005; exh[150]=0.005; eyl[150]=0.0189325; eyh[150]=0.0168006;
  x[151]=1.515; y[151]=0.90681; exl[151]=0.005; exh[151]=0.005; eyl[151]=0.0183658; eyh[151]=0.0164257;
  x[152]=1.525; y[152]=0.859438; exl[152]=0.005; exh[152]=0.005; eyl[152]=0.0229382; eyh[152]=0.0210222;
  x[153]=1.535; y[153]=0.873077; exl[153]=0.005; exh[153]=0.005; eyl[153]=0.0215628; eyh[153]=0.019657;
  x[154]=1.545; y[154]=0.925197; exl[154]=0.005; exh[154]=0.005; eyl[154]=0.0176326; eyh[154]=0.0154021;
  x[155]=1.555; y[155]=0.917603; exl[155]=0.005; exh[155]=0.005; eyl[155]=0.0178708; eyh[155]=0.0157882;
  x[156]=1.565; y[156]=0.870036; exl[156]=0.005; exh[156]=0.005; eyl[156]=0.0210579; eyh[156]=0.0192834;
  x[157]=1.575; y[157]=0.898113; exl[157]=0.005; exh[157]=0.005; eyl[157]=0.0195653; eyh[157]=0.0175676;
  x[158]=1.585; y[158]=0.880309; exl[158]=0.005; exh[158]=0.005; eyl[158]=0.0211149; eyh[158]=0.0191642;
  x[159]=1.595; y[159]=0.897338; exl[159]=0.005; exh[159]=0.005; eyl[159]=0.0197033; eyh[159]=0.0176944;
  x[160]=1.605; y[160]=0.85159; exl[160]=0.005; exh[160]=0.005; eyl[160]=0.0219154; eyh[160]=0.0202659;
  x[161]=1.615; y[161]=0.855019; exl[161]=0.005; exh[161]=0.005; eyl[161]=0.0222994; eyh[161]=0.0205473;
  x[162]=1.625; y[162]=0.872428; exl[162]=0.005; exh[162]=0.005; eyl[162]=0.0223797; eyh[162]=0.0203448;
  x[163]=1.635; y[163]=0.872881; exl[163]=0.005; exh[163]=0.005; eyl[163]=0.0226918; eyh[163]=0.0205941;
  x[164]=1.645; y[164]=0.897638; exl[164]=0.005; exh[164]=0.005; eyl[164]=0.020043; eyh[164]=0.0179616;
  x[165]=1.655; y[165]=0.9; exl[165]=0.005; exh[165]=0.005; eyl[165]=0.0204562; eyh[165]=0.0182403;
  x[166]=1.665; y[166]=0.884462; exl[166]=0.005; exh[166]=0.005; eyl[166]=0.0211664; eyh[166]=0.0191314;
  x[167]=1.675; y[167]=0.872428; exl[167]=0.005; exh[167]=0.005; eyl[167]=0.0223797; eyh[167]=0.0203448;
  x[168]=1.685; y[168]=0.89916; exl[168]=0.005; exh[168]=0.005; eyl[168]=0.0206156; eyh[168]=0.018386;
  x[169]=1.695; y[169]=0.88716; exl[169]=0.005; exh[169]=0.005; eyl[169]=0.0207118; eyh[169]=0.0187099;
  x[170]=1.705; y[170]=0.889734; exl[170]=0.005; exh[170]=0.005; eyl[170]=0.0202762; eyh[170]=0.0183066;
  x[171]=1.715; y[171]=0.900794; exl[171]=0.005; exh[171]=0.005; eyl[171]=0.019874; eyh[171]=0.0177591;
  x[172]=1.725; y[172]=0.889344; exl[172]=0.005; exh[172]=0.005; eyl[172]=0.0211174; eyh[172]=0.0189971;
  x[173]=1.735; y[173]=0.896825; exl[173]=0.005; exh[173]=0.005; eyl[173]=0.0201904; eyh[173]=0.0180969;
  x[174]=1.745; y[174]=0.844961; exl[174]=0.005; exh[174]=0.005; eyl[174]=0.0233689; eyh[174]=0.0215946;
  x[175]=1.755; y[175]=0.880658; exl[175]=0.005; exh[175]=0.005; eyl[175]=0.0218042; eyh[175]=0.0197236;
  x[176]=1.765; y[176]=0.886792; exl[176]=0.005; exh[176]=0.005; eyl[176]=0.020409; eyh[176]=0.0184692;
  x[177]=1.775; y[177]=0.875502; exl[177]=0.005; exh[177]=0.005; eyl[177]=0.0218883; eyh[177]=0.0198855;
  x[178]=1.785; y[178]=0.897321; exl[178]=0.005; exh[178]=0.005; eyl[178]=0.0214384; eyh[178]=0.0190811;
  x[179]=1.795; y[179]=0.926471; exl[179]=0.005; exh[179]=0.005; eyl[179]=0.0168814; eyh[179]=0.0147919;
  x[180]=1.805; y[180]=0.842105; exl[180]=0.005; exh[180]=0.005; eyl[180]=0.0240635; eyh[180]=0.0222261;
  x[181]=1.815; y[181]=0.912351; exl[181]=0.005; exh[181]=0.005; eyl[181]=0.0189386; eyh[181]=0.0167523;
  x[182]=1.825; y[182]=0.881679; exl[182]=0.005; exh[182]=0.005; eyl[182]=0.0208933; eyh[182]=0.0189578;
  x[183]=1.835; y[183]=0.904762; exl[183]=0.005; exh[183]=0.005; eyl[183]=0.0187424; eyh[183]=0.0167699;
  x[184]=1.845; y[184]=0.869919; exl[184]=0.005; exh[184]=0.005; eyl[184]=0.0224058; eyh[184]=0.0204093;
  x[185]=1.855; y[185]=0.875; exl[185]=0.005; exh[185]=0.005; eyl[185]=0.0212649; eyh[185]=0.0193781;
  x[186]=1.865; y[186]=0.857708; exl[186]=0.005; exh[186]=0.005; eyl[186]=0.0228557; eyh[186]=0.020979;
  x[187]=1.875; y[187]=0.866667; exl[187]=0.005; exh[187]=0.005; eyl[187]=0.0215528; eyh[187]=0.0197492;
  x[188]=1.885; y[188]=0.87395; exl[188]=0.005; exh[188]=0.005; eyl[188]=0.0225183; eyh[188]=0.0204321;
  x[189]=1.895; y[189]=0.892704; exl[189]=0.005; exh[189]=0.005; eyl[189]=0.0213707; eyh[189]=0.019131;

  TGraphAsymmErrors* graph = new TGraphAsymmErrors(190,x,y,exl,exh,eyl,eyh);
  TF1* fitResult = PhysMathHistogramUtilities::determineGoodLegendreOrder(graph, 0.0, 1.5);
  fitResult->SetLineColor(kRed);
  graph->SetTitle(Form("Fit with Legendre Polynomial (order = %d)",fitResult->GetNpar()-3));

  TCanvas* canv = new TCanvas();
  graph->Draw("AP");
  fitResult->Draw("same");


}

#endif
