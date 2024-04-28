#include <iostream>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TKey.h>
#include <TList.h>
#include <vector>
#include <TMath.h>
#include <TLegend.h>
#include <TPad.h>

using namespace std;

Double_t Landau(Double_t* x, Double_t* par)
{
  return TMath::Landau(x[0], par[0], par[1]);
}

Double_t Gaus(Double_t* x, Double_t* par)
{
  return par[3] * TMath::Gaus(x[0], par[4], par[5]);
}

Double_t Convolution(Double_t* x, Double_t* par)
{
  return Landau(x, par) + Gaus(x, par);
}

Double_t Vavilov(Double_t* x, Double_t* par)
{
  return par[4] * TMath::Vavilov(((x[0] - par[2])/par[3]), par[0], par[1]);
}

// Double_t Vavilov(Double_t* x, Double_t* par)
// {
//   return par[0] * ROOT::Math::vavilov_accurate_pdf(x[0], par[1], par[2]);
// }

// Double_t Vavilov(Double_t* x, Double_t* par)
// {
//   return TMath::VavilovAccurate(x[0], par[0]);
// }

void FitAllHistograms(const char* inputFileName,
		      const char* beginPattern,
		      const char* endPattern)
{
  // Open the ROOT file containg the histograms
  TFile* inputFile = new TFile(inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie())
    {
      cerr << "Error: Canont open ROOT file " << inputFileName << endl;
      return;
    }

  // Get the list of keys (histogram names) in the file
  TList* keyList = inputFile->GetListOfKeys();

  // Declaring a vector of vectors that we will store the fit parameters in
  vector<vector<double>> fitParamsVector;

  // Declaring vectors to hold the chi2Values
  vector<double> vectorChi2Landau;
  vector<double> vectorChi2Vavilov;

  // Save the canvs as a PNG file with a nice name
  // Taking run part of file name
  string input(inputFileName);
  string result = "";

  // Find the position of "_run"
  size_t startPos = input.find("_run");
  if (startPos != string::npos)
    {
      // Extract the substring starting from "_run" and its length
      result = input.substr(startPos);
    
      // Remove ".root" from the end
      size_t dotPos = result.find(".root");
      if (dotPos != string::npos)
	{
	  result.erase(dotPos);
	}
    }
  
  // Loop through all keys (histograms) in the file
  for (int i = 0; i < keyList->GetSize(); ++i)
    {
      // Get the key (histogram) at index i
      TKey* key = dynamic_cast<TKey*>(keyList->At(i));
      if(!key) continue;

      // Check if the key is a histogram
      if (key->ReadObj()->IsA() == TH1::Class());
	{
	  // Getting the histogram
	  TH1* hist = static_cast<TH1*>(key->ReadObj());

	  // Getting the name of the histogram
	  TString objName = hist->GetName();

	  //Check if the object meets the matching pattern that we want`
	  if(objName.BeginsWith(beginPattern) && objName.EndsWith(endPattern))
	    {
	      // Print out the histogram name (sanity check)
	      cout << objName << endl;

	      // Create canvas for the histogram
	      TCanvas* canvas = new TCanvas("canvas", "Fitted Histogram");
	      hist->SetStats(0);
	      // gPad->SetLogy();
	      hist->Draw(); // Draw the original histogram on the canvas
	      
	      // Determine range of to use for the fits
	      Double_t minX = hist->GetXaxis()->GetXmin();
	      Double_t maxX = hist->GetXaxis()->GetXmax();
	      hist->GetXaxis()->SetRangeUser(minX, 1500); // Set xMin and xMax to the desired range

	      // Define fit function Landau 
	      TF1* fitLandau = new TF1("fitLandau", "landau", 20, 1000);
	      
	      hist->Fit("fitLandau", "RQ0");
	      hist->Fit("fitLandau", "RWQ");

	      // Could use Landau fit parameters as guesses for Vavilov fit
	      double mpv = fitLandau->GetParameter(1);
	      double width = fitLandau->GetParameter(2);
	      // cout << "mpv: " << mpv << endl;
	      // cout << "width: " << width << endl;

	      // Plot the fit on histogram
	      fitLandau->SetNpx(1000000);
	      fitLandau->Draw("same");

	      // Define fit Vavilov
	      TF1* fitVavilov = new TF1("fitVavilov", Vavilov, 20, 1000, 5);

	      // Name fit parameters
	      fitVavilov->SetParNames("kappa", "beta2", "shift", "width", "scale");

	      // Set initial fit parameters
	      fitVavilov->SetParameters(0.05,0.4, mpv,50.0, 1500.0);
	      // Kappa paramater [0.01, 12.0]
	      fitVavilov->SetParLimits(0, 0.01, 12.0);

	      // Beta2 paramater [0, 1.0]
	      fitVavilov->SetParLimits(1, 0.0, 1.0);

	      // Plotting more points 
	      fitVavilov->SetNpx(1000000);
	      	      
	      // Fit the histogram using the fitting function
	      //hist->Fit("fitVavilov", "WQ");
	      hist->Fit("fitVavilov", "RWQ0");
	      hist->Fit("fitVavilov","RQ+");

	      // Plot the Vavilov fit
	      fitVavilov->SetLineColor(3);
	      fitVavilov->Draw("same");

	      
	      // // Print out parameters
	      // cout << "kappa: " << fitVavilov->GetParameter(0) << endl;
	      // cout << "beta2: " << fitVavilov->GetParameter(1) << endl;
	      // cout << "shift: " << fitVavilov->GetParameter(2) << endl;
	      // cout << "scale: " << fitVavilov->GetParameter(3) << endl;
	      // cout << "constant: " << fitVavilov->GetParameter(4) << endl;

	      // // Save the canvs as a PNG file with a nice name
	      // // Taking run part of file name
	      // string input(inputFileName);
	      // string result = "";

	      // // Find the position of "_run"
	      // size_t startPos = input.find("_run");
	      // if (startPos != string::npos)
	      // 	{
	      // 	  // Extract the substring starting from "_run" and its length
	      // 	  result = input.substr(startPos);
    
	      // 	  // Remove ".root" from the end
	      // 	  size_t dotPos = result.find(".root");
	      // 	  if (dotPos != string::npos)
	      // 	    {
	      // 	      result.erase(dotPos);
	      // 	    }
	      // 	}

	      // Calculate chi-square values for both fits
	      double chi2Landau = fitLandau->GetChisquare();
	      double chi2Vavilov = fitVavilov->GetChisquare();

	      // Get the number of bins in the histogram
	      int nBins = hist->GetNbinsX();

	      // Calculate the number of parameters in the fits
	      int nParamsLandau = fitLandau->GetNpar();
	      int nParamsVavilov = fitVavilov->GetNpar();

	      // Calculate the number of degrees of freedom
	      int ndofLandau = nBins - nParamsLandau;
	      int ndofVavilov = nBins - nParamsVavilov;

	      // Calculate reduced chi-square values
	      double reducedChi2Landau = chi2Landau/ ndofLandau;
	      double reducedChi2Vavilov = chi2Vavilov/ ndofVavilov;

	      // Store the values
	      vectorChi2Landau.push_back(reducedChi2Landau);
	      vectorChi2Vavilov.push_back(reducedChi2Vavilov);
	      // Create a new title with desired information
	      TString newTitle = string(objName) + result;
	      // Setting new title
	      canvas->SetTitle(newTitle);
	      
	      // Create a new legend with fit information and statistics
	      TLegend *newLegend = new TLegend(0.7, 0.7, 0.95, 0.95); // Adjust the position as needed
	      // newLegend->AddEntry(hist, "Histogram", "l");
	      newLegend->SetHeader(result.c_str());
	      newLegend->AddEntry((TObject*)0, Form("Entries: %.0f", hist->GetEntries()), "");
	      newLegend->AddEntry((TObject*)0, Form("Mean: %.2f", hist->GetMean()), "");
	      newLegend->AddEntry((TObject*)0, Form("Std Dev: %.2f", hist->GetStdDev()), "");
	      newLegend->AddEntry(fitLandau, Form("Landau #chi_{R}^{2}=%.2f", reducedChi2Landau), "l");
	      newLegend->AddEntry(fitVavilov, Form("Vavilov #chi_{R}^{2}=%.2f", reducedChi2Vavilov), "l");

	      // Set text size
	      // newLegend->SetTextSize(0.04);

		// Draw the new legend
	      newLegend->Draw();//"same");
	      // Update the canvas
	      canvas->Update();

	      // Saving canvas 
	      string canvasName = "WireEnd/Fitted/Fitted_" + string(objName) + result + ".png";
	      canvas->SaveAs(canvasName.c_str());
	      
	      // // Plot fit parameters as function vs wireend id

	      // // Get the number of fit parameters from the fitting function
	      // int numParams = fitVavilov->GetNpar();
	      
	      // // Store the fit parameters in a vector
	      // vector<double> params(numParams);
	      // for (int j = 0; j < numParams; ++j)
	      // 	{
	      // 	  params[j] = fitVavilov->GetParameter(j);
	      // 	}
	      // fitParamsVector.push_back(params);
	      
	      

	      // Calculate residuals, histogram, etc.
	      TH1F* residHist = new TH1F("residHist", "Residual Histogram", 20, -30, 30);
	      TGraph* residGraph = new TGraph();

	      for (int i = 2; i <= hist->GetNbinsX(); ++i)
		{
		  double binCenter = hist->GetBinCenter(i);
		  double binContent = hist->GetBinContent(i);
		  double fittedValueLandau = fitLandau->Eval(binCenter);
		  double fittedValueVavilov = fitVavilov->Eval(binCenter);
		  double residLandau = binContent - fittedValueLandau;
		  double residVavilov = binContent - fittedValueVavilov;

		  // Fill residuals histograms
		  residHist->Fill(residVavilov);
		  // Plot TGraph
		  residGraph->SetPoint(i, binCenter, residVavilov);
		}

	      // Create canvas for these plots
	      TCanvas* residCanvas = new TCanvas("residCanvas", "Residual Plots", 1600, 600);

	      // Divide the canvas
	      residCanvas->Divide(2,1);

	      // // plot the original histogram
	      // residCanvas->cd(1);
	      // // TPad* pad = new TPad("pad", "Fitted Hist", 0, 0, 1, 0.8);
	      // // pad->Draw();
	      // // pad->cd();
	      // hist->Draw();
	      // fitLandau->Draw("same");
	      // fitVavilov->Draw("same");

	      
	      // Get the pad for the histogram
	      // TPad* pad = (TPad*)residCanvas->GetPad(1);
	      // pad->SetPad(0,0.6,2,2);
	      
	      // Plot the residuals
	      residCanvas->cd(1);
	      residGraph->SetMarkerStyle(52);
	      residGraph->SetMarkerColor(kBlue);
	      
	      residGraph->Draw("AP");

	      // plot the resid hist
	      residCanvas->cd(2);
	      // gPad->SetLogy();
	      residHist->SetLineColor(kRed);
	      residHist->SetFillColorAlpha(kGreen, 0.5);
	      residHist->Draw("HIST");

	      residCanvas->Update();

	      // Save
     	      string canvasResidName = "WireEnd/Resid/" + string(objName) + result + ".png";
	      
	      // canvasResid->SaveAs(canvasResidName.c_str());

	      residCanvas->SaveAs(canvasResidName.c_str());


	      // Multi plot of all
	      TCanvas* multiCanvas = new TCanvas("multiCanvas", "Fitted Hist + Residula", 3200, 1600);
	      multiCanvas->Divide(2,1);
	      
	      // upper pad original histogram with the fitted function
	      multiCanvas->cd(1);
	      hist->Draw();
	      fitLandau->Draw("same");
	      fitVavilov->Draw("same");
	      newLegend->Draw("same");
	      // Lower pad plot of residuals
	      multiCanvas->cd(2);

	      // TPad* pad1 = new TPad("pad1", "pad1", 0.05, 0.05, 0.95, 0.5);
	      TPad* pad = new TPad("pad", "pad", 0, 0, 1, 1);
	      pad->Draw();
	      pad->cd();
	      pad->Divide(1,2);

	      //Subpad1
	      pad->cd(1);
	      	      
	      residGraph->SetMarkerStyle(52);
	      residGraph->SetMarkerColor(kBlue);
	      residGraph->Draw("AP");
	      residGraph->GetXaxis()->SetTitle("X");
	      residGraph->GetYaxis()->SetTitle("Residuals");

	      //subpad 2
	      //TPad* pad2 = new TPad("pad2", "pad2", 0.05, 0.5, 0.95, 0.95);
	      // pad2->Draw();
	      pad->cd(2);
	      residHist->SetLineColor(kRed);
	      residHist->SetFillColorAlpha(kGreen, 0.5);
	      residHist->Draw("HIST");
	      residHist->GetXaxis()->SetTitle("Residuals");
	      residHist->GetYaxis()->SetTitle("Events");

	      multiCanvas->Update();

	      string canvasMultiName = "WireEnd/Multi/" + string(objName) + result + ".png";
	      multiCanvas->SaveAs(canvasMultiName.c_str());

	      
	      // // Get parameters of the fitted function
	      // Double_t paramsVavilov[fitVavilov->GetNpar()]; // Assuming a Gaussian fit with 3 parameters
	      // fitVavilov->GetParameters(paramsVavilov);

	      // // Calculate residuals
	      // Int_t numBins = hist->GetNbinsX();
	      // Double_t residuals[numBins];
	      // Double_t binCenters[numBins];
	      // Double_t binWidth = hist->GetBinWidth(1); // Assuming uniform bin width
	      // for (Int_t i = 1; i < numBins; ++i)
	      // 	{
	      // 	  binCenters[i] = hist->GetBinCenter(i+1);
	      // 	  cout << binCenters[i] << endl;
	      // 	  residuals[i] = hist->GetBinContent(i+1) - fitVavilov->Eval(binCenters[i+1]);
	      // 	  cout << residuals[i - 1] << endl;
	      // 	}

	     
	      // // // Optionally, normalize residuals by bin uncertainties
	      // // for (Int_t i = 1; i < numBins; ++i)
	      // // 	{
	      // // 	  Double_t binError = hist->GetBinError(i);
	      // // 	  residuals[i-1] /= binError;
	      // // 	}

	      // // // Plot residuals
	      // TCanvas* canvasResid = new TCanvas("canvasResid", "Residuals", 800, 600);
	      // // TGraph* graph = new TGraph(numBins - 1, binCenters, residuals);

	      // TGraph* graph = new TGraph();
	      // for (int i = 1; i < hist->GetNbinsX(); ++i)
	      // 	{
	      // 	  Double_t binCenter = hist->GetBinCenter(i + 1);
	      // 	  Double_t binContent = hist->GetBinContent(i + 1);
	      // 	  Double_t fittedValue = fitVavilov->Eval(binCenter);
	      // 	  Double_t resid = binContent - fittedValue;

	      // 	  graph->SetPoint(i - 1, binCenter, resid);
		  
	      // 	}
	      

	      // string residName = string(objName) + " Residuals" + result;
	      // graph->SetTitle(residName.c_str());
	      // graph->GetXaxis()->SetTitle("X");
	      // graph->GetYaxis()->SetTitle("Residuals");
	      // graph->SetMarkerStyle(20); // Set marker style
	      // graph->SetMarkerSize(.25); // Set marker size
	      // graph->Draw("AP");
	      // // canvasResid->Draw();

	      // string canvasResidName = "WireEnd/Resid/" + string(objName) + result + ".png";
	      
	      // canvasResid->SaveAs(canvasResidName.c_str());

	      // // // New canvas
	      // // TCanvas *canvasResidHist = new TCanvas("canvasResidHist", "Residuals", 800, 600);

	      // // //	      canvasResidHist->cd();
	      // // // Plot residuals as a histogram
	      // // TH1F* residHist = new TH1F("residHist", "Residuals",100,-100,100);//numBins, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); // GetBinLowEdge(2), hist->GetXaxis()->GetBinUpEdge(numBins)); // Skip the first bin
	      // // for (Int_t i = 0; i < numBins - 1; ++i)
	      // // 	{ // Fill the histogram with residuals
	      // // 	  residHist->Fill(residuals[i]);
	      // // 	  // cout << residuals[i] << endl;
	      // // 	}
		    
	      // // 	    //		    cout << residHist << endl;
	      
	      // // // Customize the histogram appearance
	      // // // residHist->SetLineColor(kBlack);
	      // // // residHist->SetFillColor(kBlue);
	      // // // residHist->SetMarkerStyle(20);
	      // // // residHist->SetMarkerSize(1.5);

	      // // // Plot the histogram of residuals
	
	      // // // residHist->Draw("HIST P");
	      // // residHist->Draw();
	      // // residHist->GetXaxis()->SetTitle("X");
	      // // residHist->GetYaxis()->SetTitle("Residuals");
	      // // canvasResidHist->Draw();

	      // // string canvasResidHistName = "WireEnd/Resid/" + string(objName) + "_Hist" + result + ".png";
	      // // canvasResidHist->SaveAs(canvasResidHistName.c_str());
	      
	      // Clean up resources
	      delete fitVavilov;
	      delete fitLandau;
	      delete hist;
	      delete residHist;
	      delete canvas;
	      delete residCanvas;
	      delete multiCanvas;
	      //delete canvasResid;
	      //      delete canvasResidHist;
	    }
	}
    }
  
  // Close the input file and delete
  inputFile->Close();
  delete inputFile;

  // Create a new canvs for the chi2 vs wireEnd
  TCanvas* chi2Canvas = new TCanvas("chi2Canvas", "chi2", 800, 600);
  
  // Plot the chi-square vs wireEnd
  TGraph* graphChi2Landau = new TGraph(vectorChi2Landau.size());
  TGraph* graphChi2Vavilov = new TGraph(vectorChi2Vavilov.size());

  // loop through to plot
  for (int i = 0; i < vectorChi2Vavilov.size(); ++i)
    {
      // customize
      graphChi2Landau->SetMarkerStyle(5); // 52
      graphChi2Landau->SetMarkerColor(kBlue);
      graphChi2Landau->SetLineColor(kRed);
      graphChi2Landau->SetPoint(i, i, vectorChi2Landau[i]);
      
      graphChi2Vavilov->SetMarkerStyle(2); // 52
      graphChi2Vavilov->SetMarkerColor(kGreen);
      graphChi2Vavilov->SetLineColor(kViolet);
      graphChi2Vavilov->SetPoint(i, i, vectorChi2Vavilov[i]);
    }

  graphChi2Landau->SetTitle("WireEndID vs. Chi^{2}");
  graphChi2Landau->GetXaxis()->SetTitle("WireEnd ID");
  graphChi2Landau->GetYaxis()->SetTitle("Chi^{2}");
  //  graphChi2Landau->Draw("AP");
  //  graphChi2Landau->GetYaxis()->SetRangeUser(-100, 5000);

  graphChi2Vavilov->SetTitle("WireEndID vs. Chi^{2}");
  graphChi2Vavilov->GetXaxis()->SetTitle("WireEnd ID");
  graphChi2Vavilov->GetYaxis()->SetTitle("Chi^{2}");
  graphChi2Vavilov->GetYaxis()->SetRangeUser(-10,30);
  graphChi2Vavilov->Draw("AP");
  // graphChi2Vavilov->Draw("psame");

  TLegend* chi2Legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  // chi2Legend->AddEntry(graphChi2Landau, "Landau Chi^{2}", "lp");
  chi2Legend->AddEntry(graphChi2Vavilov, "Vavilov Chi^{2}", "lp");
  chi2Legend->Draw();

  chi2Canvas->Update();
  
  string chi2CanvasName = "WireEnd/Chi2_WireEnd" + result + ".png";
  
  chi2Canvas->SaveAs(chi2CanvasName.c_str());

  delete chi2Canvas;

  
  // // Create a new canvas for the parameter histogram
  // TCanvas* paramCanvas = new TCanvas("paramCanvas", "Parameter Histogram", 800, 600);

  // // Loop over each parameter indexed by j (Each element of the params vector)
  // for (size_t j = 0; j < fitParamsVector[0].size(); ++j)
  //   {
  //     // Create a new histogram for the parameter
  //     TH1F* paramHist = new TH1F("paramHist",
  // 				 Form("Histogram Param %zu",j),
  // 				 100,0,0);//, Form("Histogram Paramater %zu",j));//, 0, 0 ,0);

  //     // Loop through each element of the fitParamsVector indexed by i (Each element of the vector)
  //     for (size_t i = 0; i < fitParamsVector.size(); ++i)
  // 	{
  // 	  // store value in a variable that we will use to fill
  // 	  double param = fitParamsVector[i][j];

  // 	  // Fill the histogram with the fit paramter value
  // 	  paramHist->Fill(param);
  // 	}
  //     // Draw the histogram on the parameter canvas
  //     paramHist->Draw();

  //     // Save the histogram as a PNG file
  //     string canvasParam = "param_" + to_string(j) + "_histogram.png";
  //     paramCanvas->SaveAs(canvasParam.c_str());

  //     // Clear canvas and hist
  //     //paramHist->Clear();
  //     paramCanvas->Clear();
  //     delete paramHist;
  //   }

  // // Clean up resources

  // delete paramCanvas;
}


int main()
{
  // Load necessary library
  gSystem->Load("/users/PAS0654/dcalderon/Research/HELIX/helix-env_2023/helix-flight-software/01build/lib/libDCTDisplayer.so");

  // Call the function to fit histograms
  FitAllHistograms("/users/PAS0654/dcalderon/Research/HELIX/Data/DCTPresenterOut_new_run13710.root", "WireEnd", "RoIADCSum"); //, outPDFName, histogramsPerPage);

  FitAllHistograms("/users/PAS0654/dcalderon/Research/HELIX/Data/DCTPresenterOut_run13738.root", "WireEnd", "RoIADCSum"); //, outPDFName, histogramsPerPage);

  // FitAllHistograms("DCTPresenterOut_run13710.root", "WireEnd_106", "RoIADCSum"); 

  return 0;
}

      
