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


Double_t Vavilov(Double_t* x, Double_t* par)
{
  return par[4] * TMath::Vavilov(((x[0] - par[2])/par[3]), par[0], par[1]);
}


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
  for (int i = 0; i < (keyList->GetSize()); ++i)
    {
      // Get the key (histogram) at index i
      TKey* key1 = dynamic_cast<TKey*>(keyList->At(i));
      TKey* key2 = dynamic_cast<TKey*>(keyList->At(i+2160+216*2));
      // cout << i << endl;
      if(!key1 && !key2) continue;

      // Check if the key is a histogram
      if ((key1->ReadObj()->IsA() == TH1::Class()) && (key2->ReadObj()->IsA() == TH1::Class()));
	{
	  // Getting the histogram
	  TH1* hist1 = static_cast<TH1*>(key1->ReadObj());
	  TH1* hist2 = static_cast<TH1*>(key2->ReadObj());
	  
	  // Getting the name of the histogram
	  TString objName1 = hist1->GetName();
	  TString objName2 = hist2->GetName();

	  //Check if the object meets the matching pattern that we want`
	  if((objName1.BeginsWith(beginPattern) && objName1.EndsWith(endPattern)) && (objName2.BeginsWith(beginPattern) && objName2.EndsWith(endPattern)))
	    {
	      // Print out the histogram name (sanity check)
	      cout << "WireEnd pair: " << objName1 << " + " << objName2 << endl;

	      // Create canvas for the histogram
	      TCanvas* canvas = new TCanvas("canvas", "Fitted Histogram", 1600, 800);
	      canvas->Divide(2,1);
	      
	      canvas->cd(1);
	      hist1->Draw(); // Draw the original histogram on the canvas
	      hist1->SetStats(0);
	
	      // Determine range of to use for the fits
	      Double_t minX = hist1->GetXaxis()->GetXmin();
	      Double_t maxX = hist1->GetXaxis()->GetXmax();
	      hist1->GetXaxis()->SetRangeUser(minX, 1500); // Set xMin and xMax to the desired range

	      // Define fit function Landau 
	      TF1* fitLandau1 = new TF1("fitLandau1", "landau", 20, 1000);
	       
	      hist1->Fit("fitLandau1", "RQ0");
	      hist1->Fit("fitLandau1", "RWQ");

	      
	      // Could use Landau fit parameters as guesses for Vavilov fit
	      double mpv1 = fitLandau1->GetParameter(1);
	      double width1 = fitLandau1->GetParameter(2);
	      
	      // Plot the fit on histogram
	      // fitLandau->SetNpx(1000000);
	      // fitLandau->Draw("same");

	      // Define fit Vavilov
	      TF1* fitVavilov1 = new TF1("fitVavilov1", Vavilov, 20, 1000, 5);
	      
	      // Name fit parameters
	      fitVavilov1->SetParNames("kappa", "beta2", "shift", "width", "scale");

	      // Set initial fit parameters
	      fitVavilov1->SetParameters(0.05,0.4, mpv1,50.0, 1500.0);
	      // Kappa paramater [0.01, 12.0]
	      fitVavilov1->SetParLimits(0, 0.01, 12.0);
	      // Beta2 paramater [0, 1.0]
	      fitVavilov1->SetParLimits(1, 0.0, 1.0);
	      // Plotting more points 
	      fitVavilov1->SetNpx(1000000);
	      	      
	      // Fit the histogram using the fitting function
	      //hist->Fit("fitVavilov", "WQ");
	      hist1->Fit("fitVavilov1", "RWQ0");
	      hist1->Fit("fitVavilov1","RQ+");

	      // Plot the Vavilov fit
	      fitVavilov1->SetLineColor(3);
	      fitVavilov1->Draw("same");

	      // Calculate chi-squared
	      double chi2_1 = fitVavilov1->GetChisquare();
	      int ndf_1 = fitVavilov1->GetNDF();
	      double reducedChi2Vavilov_1 = (ndf_1 > 0) ? chi2_1/ndf_1 : 0.0;
	      
	      // Create a new legend with fit information and statistics
	      TLegend *newLegend1 = new TLegend(0.7, 0.7, 0.95, 0.95); // Adjust the position as needed
	      // newLegend->AddEntry(hist, "Histogram", "l");
	      newLegend1->SetHeader(result.c_str());
	      newLegend1->AddEntry((TObject*)0, Form("Entries: %.0f", hist1->GetEntries()), "");
	      newLegend1->AddEntry((TObject*)0, Form("Mean: %.2f", hist1->GetMean()), "");
	      newLegend1->AddEntry((TObject*)0, Form("Std Dev: %.2f", hist1->GetStdDev()), "");
	      // newLegend1->AddEntry(fitLandau, Form("Landau #chi_{R}^{2}=%.2f", reducedChi2Landau), "l");
	      newLegend1->AddEntry(fitVavilov1, Form("Vavilov #chi_{R}^{2}=%.2f", reducedChi2Vavilov_1), "l");

	      // Set text size
	      // newLegend1->SetTextSize(0.04);

	      // Draw the new legend
	      newLegend1->Draw();//"same");

	      
	      // second plot
	      canvas->cd(2);
	      hist2->Draw(); // Draw the original histogram on the canvas
	      hist2->GetXaxis()->SetRangeUser(minX, 1500); // Set xMin and xMax to the desired range
	      hist2->SetStats(0);

	      TF1* fitLandau2 = new TF1("fitLandau2", "landau", 20, 1000);
	      hist2->Fit("fitLandau2", "RQ0");
	      hist2->Fit("fitLandau2", "RWQ");

	      double mpv2 = fitLandau2->GetParameter(1);
	      double width2 = fitLandau2->GetParameter(2);

	      TF1* fitVavilov2 = new TF1("fitVavilov2", Vavilov, 20, 1000, 5);

	      fitVavilov2->SetParNames("kappa", "beta2", "shift", "width", "scale");

	      fitVavilov2->SetParameters(0.05,0.4, mpv2,50.0, 1500.0);
	      fitVavilov2->SetParLimits(0, 0.01, 12.0);
	      fitVavilov2->SetParLimits(1, 0.0, 1.0);
	      fitVavilov2->SetNpx(1000000);

	      hist2->Fit("fitVavilov2", "RWQ0");
	      hist2->Fit("fitVavilov2","RQ+");

	      fitVavilov2->SetLineColor(3);
	      fitVavilov2->Draw("same");

	      // Calculate chi-squared
	      double chi2_2 = fitVavilov2->GetChisquare();
	      int ndf_2 = fitVavilov2->GetNDF();
	      double reducedChi2Vavilov_2 = (ndf_2 > 0) ? chi2_2/ndf_2 : 0.0;
	      
	      // Create a new legend with fit information and statistics
	      TLegend *newLegend2 = new TLegend(0.7, 0.7, 0.95, 0.95); // Adjust the position as needed
	      // newLegend->AddEntry(hist, "Histogram", "l");
	      newLegend2->SetHeader(result.c_str());
	      newLegend2->AddEntry((TObject*)0, Form("Entries: %.0f", hist2->GetEntries()), "");
	      newLegend2->AddEntry((TObject*)0, Form("Mean: %.2f", hist2->GetMean()), "");
	      newLegend2->AddEntry((TObject*)0, Form("Std Dev: %.2f", hist2->GetStdDev()), "");
	      // newLegend2->AddEntry(fitLandau, Form("Landau #chi_{R}^{2}=%.2f", reducedChi2Landau), "l");
	      newLegend2->AddEntry(fitVavilov2, Form("Vavilov #chi_{R}^{2}=%.2f", reducedChi2Vavilov_2), "l");

	      // Set text size
	      // newLegend2->SetTextSize(0.04);

	      // Draw the new legend
	      newLegend2->Draw();//"same");

	      // Update the canvas
	      canvas->Update();

	      // Saving canvas 
	      string canvasName = "WireEnd/Pair/Fitted_" + string(objName1) + "_" + string(objName2) + result + ".png";
	      canvas->SaveAs(canvasName.c_str());
	      

	      // // Make a correlation plot of the two wire ends

	      // // normalize the histograms
	      // hist1->Scale(1.0/hist1->Integrate());
	      // hist2->Scale(1.0/hist2->Integrate());

	      // // Create a 2D histogram for the corr
	      // TH2D* correlationPlot = new TH2D("correlationPlot", "Correlation Plot", hist1->GetNbinsX(), hist1->GetXaxis()->GetXmin(), hist1->GetXaxis()->GetXmax(),hist2->GetNbinsX(), hist2->GetXaxis()->GetXmin(), hist2->GetXaxis()->GetXmax());

	      // // Fill the plot
	      // for (int i = 1; i <= hist1->GetNbinsX(); ++i)
	      // 	{
	      // 	  for (int j = 1; j <= hist2->GetNbinsX(); ++j)
	      // 	    {
	      // 	      correlationPlot->SetBinContent(i, j, hist1->GetBinContent(i) * hist2->GetBinContent(j));
	      // 	    }
	      // 	}

	      // // Draw the plot
	      // TCanvas* corrCanvas = new TCanvas("corrCanvas", "Correlation Plot", 800, 600);

	      
	      // Clean up resources
	      delete fitVavilov1;
	      delete fitLandau1;
	      delete hist1;
	      delete fitVavilov2;
	      delete fitLandau2;
	      delete hist2;
	      
	      delete canvas;
	      
	    }
	}
    }
  
  // Close the input file and delete
  inputFile->Close();
  delete inputFile;
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

      
