
//____________________________________________________________________________________________________
// outputFileName: default "ncoll_npart.root" (what doScanX_OO200.csh reads). For the Glauber systematic
// variations write ncoll_npart_<type>.root instead, so the default file is never touched while scans run.
void addNcollVsNpart(
		const Char_t* inputFileList = "./LIST/tree.default.list",
		const Char_t* outputFileName = "ncoll_npart.root"
		)
{
	gSystem->Load("St_base");
	gSystem->Load("StUtilities");
	gSystem->Load("StGlauberUtilities");
	gSystem->Load("StCentralityMaker");
	gSystem->Load("StGlauberTree");
	gSystem->Load("StGlauberAnalysisMaker");

	ifstream fin(inputFileList);
	if(!fin)
	{
		Error("addNcollVsNpart", "can't open %s", inputFileList);
		return;
	}

	TChain t("tree");
	TString file;
	while(fin>>file)
	{
		cout << "Add " << file << " into chain" << endl;
		t.Add(file);
	}

	TH2* hNcoll_Npart = new TH2D
		("hNcoll_Npart", "hNcoll_Npart;Npart;Ncoll;",
		 //      StGlauberConstUtilities::GetNpartBin(), 0, StGlauberConstUtilities::GetNpartMax(),
		 //      StGlauberConstUtilities::GetNcollBin(), 0, StGlauberConstUtilities::GetNcollMax()
		 500, 0, 500, 1800, 0, 1800
		);
	t.Draw("ncoll:npart>>hNcoll_Npart", "npart>=2&&ncoll>=1");

	// Write output ROOT file
	TFile* outputFile = TFile::Open(outputFileName, "recreate");
	hNcoll_Npart->Print();
	hNcoll_Npart->Write();
	outputFile->Close();
}
