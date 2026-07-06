// FindHistoInFile.C
//
// Recursively walks every directory in a ROOT file and prints the full path of any
// key whose name contains a_nameSubstring. Written because PicoBinner.cxx's actual
// output directory structure ("yields/<Species>_space", guessed from reading the
// source) didn't match what's really in yieldHistos_OOGeV_proton_2026_07_01.root --
// rather than guess again, just search the whole file.
//
// Usage:
//   root -l -q 'FindHistoInFile.C("/Users/aliggett/data/OO/yieldHistos_OOGeV_proton_2026_07_01.root","zTPCPlus_Proton")'

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TList.h"
#include <iostream>
#include <string>
using namespace std;

void searchDir(TDirectory* a_dir, string a_path, string a_nameSubstring, int& a_nFound){
  TList* keys = a_dir->GetListOfKeys();
  if(!keys) return;
  for(int i = 0; i < keys->GetEntries(); i++){
    TKey* key = (TKey*) keys->At(i);
    string keyName = key->GetName();
    string className = key->GetClassName();
    string fullPath = a_path + "/" + keyName;

    if(keyName.find(a_nameSubstring) != string::npos){
      cout << "  MATCH: " << fullPath << "   (class " << className << ")" << endl;
      a_nFound++;
    }

    if(className == "TDirectoryFile" || className == "TDirectory"){
      TDirectory* subDir = (TDirectory*) a_dir->Get(keyName.c_str());
      if(subDir) searchDir(subDir, fullPath, a_nameSubstring, a_nFound);
    }
  }
}

void FindHistoInFile(string a_file, string a_nameSubstring){
  TFile* f = TFile::Open(a_file.c_str());
  if(!f || f->IsZombie()){ cout << "Could not open " << a_file << endl; return; }

  cout << "Searching " << a_file << " for keys containing \"" << a_nameSubstring << "\"..." << endl;
  int nFound = 0;
  searchDir(f, "", a_nameSubstring, nFound);
  cout << nFound << " match(es) found." << endl;

  if(nFound == 0){
    cout << "No matches. Top-level directory contents:" << endl;
    TList* keys = f->GetListOfKeys();
    for(int i = 0; i < keys->GetEntries(); i++){
      TKey* key = (TKey*) keys->At(i);
      cout << "  " << key->GetName() << "   (class " << key->GetClassName() << ")" << endl;
    }
  }
}
