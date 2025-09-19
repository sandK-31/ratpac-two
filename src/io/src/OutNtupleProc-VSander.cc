#include <TFile.h>
#include <TTimeStamp.h>
#include <TTree.h>
#include <TVector3.h>

#include <RAT/DB.hh>
#include <RAT/DS/DigitPMT.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/MCSummary.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/WaveformAnalysisResult.hh>
#include <RAT/Log.hh>
#include <RAT/OutNtupleProc.hh>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "RAT/DS/ChannelStatus.hh"

namespace RAT {

OutNtupleProc::OutNtupleProc() : Processor("outntuple") {
  outputFile = nullptr;
  outputTree = nullptr;
  metaTree = nullptr;
  runBranch = new DS::Run();
  done_writing_calib = false;

  // Load options from the database
  DB* db = DB::Get();
  DBLinkPtr table = db->GetLink("IO", "NtupleProc");
  try {
    defaultFilename = table->GetS("default_output_filename");
    if (defaultFilename.find(".") == std::string::npos) {
      defaultFilename += ".ntuple.root";
    }
  } catch (DBNotFoundError& e) {
    defaultFilename = "output.ntuple.root";
  }
}

bool OutNtupleProc::OpenFile(std::string filename) {
  outputFile = TFile::Open(filename.c_str(), "RECREATE");
  // Meta Tree
  metaTree = new TTree("meta", "meta");
  metaTree->Branch("runId", &runId);
  metaTree->Branch("runType", &runType);
  metaTree->Branch("runTime", &runTime);
  metaTree->Branch("dsentries", &dsentries);
  metaTree->Branch("macro", &macro);
  metaTree->Branch("pmtType", &pmtType);
  metaTree->Branch("pmtId", &pmtId);
  metaTree->Branch("pmtChannel", &pmtChannel);
  metaTree->Branch("pmtIsOnline", &pmtIsOnline);
  metaTree->Branch("pmtCableOffset", &pmtCableOffset);
  metaTree->Branch("pmtChargeScale", &pmtChargeScale);
  metaTree->Branch("pmtPulseWidthScale", &pmtPulseWidthScale);
  metaTree->Branch("pmtX", &pmtX);
  metaTree->Branch("pmtY", &pmtY);
  metaTree->Branch("pmtZ", &pmtZ);
  metaTree->Branch("pmtU", &pmtU);
  metaTree->Branch("pmtV", &pmtV);
  metaTree->Branch("pmtW", &pmtW);
  metaTree->Branch("digitizerWindowSize", &digitizerWindowSize);
  metaTree->Branch("digitizerSampleRate_GHz", &digitizerSampleRate);
  metaTree->Branch("digitizerDynamicRange_mV", &digitizerDynamicRange);
  metaTree->Branch("digitizerResolution_mVPerADC", &digitizerVoltageResolution);

  this->AssignAdditionalMetaAddresses();
  dsentries = 0;
  // Data Tree
  outputTree = new TTree("output", "output");
  // These are the *first* particles MC positions, directions, and time

  // Save particle tracking information
  outputTree->Branch("trackPDG", &trackPDG);
  outputTree->Branch("trackID", &trackID);
  outputTree->Branch("trackParentID", &trackParentID);
  outputTree->Branch("trackPosX", &trackPosX);
  outputTree->Branch("trackPosY", &trackPosY);
  outputTree->Branch("trackPosZ", &trackPosZ);
  outputTree->Branch("trackMomX", &trackMomX);
  outputTree->Branch("trackMomY", &trackMomY);
  outputTree->Branch("trackMomZ", &trackMomZ);
  outputTree->Branch("trackEdep", &trackEdep);
  outputTree->Branch("trackKE", &trackKE);
  outputTree->Branch("trackTime", &trackTime);
  outputTree->Branch("trackProcess", &trackProcess);
  metaTree->Branch("processCodeMap", &processCodeMap);
  outputTree->Branch("trackVolume", &trackVolume);
  metaTree->Branch("volumeCodeMap", &volumeCodeMap);

  this->AssignAdditionalAddresses();

  // At the start of your simulation
  auto* volumeStore = G4PhysicalVolumeStore::GetInstance();

  for (auto* vol : *volumeStore) {
    std::string volName = vol->GetName();
    if (volumeCodeMap.find(volName) == volumeCodeMap.end()) {
      volumeCodeMap[volName] = volumeCodeMap.size();
      volumeCodeIndex.push_back(volumeCodeMap.size() - 1);
      volumeName.push_back(volName);
    }
  }

  // Define the list of processes we care about
  std::vector<std::string> predefinedProcesses = {"Attenuation",
                                                  "B+Inelastic",
                                                  "B-Inelastic",
                                                  "B0Inelastic",
                                                  "Bc+Inelastic",
                                                  "Bc-Inelastic",
                                                  "Bs0Inelastic",
                                                  "Cerenkov",
                                                  "CoulombScat",
                                                  "D+Inelastic",
                                                  "D-Inelastic",
                                                  "D0Inelastic",
                                                  "Decay",
                                                  "Ds+Inelastic",
                                                  "Ds-Inelastic",
                                                  "G4FastSimulationManagerProcess",
                                                  "GammaGeneralProc",
                                                  "He3Inelastic",
                                                  "OpBoundary",
                                                  "OpRayleigh",
                                                  "RadioactiveDecay",
                                                  "Rayl",
                                                  "Transportation",
                                                  "alphaInelastic",
                                                  "annihil",
                                                  "anti_B0Inelastic",
                                                  "anti_Bs0Inelastic",
                                                  "anti_D0Inelastic",
                                                  "anti_He3Inelastic",
                                                  "anti_alphaInelastic",
                                                  "anti_deuteronInelastic",
                                                  "anti_lambdaInelastic",
                                                  "anti_lambda_bInelastic",
                                                  "anti_lambda_c+Inelastic",
                                                  "anti_neutronInelastic",
                                                  "anti_omega-Inelastic",
                                                  "anti_omega_b-Inelastic",
                                                  "anti_omega_c0Inelastic",
                                                  "anti_protonInelastic",
                                                  "anti_sigma+Inelastic",
                                                  "anti_sigma-Inelastic",
                                                  "anti_tritonInelastic",
                                                  "anti_xi-Inelastic",
                                                  "anti_xi0Inelastic",
                                                  "anti_xi_b-Inelastic",
                                                  "anti_xi_b0Inelastic",
                                                  "anti_xi_c+Inelastic",
                                                  "anti_xi_c0Inelastic",
                                                  "compt",
                                                  "dInelastic",
                                                  "eBrem",
                                                  "eIoni",
                                                  "electronNuclear",
                                                  "hBertiniCaptureAtRest",
                                                  "hBrems",
                                                  "hFritiofCaptureAtRest",
                                                  "hIoni",
                                                  "hPairProd",
                                                  "hadElastic",
                                                  "ionElastic",
                                                  "ionInelastic",
                                                  "ionIoni",
                                                  "kaon+Inelastic",
                                                  "kaon-Inelastic",
                                                  "kaon0LInelastic",
                                                  "kaon0SInelastic",
                                                  "lambdaInelastic",
                                                  "lambda_bInelastic",
                                                  "lambda_c+Inelastic",
                                                  "msc",
                                                  "muBrems",
                                                  "muIoni",
                                                  "muMinusCaptureAtRest",
                                                  "muPairProd",
                                                  "muonNuclear",
                                                  "nCapture",
                                                  "nFission",
                                                  "neutronInelastic",
                                                  "omega-Inelastic",
                                                  "omega_b-Inelastic",
                                                  "omega_c0Inelastic",
                                                  "phot",
                                                  "pi+Inelastic",
                                                  "pi-Inelastic",
                                                  "positronNuclear",
                                                  "protonInelastic",
                                                  "sigma+Inelastic",
                                                  "sigma-Inelastic",
                                                  "start",
                                                  "tInelastic",
                                                  "xi-Inelastic",
                                                  "xi0Inelastic",
                                                  "xi_b-Inelastic",
                                                  "xi_b0Inelastic",
                                                  "xi_c+Inelastic",
                                                  "xi_c0Inelastic"};

  // Fill the process maps directly
  for (const auto& procName : predefinedProcesses) {
    if (processCodeMap.find(procName) == processCodeMap.end()) {
      processCodeMap[procName] = processCodeMap.size();
      processCodeIndex.push_back(processCodeMap.size() - 1);
      processName.push_back(procName);
    }
  }

  return true;
}

Processor::Result OutNtupleProc::DSEvent(DS::Root* ds) {
  if (!this->outputFile) {
    if (!OpenFile(this->defaultFilename.c_str())) {
      Log::Die("No output file specified");
    }
  }
  DS::MC* mc = ds->GetMC();
  runBranch = DS::RunStore::GetRun(ds);
  dsentries++;
  std::map<std::string, double> edep_per_volume;

  int nTracks = mc->GetMCTrackCount();
  // Clear previous event
  edep_per_volume.clear();
  trackPDG.clear();
  trackID.clear();
  trackParentID.clear();
  trackPosX.clear();
  trackPosY.clear();
  trackPosZ.clear();
  trackMomX.clear();
  trackMomY.clear();
  trackMomZ.clear();
  trackKE.clear();
  trackEdep.clear();
  trackEdep.clear();
  trackTime.clear();
  trackProcess.clear();
  trackVolume.clear();
  std::vector<double> xtrack, ytrack, ztrack;
  std::vector<double> pxtrack, pytrack, pztrack;
  std::vector<double> kinetic, globaltime, energy_deposited;
  std::vector<int> processMapID;
  std::vector<int> volumeMapID;
  double edep;

  for (int trk = 0; trk < nTracks; trk++) {
    DS::MCTrack* track = mc->GetMCTrack(trk);
    trackPDG.push_back(track->GetPDGCode());
    trackID.push_back(track->GetID());
    trackParentID.push_back(track->GetParentID());
    xtrack.clear();
    ytrack.clear();
    ztrack.clear();
    pxtrack.clear();
    pytrack.clear();
    pztrack.clear();
    kinetic.clear();
    energy_deposited.clear();
    globaltime.clear();
    processMapID.clear();
    volumeMapID.clear();
    int nSteps = track->GetMCTrackStepCount();
    // You can now get the total energy deposited in, e.g., "detector" as:

    for (int stp = 0; stp < nSteps; stp++) {
      DS::MCTrackStep* step = track->GetMCTrackStep(stp);
      // Process
      std::string proc = step->GetProcess();
      // Volume
      std::string vol = step->GetVolume();

      edep = step->GetDepositedEnergy();
      volumeMapID.push_back(volumeCodeMap.at(vol));
      processMapID.push_back(processCodeMap.at(proc));
      TVector3 tv = step->GetEndpoint();
      TVector3 momentum = step->GetMomentum();
      energy_deposited.push_back(edep);
      kinetic.push_back(step->GetKE());
      globaltime.push_back(step->GetGlobalTime());
      xtrack.push_back(tv.X());
      ytrack.push_back(tv.Y());
      ztrack.push_back(tv.Z());
      pxtrack.push_back(momentum.X());
      pytrack.push_back(momentum.Y());
      pztrack.push_back(momentum.Z());
      // or step->GetTotalEnergyDeposit(
      edep_per_volume[vol] += edep;
    }

    trackKE.push_back(kinetic);
    trackEdep.push_back(energy_deposited);
    trackTime.push_back(globaltime);
    trackPosX.push_back(xtrack);
    trackPosY.push_back(ytrack);
    trackPosZ.push_back(ztrack);
    trackMomX.push_back(pxtrack);
    trackMomY.push_back(pytrack);
    trackMomZ.push_back(pztrack);
    trackProcess.push_back(processMapID);
    trackVolume.push_back(volumeMapID);
  }

  /*
  if (edep_per_volume["cebr"] != 0 || edep_per_volume["scintillator"] != 0) {
    outputTree->Fill();
  }
  */
  outputTree->Fill();

  return Processor::OK;
}

void OutNtupleProc::EndOfRun(DS::Run* run) {
  if (outputFile) {
    G4cout << "Closing output file " << outputFile->GetName() << std::endl;
    outputFile->cd();
    runId = runBranch->GetID();
    runType = runBranch->GetType();
    // Converting to unix time
    TTimeStamp rootTime = runBranch->GetStartTime();
    runTime = TTimeStamp_to_UnixTime(rootTime.GetSec());
    macro = Log::GetMacro();
    FillMeta();
    metaTree->Fill();
    metaTree->Write();
    outputTree->Write();

    /*
    TMap* dbtrace = Log::GetDBTraceMap();
    dbtrace->Write("db", TObject::kSingleKey);
    */
    // outputFile->Write(0, TObject::kOverwrite);
    outputFile->Close();
    delete outputFile;
  }
}

void OutNtupleProc::SetS(std::string param, std::string value) {
  if (param == "file") {
    this->defaultFilename = value;
  }
}

void OutNtupleProc::SetI(std::string param, int value) {
  if (param == "include_tracking") {
    options.tracking = value ? true : false;
  }
  if (param == "include_mcparticles") {
    options.mcparticles = value ? true : false;
  }
  if (param == "include_pmthits") {
    options.pmthits = value ? true : false;
  }
  if (param == "include_nestedtubehits") {
    options.nthits = value ? true : false;
  }
  if (param == "include_untriggered_events") {
    options.untriggered = value ? true : false;
  }
  if (param == "include_mchits") {
    options.mchits = value ? true : false;
  }
  if (param == "include_digitizerwaveforms") {
    options.digitizerwaveforms = value ? true : false;
  }
  if (param == "include_digitizerhits") {
    options.digitizerhits = value ? true : false;
  }
  if (param == "include_digitizerfits") {
    options.digitizerfits = value ? true : false;
  }
}

}  // namespace RAT
