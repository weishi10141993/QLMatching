////////////////////////////////////////////////////////////////////////////////////
// Class:       QLMatchAna                                                        //
// Module Type: analyzer                                                          //
// File:        QLMatchAna_module.cc                                              //
//                                                                                //
// Written by Wei Shi                                                             //
////////////////////////////////////////////////////////////////////////////////////

// C++ includes
#ifndef QLMatchAna_h
#define QLMatchAna_h

// ROOT includes
#include <TApplication.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TRandom.h"
#include <fcntl.h>

// Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SupernovaTruth.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Art includes and others
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

namespace solar
{
  class QLMatchAna : public art::EDAnalyzer
  {
  public:
    // --- Standard constructor and destructor for an ART module.
    explicit QLMatchAna(fhicl::ParameterSet const &p);
    QLMatchAna(QLMatchAna const &) = delete;
    QLMatchAna(QLMatchAna &&) = delete;
    QLMatchAna &operator=(QLMatchAna const &) = delete;
    QLMatchAna &operator=(QLMatchAna &&) = delete;
    void analyze(art::Event const &evt) override;
    void reconfigure(fhicl::ParameterSet const &p);
    void beginJob() override;

  private:
    // --- Some of our own functions.
    void ResetVariables();

    // --- Our fcl parameter labels for the modules that made the data products
    std::string fSignalLabel, fHitLabel, fOpHitLabel, fGEANT4Label, fIonAndScintLabel;
    std::vector<std::string> fLabels, fBackgroundLabels;
    std::string marleynueInteraction;
    int fDetectorSizeX, fDetectorSizeY, fDetectorSizeZ, fExampleevt, NSignalParticles;

    TTree *fMCTruthTree;
    std::vector<std::map<int, simb::MCParticle>> GeneratorParticles = {};
    int Event, Flag, MNHit, MGen, MTPC, MInd0TPC, MInd1TPC, MInd0NHits, MInd1NHits, MMainID, MMainPDG, MMainParentPDG, TrackNum, OpHitNum, OpFlashNum, MTrackNPoints, MAdjClNum, MSignalAdjClNum, marleynuePDG;
    float marleynueE, marleynueP, marleynueK, marleynueX, marleynueY, marleynueZ, marleynueTime, MTime, MCharge, MMaxCharge, MInd0Charge, MInd1Charge, MInd0MaxCharge, MInd1MaxCharge;
    float MInd0dTime, MInd1dTime, MInd0RecoY, MInd1RecoY, MRecX, MRecY, MRecZ, MPur, MInd0Pur, MInd1Pur, MGenPur, MMainE, MMainP, MMainK, MMainTime, MMainParentE, MMainParentP, MMainParentK, MMainParentTime, MTrackChi2;
    std::vector<int> TPart, SignalPDGList, SignalIDList, SignalMotherList, HitNum;
    std::vector<float> OpHitPE, OpHitX, OpHitY, OpHitZ, OpHitTime, OpHitChannel;
    std::vector<float> SignalEList, SignalPList, SignalKList, SignalTimeList, SignalEndXList, SignalEndYList, SignalEndZList;
    std::vector<double> MMainVertex, MEndVertex, MMainParentVertex;
    std::vector<double> MTrackStart, MTrackEnd;
    bool MPrimary;
    unsigned int Nstep;
    std::vector<float> edepe;
    std::vector<float> edepx;
    std::vector<float> edepy;
    std::vector<float> edepz;
    std::vector<float> edept;
    std::vector<float> num_photons;
    std::vector<float> num_electrons;

    // --- Maps to hold the geo::TPCID object for each TPCid
    std::map<unsigned int, geo::TPCID> TPCIDMap; // Key is the TPC index, value is the TPCID object
    std::map<unsigned int, float> TPCIDdriftLength; // Key is the TPC index, value is the drift length in cm
    std::map<unsigned int, float> TPCIDdriftTime; // Key is the TPC index, value is the drift time in us

    // --- Histograms to fill about collection plane hits
    TH2F *hedeps_zy_noEthres;
    TH2F *hedeps_yx_noEthres;
    TH2F *hedeps_zx_noEthres;
    TH2F *hedeps_zy_100keV;
    TH2F *hedeps_yx_100keV;
    TH2F *hedeps_zx_100keV;
    TH2F *hedeps_zy_200keV;
    TH2F *hedeps_yx_200keV;
    TH2F *hedeps_zx_200keV;
    TH2F *hnuvtx_zy;
    TH2F *hnuvtx_zx;
    TH2F *hnuvtx_yx;
    TH3F *hophit_zyt;
    TH2F *hophit_zy_t0_20us;
    TH2F *hophit_zy_t0_40us;
    TH2F *hophit_yx_t0_20us_posz;
    TH2F *hophit_yx_t0_20us_negz;
    TH2F *hophit_zx_t0_20us_posy;
    TH2F *hophit_zx_t0_20us_negy;

    // --- Declare our services
    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  };
#endif

  //......................................................
  QLMatchAna::QLMatchAna(fhicl::ParameterSet const &p)
      : EDAnalyzer(p)
  {
    this->reconfigure(p);
  }

  //......................................................
  void QLMatchAna::reconfigure(fhicl::ParameterSet const &p)
  {
    fSignalLabel = p.get<std::string>("SignalLabel");
    fBackgroundLabels = p.get<std::vector<std::string>>("BackgroundLabelVector");
    fHitLabel = p.get<std::string>("HitLabel");
    fOpHitLabel = p.get<std::string>("OpHitLabel");
    fGEANT4Label = p.get<std::string>("GEANT4Label");
    fIonAndScintLabel = p.get<std::string>("IonAndScintLabel");
    fDetectorSizeY = p.get<int>("DetectorSizeX");
    fDetectorSizeY = p.get<int>("DetectorSizeY");
    fDetectorSizeZ = p.get<int>("DetectorSizeZ");
    fExampleevt = p.get<int>("Exampleevt");

    // Generate the list of labels to be used in the analysis
    fLabels.push_back(fSignalLabel);
    for (auto const &label : fBackgroundLabels)
    {
      fLabels.push_back(label);
    }
  } // Reconfigure

  //......................................................
  void QLMatchAna::beginJob()
  {

    // --- Make our handle to the TFileService
    art::ServiceHandle<art::TFileService> tfs;
    // Histograms...
    hedeps_zy_noEthres = tfs->make<TH2F>("hedeps_zy_noEthres", "edeps zy view; z [cm]; y [cm]", 2100, 0, 2100, 1400, -700, 700);
    hedeps_yx_noEthres = tfs->make<TH2F>("hedeps_yx_noEthres", "edeps yx view; y [cm]; x [cm]", 1400, -700, 700, 500, -250, 250);
    hedeps_zx_noEthres = tfs->make<TH2F>("hedeps_zx_noEthres", "edeps zx view; z [cm]; x [cm]", 2100, 0, 2100, 500, -250, 250);
    hedeps_zy_100keV   = tfs->make<TH2F>("hedeps_zy_100keV",   "edeps zy view; z [cm]; y [cm]", 2100, 0, 2100, 1400, -700, 700);
    hedeps_yx_100keV   = tfs->make<TH2F>("hedeps_yx_100keV",   "edeps yx view; y [cm]; x [cm]", 1400, -700, 700, 500, -250, 250);
    hedeps_zx_100keV   = tfs->make<TH2F>("hedeps_zx_100keV",   "edeps zx view; z [cm]; x [cm]", 2100, 0, 2100, 500, -250, 250);
    hedeps_zy_200keV   = tfs->make<TH2F>("hedeps_zy_200keV",   "edeps zy view; z [cm]; y [cm]", 2100, 0, 2100, 1400, -700, 700);
    hedeps_yx_200keV   = tfs->make<TH2F>("hedeps_yx_200keV",   "edeps yx view; y [cm]; x [cm]", 1400, -700, 700, 500, -250, 250);
    hedeps_zx_200keV   = tfs->make<TH2F>("hedeps_zx_200keV",   "edeps zx view; z [cm]; x [cm]", 2100, 0, 2100, 500, -250, 250);
    hnuvtx_zy          = tfs->make<TH2F>("hnuvtx_zy",          "edeps zy view; z [cm]; y [cm]", 210, 0, 2100, 140, -700, 700);
    hnuvtx_zx          = tfs->make<TH2F>("hnuvtx_zx",          "edeps zx view; z [cm]; x [cm]", 210, 0, 2100, 50, -250, 250);
    hnuvtx_yx          = tfs->make<TH2F>("hnuvtx_yx",          "edeps yx view; y [cm]; x [cm]", 140, -700, 700, 50, -250, 250);
    hophit_zyt         = tfs->make<TH3F>("hophit_zyt",         "ophit zyt view; z [cm]; y [cm]; t[us]", 42, 0, 2100, 28, -700, 700, 430, -4300, 4300);
    hophit_zy_t0_40us  = tfs->make<TH2F>("hophit_zy_t0_40us",  "opchs zy t0 - t0+40us view; z [cm]; y [cm];", 42, 0, 2100, 28, -700, 700);
    hophit_zy_t0_20us  = tfs->make<TH2F>("hophit_zy_t0_20us",  "opchs zy t0 - t0+20us view; z [cm]; y [cm];", 42, 0, 2100, 28, -700, 700);
    hophit_zx_t0_20us_posy  = tfs->make<TH2F>("hophit_zx_t0_20us_posy",  "opchs zx t0 - t0+20us view; z [cm]; x [cm];", 42, 0, 2100,   5, 0, 250); // only interested in membrane PD in upper vol.
    hophit_zx_t0_20us_negy  = tfs->make<TH2F>("hophit_zx_t0_20us_negy",  "opchs zx t0 - t0+20us view; z [cm]; x [cm];", 42, 0, 2100,   5, 0, 250);
    hophit_yx_t0_20us_posz  = tfs->make<TH2F>("hophit_yx_t0_20us_posz",  "opchs yx t0 - t0+20us view; y [cm]; x [cm];", 28, -700, 700, 5, 0, 250);
    hophit_yx_t0_20us_negz  = tfs->make<TH2F>("hophit_yx_t0_20us_negz",  "opchs yx t0 - t0+20us view; y [cm]; x [cm];", 28, -700, 700, 5, 0, 250);


    fMCTruthTree = tfs->make<TTree>("MCTruthTree", "MC Truth Tree");

    // MC Truth info.
    fMCTruthTree->Branch("Event",                &Event,                 "Event/I");               // Event number
    fMCTruthTree->Branch("Flag",                 &Flag,                  "Flag/I");                // Flag used to match truth with reco tree entries
    fMCTruthTree->Branch("NSignalParticles",     &NSignalParticles,      "NSignalParticles/I");    // Number particles per generator
    fMCTruthTree->Branch("marleynueInteraction", &marleynueInteraction);                           // True signal interaction process
    fMCTruthTree->Branch("marleynueE",           &marleynueE,            "marleynueE/F");          // True signal energy [MeV]
    fMCTruthTree->Branch("marleynueP",           &marleynueP,            "marleynueP/F");          // True signal momentum [MeV]
    fMCTruthTree->Branch("marleynueK",           &marleynueK,            "marleynueK/F");          // True signal K.E. [MeV]
    fMCTruthTree->Branch("marleynueX",           &marleynueX,            "marleynueX/F");          // True signal X [cm]
    fMCTruthTree->Branch("marleynueY",           &marleynueY,            "marleynueY/F");          // True signal Y [cm]
    fMCTruthTree->Branch("marleynueZ",           &marleynueZ,            "marleynueZ/F");          // True signal Z [cm]
    fMCTruthTree->Branch("marleynuePDG",         &marleynuePDG,          "marleynuePDG/I");        // True signal PDG
    fMCTruthTree->Branch("marleynueTime",        &marleynueTime,         "marleynueTime/F");       // True signal time [tick]

    // Save Signal Daughters. (Only makes sense for marley)
    fMCTruthTree->Branch("TSignalPDG",    &SignalPDGList);         // PDG of Signal marticles
    fMCTruthTree->Branch("TSignalE",      &SignalEList);             // Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalP",      &SignalPList);             // Energy of Signal momentum [MeV]
    fMCTruthTree->Branch("TSignalK",      &SignalKList);             // Kinetik Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalT",      &SignalTimeList);          // Time of Signal particles [ticks]
    fMCTruthTree->Branch("TSignalEndX",   &SignalEndXList);       // X of Signal particles [cm]
    fMCTruthTree->Branch("TSignalEndY",   &SignalEndYList);       // Y of Signal particles [cm]
    fMCTruthTree->Branch("TSignalEndZ",   &SignalEndZList);       // Z of Signal particles [cm]
    fMCTruthTree->Branch("TSignalID",     &SignalIDList);           // TrackID of Signal particles
    fMCTruthTree->Branch("TSignalMother", &SignalMotherList);   // TrackID of Signal mother

    // Save OpHits. (Can be very heavy for background productions)
    fMCTruthTree->Branch("OpHitNum",     &OpHitNum, "OpHitNum/I");                               // Number of OpHits
    fMCTruthTree->Branch("OpHitPE",      &OpHitPE);           // OpHit PE
    fMCTruthTree->Branch("OpHitX",       &OpHitX);             // OpHit X
    fMCTruthTree->Branch("OpHitY",       &OpHitY);             // OpHit Y
    fMCTruthTree->Branch("OpHitZ",       &OpHitZ);             // OpHit Z
    fMCTruthTree->Branch("OpHitTime",    &OpHitTime);       // OpHit Time
    fMCTruthTree->Branch("OpHitChannel", &OpHitChannel); // OpHit Channel

    fMCTruthTree->Branch("HitNum", &HitNum);                                                 // Number of hits in each TPC plane


  } // BeginJob

  //......................................................
  void QLMatchAna::analyze(art::Event const &evt)
  {
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------- Prepare everything for new event ----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::vector<std::set<int>> trackids = {};
    std::map<int, simb::MCParticle> ThisGeneratorParts;
    std::vector<recob::Hit> ColHits0, ColHits1, ColHits2, ColHits3;
    std::vector<std::vector<recob::Hit>> ColHits = {ColHits0, ColHits1, ColHits2, ColHits3};
    std::vector<std::vector<recob::Hit>> Clusters0, Clusters1, Clusters2, Clusters3;

    // --- We want to reset all of our previous run and TTree variables ---
    ResetVariables();
    ThisGeneratorParts.clear();
    Event = evt.event();

    //GET INFORMATION ABOUT THE DETECTOR'S GEOMETRY.
    auto const* geo = lar::providerFrom<geo::Geometry>();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
    Flag = rand() % 10000000000;

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---
    std::string sMcTruth = "";
    sMcTruth = sMcTruth + "\nThere are a total of ?? " + " Particles in the event\n";

    // Loop over all signal+bkg handles and collect track IDs
    for (size_t i = 0; i < fLabels.size(); i++)
    {
      GeneratorParticles.push_back(ThisGeneratorParts); // For each label insert empty list

      art::Handle<std::vector<simb::MCTruth>> ThisHandle;
      evt.getByLabel(fLabels[i], ThisHandle);

      if (ThisHandle)
      {
        auto ThisValidHanlde = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[i]); // Get generator handles
        art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde, evt, fGEANT4Label);          // Assign labels to MCPArticles
        //producer->FillMyMaps(GeneratorParticles[i], Assn, ThisValidHanlde);                          // Fill empty list with previously assigned particles
        if (GeneratorParticles[i].size() < 1000)
        {
          sMcTruth = sMcTruth + "\n# of particles ?? "+ "\tfrom gen ?" + " " + fLabels[i];
        }
        else
        {
          sMcTruth = sMcTruth + "\n# of particles ??" + "\tfrom gen ?" + " " + fLabels[i];
        }
        TPart.push_back(GeneratorParticles[i].size());
        if (GeneratorParticles[i].size() > 0)
        {
          for (std::map<int, simb::MCParticle>::iterator iter = GeneratorParticles[i].begin(); iter != GeneratorParticles[i].end(); iter++)
          {
            std::set<int> ThisGeneratorIDs = {};
            trackids.push_back(ThisGeneratorIDs);
            trackids[i].insert(iter->first);
          }
        }
        else
        {
          std::set<int> ThisGeneratorIDs = {};
          trackids.push_back(ThisGeneratorIDs);
        }
      }
      else
      {
        sMcTruth = sMcTruth + "\n# of particles ??" + "\tfrom gen " + " " + fLabels[i] + " *not generated!";
        TPart.push_back(0);
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
      }
    }

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Some MC Truth information -------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::set<int> SignalTrackIDs;                                    // Signal TrackIDs to be used in OpFlash matching
    std::string sSignalTruth = "";
    std::vector<std::vector<int>> ClPartTrackIDs = {{}, {}, {}, {}}; // Track IDs corresponding to each kind of MCTruth particle  {11,22,2112,else}
    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    evt.getByLabel(fLabels[0], ThisHandle);
    if (ThisHandle)
    {
      auto Signal = evt.getValidHandle<std::vector<simb::MCTruth>>(fLabels[0]); // Get handle for SIGNAL MCTruths
      // --- Loop over all neutrinos in the event ---
      for (auto const &SignalTruth : *Signal)
      {
        NSignalParticles = SignalTruth.NParticles();

        if (fLabels[0] == "marley")
        {
          const simb::MCNeutrino &nue = SignalTruth.GetNeutrino();
          marleynueInteraction = std::to_string(nue.InteractionType());
          marleynueE = 1e3 * nue.Nu().E(); // MeV
          marleynueP = 1e3 * nue.Nu().P();
          marleynueK = 1e3 * nue.Nu().E() - 1e3 * nue.Nu().Mass();
          marleynueX = nue.Nu().Vx(); // cm
          marleynueY = nue.Nu().Vy();
          marleynueZ = nue.Nu().Vz();
          marleynuePDG = nue.Nu().PdgCode();
          marleynueTime = nue.Nu().T();
        }
        else {
          mf::LogError("QLMatchAna") << "Not a marley generated signal";
        } // end marley

      } //end SignalTruth

      // Find all secondary particles in g4 (simb::MCParticle) associated with the mc truth signal
      art::FindManyP<simb::MCParticle> SignalAssn(Signal, evt, fGEANT4Label);

      for (size_t i = 0; i < SignalAssn.size(); i++)
      {
        auto SignalParticles = SignalAssn.at(i);
        for (auto SignalParticle = SignalParticles.begin(); SignalParticle != SignalParticles.end(); SignalParticle++)
        {

          SignalPDGList.push_back((*SignalParticle)->PdgCode());
          SignalEList.push_back(1e3 * (*SignalParticle)->E()); // MeV
          SignalPList.push_back(1e3 * (*SignalParticle)->P());
          SignalKList.push_back(1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass());
          SignalTimeList.push_back((*SignalParticle)->T());
          SignalEndXList.push_back((*SignalParticle)->EndX());
          SignalEndYList.push_back((*SignalParticle)->EndY());
          SignalEndZList.push_back((*SignalParticle)->EndZ());
          SignalIDList.push_back((*SignalParticle)->TrackId());
          SignalMotherList.push_back((*SignalParticle)->Mother());

          if (Event == fExampleevt) std::cout << "G4 assn with true signal i=: " << i << ", pdg: "<< (*SignalParticle)->PdgCode() << ", e: " << (*SignalParticle)->E() << ", trkid:" << (*SignalParticle)->TrackId() << ", mother trkid: " << (*SignalParticle)->Mother() << std::endl;

        } // end daughters loop
      } // end signal assn
    } // end handle
    else
    {
      mf::LogWarning("QLMatchAna") << "No SIGNAL MCTruths found.";
    }


    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------- Access energy deposits at G4 stage 1 and make a plot -----------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    auto simedepHandle = evt.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintLabel);
    std::vector<art::Ptr<sim::SimEnergyDeposit>> simedeplist;

    if (simedepHandle.isValid())
      art::fill_ptr_vector(simedeplist, simedepHandle);

    Nstep = 0;
    edepe.clear();
    edepx.clear();
    edepy.clear();
    edepz.clear();
    edept.clear();
    num_photons.clear();
    num_electrons.clear();

    Nstep = simedeplist.size();
    std::cout << "# of edeps: " << Nstep << std::endl;

    edepe.resize(Nstep);
    edepx.resize(Nstep);
    edepy.resize(Nstep);
    edepz.resize(Nstep);
    edept.resize(Nstep);
    num_photons.resize(Nstep);
    num_electrons.resize(Nstep);

    if (Event == fExampleevt) {

      for (size_t i=0; i<Nstep; i++) {
        edepe[i] = simedeplist[i]->Energy(); // MeV
        edepx[i] = simedeplist[i]->X(); // cm, midpoint
        edepy[i] = simedeplist[i]->Y();
        edepz[i] = simedeplist[i]->Z();
        edept[i] = simedeplist[i]->T(); // ns
        num_photons[i]   = simedeplist[i]->NumPhotons();
        num_electrons[i] = simedeplist[i]->NumElectrons();

        hedeps_zy_noEthres->Fill(edepz[i], edepy[i], edepe[i]); // x is drift
        hedeps_yx_noEthres->Fill(edepy[i], edepx[i], edepe[i]);
        hedeps_zx_noEthres->Fill(edepz[i], edepx[i], edepe[i]);
        if (edepe[i] > 0.1) {
          hedeps_zy_100keV->Fill(edepz[i], edepy[i], edepe[i]); // x is drift
          hedeps_yx_100keV->Fill(edepy[i], edepx[i], edepe[i]);
          hedeps_zx_100keV->Fill(edepz[i], edepx[i], edepe[i]);
        } // 100 keV thres
        if (edepe[i] > 0.2) {
          hedeps_zy_200keV->Fill(edepz[i], edepy[i], edepe[i]); // x is drift
          hedeps_yx_200keV->Fill(edepy[i], edepx[i], edepe[i]);
          hedeps_zx_200keV->Fill(edepz[i], edepx[i], edepe[i]);
        } // 200 keV thres
      } // end loop edeps

      // Fill the neutrino true position
      hnuvtx_zy->Fill(marleynueZ, marleynueY, 1000); // put high weight so its visible
      hnuvtx_zx->Fill(marleynueZ, marleynueX, 1000);
      hnuvtx_yx->Fill(marleynueY, marleynueX, 1000);

    } // end example evt


    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------- OpHit analysis -------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//

    std::vector<art::Ptr<recob::OpHit>> OpHitList;
    art::Handle<std::vector<recob::OpHit>> OpHitHandle;
    if (evt.getByLabel(fOpHitLabel, OpHitHandle))
    {
      art::fill_ptr_vector(OpHitList, OpHitHandle);
    }
    OpHitNum = int(OpHitList.size());
    int ophitch = 0;
    double ophittime = 0;
    double ophitpe = 0;
    double ophitx = 0;
    double ophity = 0;
    double ophitz = 0;
    std::cout << "OpHitNum: " << OpHitNum << std::endl;
    for (int j = 0; j < OpHitNum; j++){
        recob::OpHit OpHit = *OpHitList[j];
        ophitch = OpHit.OpChannel();
        ophittime = OpHit.PeakTime();
        ophitpe = OpHit.PE();
        auto const OpHitXYZ = geo->OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
        ophitx = OpHitXYZ.X();
        ophity = OpHitXYZ.Y();
        ophitz = OpHitXYZ.Z();
        OpHitChannel.push_back(ophitch);
        OpHitTime.push_back(ophittime);
        OpHitPE.push_back(ophitpe);
        OpHitX.push_back(ophitx);
        OpHitY.push_back(ophity);
        OpHitZ.push_back(ophitz);
        if (Event == fExampleevt) {
          hophit_zyt->Fill(ophitz, ophity, ophittime, ophitpe); // weight by pe
          if (ophittime >= 0 && ophittime < 40) hophit_zy_t0_40us->Fill(ophitz, ophity, ophitpe);
          if (ophittime >= 0 && ophittime < 20) {
            hophit_zy_t0_20us->Fill(ophitz, ophity, ophitpe);
            if (ophitz > 0) hophit_yx_t0_20us_posz->Fill(ophity, ophitx, ophitpe);
            if (ophitz < 0) hophit_yx_t0_20us_negz->Fill(ophity, ophitx, ophitpe);
            if (ophity > 0) hophit_zx_t0_20us_posy->Fill(ophitz, ophitx, ophitpe);
            if (ophity < 0) hophit_zx_t0_20us_negy->Fill(ophitz, ophitx, ophitpe);
          }
          if (j%1000 == 0) std::cout << "finished filling: " << j << "/"<< OpHitNum << std::endl;
          //std::cout << "OpHit Num #" << j << ": Channel="<< ophitch << ", PeakTime = " << ophittime << ", PE: " << ophitpe << " X = " << ophitx << ", Y = "<< ophity << ", Z = "<< ophitz << std::endl;
        }
    }
    std::cout << "vtx x [cm]: " << marleynueX << ", vtx y [cm]: " << marleynueY << ", vtx z [cm]: " << marleynueZ << std::endl;
    std::cout << "Event: "<< Event << ", NSignalParticles: " << NSignalParticles << ", marleynueE: " << marleynueE << std::endl;


    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------- Hit collection and assignment ----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Lift out the reco hits:
    auto RecoHits = evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    int NTotHits = RecoHits->size();

    for (int i = 0; i < NTotHits; ++i)
    {
      // --- Loop over the reconstructed hits to separate them among tpc planes according to view and signal type
      recob::Hit const &ThisHit = RecoHits->at(i);
      if (ThisHit.PeakTime() < 0) mf::LogWarning("QLMatchAna") << "Negtive hit time";

      mf::LogDebug("QLMatchAna") << "Hit " << i << " has view " << ThisHit.View() << " and signal type " << ThisHit.SignalType();

      if (ThisHit.SignalType() == 0 && ThisHit.View() == 0)
      {
        ColHits0.push_back(ThisHit);
      } // SignalType = 0
      else if (ThisHit.SignalType() == 0 && ThisHit.View() == 1)
      {
        ColHits1.push_back(ThisHit);
      } // SignalType = 0
      else if (ThisHit.SignalType() == 1)
      {
        ColHits2.push_back(ThisHit);
      } // SignalType = 1
      else
      {
        ColHits3.push_back(ThisHit);
        mf::LogError("QLMatchAna") << "Hit was found with view out of scope";
      }
    }

    fMCTruthTree->Fill();
  } // end analyze

  //......................................................
  // Reset variables for each event
  void QLMatchAna::ResetVariables()
  {
    marleynueE = 0;
    marleynueP = 0;
    marleynueK = 0;
    marleynueX = 0;
    marleynueY = 0;
    marleynueZ = 0;
    marleynuePDG = 0;
    NSignalParticles = 0;
    marleynueTime = 0;
    OpHitNum = 0;
    OpHitChannel = {}, OpHitPE = {}, OpHitX = {}, OpHitY = {}, OpHitZ = {}, OpHitTime = {};
    HitNum = {};
    SignalPDGList = {};
    SignalIDList = {};
    SignalMotherList = {};
    SignalEList = {}, SignalPList = {}, SignalKList = {}, SignalTimeList = {}, SignalEndXList = {}, SignalEndYList = {}, SignalEndZList = {};
    TPart = {}, GeneratorParticles = {};
    edepe = {}, edepx = {}, edepy = {}, edepz = {}, edept = {}, num_photons = {}, num_electrons = {};
  }
} // namespace solar
DEFINE_ART_MODULE(solar::QLMatchAna)
