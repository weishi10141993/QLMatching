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
#include "TTree.h"
#include "TVector3.h"
#include "TRandom.h"
#include <fcntl.h>

// Larsoft includes (not all might be necessary)
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/OpDetGeo.h"
#include "lardataobj/RecoBase/Hit.h"
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
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
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
    std::string fSignalLabel, fHitLabel, fOpHitLabel, fGEANTLabel, TNuInteraction;
    std::vector<std::string> fLabels, fBackgroundLabels;
    int fDetectorSizeX, fDetectorSizeY, fDetectorSizeZ;

    TTree *fMCTruthTree;
    std::vector<std::map<int, simb::MCParticle>> GeneratorParticles = {};
    int Event, Flag, MNHit, MGen, MTPC, MInd0TPC, MInd1TPC, MInd0NHits, MInd1NHits, MMainID, MMainPDG, MMainParentPDG, TrackNum, OpHitNum, OpFlashNum, MTrackNPoints, MAdjClNum, MSignalAdjClNum, SignalParticlePDG;
    float SignalParticleE, SignalParticleP, SignalParticleK, SignalParticleX, SignalParticleY, SignalParticleZ, SignalParticleTime, MTime, MCharge, MMaxCharge, MInd0Charge, MInd1Charge, MInd0MaxCharge, MInd1MaxCharge;
    float MInd0dTime, MInd1dTime, MInd0RecoY, MInd1RecoY, MRecX, MRecY, MRecZ, MPur, MInd0Pur, MInd1Pur, MGenPur, MMainE, MMainP, MMainK, MMainTime, MMainParentE, MMainParentP, MMainParentK, MMainParentTime, MTrackChi2;
    std::vector<int> TPart, SignalPDGList, SignalPDGDepList, SignalIDList, SignalMotherList, SignalIDDepList, HitNum, SignalElectronDepList;
    std::vector<float> SignalEDepList, SignalXDepList, SignalYDepList, SignalZDepList;
    std::vector<float> OpHitPE, OpHitX, OpHitY, OpHitZ, OpHitTime, OpHitChannel;
    std::vector<float> SignalEList, SignalPList, SignalKList, SignalTimeList, SignalEndXList, SignalEndYList, SignalEndZList, SignalMaxEDepList, SignalMaxEDepXList, SignalMaxEDepYList, SignalMaxEDepZList;
    std::vector<double> MMainVertex, MEndVertex, MMainParentVertex;
    std::vector<double> MTrackStart, MTrackEnd;
    bool MPrimary;

    // --- Maps to hold the geo::TPCID object for each TPCid
    std::map<unsigned int, geo::TPCID> TPCIDMap; // Key is the TPC index, value is the TPCID object
    std::map<unsigned int, float> TPCIDdriftLength; // Key is the TPC index, value is the drift length in cm
    std::map<unsigned int, float> TPCIDdriftTime; // Key is the TPC index, value is the drift time in us

    // --- Histograms to fill about collection plane hits
    float MainElectronEndPointX;
    TH2F *hXTruth;
    TH2F *hYTruth;
    TH2F *hZTruth;
    TH1I *hAdjHits;
    TH1F *hAdjHitsADCInt;
    TH2F *hDriftTime;

    // --- Declare our services
    //art::ServiceHandle<geo::GeometryCore> geom;
    //geo::GeometryCore const* fGeometryService; // pointer to Geometry provider
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
    //fGeometryService = lar::providerFrom<geo::Geometry>();
  }

  //......................................................
  void QLMatchAna::reconfigure(fhicl::ParameterSet const &p)
  {
    fSignalLabel = p.get<std::string>("SignalLabel");
    fBackgroundLabels = p.get<std::vector<std::string>>("BackgroundLabelVector");
    fHitLabel = p.get<std::string>("HitLabel");
    fOpHitLabel = p.get<std::string>("OpHitLabel");
    fGEANTLabel = p.get<std::string>("GEANT4Label");
    fDetectorSizeY = p.get<int>("DetectorSizeX");
    fDetectorSizeY = p.get<int>("DetectorSizeY");
    fDetectorSizeZ = p.get<int>("DetectorSizeZ");

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
    fMCTruthTree = tfs->make<TTree>("MCTruthTree", "MC Truth Tree");

    // MC Truth info.
    fMCTruthTree->Branch("Event", &Event, "Event/I");                                        // Event number
    fMCTruthTree->Branch("Flag", &Flag, "Flag/I");                                           // Flag used to match truth with reco tree entries
    fMCTruthTree->Branch("TruthPart", &TPart);                                               // Number particles per generator
    fMCTruthTree->Branch("Interaction", &TNuInteraction);                                    // True signal interaction process
    fMCTruthTree->Branch("SignalParticleE", &SignalParticleE, "SignalParticleE/F");          // True signal energy [MeV]
    fMCTruthTree->Branch("SignalParticleP", &SignalParticleP, "SignalParticleP/F");          // True signal momentum [MeV]
    fMCTruthTree->Branch("SignalParticleK", &SignalParticleK, "SignalParticleK/F");          // True signal K.E. [MeV]
    fMCTruthTree->Branch("SignalParticleX", &SignalParticleX, "SignalParticleX/F");          // True signal X [cm]
    fMCTruthTree->Branch("SignalParticleY", &SignalParticleY, "SignalParticleY/F");          // True signal Y [cm]
    fMCTruthTree->Branch("SignalParticleZ", &SignalParticleZ, "SignalParticleZ/F");          // True signal Z [cm]
    fMCTruthTree->Branch("SignalParticlePDG", &SignalParticlePDG, "SignalParticlePDG/I");    // True signal PDG
    fMCTruthTree->Branch("SignalParticleTime", &SignalParticleTime, "SignalParticleTime/F"); // True signal time [tick]

    // Save OpHits. (Can be very heavy for background productions)
    fMCTruthTree->Branch("OpHitNum", &OpHitNum, "OpHitNum/I");                               // Number of OpHits
    fMCTruthTree->Branch("OpHitPE", &OpHitPE);           // OpHit PE
    fMCTruthTree->Branch("OpHitX", &OpHitX);             // OpHit X
    fMCTruthTree->Branch("OpHitY", &OpHitY);             // OpHit Y
    fMCTruthTree->Branch("OpHitZ", &OpHitZ);             // OpHit Z
    fMCTruthTree->Branch("OpHitTime", &OpHitTime);       // OpHit Time
    fMCTruthTree->Branch("OpHitChannel", &OpHitChannel); // OpHit Channel

    fMCTruthTree->Branch("HitNum", &HitNum);                                                 // Number of hits in each TPC plane

    // Save Signal Daughters. (Only makes sense for marley)
    fMCTruthTree->Branch("TSignalPDG", &SignalPDGList);         // PDG of Signal marticles
    fMCTruthTree->Branch("TSignalE", &SignalEList);             // Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalP", &SignalPList);             // Energy of Signal momentum [MeV]
    fMCTruthTree->Branch("TSignalK", &SignalKList);             // Kinetik Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalT", &SignalTimeList);          // Time of Signal particles [ticks]
    fMCTruthTree->Branch("TSignalEndX", &SignalEndXList);       // X of Signal particles [cm]
    fMCTruthTree->Branch("TSignalEndY", &SignalEndYList);       // Y of Signal particles [cm]
    fMCTruthTree->Branch("TSignalEndZ", &SignalEndZList);       // Z of Signal particles [cm]
    fMCTruthTree->Branch("TSignalMaxEDep", &SignalMaxEDepList); // Energy of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalX", &SignalMaxEDepXList);      // X of Signal particles [cm]
    fMCTruthTree->Branch("TSignalY", &SignalMaxEDepYList);      // Y of Signal particles [cm]
    fMCTruthTree->Branch("TSignalZ", &SignalMaxEDepZList);      // Z of Signal particles [cm]
    fMCTruthTree->Branch("TSignalID", &SignalIDList);           // TrackID of Signal particles
    fMCTruthTree->Branch("TSignalMother", &SignalMotherList);   // TrackID of Signal mother

    fMCTruthTree->Branch("TSignalPDGDepList", &SignalPDGDepList);           // PDG for Energy deposited of Signal particles
    fMCTruthTree->Branch("TSignalEDepList", &SignalEDepList);               // Energy deposited of Signal particles [MeV]
    fMCTruthTree->Branch("TSignalXDepList", &SignalXDepList);               // X deposited of Signal particles [cm]
    fMCTruthTree->Branch("TSignalYDepList", &SignalYDepList);               // Y deposited of Signal particles [cm]
    fMCTruthTree->Branch("TSignalZDepList", &SignalZDepList);               // Z deposited of Signal particles [cm]
    fMCTruthTree->Branch("TSignalIDDepList", &SignalIDDepList);             // ParentID of Signal particles
    fMCTruthTree->Branch("TSignalElectronDepList", &SignalElectronDepList); // Number of electrons in the Signal particles



    // --- Our Histograms...
    hDriftTime = tfs->make<TH2F>("hDriftTime", "hDriftTime", 100, -400., 400., 100, 0., 10000.);
    hXTruth = tfs->make<TH2F>("hXTruth", "Missmatch in X distance; Distance [cm]; True X position [cm]", 100, -600, 600, 100, -600, 600);
    hYTruth = tfs->make<TH2F>("hYTruth", "Missmatch in Y distance; Distance [cm]; True Y position [cm]", 100, -600, 600, 100, -600, 600);
    hZTruth = tfs->make<TH2F>("hZTruth", "Missmatch in Z distance; Distance [cm]; True Z position [cm]", 100, -600, 600, 100, 0, 1600);
    hAdjHits = tfs->make<TH1I>("hAdjHits", "Number of adjacent collection plane hits; Number of adjacent collection plane hits; Number of events", 21, -0.5, 20.5);
    hAdjHitsADCInt = tfs->make<TH1F>("hAdjHitsADCInt", "Total summed ADC Integrals for clusters; Total summed ADC Integrals for clusters; Number of events", 1000, 0, 10000);
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
        art::FindManyP<simb::MCParticle> Assn(ThisValidHanlde, evt, fGEANTLabel);          // Assign labels to MCPArticles
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
        int NSignalParticles = SignalTruth.NParticles();

        if (fLabels[0] == "marley")
        {
          const simb::MCNeutrino &nue = SignalTruth.GetNeutrino();
          TNuInteraction = std::to_string(nue.InteractionType());
          SignalParticleE = 1e3 * nue.Nu().E();
          SignalParticleP = 1e3 * nue.Nu().P();
          SignalParticleK = 1e3 * nue.Nu().E() - 1e3 * nue.Nu().Mass();
          SignalParticleX = nue.Nu().Vx();
          SignalParticleY = nue.Nu().Vy();
          SignalParticleZ = nue.Nu().Vz();
          SignalParticlePDG = nue.Nu().PdgCode();
          SignalParticleTime = nue.Nu().T();
          sSignalTruth = sSignalTruth + "\nNeutrino Interaction: " + TNuInteraction;
          sSignalTruth = sSignalTruth + "\nNeutrino Energy: " + std::to_string(SignalParticleE) + " MeV";
          sSignalTruth = sSignalTruth + "\nPosition (" + std::to_string(SignalParticleX) + ", " + std::to_string(SignalParticleY) + ", " + std::to_string(SignalParticleZ) + ") cm";
        }
        if (fLabels[0] == "generator")
        {
          sSignalTruth = sSignalTruth + "\nFound generator label: " + fLabels[0] + ". Using single generator config.\n";
          if (NSignalParticles > 1)
          {
            sSignalTruth = sSignalTruth + "\n[WARNING] Multiple particles found in the Signal MCTruth. Using the first one.\n";
          }
          const simb::MCParticle &SignalParticle = SignalTruth.GetParticle(0);
          SignalParticleE = 1e3 * SignalParticle.E();
          SignalParticleP = 1e3 * SignalParticle.P();
          SignalParticleK = 1e3 * SignalParticle.E() - 1e3 * SignalParticle.Mass();
          SignalParticleX = SignalParticle.Vx();
          SignalParticleY = SignalParticle.Vy();
          SignalParticleZ = SignalParticle.Vz();
          SignalParticlePDG = SignalParticle.PdgCode();
          SignalParticleTime = SignalParticle.T();
          std::string sSignalParticle = "";
          if (abs(SignalParticle.PdgCode()) == 12)
          {
            sSignalParticle = "Neutrino";
          }
          else if (abs(SignalParticle.PdgCode()) == 11)
          {
            sSignalParticle = "Electron";
          }
          else if (abs(SignalParticle.PdgCode()) == 22)
          {
            sSignalParticle = "Photon";
          }
          else if (abs(SignalParticle.PdgCode()) == 2112)
          {
            sSignalParticle = "Neutron";
          }
          else
          {
            sSignalParticle = "Other";
          }
          TNuInteraction = "Single " + sSignalParticle;
          sSignalTruth = sSignalTruth + "\n" +  sSignalParticle + " Energy: " + std::to_string(SignalParticleE) + " MeV";
          sSignalTruth = sSignalTruth + "\nPosition (" + std::to_string(SignalParticleX) + ", " + std::to_string(SignalParticleY) + ", " + std::to_string(SignalParticleZ) + ") cm\n";
        }
      }
      art::FindManyP<simb::MCParticle> SignalAssn(Signal, evt, fGEANTLabel);
      sSignalTruth = sSignalTruth + "\nGen.\tPdgCode\t\tEnergy\t\tEndPosition\t\tMother";
      sSignalTruth = sSignalTruth + "\n------------------------------------------------------------------------";

      for (size_t i = 0; i < SignalAssn.size(); i++)
      {
        auto SignalParticles = SignalAssn.at(i);
        for (auto SignalParticle = SignalParticles.begin(); SignalParticle != SignalParticles.end(); SignalParticle++)
        {

          SignalPDGList.push_back((*SignalParticle)->PdgCode());
          SignalEList.push_back(1e3 * (*SignalParticle)->E());
          SignalPList.push_back(1e3 * (*SignalParticle)->P());
          SignalKList.push_back(1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass());
          SignalTimeList.push_back((*SignalParticle)->T());
          SignalEndXList.push_back((*SignalParticle)->EndX());
          SignalEndYList.push_back((*SignalParticle)->EndY());
          SignalEndZList.push_back((*SignalParticle)->EndZ());
          SignalIDList.push_back((*SignalParticle)->TrackId());
          SignalMotherList.push_back((*SignalParticle)->Mother());

          std::map<int, float> SignalMaxEDepMap, SignalMaxEDepXMap, SignalMaxEDepYMap, SignalMaxEDepZMap;
          std::vector<const sim::IDE *> ides = bt_serv->TrackIdToSimIDEs_Ps((*SignalParticle)->TrackId());
          for (auto const &ide : ides)
          {
            if (ide->numElectrons < 1 || ide->energy < 1e-6 || abs(ide->x) > TPCIDdriftLength[0] || abs(ide->y) > fDetectorSizeY || abs(ide->z) > fDetectorSizeZ)
            {
              continue;
            }

            if (ide->energy > SignalMaxEDepMap[(*SignalParticle)->TrackId()])
            {
              SignalMaxEDepMap[(*SignalParticle)->TrackId()] = ide->energy;
              SignalMaxEDepXMap[(*SignalParticle)->TrackId()] = ide->x;
              SignalMaxEDepYMap[(*SignalParticle)->TrackId()] = ide->y;
              SignalMaxEDepZMap[(*SignalParticle)->TrackId()] = ide->z;
            }
            if (abs((*SignalParticle)->PdgCode()) == 11 || abs((*SignalParticle)->PdgCode()) == 22 || abs((*SignalParticle)->PdgCode()) == 2112)
            {
              SignalIDDepList.push_back((*SignalParticle)->TrackId());
              SignalEDepList.push_back(ide->energy);
              SignalPDGDepList.push_back((*SignalParticle)->PdgCode());
              SignalXDepList.push_back(ide->x);
              SignalYDepList.push_back(ide->y);
              SignalZDepList.push_back(ide->z);
              SignalElectronDepList.push_back(ide->numElectrons);
            }
          }
          SignalMaxEDepList.push_back(SignalMaxEDepMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepXList.push_back(SignalMaxEDepXMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepYList.push_back(SignalMaxEDepYMap[(*SignalParticle)->TrackId()]);
          SignalMaxEDepZList.push_back(SignalMaxEDepZMap[(*SignalParticle)->TrackId()]);
          SignalTrackIDs.emplace((*SignalParticle)->TrackId());

          if ((*SignalParticle)->PdgCode() < 1000000)
          {
            sSignalTruth = sSignalTruth + "\n" + fLabels[0] + "\t" + std::to_string((*SignalParticle)->PdgCode()) + "\t\t" + std::to_string(1e3 * (*SignalParticle)->E()) + "\t (" + std::to_string((*SignalParticle)->EndX()) + ", " + std::to_string((*SignalParticle)->EndY()) + ", " + std::to_string((*SignalParticle)->EndZ()) + ")\t" + std::to_string((*SignalParticle)->Mother());
          }
          else
          {
            sSignalTruth = sSignalTruth + "\n" + fLabels[0] + "\t" + std::to_string((*SignalParticle)->PdgCode()) + "\t" + std::to_string(1e3 * (*SignalParticle)->E()) + " (" + std::to_string((*SignalParticle)->EndX()) + ", " + std::to_string((*SignalParticle)->EndY()) + ", " + std::to_string((*SignalParticle)->EndZ()) + ")\t" + std::to_string((*SignalParticle)->Mother());
          }

          if ((*SignalParticle)->PdgCode() == 11) // Electrons
          {
            const TLorentzVector &MainElectronEndPoint = (*SignalParticle)->EndPosition();
            MainElectronEndPointX = MainElectronEndPoint.X();
            ClPartTrackIDs[0].push_back((*SignalParticle)->TrackId());
            mf::LogDebug("QLMatchAna") << "\nMC Electron truth position x = " << MainElectronEndPoint.X() << ", y = " << MainElectronEndPoint.Y() << ", z = " << MainElectronEndPoint.Z();
            mf::LogDebug("QLMatchAna") << "Initial KE " << 1e3 * (*SignalParticle)->E() - 1e3 * (*SignalParticle)->Mass();
          }
          if ((*SignalParticle)->PdgCode() == 22) // Gammas
          {
            ClPartTrackIDs[1].push_back((*SignalParticle)->TrackId());
          }
          if ((*SignalParticle)->PdgCode() == 2112) // Neutrons
          {
            ClPartTrackIDs[2].push_back((*SignalParticle)->TrackId());
          }
          if ((*SignalParticle)->PdgCode() != 11 && (*SignalParticle)->PdgCode() != 22 && (*SignalParticle)->PdgCode() != 2112) // Others
          {
            ClPartTrackIDs[3].push_back((*SignalParticle)->TrackId());
          }
        }
      }
    }
    else
    {
      mf::LogWarning("QLMatchAna") << "No SIGNAL MCTruths found.";
    }

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
    for (int j = 0; j < OpHitNum; j++){
        recob::OpHit OpHit = *OpHitList[j];
        //auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadoutGeom const>()->Get();
        //auto OpHitXYZ = wireReadoutGeom->OpDetGeoFromOpChannel(OpHit.OpChannel()).GetCenter();
        OpHitChannel.push_back(OpHit.OpChannel());
        OpHitTime.push_back(OpHit.PeakTime());
        OpHitPE.push_back(OpHit.PE());
        //OpHitX.push_back(OpHitXYZ.X());
        //OpHitY.push_back(OpHitXYZ.Y());
        //OpHitZ.push_back(OpHitXYZ.Z());
    }


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
  } // end analyze

  //......................................................
  // Reset variables for each event
  void QLMatchAna::ResetVariables()
  {
    SignalParticleE = 0;
    SignalParticleP = 0;
    SignalParticleK = 0;
    SignalParticleX = 0;
    SignalParticleY = 0;
    SignalParticleZ = 0;
    SignalParticlePDG = 0;
    SignalParticleTime = 0;
    OpHitNum = 0;
    OpHitChannel = {}, OpHitPE = {}, OpHitX = {}, OpHitY = {}, OpHitZ = {}, OpHitTime = {};
    HitNum = {};
    SignalElectronDepList = {};
    SignalPDGList = {};
    SignalPDGDepList = {};
    SignalIDList = {}, SignalIDDepList = {};
    SignalMotherList = {};
    SignalEList = {}, SignalPList = {}, SignalKList = {}, SignalTimeList = {}, SignalEndXList = {}, SignalEndYList = {}, SignalEndZList = {};
    SignalEDepList = {}, SignalXDepList = {}, SignalYDepList = {}, SignalZDepList = {};
    SignalMaxEDepList = {}, SignalMaxEDepXList = {}, SignalMaxEDepYList = {}, SignalMaxEDepZList = {};
    TPart = {}, GeneratorParticles = {};
  }
} // namespace solar
DEFINE_ART_MODULE(solar::QLMatchAna)
