#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include "AthenaBaseComps/AthAlgorithm.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"
#include "xAODTrigMuon/L2CombinedMuonContainer.h"
#include "xAODTrigger/MuonRoIContainer.h"
#include "xAODEventInfo/EventInfo.h"
#include "TruthClassification/TruthClassificationTool.h"
#include "MuonAnalysisInterfaces/IMuonSelectionTool.h"
#include <TrigDecisionTool/TrigDecisionTool.h>
#include <TrkExInterfaces/IExtrapolator.h>
#include "StoreGate/ReadHandleKey.h"
#include "GaudiKernel/ITHistSvc.h"

class MyxAODAnalysis : public AthAlgorithm
{
 public:
  // this is a standard algorithm constructor
  MyxAODAnalysis (const std::string& name, ISvcLocator* pSvcLocator);

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;

 private:
  ServiceHandle <ITHistSvc> m_tHistSvc;

  ToolHandle<Trig::TrigDecisionTool> m_trigDecisionTool;

  // extrapolation tool
  ToolHandle<Trk::IExtrapolator> m_extrapolator;

  ToolHandle<CP::IMuonSelectionTool> m_selTool;
  ToolHandle<CP::IClassificationTool> m_truthClassificationTool;

  SG::ReadHandleKey<xAOD::EventInfo> m_eventInfo_key{this, "EventInfoName", "EventInfo", "event info key"};
  SG::ReadHandleKey<xAOD::MuonContainer> m_muon_key{this, "MuonContName", "Muons", "muon key"};

  SG::ReadHandleKey<xAOD::L2StandAloneMuonContainer> m_L2SAKey{this, "xAODL2SAContainer", "HLT_MuonL2SAInfo", "Name of L2SAContainer object"};
  SG::ReadHandleKey<xAOD::L2StandAloneMuonContainer> m_L2SAIOKey{this, "xAODL2SAContainer_IO", "HLT_MuonL2SAInfoIOmode", "Name of L2SAIOContainer object"};
  SG::ReadHandleKey<xAOD::L2CombinedMuonContainer> m_L2CBKey{this, "xAODL2CBContainer", "HLT_MuonL2CBInfo", "Name of L2CBContainer object"};
  SG::ReadHandleKey<xAOD::L2CombinedMuonContainer> m_L2CBIOKey{this, "xAODL2CBContainer_IO", "HLT_MuonL2CBInfoIOmode", "Name of L2CBIOContainer object"};

  TTree* mytree;

  unsigned int m_runNumber = 0;
  unsigned int m_lumiBlock = 0;
  unsigned long long m_eventNumber = 0;
  float m_mcEventWeight = 0;
  float m_actualInteractionsPerCrossing = 0;
  float m_averageInteractionsPerCrossing = 0;

  std::vector<float> m_extZposition;
  std::vector<float> m_b_extPosition;

  std::unique_ptr<std::vector<float>> m_muon_e;
  std::unique_ptr<std::vector<float>> m_muon_pt;
  std::unique_ptr<std::vector<float>> m_muon_eta;
  std::unique_ptr<std::vector<float>> m_muon_phi;
  std::unique_ptr<std::vector<float>> m_muon_charge;
  std::unique_ptr<std::vector<int>> m_muon_quality;
  std::unique_ptr<std::vector<bool>> m_muon_isBadMuon_other;
  std::unique_ptr<std::vector<int>> m_muon_truthType;
  std::unique_ptr<std::vector<int>> m_muon_truthTypeMCTC;
  std::unique_ptr<std::vector<int>> m_muon_truthOriginMCTC;
  std::unique_ptr<std::vector<int>> m_muon_IFFtruthType;
  std::unique_ptr<std::vector<int>> m_muon_muonType;
  std::unique_ptr<std::vector<int>> m_muon_author;

  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_targetPlaneVec;
  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_targetEtaVec;
  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_targetPhiVec;
  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_targetPxVec;
  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_targetPyVec;
  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_targetPzVec;

  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_b_targetPlaneVec;
  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_b_targetEtaVec;
  std::unique_ptr<std::vector< std::vector< float>>> m_mu_ext_b_targetPhiVec;

  std::unique_ptr<std::vector<float>> m_truthmuon_e;
  std::unique_ptr<std::vector<float>> m_truthmuon_pt;
  std::unique_ptr<std::vector<float>> m_truthmuon_eta;
  std::unique_ptr<std::vector<float>> m_truthmuon_phi;
  std::unique_ptr<std::vector<int>> m_truthmuon_pdgId;
  std::unique_ptr<std::vector<float>> m_muon_mePt;
  std::unique_ptr<std::vector<float>> m_muon_idPt;
  std::unique_ptr<std::vector<float>> m_muon_cbPt;
  std::unique_ptr<std::vector<float>> m_muon_meP;
  std::unique_ptr<std::vector<float>> m_muon_idP;
  std::unique_ptr<std::vector<float>> m_muon_etaMS;
  std::unique_ptr<std::vector<float>> m_muon_phiMS;

  std::unique_ptr<std::vector<int>> m_muon_innerSmallHits;
  std::unique_ptr<std::vector<int>> m_muon_innerLargeHits;
  std::unique_ptr<std::vector<int>> m_muon_middleSmallHits;
  std::unique_ptr<std::vector<int>> m_muon_middleLargeHits;
  std::unique_ptr<std::vector<int>> m_muon_outerSmallHits;
  std::unique_ptr<std::vector<int>> m_muon_outerLargeHits;
  std::unique_ptr<std::vector<int>> m_muon_extendedSmallHits;
  std::unique_ptr<std::vector<int>> m_muon_extendedLargeHits;

  // sTGC
  std::unique_ptr<std::vector< uint8_t>> m_muon_phiLayer1STGCHits;
  std::unique_ptr<std::vector< uint8_t>> m_muon_phiLayer2STGCHits;
  std::unique_ptr<std::vector< uint8_t>> m_muon_etaLayer1STGCHits;
  std::unique_ptr<std::vector< uint8_t>> m_muon_etaLayer2STGCHits;
  std::unique_ptr<std::vector< uint8_t>> m_muon_phiLayer1STGCHoles;
  std::unique_ptr<std::vector< uint8_t>> m_muon_phiLayer2STGCHoles;
  std::unique_ptr<std::vector< uint8_t>> m_muon_etaLayer1STGCHoles;
  std::unique_ptr<std::vector< uint8_t>> m_muon_etaLayer2STGCHoles;

  // MM
  std::unique_ptr<std::vector< uint8_t>> m_muon_MMHits;
  std::unique_ptr<std::vector< uint8_t>> m_muon_MMHoles;

  // segment
  std::unique_ptr<std::vector< size_t>> m_muon_nSegments;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_chiSquared;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_numberDoF;
  std::unique_ptr<std::vector< std::vector< int>>> m_muon_seg_sector;
  std::unique_ptr<std::vector< std::vector< int>>> m_muon_seg_chamberIndex;
  std::unique_ptr<std::vector< std::vector< int>>> m_muon_seg_nPrecisionHits;
  std::unique_ptr<std::vector< std::vector< int>>> m_muon_seg_nPhiLayers;
  std::unique_ptr<std::vector< std::vector< int>>> m_muon_seg_nTrigEtaLayers;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_x;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_y;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_z;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_px;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_py;
  std::unique_ptr<std::vector< std::vector< float>>> m_muon_seg_pz;

  std::unique_ptr<std::vector<float>> m_l1_eta;
  std::unique_ptr<std::vector<float>> m_l1_phi;
  std::unique_ptr<std::vector<float>> m_l1_dRoff;
  std::unique_ptr<std::vector<int>>   m_l1_BOM;
  std::unique_ptr<std::vector<int>>   m_l1_thrNum;
  std::unique_ptr<std::vector<int>>   m_l1_roiNum;

  std::unique_ptr<std::vector<float>> m_l2sa_e;
  std::unique_ptr<std::vector<float>> m_l2sa_pt;
  std::unique_ptr<std::vector<float>> m_l2sa_eta;
  std::unique_ptr<std::vector<float>> m_l2sa_phi;
  std::unique_ptr<std::vector<float>> m_l2sa_etaMS;
  std::unique_ptr<std::vector<float>> m_l2sa_phiMS;
  std::unique_ptr<std::vector<bool>> m_l2sa_pass;

  std::unique_ptr<std::vector< std::vector< float >>> m_l2sa_superPointR;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2sa_superPointZ;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2sa_superPointSlope;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2sa_superPointIntercept;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2sa_superPointChi2;

  // sTGC clusters
  std::unique_ptr<std::vector< std::vector< unsigned int >>> m_l2sa_stgcClusterLayer;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_stgcClusterIsOutlier;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_stgcClusterType;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_stgcClusterEta;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_stgcClusterPhi;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_stgcClusterR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_stgcClusterZ;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_stgcClusterResidualR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_stgcClusterResidualPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_stgcClusterStationEta;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_stgcClusterStationPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_stgcClusterStationName;

  // MM clusters
  std::unique_ptr<std::vector< std::vector< unsigned int >>> m_l2sa_mmClusterLayer;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_mmClusterIsOutlier;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_mmClusterEta;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_mmClusterPhi;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_mmClusterR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_mmClusterZ;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_mmClusterResidualR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2sa_mmClusterResidualPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_mmClusterStationEta;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_mmClusterStationPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2sa_mmClusterStationName;

  std::unique_ptr<std::vector<float>> m_l2sa_ptEndcapAlpha;
  std::unique_ptr<std::vector<float>> m_l2sa_ptEndcapBeta;

  std::unique_ptr<std::vector<float>> m_l2cb_e;
  std::unique_ptr<std::vector<float>> m_l2cb_pt;
  std::unique_ptr<std::vector<float>> m_l2cb_eta;
  std::unique_ptr<std::vector<float>> m_l2cb_phi;
  std::unique_ptr<std::vector<bool>> m_l2cb_pass;

  // for l2mt
  std::unique_ptr<std::vector<float>> m_l1_eta_mt;
  std::unique_ptr<std::vector<float>> m_l1_phi_mt;
  //std::unique_ptr<std::vector<std::vector<float>>> m_l1_dRoff;
  std::unique_ptr<std::vector<int>>   m_l1_BOM_mt;
  std::unique_ptr<std::vector<int>>   m_l1_thrNum_mt;
  std::unique_ptr<std::vector<int>>   m_l1_roiNum_mt;

  std::unique_ptr<std::vector<float>> m_l2mt_e;
  std::unique_ptr<std::vector<float>> m_l2mt_pt;
  std::unique_ptr<std::vector<float>> m_l2mt_eta;
  std::unique_ptr<std::vector<float>> m_l2mt_phi;
  std::unique_ptr<std::vector<float>> m_l2mt_etaMS;
  std::unique_ptr<std::vector<float>> m_l2mt_phiMS;
  std::unique_ptr<std::vector<bool>> m_l2mt_pass;

  std::unique_ptr<std::vector< std::vector< float >>> m_l2mt_superPointR;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2mt_superPointZ;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2mt_superPointSlope;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2mt_superPointIntercept;
  std::unique_ptr<std::vector< std::vector< float >>> m_l2mt_superPointChi2;

  // sTGC clusters
  std::unique_ptr<std::vector< std::vector< unsigned int >>> m_l2mt_stgcClusterLayer;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_stgcClusterIsOutlier;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_stgcClusterType;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_stgcClusterEta;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_stgcClusterPhi;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_stgcClusterR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_stgcClusterZ;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_stgcClusterResidualR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_stgcClusterResidualPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_stgcClusterStationEta;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_stgcClusterStationPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_stgcClusterStationName;

  // MM clusters
  std::unique_ptr<std::vector< std::vector< unsigned int >>> m_l2mt_mmClusterLayer;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_mmClusterIsOutlier;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_mmClusterEta;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_mmClusterPhi;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_mmClusterR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_mmClusterZ;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_mmClusterResidualR;
  std::unique_ptr<std::vector< std::vector< float >>>        m_l2mt_mmClusterResidualPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_mmClusterStationEta;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_mmClusterStationPhi;
  std::unique_ptr<std::vector< std::vector< int >>>          m_l2mt_mmClusterStationName;

  std::unique_ptr<std::vector<float>> m_road_Aw;
  std::unique_ptr<std::vector<float>> m_road_Bw;
  std::unique_ptr<std::vector<float>> m_road_R;
  std::unique_ptr<std::vector<float>> m_road_Z;
  std::unique_ptr<std::vector<float>> m_road_eta;
  std::unique_ptr<std::vector<float>> m_road_phi;

  std::unique_ptr<std::vector<float>> m_tgcHitEta;
  std::unique_ptr<std::vector<float>> m_tgcHitPhi;
  std::unique_ptr<std::vector<float>> m_tgcHitR;
  std::unique_ptr<std::vector<float>> m_tgcHitZ;
  std::unique_ptr<std::vector<int>> m_tgcHitStationNum;
  std::unique_ptr<std::vector<bool>> m_tgcHitIsStrip;
  std::unique_ptr<std::vector<int>> m_tgcHitBCTag;
  std::unique_ptr<std::vector<bool>> m_tgcHitInRoad;

  std::unique_ptr<std::vector<std::string>>                  m_chain;
  std::unique_ptr<std::vector< std::vector<int>>>            m_typeVec;
  std::unique_ptr<std::vector< std::vector<float>>>          m_ptVec;
  std::unique_ptr<std::vector< std::vector<float>>>          m_etaVec;
  std::unique_ptr<std::vector< std::vector<float>>>          m_phiVec;

  // parameters
  const double m_endcapPivotPlaneMinimumRadius{0.};      // minimum radius of pivot plane in endcap region
  const double m_endcapPivotPlaneMaximumRadius{11977.};  // maximum radius of pivot plane in endcap region

};

#endif
