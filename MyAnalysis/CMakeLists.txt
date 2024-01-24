# The name of the package:
atlas_subdir (MyAnalysis)

# Add the shared library:
atlas_add_library (MyAnalysisLib
  MyAnalysis/*.h Root/*.cxx
  PUBLIC_HEADERS MyAnalysis
  LINK_LIBRARIES AnaAlgorithmLib TrigConfInterfaces TrigDecisionToolLib TrkExInterfaces MuonAnalysisInterfacesLib TruthClassificationLib xAODMuon xAODTrigMuon xAODEventInfo)

atlas_add_component (MyAnalysis
  src/components/*.cxx
  LINK_LIBRARIES MyAnalysisLib)

# Install files from the package:
atlas_install_python_modules( python/*.py )
