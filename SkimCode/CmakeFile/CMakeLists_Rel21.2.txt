## automatically generated CMakeLists.txt file

# Declare the package
atlas_subdir( HFtrigValidation )

# Declare external dependencies ... default here is to include ROOT
find_package( ROOT COMPONENTS MathCore RIO Core Tree Hist )

# Declare public and private dependencies
# Public dependencies are anything that appears in the headers in public include dir
# Private is anything else

# An example is included
atlas_depends_on_subdirs(
    PUBLIC
    DataQuality/GoodRunsLists
    Event/xAOD/xAODMuon
    Event/xAOD/xAODTracking
    Event/xAOD/xAODTruth
    External/AtlasROOT
    InnerDetector/InDetRecTools/InDetTrackSelectionTool
    PhysicsAnalysis/MuonID/MuonSelectorTools
    PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonMomentumCorrections
    Trigger/TrigAnalysis/TrigDecisionTool
    Trigger/TrigAnalysis/TriggerMatchingTool
    ForwardDetectors/ZDC/ZdcAnalysis
    PhysicsAnalysis/HeavyIonPhys/HIEventUtils

    PRIVATE
    Control/AthAnalysisBaseComps
    Event/xAOD/xAODEventInfo
    Event/xAOD/xAODTrigger
    Event/xAOD/xAODForward
    Event/xAOD/xAODMissingET
    Event/xAOD/xAODHIEvent
    Event/xAOD/ZdcModuleContainer
)

# Declare package as a library
# Note the convention that library names get "Lib" suffix
# Any package you add to dependencies above, you should add
# to LINK_LIBRARIES line below (see the example)
atlas_add_library( HFtrigValidationLib src/*.cxx
                   PUBLIC_HEADERS HFtrigValidation
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES ${ROOT_LIBRARIES}
                   AthAnalysisBaseCompsLib 
                   GoodRunsLists 
                   xAODMuon
                   xAODTracking 
                   xAODTruth 
                   InDetTrackSelectionTool 
                   MuonSelectorTools
                   MuonMomentumCorrections
                   TrigDecisionTool 
                   TriggerMatchingTool
                   xAODEventInfo 
                   xAODTrigger
                   xAODForward
                   xAODMissingET 
                   xAODHIEvent 
                   ZdcAnalysis
                   HIEventUtils
                   GeneratorObjects
)

# if you add components (tools, algorithms) to this package
# these lines are needed so you can configure them in joboptions
atlas_add_component( HFtrigValidation src/components/*.cxx
                      LINK_LIBRARIES HFtrigValidationLib 
)

# if you add an application (exe) to this package
# declare it like this (note convention that apps go in the util dir)
# atlas_add_executable( MyApp util/myApp.cxx
#                       LINK_LIBRARIES HFtrigValidationLib
# )

# Install python modules, joboptions, and share content
atlas_install_python_modules( python/*.py )
atlas_install_joboptions( share/*.py )
atlas_install_data( data/* )
# You can access your data from code using path resolver, e.g.
# PathResolverFindCalibFile("HFtrigValidation/file.txt")



