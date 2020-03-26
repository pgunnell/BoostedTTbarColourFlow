#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/BoostedTTbarColourFlow/include/BoostedTTbarColourFlowSelections.h"
#include "UHH2/BoostedTTbarColourFlow/include/BoostedTTbarColourFlowHists.h"
#include "UHH2/BoostedTTbarColourFlow/include/BoostedTTbarColourFlowGenHists.h"
#include "UHH2/BoostedTTbarColourFlow/include/imagefunctions.h"
#include "UHH2/BoostedTTbarColourFlow/include/RecoSelections.h"
#include <UHH2/common/include/MuonIds.h>

using namespace std;
using namespace uhh2;

namespace uhh2examples {

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class BoostedTTbarColourFlowModule: public AnalysisModule {
public:
    
    explicit BoostedTTbarColourFlowModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
  std::unique_ptr<CommonModules> common;
  
  std::unique_ptr<JetCleaner> jetcleaner;
  
  // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.
  std::unique_ptr<Selection> njet_sel;
  std::unique_ptr<uhh2::Selection> njet_had;
  std::unique_ptr<uhh2::Selection> njet_lep;
  std::unique_ptr<uhh2::Selection> pt_sel;
  std::unique_ptr<uhh2::Selection> pt2_sel;
  std::unique_ptr<uhh2::Selection> eta_sel;
  std::unique_ptr<uhh2::Selection> mass_sel;
  std::unique_ptr<uhh2::Selection> subjet_quality;
  std::unique_ptr<uhh2::Selection> lepton_sel;
  
  // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_nocuts, h_njet, h_ele, h_ngenjet;

  Event::Handle<bool> h_recsel_2;
  Event::Handle<std::vector<TopJet>>h_jets;

  uhh2::Event::Handle<vector<double>> h_ptparticlesbjet;
  uhh2::Event::Handle<vector<double>> h_etaparticlesbjet;
  uhh2::Event::Handle<vector<double>> h_phiparticlesbjet;

  uhh2::Event::Handle<double> h_bjet_phi;
  uhh2::Event::Handle<double> h_bjet_eta;
  uhh2::Event::Handle<double> h_bjet_pt;

  uhh2::Event::Handle<double> h_lightjet_phi;
  uhh2::Event::Handle<double> h_lightjet_eta;
  uhh2::Event::Handle<double> h_lightjet_pt;

  uhh2::Event::Handle<vector<double>> h_ptparticleslight;
  uhh2::Event::Handle<vector<double>> h_etaparticleslight;
  uhh2::Event::Handle<vector<double>> h_phiparticleslight;
  
};


BoostedTTbarColourFlowModule::BoostedTTbarColourFlowModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from BoostedTTbarColourFlowModule!" << endl;
    
    // If needed, access the configuration of the module here, e.g.:
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    // 1. setup other modules. CommonModules and the JetCleaner:
    common.reset(new CommonModules());
    // TODO: configure common here, e.g. by 
    // calling common->set_*_id or common->disable_*
    common->switch_jetlepcleaner();
    common->switch_jetPtSorter();

    JetId btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
    common->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_PUPPI),PtEtaCut(200.0, 2.4)));  

    //common->init(ctx);
    
    // note that the JetCleaner is only kept for the sake of example;
    // instead of constructing a jetcleaner explicitly,
    // the cleaning can also be achieved with less code via CommonModules with:
    // common->set_jet_id(PtEtaCut(30.0, 2.4));
    // before the 'common->init(ctx)' line.
    
    // 2. set up selections
    njet_sel.reset(new NJetSelection(2)); // see common/include/NSelections.h

    // 3. Set up Hists classes:
    h_nocuts.reset(new BoostedTTbarColourFlowHists(ctx, "NoCuts"));
    h_njet.reset(new BoostedTTbarColourFlowHists(ctx, "PullDetectorLevel"));
    h_ngenjet.reset(new BoostedTTbarColourFlowGenHists(ctx, "PullGeneratorLevel"));
    h_ele.reset(new ElectronHists(ctx, "ele_nocuts"));

    h_recsel_2 = ctx.get_handle<bool>("passed_recsel_2");

    h_jets = ctx.get_handle<std::vector<TopJet> > ("XCone33_lep_Combined_Corrected");

    MuonId muid = AndId<Muon>(MuonID(Muon::Tight), PtEtaCut(60., 2.4));
    lepton_sel.reset(new NMuonSelection(1, -1, muid));
    
    h_ptparticlesbjet = ctx.declare_event_output<vector<double> >("h_ptparticlesbjet");
    h_etaparticlesbjet = ctx.declare_event_output<vector<double> >("h_etaparticlesbjet");
    h_phiparticlesbjet = ctx.declare_event_output<vector<double> >("h_phiparticlesbjet");

    h_ptparticleslight = ctx.declare_event_output<vector<double> >("h_ptparticleslight");
    h_etaparticleslight = ctx.declare_event_output<vector<double> >("h_etaparticleslight");
    h_phiparticleslight = ctx.declare_event_output<vector<double> >("h_phiparticleslight");

    h_bjet_eta = ctx.declare_event_output<double> ("h_bjet_eta");
    h_bjet_pt = ctx.declare_event_output<double> ("h_bjet_pt");
    h_bjet_phi = ctx.declare_event_output<double> ("h_bjet_phi");

    h_lightjet_eta = ctx.declare_event_output<double> ("h_lightjet_eta");
    h_lightjet_pt = ctx.declare_event_output<double> ("h_lightjet_pt");
    h_lightjet_phi = ctx.declare_event_output<double> ("h_lightjet_phi");
}


bool BoostedTTbarColourFlowModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.

    // 1. run all modules other modules.
  //common->process(event); //I don't apply common modules because I use Ntuples where common modules have been already applied
  sort_by_pt<PFParticle>(*event.pfparticles);
  
  // 2. test selections and fill histograms
  h_ele->fill(event);
  
  bool njet_selection = false;

  h_ngenjet->fill(event);

  assert(event.topjets); 
  if(event.topjets->size() < 1) return false;

  bool passed_recsel = false;
  passed_recsel = event.get(h_recsel_2);

  std::vector<TopJet> xconecorrectedjets_lep = event.get(h_jets);

  if(xconecorrectedjets_lep.size() < 1) return false;

  // the leptonic jet does not include the muon, so we need to re-include it
  int Nmuons = event.muons->size();
  if(Nmuons<1) return false;
  auto leading_muon = event.muons->at(0);

  if((xconecorrectedjets_lep.at(0).v4()+leading_muon.v4()).mass() > event.topjets->at(0).v4().mass()) return false;
  if( xconecorrectedjets_lep.at(0).pt() < 10) return false;

  //if(passed_recsel)
  if(passed_recsel && event.topjets->at(0).pt()>400 && fabs(event.topjets->at(0).eta())<2.5 && lepton_sel->passes(event)){// && pt_sel->passes(event) && pt2_sel->passes(event) && eta_sel->passes(event) && mass_sel->passes(event) && subjet_quality->passes(event) && lepton_sel->passes(event)) {
    
    double bjeteta = -10;
    double bjetphi = -10;
    double bjetpt = -10;

    double lightjeteta = -10;
    double lightjetphi = -10;
    double lightjetpt = -10;

    vector<double> ptparticlesbjet;
    vector<double> etaparticlesbjet;
    vector<double> phiparticlesbjet;
    
    vector<double> ptparticleslight;
    vector<double> etaparticleslight;
    vector<double> phiparticleslight;

    Jet bjet=jets_for_images_withWMass(event,event.topjets->at(0)).at(0);
    Jet lightjet1=jets_for_images_withWMass(event,event.topjets->at(0)).at(1);
    Jet lightjet2=jets_for_images_withWMass(event,event.topjets->at(0)).at(2);

    bjeteta = bjet.eta();
    bjetphi = bjet.phi();
    bjetpt = bjet.pt();

    lightjeteta = lightjet1.eta();
    lightjetphi = lightjet1.phi();
    lightjetpt = lightjet1.pt();
 
    event.set(h_lightjet_pt, lightjetpt);
    event.set(h_lightjet_eta, lightjeteta);
    event.set(h_lightjet_phi, lightjetphi);

    event.set(h_bjet_pt, bjetpt);
    event.set(h_bjet_eta, bjeteta);
    event.set(h_bjet_phi, bjetphi);
    
    vector<double> infoparticlesbjet_fill = info_particles(event,bjet);
    
    vector<double> infoparticleslight_fill = info_particles(event,lightjet1);
    
    if(infoparticlesbjet_fill.size()==0){
      std::cout<<"No particles selected,exiting"<<std::endl;
      return false;
    }
    
    for(unsigned int i=0; i<infoparticlesbjet_fill.size();i=i+3){
      ptparticlesbjet.push_back(infoparticlesbjet_fill.at(i));
      etaparticlesbjet.push_back(infoparticlesbjet_fill.at(i+1));
      phiparticlesbjet.push_back(infoparticlesbjet_fill.at(i+2));
    }
    
    event.set(h_ptparticlesbjet, ptparticlesbjet);
    event.set(h_etaparticlesbjet, etaparticlesbjet);
    event.set(h_phiparticlesbjet, phiparticlesbjet);
    
    for(unsigned int i=0; i<infoparticleslight_fill.size();i=i+3){
      ptparticleslight.push_back(infoparticleslight_fill.at(i));
      etaparticleslight.push_back(infoparticleslight_fill.at(i+1));
      phiparticleslight.push_back(infoparticleslight_fill.at(i+2));
    }
    
    event.set(h_ptparticleslight, ptparticleslight);
    event.set(h_etaparticleslight, etaparticleslight);
    event.set(h_phiparticleslight, phiparticleslight);
    
    infoparticlesbjet_fill.clear();
    infoparticleslight_fill.clear();
    ptparticlesbjet.clear();
    etaparticlesbjet.clear();
    phiparticlesbjet.clear();
    ptparticleslight.clear();
    etaparticleslight.clear();
    phiparticleslight.clear();

    h_njet->fill(event);
    njet_selection = true;
    
  }
  // 3. decide whether or not to keep the current event in the output:
  return njet_selection;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the BoostedTTbarColourFlowModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(BoostedTTbarColourFlowModule)

}
