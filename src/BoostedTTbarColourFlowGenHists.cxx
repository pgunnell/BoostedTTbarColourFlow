#include "UHH2/BoostedTTbarColourFlow/include/BoostedTTbarColourFlowGenHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/BoostedTTbarColourFlow/include/imagefunctions.h"
#include "UHH2/BoostedTTbarColourFlow/include/AnalysisFunctions.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

//TRY TO MATCH THE QUARKS TO THE JETS AND SEE WHAT HAPPENS

BoostedTTbarColourFlowGenHists::BoostedTTbarColourFlowGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname) {
  // book all histograms here
  // jets
  book<TH1F>("N_genjets", "N_{gen-jets}", 20, 0, 20);  

  book<TH2D>("WsubjetsImage","Det. eta vs phi; eta; phi",400,-3.0,3.0,400,-3.1415,3.1415);   
  book<TH2D>("BjetOnesubjetImage","Det. eta vs phi; eta; phi",400,-3.0,3.0,400,-3.1415,3.1415);   

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  // primary vertices
  book<TH1F>("PullParallel_LightRivet", "Parallel Pull", 50, 0., 0.5);
  book<TH1F>("PullPerpendicular_LightRivet", "Perpendicular Pull", 50, 0., 0.5);
  book<TH1F>("PullParallelRivet", "Parallel Pull", 10, 0., 0.02);
  book<TH1F>("PullPerpendicularRivet", "Perpendicular Pull", 10, 0., 0.02);

  book<TH1F>("PullModuleRivet", "PullModuleRivet", 50, 0., 0.5);
  book<TH1F>("PullHistoRivet", "ThetaPull", 10, 0., 1.);
  book<TH1F>("PullModuleRivetSmallBins", "PullModuleRivet", 10, 0., 0.02);

  book<TH1F>("PullModuleRivet_Light", "PullModuleRivet", 50, 0., 0.5);
  book<TH1F>("PullHistoRivet_Light", "ThetaPull", 10, 0., 1.);
  book<TH1F>("PullModuleRivet_LightSmallBins", "PullModuleRivet", 10, 0., 0.02);

  book<TH1F>("PullHisto_Light", "ThetaPull", 50, 0., 1.);
  book<TH1F>("PullModule_Light", "PullModule", 50, 0., 0.5);

  book<TH1F>("RecoTopMass", "Reconstructed three-jet mass", 50, 0, 400);
  book<TH1F>("RecoWMass", "Reconstructed two-jet mass", 40, 0, 200);
  book<TH1F>("RecoWMass_best", "Reconstructed two-jet mass", 40, 0, 200);

  book<TH1F>("RecoTopMassXConeJet", "Reconstructed jet mass", 50, 0, 400);

  book<TH1F>("NumberOfSubJets", "Number of XCone subjets", 6, -0.5, 5.5);

  h_genjets_had = ctx.declare_event_input<std::vector<GenTopJet> > ("GEN_XCone33_had_Combined");
  h_gensel_2 = ctx.declare_event_input<bool> ("passed_gensel_2");

}

void BoostedTTbarColourFlowGenHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronGenHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  double pi = 3.1415926535;

  std::vector<GenTopJet>* genxconelepjets = event.gentopjets;

  std::vector<GenParticle>* genparticles = event.genparticles;
  
  std::vector<GenTopJet> genxconecorrectedjets = event.get(h_genjets_had);

  bool passed_gensel = event.get(h_gensel_2);

  int Ngenjets = genxconecorrectedjets.size();
  int Ngenlepjets = genxconelepjets->size();
  hist("N_genjets")->Fill(Ngenjets, weight);
  
  int Ngenmuons = 0;

  vector <GenParticle> genmuon;

  for(const auto & genparticle: *genparticles){
    if(genparticle.status()==1 && fabs(genparticle.pdgId())==13){
      Ngenmuons++;
      genmuon.push_back(genparticle);
    }
  }
  
  if(Ngenmuons<1) return;
  if(Ngenjets<1) return;
  if(Ngenlepjets<1) return;

  auto & jet0 = genxconecorrectedjets.at(0); 
  
  if(!passed_gensel) return;
  if(genxconecorrectedjets.at(0).pt() < 400 || fabs(genxconecorrectedjets.at(0).eta()) > 2.4) return;
  if(genmuon.at(0).pt() < 60 || fabs(genmuon.at(0).eta()) > 2.4) return;
  if(genxconelepjets->at(0).pt()<10) return;
  if((genxconelepjets->at(0).v4()+genmuon.at(0).v4()).mass() > jet0.v4().mass()) return;
  
  GenJet bjet;  GenJet lightjet1;  GenJet lightjet2; bool matchbjetgen=false; bool matchlightjetgen1=false; bool matchlightjetgen2=false;

  if(jet0.subjets().size()<3) return;  

  for(const auto &subjet: jet0.subjets()){
    for(const auto & genparticle: *event.genparticles){
      if(fabs(genparticle.pdgId())==5){
	//cout<<genparticle.pdgId()<<" "<<genparticle.status()<<" "<<genparticle.phi()<<" "<<genparticle.eta()<<endl;
	//cout<<deltaR(genparticle,bjet)<<endl;
	if(deltaR(genparticle,subjet)<0.3) {
	  bjet=subjet;
	  matchbjetgen=true;
	}
      }
    }
  } 

  //cout<<"New Event"<<endl;

  if(matchbjetgen){
    for(const auto &subjet: jet0.subjets()){
      if(subjet==bjet) continue;
      if(lightjet1.pt()<0.1) { lightjet1=subjet; matchlightjetgen1=true; continue;}
      if(lightjet2.pt()<0.1) { lightjet2=subjet; matchlightjetgen2=true; continue;}
    }
  }

  if(!matchlightjetgen2 || !matchlightjetgen1) return;

  if(bjet.pt()<30 || fabs(bjet.eta())>2.4) return;
  if(lightjet1.pt()<30 || fabs(lightjet1.eta())>2.4) return;
  if(lightjet2.pt()<30 || fabs(lightjet2.eta())>2.4) return;

  /*bool genmatching = false;
  bool matchbjet = true;
  bool matchlightjet1 = true;
  bool matchlightjet2 = true;

  bool debug=false;

  if(debug)
    cout<<"New event"<<endl;
 
  if(genmatching) {

    for(const auto & genparticle: *event.genparticles){

      if(genparticle.status()==23 && fabs(genparticle.pdgId())==5){
	if(deltaR(genparticle,bjet)<0.4) matchbjet=false;
      }
      if(genparticle.status()==23 && fabs(genparticle.pdgId())!=5 && fabs(genparticle.pdgId())!=6){
	if(deltaR(genparticle,lightjet1)<0.4) matchlightjet1=false;
	if(deltaR(genparticle,lightjet2)<0.4) matchlightjet2=false;
      }
      
    }
    }*/

  double pull_angle = fabs(CalculatePullGenAngle(bjet, lightjet1, event));
  hist("PullHistoRivet")->Fill(pull_angle / pi, weight);

  TVector3 pull_module_threevector = CalculatePullGenRivet(bjet, event);
  double pull_module = sqrt(pull_module_threevector.X()*pull_module_threevector.X()+pull_module_threevector.Y()*pull_module_threevector.Y()+pull_module_threevector.Z()*pull_module_threevector.Z());
  hist("PullModuleRivet")->Fill(pull_module, weight);
  if(pull_module>=0.2) pull_module=0.0199999;
  hist("PullModuleRivetSmallBins")->Fill(pull_module, weight);

  hist("PullParallelRivet")->Fill(pull_module*cos(pull_angle),weight);
  hist("PullPerpendicularRivet")->Fill(pull_module*sin(pull_angle),weight);

  double pull_angle_light = fabs(CalculatePullGenAngle(lightjet1, lightjet2, event));
  hist("PullHistoRivet_Light")->Fill(pull_angle_light / pi, weight);

  TVector3 pull_module_threevector_light = CalculatePullGenRivet(lightjet1, event);
  double pull_module_light = sqrt(pull_module_threevector_light.X()*pull_module_threevector_light.X()+pull_module_threevector_light.Y()*pull_module_threevector_light.Y()+pull_module_threevector_light.Z()*pull_module_threevector_light.Z());
  hist("PullModuleRivet_Light")->Fill(pull_module_light, weight);
  if(pull_module_light>=0.2) pull_module_light=0.0199999;
  hist("PullModuleRivet_LightSmallBins")->Fill(pull_module_light, weight);

  hist("PullParallel_LightRivet")->Fill(pull_module_light*cos(pull_angle_light),weight);
  hist("PullPerpendicular_LightRivet")->Fill(pull_module_light*sin(pull_angle_light),weight);

  ////////////////////////////////////////////////////

  double rapidityleadlightjet1=lightjet1.v4().Rapidity();
  double phileadlightjet1=lightjet1.phi();

  double rapidityleadlightjet2=lightjet2.v4().Rapidity();
  double phileadlightjet2=lightjet2.phi();

  double rotationangle_light=atan(DeltaPhiGen(phileadlightjet2,phileadlightjet1)/(rapidityleadlightjet2-rapidityleadlightjet1));
  if((rapidityleadlightjet2-rapidityleadlightjet1)<0 && DeltaPhiGen(phileadlightjet2,phileadlightjet1)>0) rotationangle_light=pi+rotationangle_light;
  if((rapidityleadlightjet2-rapidityleadlightjet1)<0 && DeltaPhiGen(phileadlightjet2,phileadlightjet1)<0) rotationangle_light=pi+rotationangle_light;

  double NormFactor_light = sqrt(pow(rapidityleadlightjet1-rapidityleadlightjet2,2)+pow(DeltaPhiGen(phileadlightjet1,phileadlightjet2),2));
  
  for (const auto candIndex : lightjet1.genparticles_indices()) {
 
    vector<double> histfill_light = FillRapidityPhiPt_Gen(event.genparticles->at(candIndex),lightjet1,lightjet2,weight,rapidityleadlightjet1,phileadlightjet1,rotationangle_light,NormFactor_light);
    ((TH2D*) hist("WsubjetsImage"))->Fill(histfill_light.at(0),histfill_light.at(1),histfill_light.at(2));
    histfill_light.clear();
  }

  double rapidityleadsubjet=bjet.v4().Rapidity();
  double phileadsubjet=bjet.phi();

  double rapiditysubleadsubjet=lightjet1.v4().Rapidity();
  double phisubleadsubjet=lightjet1.phi();

  double rotationangle=atan((DeltaPhiGen(phisubleadsubjet,phileadsubjet))/(rapiditysubleadsubjet-rapidityleadsubjet));
  if((rapiditysubleadsubjet-rapidityleadsubjet)<0 && (DeltaPhiGen(phisubleadsubjet,phileadsubjet))>0) rotationangle=pi+rotationangle;
  if((rapiditysubleadsubjet-rapidityleadsubjet)<0 && (DeltaPhiGen(phisubleadsubjet,phileadsubjet))<0) rotationangle=pi+rotationangle;

  //we want to normalize to deltaEta between the jets
  double NormFactor = sqrt(pow(rapidityleadsubjet-rapiditysubleadsubjet,2)+pow(DeltaPhiGen(phileadsubjet,phisubleadsubjet),2));

  for (const auto candIndex : bjet.genparticles_indices()) {
    
    vector<double> histfill = FillRapidityPhiPt_Gen(event.genparticles->at(candIndex),bjet,lightjet1,weight,rapidityleadsubjet,phileadsubjet,rotationangle,NormFactor);
    ((TH2D*) hist("BjetOnesubjetImage"))->Fill(histfill.at(0),histfill.at(1),histfill.at(2));
    histfill.clear();
  }

  //I am able to reproduce what is in the original paper (for the background)

  LorentzVector recoTop=bjet.v4()+lightjet1.v4()+lightjet2.v4();
  LorentzVector recoW=lightjet1.v4()+lightjet2.v4();

  double massW=-10;
  double massWdiff = 100;

  double massW_1 = (lightjet1.v4()+lightjet2.v4()).mass();
  double massW_2 = (lightjet1.v4()+bjet.v4()).mass();
  double massW_3 = (lightjet2.v4()+bjet.v4()).mass();

  if(fabs(massW_1-80.3) < massWdiff) { massW = massW_1; massWdiff = massW_1=80.3; }
  if(fabs(massW_2-80.3) < massWdiff) { massW = massW_2; massWdiff =massW_1=80.3; }
  if(fabs(massW_3-80.3) < massWdiff) { massW = massW_3; massWdiff =massW_1=80.3; }

  hist("RecoTopMass")->Fill(recoTop.mass(),weight);
  hist("RecoWMass")->Fill(recoW.mass(),weight);
  hist("RecoWMass_best")->Fill(massW,weight);
  hist("RecoTopMassXConeJet")->Fill(jet0.v4().mass(),weight);
    
}

BoostedTTbarColourFlowGenHists::~BoostedTTbarColourFlowGenHists(){}
