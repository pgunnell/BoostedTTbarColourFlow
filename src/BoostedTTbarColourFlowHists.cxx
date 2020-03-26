#include "UHH2/BoostedTTbarColourFlow/include/BoostedTTbarColourFlowHists.h"
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

BoostedTTbarColourFlowHists::BoostedTTbarColourFlowHists(Context & ctx, const string & dirname): Hists(ctx, dirname) {
  // book all histograms here
  // jets
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);  
  book<TH1F>("N_PU", "N_{PU}", 100, 0, 100);  

  book<TH2D>("WsubjetsImage","Det. eta vs phi; eta; phi",400,-3.0,3.0,400,-3.1415,3.1415);   
  book<TH2D>("BjetOnesubjetImage","Det. eta vs phi; eta; phi",400,-3.0,3.0,400,-3.1415,3.1415);   

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);

  book<TH1F>("PullHisto", "ThetaPull", 50, 0., 1.);
  book<TH1F>("PuppiWeights", "PuppiWeight", 50, 0., 1.);

  book<TH1F>("PullParallel_LightRivet", "Parallel Pull", 50, 0., 0.5);
  book<TH1F>("PullPerpendicular_LightRivet", "Perpendicular Pull", 50, 0., 0.5);

  book<TH1F>("PullParallelRivet", "Parallel Pull", 10, 0., 0.02);
  book<TH1F>("PullPerpendicularRivet", "Perpendicular Pull", 10, 0., 0.02);

  book<TH1F>("PullModule", "PullModule", 50, 0., 0.5);

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

}

void BoostedTTbarColourFlowHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  double pi = 3.1415926535;

  std::vector<TopJet>* xconecorrectedjets = event.topjets;

  std::vector<PFParticle>* pfparticles = event.pfparticles;
  if (pfparticles == nullptr) { throw std::runtime_error("pfparticles is nullptr"); }

  int Njets = xconecorrectedjets->size();
  hist("N_jets")->Fill(Njets, weight);
  if(!event.isRealData)  hist("N_PU")->Fill(event.genInfo->pileup_TrueNumInteractions(), weight);

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);

  double deltaRMuonMin=0;

  auto & jet0 = xconecorrectedjets->at(0); 
  for(const auto & topjet: *xconecorrectedjets){
    if(deltaR(topjet,event.muons->at(0))>deltaRMuonMin) {
      jet0 = topjet;
      deltaRMuonMin=deltaR(topjet,event.muons->at(0));
    }
  }

  //const auto & smalljet0 = event.jets->at(0);
  //const auto & genjet0 = event.genjets->at(0);

  Jet bjet; bjet=jets_for_images_withWMass(event,jet0).at(0);
  Jet lightjet1; lightjet1=jets_for_images_withWMass(event,jet0).at(1);
  Jet lightjet2; lightjet2=jets_for_images_withWMass(event,jet0).at(2);

  //Kinematic cuts are already included in the jets_for_images function
  if(bjet.pt()<30 || fabs(bjet.eta())>2.4) return;
  if(lightjet1.pt()<30 || fabs(lightjet1.eta())>2.4) return;
  if(lightjet2.pt()<30 || fabs(lightjet2.eta())>2.4) return;

  bool genmatching = false;
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
  }

  if(!matchbjet || !matchlightjet1 || !matchlightjet2) return;

  vector<double> Pull;
  Pull.push_back(0.);
  Pull.push_back(0.);

  //try to build colour flow variables

  double rapidityleadsubjet=bjet.v4().Rapidity();
  double phileadsubjet=bjet.phi();

  double rapiditysubleadsubjet=lightjet1.v4().Rapidity();
  double phisubleadsubjet=lightjet1.phi();

  double rotationangle=atan((DeltaPhi(phisubleadsubjet,phileadsubjet))/(rapiditysubleadsubjet-rapidityleadsubjet));
  if((rapiditysubleadsubjet-rapidityleadsubjet)<0 && (DeltaPhi(phisubleadsubjet,phileadsubjet))>0) rotationangle=pi+rotationangle;
  if((rapiditysubleadsubjet-rapidityleadsubjet)<0 && (DeltaPhi(phisubleadsubjet,phileadsubjet))<0) rotationangle=pi+rotationangle;

  //we want to normalize to deltaEta between the jets
  double NormFactor = sqrt(pow(rapidityleadsubjet-rapiditysubleadsubjet,2)+pow(DeltaPhi(phileadsubjet,phisubleadsubjet),2));

  ///////////////////////////////////////////

  if(debug)
    cout<<" Large jet Area "<<jet0.jetArea()<<" "<<jet0.pt()<<endl;
  LorentzVector consistSum_big;

  for (const auto candIndex : jet0.pfcand_indexs()) {

    consistSum_big += pfparticles->at(candIndex).v4();
  }

  if(debug)
    cout<<consistSum_big.pt()<<" with no subjets "<<endl;

  for (auto subjet : jet0.subjets() ) {
    consistSum_big += subjet.v4();
  }

  if(debug)
    cout<<consistSum_big.pt()<<" with subjets "<<endl;

  /*if(debug)
    cout<<" AK4 jet Area "<<smalljet0.jetArea()<<" "<<smalljet0.pt()<<endl;
  LorentzVector consistSum_small;

  //for (const auto candIndex : smalljet0.pfcand_indexs()) { 
  for (const PFParticle & candInd : *event.pfparticles){
  
    if(deltaR(candInd,smalljet0)<0.4)
      consistSum_small += candInd.v4();
    //consistSum_small += pfparticles->at(candIndex).v4();
  }

  if(debug)
    cout<<consistSum_small.pt()<<endl;

  if(debug)
    cout<<" Gen jet Area Not Defined "<<genjet0.pt()<<endl;

  LorentzVector consistSum_gen;
  for (const auto candInd : genjet0.genparticles_indices()) {
    consistSum_gen += genparticles->at(candInd).v4();
  }

  if(debug)
  cout<<consistSum_gen.pt()<<endl;*/

  if(debug)
    cout<<" Subjet Area "<<bjet.jetArea()<<" "<<bjet.pt()<<" number of daughters "<<bjet.numberOfDaughters()<<endl;
  LorentzVector consistSum;

  double radiusbjet = CalculateRadiusJet(bjet);
    
  for (const PFParticle & candInd : *event.pfparticles){
    
    if(deltaR(candInd,bjet)>radiusbjet) continue;
    if(candInd.puppiWeight()<0.4) continue;
    //cout<<fabs(pfparticles->at(candIndex).eta()-bjet.eta())<<" delta Eta bjet "<<pfparticles->at(candIndex).pt()<<endl;
    //cout<<pfparticles->at(candIndex).pt()<<" "<<i<<endl;
    vector<double> histfill = FillRapidityPhiPt(candInd,bjet,lightjet1,weight,rapidityleadsubjet,phileadsubjet,rotationangle,NormFactor);
    consistSum += candInd.v4();

    ((TH2D*) hist("BjetOnesubjetImage"))->Fill(histfill.at(0),histfill.at(1),histfill.at(2));
    vector<double> OldPull = Pull;
    Pull = GetPull(candInd, OldPull, bjet.pt(), bjet.v4().Rapidity(), bjet.phi());

    OldPull.clear();
    histfill.clear();
  }

  if(debug)
    cout<<consistSum.pt()<<endl;

  LorentzVector consistSumCand;

  for (const PFParticle & candInd : *event.pfparticles){
    
    if(deltaR(candInd,bjet)<0.4)
      consistSumCand += candInd.v4();
  }

  if(debug)
    cout<<consistSumCand.pt()<<endl;

  double pull_angle = fabs(CalculatePullAngle(bjet, lightjet1, event));
  hist("PullHistoRivet")->Fill(pull_angle / pi, weight);

  TVector3 pull_module_threevector = CalculatePullRivet(bjet, event);
  double pull_module = sqrt(pull_module_threevector.X()*pull_module_threevector.X()+pull_module_threevector.Y()*pull_module_threevector.Y()+pull_module_threevector.Z()*pull_module_threevector.Z());
  hist("PullModuleRivet")->Fill(pull_module, weight);
  if(pull_module>=0.2) pull_module=0.0199999;
  hist("PullModuleRivetSmallBins")->Fill(pull_module, weight);

  hist("PullParallelRivet")->Fill(pull_module*cos(pull_angle),weight);
  hist("PullPerpendicularRivet")->Fill(pull_module*sin(pull_angle),weight);

  double pull_angle_light = fabs(CalculatePullAngle(lightjet1, lightjet2, event));
  hist("PullHistoRivet_Light")->Fill(pull_angle_light / pi, weight);

  TVector3 pull_module_threevector_light = CalculatePullRivet(lightjet1, event);
  double pull_module_light = sqrt(pull_module_threevector_light.X()*pull_module_threevector_light.X()+pull_module_threevector_light.Y()*pull_module_threevector_light.Y()+pull_module_threevector_light.Z()*pull_module_threevector_light.Z());
  hist("PullModuleRivet_Light")->Fill(pull_module_light, weight);
  if(pull_module_light>=0.2) pull_module_light=0.0199999;
  hist("PullModuleRivet_LightSmallBins")->Fill(pull_module_light, weight);

  hist("PullParallel_LightRivet")->Fill(pull_module_light*cos(pull_angle_light),weight);
  hist("PullPerpendicular_LightRivet")->Fill(pull_module_light*sin(pull_angle_light),weight);

  ////////////////////////////////////////////////////

  vector<double> Pull_light;
  Pull_light.push_back(0.);
  Pull_light.push_back(0.);

  double rapidityleadlightjet1=lightjet1.v4().Rapidity();
  double phileadlightjet1=lightjet1.phi();

  double rapidityleadlightjet2=lightjet2.v4().Rapidity();
  double phileadlightjet2=lightjet2.phi();

  double rotationangle_light=atan(DeltaPhi(phileadlightjet2,phileadlightjet1)/(rapidityleadlightjet2-rapidityleadlightjet1));
  if((rapidityleadlightjet2-rapidityleadlightjet1)<0 && DeltaPhi(phileadlightjet2,phileadlightjet1)>0) rotationangle_light=pi+rotationangle_light;
  if((rapidityleadlightjet2-rapidityleadlightjet1)<0 && DeltaPhi(phileadlightjet2,phileadlightjet1)<0) rotationangle_light=pi+rotationangle_light;

  //we want to normalize to deltaEta between the jets
  double NormFactor_light = sqrt(pow(rapidityleadlightjet1-rapidityleadlightjet2,2)+pow(DeltaPhi(phileadlightjet1,phileadlightjet2),2));

  double radiusljet = CalculateRadiusJet(lightjet1);
  //for (const auto candIndex : lightjet1.pfcand_indexs()) {
  for (const PFParticle & candInd : *event.pfparticles){

    if(deltaR(candInd,lightjet1)>radiusljet) continue;

    hist("PuppiWeights")->Fill(candInd.puppiWeight(),weight);
    //vector<double> histfill_light = FillRapidityPhiPt(pfparticles->at(candInd),lightjet1,lightjet2,weight,rapidityleadlightjet1,phileadlightjet1,rotationangle_light,NormFactor_light);
    vector<double> histfill_light = FillRapidityPhiPt(candInd,lightjet1,lightjet2,weight,rapidityleadlightjet1,phileadlightjet1,rotationangle_light,NormFactor_light);

    ((TH2D*) hist("WsubjetsImage"))->Fill(histfill_light.at(0),histfill_light.at(1),histfill_light.at(2));
    
    vector<double> OldPull_light=Pull_light;
    Pull_light = GetPull(candInd, OldPull_light, lightjet1.pt(), lightjet1.v4().Rapidity(), lightjet1.phi());

    OldPull_light.clear();
    histfill_light.clear();
  }

  //I am able to reproduce what is in the original paper (for the background)
  double AxisPhi = DeltaPhi(phileadsubjet,phisubleadsubjet);
  double AxisRapidity = rapidityleadsubjet-rapiditysubleadsubjet;

  double theta = atan((AxisPhi-Pull.at(1))/(AxisRapidity-Pull.at(0)));

  //this is to bring everything in the range 0-2Pi
  if((Pull.at(0)-AxisRapidity)<0 && (Pull.at(1)-AxisPhi)>0) theta=pi+theta;
  if((Pull.at(0)-AxisRapidity)<0 && (Pull.at(1)-AxisPhi)<0) theta=pi-theta;
  if((Pull.at(0)-AxisRapidity)>0 && (Pull.at(1)-AxisPhi)<0) theta=fabs(theta);
    
  //this is to bring everything in the range 0-1
  theta = (theta)/pi;

  hist("PullHisto")->Fill(theta,weight);
  double pullmodule=sqrt(pow(Pull.at(1),2)+pow(Pull.at(0),2));
  hist("PullModule")->Fill(pullmodule,weight);

  double AxisPhi_light = DeltaPhi(phileadlightjet1,phileadlightjet2);
  double AxisRapidity_light = rapidityleadlightjet1-rapidityleadlightjet2;

  double theta_light = atan((AxisPhi_light-Pull_light.at(1))/(AxisRapidity_light-Pull_light.at(0)));

  //this is to bring everything in the range 0-2Pi
  if((Pull_light.at(0)-AxisRapidity_light)<0 && (Pull_light.at(1)-AxisPhi_light)>0) theta_light=pi+theta_light;
  if((Pull_light.at(0)-AxisRapidity_light)<0 && (Pull_light.at(1)-AxisPhi_light)<0) theta_light=pi-theta_light;
  if((Pull_light.at(0)-AxisRapidity_light)>0 && (Pull_light.at(1)-AxisPhi_light)<0) theta_light=fabs(theta_light);

  //this is to bring everything in the range 0-2
  theta_light = theta_light/pi;

  hist("PullHisto_Light")->Fill(theta_light,weight);
  double pullmodule_light=sqrt(pow(Pull_light.at(1),2)+pow(Pull_light.at(0),2));
  hist("PullModule_Light")->Fill(pullmodule_light,weight);

  LorentzVector recoTop=bjet.v4()+lightjet1.v4()+lightjet2.v4();

  double massW=-10;
  double massWdiff = 100;

  double massW_1 = (lightjet1.v4()+lightjet2.v4()).mass();
  double massW_2 = (lightjet1.v4()+bjet.v4()).mass();
  double massW_3 = (lightjet2.v4()+bjet.v4()).mass();

  if(fabs(massW_1-80.3) < massWdiff) { massW = massW_1; massWdiff = fabs(massW_1-80.3); }
  if(fabs(massW_2-80.3) < massWdiff) { massW = massW_2; massWdiff = fabs(massW_2-80.3); }
  if(fabs(massW_3-80.3) < massWdiff) { massW = massW_3; massWdiff = fabs(massW_3-80.3); }

  LorentzVector recoW=lightjet1.v4()+lightjet2.v4();

  hist("RecoTopMass")->Fill(recoTop.mass(),weight);
  hist("RecoWMass_best")->Fill(massW,weight);
  hist("RecoWMass")->Fill(recoW.mass(),weight);
  hist("RecoTopMassXConeJet")->Fill(jet0.v4().mass(),weight);
  hist("NumberOfSubJets")->Fill(jet0.subjets().size(),weight);
    
}

BoostedTTbarColourFlowHists::~BoostedTTbarColourFlowHists(){}
