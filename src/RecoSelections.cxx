#include <UHH2/MTopJet/include/RecoSelections.h>

uhh2::ElectronEtaVeto::ElectronEtaVeto(double lower_, double upper_):
lower(lower_),
upper(upper_) {}

bool uhh2::ElectronEtaVeto::passes(const uhh2::Event& event){
  bool pass_veto = false;
  if(event.electrons->size() != 0){
    double eta = fabs(event.electrons->at(0).eta());
    if(eta < lower || eta > upper) pass_veto = true;
  }
  return pass_veto;
}
////////////////////////////////////////////////////////

uhh2::NJetXCone::NJetXCone(uhh2::Context& ctx, const std::string & name, unsigned int njet):
h_jets(ctx.get_handle<std::vector<TopJet>>(name)),
njet_(njet) {}

bool uhh2::NJetXCone::passes(const uhh2::Event& event){
  bool pass_njet = false;
  std::vector<TopJet> jets = event.get(h_jets);
  if(jets.size() >= njet_) pass_njet = true;
  return pass_njet;
}
////////////////////////////////////////////////////////

uhh2::LeadingRecoJetETA::LeadingRecoJetETA(uhh2::Context& ctx, const std::string & name, float etacut):
h_jets(ctx.get_handle<std::vector<TopJet>>(name)),
etacut_(etacut) {}

bool uhh2::LeadingRecoJetETA::passes(const uhh2::Event& event){
  bool pass_jeteta = false;
  std::vector<TopJet> jets = event.get(h_jets);
  TopJet jet1;
  if(jets.size()>0){
    jet1 = jets.at(0);
    float eta = sqrt(jet1.eta()*jet1.eta());
    if(eta < etacut_) pass_jeteta = true;
  }
  return pass_jeteta;
}
////////////////////////////////////////////////////////

uhh2::LeadingRecoJetPT::LeadingRecoJetPT(uhh2::Context& ctx, const std::string & name, float ptcut):
h_jets(ctx.get_handle<std::vector<TopJet>>(name)),
ptcut_(ptcut) {}

bool uhh2::LeadingRecoJetPT::passes(const uhh2::Event& event){
  bool pass_jetpt = false;
  std::vector<TopJet> jets = event.get(h_jets);
  TopJet jet1;
  if(jets.size()>0){
    jet1 = jets.at(0);
    float pt = jet1.pt();
    if(pt > ptcut_) pass_jetpt = true;
  }
  return pass_jetpt;
}
////////////////////////////////////////////////////////

uhh2::SubjetQuality::SubjetQuality(uhh2::Context& ctx, const std::string & name, float ptmin_, float etamax_):
h_jets(ctx.get_handle<std::vector<TopJet>>(name)),
ptmin(ptmin_),
etamax(etamax_){}

bool uhh2::SubjetQuality::passes(const uhh2::Event& event){
  bool pass = true;
  std::vector<TopJet> jets = event.get(h_jets);
  if(jets.size() == 0) return false;
  std::vector<Jet> subjets = jets[0].subjets();
  if(subjets.size() != 3) pass = false;
  for(auto subjet: subjets){
    if(subjet.pt() < ptmin) pass = false;
    if(fabs(subjet.eta()) > etamax) pass = false;
  }
  return pass;
}
////////////////////////////////////////////////////////

uhh2::MassCutReco::MassCutReco(uhh2::Context& ctx, const std::string & name):
h_jets(ctx.get_handle<std::vector<Jet>>(name)) {}

bool uhh2::MassCutReco::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }
  TLorentzVector jet1_v4, jet2_v4, lepton1_v4, jet2_lep_v4;
  bool pass_masscut = false;

  if(jets.size()>1){

    jet1_v4.SetPtEtaPhiE(jets.at(0).pt(),jets.at(0).eta(),jets.at(0).phi(),jets.at(0).energy()); //v4 of first jet
    jet2_v4.SetPtEtaPhiE(jets.at(1).pt(),jets.at(1).eta(),jets.at(1).phi(),jets.at(1).energy()); //v4 of first jet
    lepton1_v4.SetPtEtaPhiE(lepton.pt(),lepton.eta(),lepton.phi(),lepton.energy()); //v4 of lepton
    jet2_lep_v4 = jet2_v4 + lepton1_v4; // v4 of (lepton+2nd jet)

    if(jet1_v4.M() > jet2_lep_v4.M()) pass_masscut = true;
  }
  return pass_masscut;
}
////////////////////////////////////////////////////////

uhh2::MassCutXCone::MassCutXCone(uhh2::Context& ctx, const std::string & hadname, const std::string & lepname):
h_hadjets(ctx.get_handle<std::vector<TopJet>>(hadname)),
h_lepjets(ctx.get_handle<std::vector<TopJet>>(lepname)){}

bool uhh2::MassCutXCone::passes(const uhh2::Event& event){
  std::vector<TopJet> hadjets = event.get(h_hadjets);
  std::vector<TopJet> lepjets = event.get(h_lepjets);

  Particle lepton;
  bool found_lep = false;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
    found_lep = true;
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
    found_lep = true;
  }

  bool pass_masscut = false;
  LorentzVector jet1_v4, jet2_v4;
  if(hadjets.size()>0 && lepjets.size()>0){
    jet1_v4 = hadjets.at(0).v4();
    if(found_lep) jet2_v4 = lepjets.at(0).v4() + lepton.v4();
    else          jet2_v4 = lepjets.at(0).v4();
    if(jet1_v4.M() > jet2_v4.M()) pass_masscut = true;
  }

  return pass_masscut;
}

////////////////////////////////////////////////////////

uhh2::DeltaRCutReco::DeltaRCutReco(uhh2::Context& ctx, const std::string & name, float jetradius):
h_jets(ctx.get_handle<std::vector<TopJet>>(name)),
jetradius_(jetradius) {}

bool uhh2::DeltaRCutReco::passes(const uhh2::Event& event){
  std::vector<TopJet> jets = event.get(h_jets);
  TopJet jet2;
  bool pass_deltaR = false;
  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }
  if(jets.size()>1){
    jet2 = jets.at(1);
    if(deltaR(lepton, jet2) < jetradius_) pass_deltaR = true;
  }
  return pass_deltaR;
}
////////////////////////////////////////////////////////

uhh2::NRecoJets::NRecoJets(uhh2::Context& ctx, const std::string & name, float min_pt, float min, float max):
h_jets(ctx.get_handle<std::vector<Jet>>(name)),
min_pt_(min_pt),
min_(min),
max_(max) {}

bool uhh2::NRecoJets::passes(const uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  bool pass_NGen = false;
  int counter = 0;
  for(unsigned int i = 0; i<jets.size(); ++i){
    if(jets.at(i).pt() >= min_pt_){
      ++counter;
    }
  }
  if(counter >= min_ && counter <= max_) pass_NGen = true;
  return pass_NGen;
}

////////////////////////////////////////////////////////

RecoJetLeptonCleaner::RecoJetLeptonCleaner(uhh2::Context& ctx, const std::string & name, float jetradius):
h_jets(ctx.get_handle<std::vector<Jet>>(name)),
jetradius_(jetradius) {}

bool RecoJetLeptonCleaner::process(uhh2::Event& event){
  std::vector<Jet> jets = event.get(h_jets);
  std::vector<Jet> cleaned_jets;
  Jet jet;

  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }


  // perform cleaning
  TLorentzVector jet_v4, lepton_v4, jetlep_v4;
  for(unsigned int i = 0; i < jets.size(); ++i){
    jet = jets.at(i);

    if(deltaR(lepton, jet) < jetradius_){

      // std::cout<<"pt before: "<< jet.pt()<<std::endl;

      jet_v4.SetPxPyPzE(jet.v4().Px(), jet.v4().Py(), jet.v4().Pz(), jet.v4().E());
      lepton_v4.SetPxPyPzE(lepton.v4().Px(), lepton.v4().Py(), lepton.v4().Pz(), lepton.v4().E());
      jetlep_v4 = jet_v4 - lepton_v4;

      jet.set_pt(jetlep_v4.Pt());
      jet.set_eta(jetlep_v4.Eta());
      jet.set_phi(jetlep_v4.Phi());
      jet.set_energy(jetlep_v4.E());

      // std::cout<<"pt after: "<< jet.pt()<<std::endl;

    }
    cleaned_jets.push_back(jet);

  }
  sort_by_pt<Jet>(cleaned_jets); // Sort Jets by pT
  event.set(h_jets, cleaned_jets);

  return true;
}

////////////////////////////////////////////////////////

RecoTopJetLeptonCleaner::RecoTopJetLeptonCleaner(uhh2::Context& ctx, float jetradius):
h_topjets(ctx.get_handle<std::vector<TopJet>>("topjets")),
jetradius_(jetradius) {}

bool RecoTopJetLeptonCleaner::process(uhh2::Event& event){
  std::vector<TopJet> jets = event.get(h_topjets);
  std::vector<TopJet> cleaned_jets;
  TopJet jet;

  Particle lepton;
  if(event.muons->size() > 0){
    lepton = event.muons->at(0);
  }
  else if(event.electrons->size() > 0){
    lepton = event.electrons->at(0);
  }


  // perform cleaning
  TLorentzVector jet_v4, lepton_v4, jetlep_v4;
  for(unsigned int i = 0; i < jets.size(); ++i){
    jet = jets.at(i);

    if(deltaR(lepton, jet) < jetradius_){

      // std::cout<<"pt before: "<< jet.pt()<<std::endl;

      jet_v4.SetPxPyPzE(jet.v4().Px(), jet.v4().Py(), jet.v4().Pz(), jet.v4().E());
      lepton_v4.SetPxPyPzE(lepton.v4().Px(), lepton.v4().Py(), lepton.v4().Pz(), lepton.v4().E());
      jetlep_v4 = jet_v4 - lepton_v4;

      jet.set_pt(jetlep_v4.Pt());
      jet.set_eta(jetlep_v4.Eta());
      jet.set_phi(jetlep_v4.Phi());
      jet.set_energy(jetlep_v4.E());

      // std::cout<<"pt after: "<< jet.pt()<<std::endl;

    }
    cleaned_jets.push_back(jet);

  }
  sort_by_pt<TopJet>(cleaned_jets); // Sort Jets by pT
  event.set(h_topjets, cleaned_jets);

  return true;
}

////////////////////////////////////////////////////////
// uhh2::TopJetMassCut::TopJetMassCut(const):
//   {}

bool uhh2::TopJetMassCut::passes(const uhh2::Event& event){

  assert(event.jets);
  if(event.topjets->size()<2) return false;
  float mass1, mass2;
  const Particle* Top1 = &event.topjets->at(0);
  const Particle* Top2 = &event.topjets->at(1);
  mass1 = Top1->v4().M();
  mass2 = Top2->v4().M();

  return (mass1 > mass2);
}

////////////////////////////////////////////////////////
uhh2::deltaRCut::deltaRCut(float max_deltaR):
max_deltaR_(max_deltaR) {}

bool uhh2::deltaRCut::passes(const uhh2::Event& event){

  assert(event.muons && event.electrons && event.jets);
  if(event.topjets->size()<2) return false;
  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- deltaRCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  float dR=100;
  if(event.muons->size()!=0){
    dR = deltaR(event.muons->at(0), event.topjets->at(1));
  }
  else dR = deltaR(event.electrons->at(0), event.topjets->at(1));

  return (dR < max_deltaR_);
}

////////////////////////////////////////////////////////
uhh2::deltaRCutAK4::deltaRCutAK4(float max_deltaR_ak4):
max_deltaR_ak4_(max_deltaR_ak4) {}

bool uhh2::deltaRCutAK4::passes(const uhh2::Event& event){

  assert(event.muons && event.electrons && event.jets);
  if(event.jets->size()<2) return false;
  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- deltaRCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  float dR1=100;
  float dR2=100;
  float dR=100;
  if(event.muons->size()!=0){
    dR1 = deltaR(event.muons->at(0), event.jets->at(0));
    dR2 = deltaR(event.muons->at(0), event.jets->at(1));
    if(dR1<dR2){
      dR=dR1;
    }
    else dR=dR2;
  }
  else{
    dR1 = deltaR(event.electrons->at(0), event.jets->at(0));
    dR2 = deltaR(event.electrons->at(0), event.jets->at(1));
    if(dR1<dR2){
      dR=dR1;
    }
    else dR=dR2;
  }
  return (dR < max_deltaR_ak4_);
}

////////////////////////////////////////////////////////
uhh2::Jet2Cut::Jet2Cut(float min_pt):
min_pt_(min_pt) {}

bool uhh2::Jet2Cut::passes(const uhh2::Event& event){
  const Particle* jet2 = &event.jets->at(1);
  float jet_pt = jet2->pt();

  return jet_pt > min_pt_;

}


////////////////////////////////////////////////////////

uhh2::HTlepCut::HTlepCut(float min_htlep, float max_htlep):
min_htlep_(min_htlep), max_htlep_(max_htlep) {}

bool uhh2::HTlepCut::passes(const uhh2::Event& event){

  auto met = event.met->pt();
  double htlep = 0.0;
  double ht_elec = 0.0;
  double ht_muon = 0.0;
  for(const auto & electron : *event.electrons){
    ht_elec += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_muon += muon.pt();
  }
  htlep = ht_elec + ht_muon + met;

  return (htlep > min_htlep_) && (htlep < max_htlep_);
}
////////////////////////////////////////////////////////

uhh2::METCut::METCut(float min_met, float max_met):
min_met_(min_met), max_met_(max_met) {}

bool uhh2::METCut::passes(const uhh2::Event& event){

  assert(event.met);

  float MET = event.met->pt();
  return (MET > min_met_) && (MET < max_met_);
}
////////////////////////////////////////////////////////

//!! uhh2::TTbarNJetSelection::TTbarNJetSelection(const float jet1_pt, const float jet2_pt, const float tjet1_pt):
//!!   jet1_pt_(jet1_pt), jet2_pt_(jet2_pt), tjet1_pt_(tjet1_pt) {
//!!
//!!   h_jets_    = ctx.get_handle<std::vector<Jet>   >("jets");
//!!   h_topjets_ = ctx.get_handle<std::vector<TopJet>>("topjets");
//!!
//!! }
//!!
//!! bool uhh2::TTbarNJetSelection::passes(const uhh2::Event& event){
//!!
//!!   const std::vector<Jet>&       jets = event.get(h_jets_);
//!!   const std::vector<TopJet>& topjets = event.get(h_topjets_);
//!!
//!!   int jetN_1(0), jetN_2(0);
//!!   for(const auto& j : jets){
//!!
//!!     if(j.pt() > jet1_pt_) ++jetN_1;
//!!     if(j.pt() > jet2_pt_) ++jetN_2;
//!!   }
//!!
//!!   int tjetN_1(0);
//!!   for(const auto& tj : topjets){
//!!
//!!     if(tj.pt() > tjet1_pt_) ++tjetN_1;
//!!   }
//!!
//!!   return (njet >= nmin) && (njet <= nmax);
//!! }
//!! ////////////////////////////////////////////////////////

bool uhh2::TwoDCut1::passes(const uhh2::Event& event){

  if(event.muons->size() != 0 || event.electrons->size() != 0){
    assert(event.muons && event.electrons && event.jets);

    const Particle* lepton = &event.muons->at(0);
    //if(event.muons->size()!=0) lepton = event.muons->at(0);
      //leading_lepton(event);

    float drmin, ptrel;
    std::tie(drmin, ptrel) = drmin_pTrel(*lepton, *event.jets);

    return ((drmin > min_deltaR_) || (ptrel > min_pTrel_));
  }
  else return true;
}
////////////////////////////////////////////////////////

bool uhh2::TwoDCutALL::passes(const uhh2::Event& event){

  assert(event.muons && event.electrons && event.jets);

  for(const auto& muo : *event.muons){

    float drmin, ptrel;
    std::tie(drmin, ptrel) = drmin_pTrel(muo, *event.jets);

    const bool pass = (drmin > min_deltaR_) || (ptrel > min_pTrel_);
    if(!pass) return false;
  }

  for(const auto& ele : *event.electrons){

    float drmin, ptrel;
    std::tie(drmin, ptrel) = drmin_pTrel(ele, *event.jets);

    const bool pass = (drmin > min_deltaR_) || (ptrel > min_pTrel_);
    if(!pass) return false;
  }

  return true;
}
////////////////////////////////////////////////////////

bool uhh2::TwoDCut::passes(const uhh2::Event& event){

  assert(event.muons && event.electrons && event.jets);
  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- TwoDCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  float drmin, ptrel;
  if(event.muons->size()) std::tie(drmin, ptrel) = drmin_pTrel(event.muons->at(0), *event.jets);
  else std::tie(drmin, ptrel) = drmin_pTrel(event.electrons->at(0), *event.jets);

  return (drmin > min_deltaR_) || (ptrel > min_pTrel_);
}
////////////////////////////////////////////////////////

bool uhh2::TriangularCuts::passes(const uhh2::Event& event){

  assert(event.muons || event.electrons);
  assert(event.jets && event.met);

  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- TriangularCuts::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  if(!event.jets->size()){
    std::cout << "\n @@@ WARNING -- TriangularCuts::passes -- unexpected number of jets in the event (==0). returning 'false'\n";
    return false;
  }

  // pt-leading charged lepton
  const Particle* lep1 = &event.muons->at(0);//leading_lepton(event);

  // 1st entry in jet collection (should be the pt-leading jet)
  const Particle* jet1 = &event.jets->at(0);

  // MET-lepton triangular cut
  bool pass_tc_lep = fabs(fabs(deltaPhi(*event.met, *lep1)) - a_) < b_ * event.met->pt();

  // MET-jet triangular cut
  bool pass_tc_jet = fabs(fabs(deltaPhi(*event.met, *jet1)) - a_) < b_ * event.met->pt();

  return pass_tc_lep && pass_tc_jet;
}
////////////////////////////////////////////////////////

bool uhh2::TriangularCutsELE::passes(const uhh2::Event& event){

  assert(event.electrons);
  assert(event.jets && event.met);

  if(event.electrons->size() != 1) std::runtime_error("TriangularCutsELE::passes -- unexpected number of electrons in the event (!=1)");
  if(event.jets     ->size() == 0) std::runtime_error("TriangularCutsELE::passes -- unexpected number of jets in the event (==0)");

  // pt-leading charged lepton
  const Particle* lep1 = &event.electrons->at(0);

  // 1st entry in jet collection (should be the pt-leading jet)
  const Particle* jet1 = &event.jets->at(0);

  // MET-lepton triangular cut
  bool pass_tc_lep = fabs(fabs(deltaPhi(*event.met, *lep1)) - a_) < b_ * event.met->pt();

  // MET-jet triangular cut
  bool pass_tc_jet = fabs(fabs(deltaPhi(*event.met, *jet1)) - a_) < b_ * event.met->pt();

  return pass_tc_lep && pass_tc_jet;
}
////////////////////////////////////////////////////////

uhh2::DiLeptonSelection::DiLeptonSelection(const std::string& channel, const bool opposite_charge, const bool veto_other_flavor):
channel_(channel), opposite_charge_(opposite_charge), veto_other_flavor_(veto_other_flavor) {}

bool uhh2::DiLeptonSelection::passes(const uhh2::Event& event){

  bool pass(false);

  assert(event.muons && event.electrons);

  if(channel_ == "muon"){

    pass = (event.muons->size() == 2);

    if(pass && opposite_charge_)   pass &= ((event.muons->at(0).charge() * event.muons->at(1).charge()) == -1);
    if(pass && veto_other_flavor_) pass &= (event.electrons->size() == 0);
  }
  else if(channel_ == "elec"){

    pass = (event.electrons->size() == 2);

    if(pass && opposite_charge_)   pass &= ((event.electrons->at(0).charge() * event.electrons->at(1).charge()) == -1);
    if(pass && veto_other_flavor_) pass &= (event.muons->size() == 0);
  }
  else throw std::runtime_error("DiLeptonSelection::passes -- undefined key for lepton channel: "+channel_);

  return pass;
}
////////////////////////////////////////////////////////

uhh2::TopJetPlusJetEventSelection::TopJetPlusJetEventSelection(const float topjet_minDR_jet, const float jet_min_pt):
topjet_minDR_jet_(topjet_minDR_jet), jet_min_pt_(jet_min_pt) {}

bool uhh2::TopJetPlusJetEventSelection::passes(const uhh2::Event& event){

  if(event.topjets->size() != 1) return false;

  for(const auto& topjet : *event.topjets){

    for(const auto& jet : *event.jets){
      if(uhh2::deltaR(topjet, jet) > topjet_minDR_jet_ && jet.pt() > jet_min_pt_) return true;
    }
  }

  return false;
}
////////////////////////////////////////////////////////

uhh2::TopTagEventSelection::TopTagEventSelection(const TopJetId& tjetID, const float minDR_jet_ttag):
topjetID_(tjetID), minDR_jet_toptag_(minDR_jet_ttag) {

  topjet1_sel_.reset(new NTopJetSelection(1, -1, topjetID_));
}

bool uhh2::TopTagEventSelection::passes(const uhh2::Event& event){

  if(!topjet1_sel_->passes(event)) return false;

  for(const auto& topjet : *event.topjets){

    if(!topjetID_(topjet, event)) continue;

    for(const auto& jet : *event.jets){
      if(deltaR(jet, topjet) > minDR_jet_toptag_) return true;
    }
  }

  return false;
}
////////////////////////////////////////////////////////

uhh2::LeptonicTopPtCut::LeptonicTopPtCut(uhh2::Context& ctx, float pt_min, float pt_max, const std::string& hyps_name, const std::string& disc_name):
tlep_pt_min_(pt_min), tlep_pt_max_(pt_max), h_hyps_(ctx.get_handle<std::vector<ReconstructionHypothesis>>(hyps_name)), disc_name_(disc_name) {}

bool uhh2::LeptonicTopPtCut::passes(const uhh2::Event& event){

  const std::vector<ReconstructionHypothesis>& hyps = event.get(h_hyps_);

  const ReconstructionHypothesis* hyp = get_best_hypothesis(hyps, disc_name_);
  if(!hyp) std::runtime_error("LeptonicTopPtCut -- best hypothesis not found (discriminator="+disc_name_+")");

  const float tlep_pt = hyp->toplep_v4().Pt();

  return (tlep_pt > tlep_pt_min_) && (tlep_pt < tlep_pt_max_);
}
////////////////////////////////////////////////////////

uhh2::HypothesisDiscriminatorCut::HypothesisDiscriminatorCut(uhh2::Context& ctx, float disc_min, float disc_max, const std::string& hyps_name, const std::string& disc_bhyp, const std::string& disc_cut):
disc_min_(disc_min), disc_max_(disc_max), h_hyps_(ctx.get_handle<std::vector<ReconstructionHypothesis>>(hyps_name)), disc_bhyp_(disc_bhyp), disc_cut_(disc_cut) {}

bool uhh2::HypothesisDiscriminatorCut::passes(const uhh2::Event& event){

  const std::vector<ReconstructionHypothesis>& hyps = event.get(h_hyps_);

  const ReconstructionHypothesis* hyp = get_best_hypothesis(hyps, disc_bhyp_);
  if(!hyp) std::runtime_error("HypothesisDiscriminatorCut -- best hypothesis not found (discriminator="+disc_bhyp_+")");

  const float disc_val = hyp->discriminator(disc_cut_);

  return (disc_val > disc_min_) && (disc_val < disc_max_);
}
////////////////////////////////////////////////////////

uhh2::GenFlavorSelection::GenFlavorSelection(const std::string& flav_key){

  flavor_key_ = flav_key;

  if(flavor_key_ != "l" && flavor_key_ != "c" && flavor_key_ != "b")
  throw std::runtime_error("GenFlavorSelection::GenFlavorSelection -- undefined key for parton flavor (must be 'l', 'c' or 'b'): "+flavor_key_);
}

bool uhh2::GenFlavorSelection::passes(const uhh2::Event& event){

  bool pass(false);

  assert(event.genparticles);

  int bottomN(0), charmN(0);
  for(const auto& genp : *event.genparticles){

    if(!(20 <= genp.status() && genp.status() <= 30)) continue;
    if(genp.mother1() == (unsigned short)(-1)) continue;
    if(genp.mother2() == (unsigned short)(-1)) continue;

    const int id = genp.pdgId();

    if(std::abs(id) == 5) ++bottomN;
    if(std::abs(id) == 4) ++charmN;
  }

  if     (flavor_key_ == "b") pass = (bottomN >= 1);
  else if(flavor_key_ == "c") pass = (bottomN == 0 && charmN >= 1);
  else if(flavor_key_ == "l") pass = (bottomN == 0 && charmN == 0);
  else throw std::runtime_error("GenFlavorSelection::GenFlavorSelection -- undefined key for parton flavor (must be 'l', 'c' or 'b'): "+flavor_key_);

  return pass;
}
////////////////////////////////////////////////////////

uhh2::JetFlavorSelection::JetFlavorSelection(const std::string& flav_key){

  flavor_key_ = flav_key;

  if(flavor_key_ != "l" && flavor_key_ != "c" && flavor_key_ != "b")
  throw std::runtime_error("JetFlavorSelection::JetFlavorSelection -- undefined key for jet flavor (must be 'l', 'c' or 'b'): "+flavor_key_);
}

bool uhh2::JetFlavorSelection::passes(const uhh2::Event& event){

  bool pass(false);

  assert(event.jets);

  int bottomN(0), charmN(0);
  for(const auto& j : *event.jets){

    const int id = j.hadronFlavour();

    if(std::abs(id) == 5) ++bottomN;
    if(std::abs(id) == 4) ++charmN;
  }

  if     (flavor_key_ == "b") pass = (bottomN >= 1);
  else if(flavor_key_ == "c") pass = (bottomN == 0 && charmN >= 1);
  else if(flavor_key_ == "l") pass = (bottomN == 0 && charmN == 0);
  else throw std::runtime_error("JetFlavorSelection::JetFlavorSelection -- undefined key for jet flavor (must be 'l', 'c' or 'b'): "+flavor_key_);

  return pass;
}
////////////////////////////////////////////////////////

uhh2::GenHTCut::GenHTCut(uhh2::Context& ctx, const float min, const float max, const std::string& meps_name):
genHT_min_(min), genHT_max_(max), h_meps_(ctx.get_handle<std::vector<GenParticle> >(meps_name)) {}

bool uhh2::GenHTCut::passes(const uhh2::Event& event){

  const std::vector<GenParticle>& me_partons = event.get(h_meps_);

  float genHT(0.);
  for(const auto& p : me_partons) genHT += p.pt();

  return (genHT_min_ <= genHT) && (genHT < genHT_max_);
}
////////////////////////////////////////////////////////

