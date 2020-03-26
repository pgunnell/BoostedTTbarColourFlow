#pragma once

#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Selection.h>
#include <UHH2/core/include/Utils.h>
#include <UHH2/core/include/LorentzVector.h>

#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/ReconstructionHypothesis.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/TTbarGen.h>

#include <UHH2/MTopJet/include/MTopJetUtils.h>
#include <UHH2/MTopJet/include/utils.h>
#include <UHH2/MTopJet/include/GenSelections.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>



namespace uhh2 {

  ////////////////////////////////////////////////////////////////
  class ElectronEtaVeto : public Selection {

  public:
    explicit ElectronEtaVeto(double, double);
    virtual bool passes(const Event&) override;

  private:
    double lower, upper;
  };

  ////////////////////////////////////////////////////////////////
  class NJetXCone : public Selection {

  public:
    explicit NJetXCone(Context&, const std::string &, unsigned int);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_jets;
    unsigned int njet_;
  };
  ////////////////////////////////////////////////////////////////
  class LeadingRecoJetPT : public Selection {

  public:
    explicit LeadingRecoJetPT(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_jets;
    float ptcut_;
  };
  ////////////////////////////////////////////////////////////////
  class LeadingRecoJetETA : public Selection {

  public:
    explicit LeadingRecoJetETA(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_jets;
    float etacut_;
  };
  ////////////////////////////////////////////////////////////////
  class SubjetQuality : public Selection {

  public:
    explicit SubjetQuality(Context&, const std::string &, float, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_jets;
    float ptmin, etamax;
  };
  ////////////////////////////////////////////////////////////////
  class MassCutReco : public Selection {

  public:
    explicit MassCutReco(Context&, const std::string &);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
  };
  ////////////////////////////////////////////////////////////////
  class MassCutXCone : public Selection {

  public:
    explicit MassCutXCone(Context&, const std::string &, const std::string &);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_hadjets;
    uhh2::Event::Handle<std::vector<TopJet>> h_lepjets;
  };
  ////////////////////////////////////////////////////////////////
  class DeltaRCutReco : public Selection {

  public:
    explicit DeltaRCutReco(Context&, const std::string &, float);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_jets;
    float jetradius_;
  };
  ////////////////////////////////////////////////////////////////
  class NRecoJets : public Selection {

  public:
    explicit NRecoJets(Context&, const std::string &, float min_pt = 0, float min = 0, float max = 9999);
    virtual bool passes(const Event&) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    float min_pt_;
    float min_;
    float max_;
  };
  ////////////////////////////////////////////////////////////////
  class RecoJetLeptonCleaner : public uhh2::AnalysisModule {

  public:
    explicit RecoJetLeptonCleaner(Context&, const std::string &, float);
    virtual bool process(Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<Jet>> h_jets;
    float jetradius_;
  };
  ////////////////////////////////////////////////////////////////
  class RecoTopJetLeptonCleaner : public uhh2::AnalysisModule {

  public:
    explicit RecoTopJetLeptonCleaner(Context&, float);
    virtual bool process(Event& ) override;

  private:
    uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
    float jetradius_;
  };
  ////////////////////////////////////////////////////////////////
  /* class TTbarGenSemilep : public Selection { */

  /* public: */
  /*   virtual bool passes (const Event&) override; */

  /* }; */
  ////////////////////////////////////////////////////////////////
  class TopJetMassCut : public Selection {

  public:
    virtual bool passes(const Event&) override;

  };
  ////////////////////////////////////////////////////////////////
  class deltaRCut : public Selection {

  public:
    explicit deltaRCut(float max_deltaR=100);
    virtual bool passes(const Event&) override;

  private:
    float max_deltaR_;
  };
  ////////////////////////////////////////////////////////////////
  class deltaRCutAK4 : public Selection {

  public:
    explicit deltaRCutAK4(float max_deltaR_ak4);
    virtual bool passes(const Event&) override;

  private:
    float max_deltaR_ak4_;
  };
  ////////////////////////////////////////////////////////////////
  class Jet2Cut : public Selection {
  public:
    explicit Jet2Cut(float min_pt=0);
    virtual bool passes(const Event&) override;

  private:
    float min_pt_;
  };
  ////////////////////////////////////////////////////////////////
  class HTlepCut : public Selection {

  public:
    explicit HTlepCut(float, float max_htlep=infinity);
    virtual bool passes(const Event&) override;

  private:
    float min_htlep_, max_htlep_;
  };
  ////////////////////////////////////////////////////////////////
  class METCut : public Selection {

  public:
    explicit METCut(float, float max_met=infinity);
    virtual bool passes(const Event&) override;

  private:
    float min_met_, max_met_;
  };
  ////////////////////////////////////////////////////////////////
  class TwoDCut : public Selection {

  public:
    explicit TwoDCut(float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {}
    virtual bool passes(const Event&) override;

  private:
    float min_deltaR_, min_pTrel_;
  };
  ////////////////////////////////////////////////////////////////
  class TwoDCut1 : public Selection {

  public:
    explicit TwoDCut1(float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {}
    virtual bool passes(const Event&) override;

  private:
    float min_deltaR_, min_pTrel_;
  };
  ////////////////////////////////////////////////////////////////
  class TwoDCutALL : public Selection {

  public:
    explicit TwoDCutALL(float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {}
    virtual bool passes(const Event&) override;

  private:
    float min_deltaR_, min_pTrel_;
  };
  ////////////////////////////////////////////////////////////////
  class TriangularCuts : public Selection {

  public:
    explicit TriangularCuts(const float a, const float b): a_(a), b_(b) {}
    virtual bool passes(const Event&) override;

  private:
    float a_, b_;
  };
  ////////////////////////////////////////////////////////////////
  class TriangularCutsELE : public Selection {

  public:
    explicit TriangularCutsELE(const float a, const float b): a_(a), b_(b) {}
    virtual bool passes(const Event&) override;

  private:
    float a_, b_;
  };
  ////////////////////////////////////////////////////////////////
  class DiLeptonSelection: public Selection {

  public:
    explicit DiLeptonSelection(const std::string&, const bool, const bool);
    virtual bool passes(const Event&) override;

  private:
    std::string channel_;
    bool opposite_charge_;
    bool veto_other_flavor_;
  };
  ////////////////////////////////////////////////////////////////
  class TopJetPlusJetEventSelection: public Selection {

  public:
    explicit TopJetPlusJetEventSelection(const float, const float);
    virtual bool passes(const Event&) override;

  private:
    float topjet_minDR_jet_;
    float jet_min_pt_;
  };
  ////////////////////////////////////////////////////////////////
  class TopTagEventSelection: public Selection {

  public:
    explicit TopTagEventSelection(const TopJetId& tjet_id=CMSTopTag(), const float minDR_jet_ttag=1.2);
    virtual bool passes(const Event&) override;

  private:
    std::unique_ptr<Selection> topjet1_sel_;
    TopJetId topjetID_;
    float minDR_jet_toptag_;
  };
  ////////////////////////////////////////////////////////////////
  class LeptonicTopPtCut: public Selection {

  public:
    explicit LeptonicTopPtCut(Context&, float, float, const std::string& hyps="TTbarReconstruction", const std::string& disc="Chi2");
    virtual bool passes(const Event&) override;

  private:
    float tlep_pt_min_, tlep_pt_max_;
    Event::Handle<std::vector<ReconstructionHypothesis>> h_hyps_;
    std::string disc_name_;
  };
  ////////////////////////////////////////////////////////////////
  class HypothesisDiscriminatorCut: public Selection {

  public:
    explicit HypothesisDiscriminatorCut(Context&, float, float, const std::string& hyps="TTbarReconstruction", const std::string& disc_bhyp="Chi2", const std::string& disc_cut="Chi2");
    virtual bool passes(const Event&) override;

  private:
    float disc_min_, disc_max_;
    Event::Handle<std::vector<ReconstructionHypothesis>> h_hyps_;
    std::string disc_bhyp_;
    std::string disc_cut_;
  };
  ////////////////////////////////////////////////////////////////
  class GenFlavorSelection: public Selection {

  public:
    explicit GenFlavorSelection(const std::string&);
    virtual bool passes(const Event&) override;

  private:
    std::string flavor_key_;
  };
  ////////////////////////////////////////////////////////////////
  class JetFlavorSelection: public Selection {

  public:
    explicit JetFlavorSelection(const std::string&);
    virtual bool passes(const Event&) override;

  private:
    std::string flavor_key_;
  };
  ////////////////////////////////////////////////////////////////
  class GenHTCut : public Selection {

  public:
    explicit GenHTCut(Context&, const float, const float, const std::string&);
    virtual bool passes(const Event&) override;

  protected:
    float genHT_min_, genHT_max_;
    Event::Handle<std::vector<GenParticle> > h_meps_;
  };
  ////////////////////////////////////////////////////////////////
  class RunLumiEventSelection : public Selection {

  public:
    explicit RunLumiEventSelection(const std::string&, const std::string& sep=":");
    virtual ~RunLumiEventSelection() {}

    virtual bool passes(const Event&) override;
    virtual bool found (const Event&);

  protected:
    std::unordered_map<unsigned long int, std::unordered_map<unsigned long int, std::vector<unsigned long int> > > rle_map_;
  };
  ////////////////////////////////////////////////////////////////
}
