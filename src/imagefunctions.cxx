#include <iostream>
#include <memory>

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/BoostedTTbarColourFlow/include/imagefunctions.h"
#include "UHH2/common/include/JetIds.h"

using namespace std;
using namespace uhh2;

double CalculateRadius(Jet jet1){
  double pi = 3.1415926535;
  double areajet = jet1.jetArea();
  double radiusjet = sqrt(areajet/pi);

  return radiusjet;

}

vector<Jet> jets_for_images(const Event & event, const TopJet HadronicTopJet){

  JetId btag_tight = CSVBTag(CSVBTag::WP_TIGHT);

  vector<Jet> jets_to_return;
  std::vector<Muon>* muons = event.muons;

  Jet bjet;
  Jet lightjet1;
  Jet lightjet2;

  if(event.muons->size()==0)
    return jets_to_return;

  for (const Jet & thisjet : HadronicTopJet.subjets()){
    //checking if the associated AK4 jet is tight b-tagged
    bool bjet_match=false;
    for(const auto & thisAK4jet: *event.jets){
      if(deltaR(thisjet,thisAK4jet)<0.2){
	if(!btag_tight(thisjet, event)) bjet_match=true;
      }
    }

    if(!bjet_match) continue;
    if(deltaR(thisjet, muons->at(0))<1.5) continue;
    if(thisjet.pt()<30) continue;
    if(fabs(thisjet.eta())>2.4) continue;
    bjet = thisjet;
  }

  double deltaRmin=10;

  for (const Jet & thisjet : HadronicTopJet.subjets()){
    if(thisjet.pt()<30) continue;
    if(fabs(thisjet.eta())>2.4) continue;
    if(deltaR(thisjet, muons->at(0))<1.5) continue;
    if(deltaR(thisjet, bjet)<deltaRmin && deltaR(thisjet, bjet)>0.4 && deltaR(thisjet, bjet)<1.0){
      deltaRmin=deltaR(thisjet, bjet);
      lightjet1=thisjet;
    }
  }

  deltaRmin=10;

  for (const Jet & thisjet : HadronicTopJet.subjets()){
    if(thisjet.pt()<30) continue;
    if(deltaR(thisjet, muons->at(0))<1.5) continue;
    if(fabs(thisjet.eta())>2.4) continue;
    if(deltaR(thisjet, bjet)>0.4 && deltaR(thisjet, lightjet1)>0.4 && deltaR(thisjet, lightjet1)<1.0 && deltaR(thisjet, lightjet1)<deltaRmin){
      deltaRmin=deltaR(thisjet, lightjet1);
      lightjet2=thisjet;
    }
  }
  
  jets_to_return.push_back(bjet);
  jets_to_return.push_back(lightjet1);
  jets_to_return.push_back(lightjet2);
  return jets_to_return;
}

vector<Jet> jets_for_images_withWMass(const Event & event, const TopJet HadronicTopJet){

  vector<Jet> jets_to_return;
  std::vector<Muon>* muons = event.muons;

  Jet bjet;
  Jet lightjet1;
  Jet lightjet2;

  if(event.muons->size()==0)
    return jets_to_return;  
  
  double massW=-10;
  double massWdiff = 100;
  
  double massW_1 = (HadronicTopJet.subjets().at(0).v4()+HadronicTopJet.subjets().at(1).v4()).mass();
  double massW_2 = (HadronicTopJet.subjets().at(0).v4()+HadronicTopJet.subjets().at(2).v4()).mass();
  double massW_3 = (HadronicTopJet.subjets().at(1).v4()+HadronicTopJet.subjets().at(2).v4()).mass();
  
  if(fabs(massW_1-80.3) < massWdiff) { massW = massW_1; massWdiff = fabs(massW_1-80.3); }
  if(fabs(massW_2-80.3) < massWdiff) { massW = massW_2; massWdiff = fabs(massW_2-80.3); }
  if(fabs(massW_3-80.3) < massWdiff) { massW = massW_3; massWdiff = fabs(massW_3-80.3); }
  
  if(massW == massW_1){
    bjet = HadronicTopJet.subjets().at(2);
    lightjet1 = HadronicTopJet.subjets().at(1);
    lightjet2 = HadronicTopJet.subjets().at(0);
  }
  
  else if(massW == massW_1){
    bjet = HadronicTopJet.subjets().at(1);
    lightjet1 =HadronicTopJet.subjets().at(0);
    lightjet2 = HadronicTopJet.subjets().at(2);
  }
  
  else if(massW == massW_1){
    bjet = HadronicTopJet.subjets().at(0);
    lightjet1 =HadronicTopJet.subjets().at(2);
    lightjet2 = HadronicTopJet.subjets().at(1);
  }

  jets_to_return.push_back(bjet);
  jets_to_return.push_back(lightjet1);
  jets_to_return.push_back(lightjet2);

  return jets_to_return;

}

vector<double> info_particles(const Event & event, const Jet jet_to_cluster){

  if (event.pfparticles == nullptr) { throw std::runtime_error("particles is nullptr"); }

  vector<double> infoparticles_to_return;
  int particles_total=0;

  //for (const PFParticle & candInd : *event.pfparticles){

    //here one needs to run on the indices of the jet
  double radius = CalculateRadius(jet_to_cluster);
    
  for (const PFParticle & candInd : *event.pfparticles){
    
    if(deltaR(candInd,jet_to_cluster)>radius) continue;

    float puppiWeight = candInd.puppiWeight();

    if(candInd.pt()<2) continue;
    
    infoparticles_to_return.push_back(candInd.pt()*puppiWeight/jet_to_cluster.pt());
    infoparticles_to_return.push_back(candInd.eta()-jet_to_cluster.eta());
    infoparticles_to_return.push_back(deltaPhi(candInd,jet_to_cluster));

    particles_total++;
    if(particles_total>=25) break;
  }

  if(infoparticles_to_return.size()<75){
    for(int i=infoparticles_to_return.size();i<75;){
      infoparticles_to_return.push_back(-10);
      infoparticles_to_return.push_back(-10);
      infoparticles_to_return.push_back(-10);
      i=i+3;
    }
  }

  return infoparticles_to_return;
}
