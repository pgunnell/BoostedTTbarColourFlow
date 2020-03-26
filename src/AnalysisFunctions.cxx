#include <iostream>
#include <memory>

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/BoostedTTbarColourFlow/include/AnalysisFunctions.h"

using namespace std;
using namespace uhh2;

double CalculateRadiusJet(Jet jet1){
  double pi = 3.1415926535;
  double areajet = jet1.jetArea();
  double radiusjet = sqrt(areajet/pi);

  return radiusjet;

}

double DeltaPhi(double p1, double p2){
  double pi = 3.1415926535;
  double deltaphi = fabs(p1 - p2);
  if(deltaphi > pi) deltaphi = 2* pi - deltaphi;
  return deltaphi;
}

vector<double> FillRapidityPhiPt(PFParticle candInd, Jet jet1, Jet jet2, double weight, double rapidityleadsubjet, double phileadsubjet, double rotationangle, double NormFactor){

  //we set the reference to the center of the leading subjet (among the two we have)                                                                                                                                                          
  double rapiditytransl=(candInd.v4().Rapidity()-rapidityleadsubjet)/NormFactor;
  double phitransl=DeltaPhi(candInd.phi(),phileadsubjet)/NormFactor;

  //we apply a rotation such that the two jets are on the horizontal axis                                                                                                                                                                     
  double rapidityrot=rapiditytransl*cos(rotationangle)+phitransl*sin(rotationangle);
  double phirot=phitransl*cos(rotationangle)-rapiditytransl*sin(rotationangle);

  //we normalize the candidate pTs to the pT of the subjet, which they belong to                                                                                                                                                              
  double weightpt=1.;
  //weightpt=candInd.pt()*candInd.puppiWeight()/jet1.pt()*weight;                                                                                                                                                                             
  weightpt=candInd.pt()/jet1.pt()*weight;

  vector<double> outputvector;
  outputvector.push_back(rapidityrot);
  outputvector.push_back(phirot);
  outputvector.push_back(weightpt);

  return outputvector;
}

vector<double> GetPull(PFParticle candInd, vector<double> Pull, double jetpt, double jetrapidity, double jetphi){

  double r_module = sqrt(pow(candInd.v4().Rapidity()-jetrapidity,2)+pow(DeltaPhi(candInd.phi(),jetphi),2));
  //double Pull_after_rapidity = Pull.at(0) + candInd.pt()*candInd.puppiWeight()*(r_module)/jetpt*(candInd.v4().Rapidity()-jetrapidity);                                                                                                      
  //double Pull_after_phi = Pull.at(1) + candInd.pt()*candInd.puppiWeight()*(r_module)/jetpt*DeltaPhi(candInd.phi(),jetphi);                                                                                                                  
  double Pull_after_rapidity = Pull.at(0) + candInd.pt()*(sqrt(r_module))/jetpt*(candInd.v4().Rapidity()-jetrapidity);
  double Pull_after_phi = Pull.at(1) + candInd.pt()*(sqrt(r_module))/jetpt*DeltaPhi(candInd.phi(),jetphi);

  vector<double> Pull_after;
  Pull_after.push_back(Pull_after_rapidity);
  Pull_after.push_back(Pull_after_phi);

  return Pull_after;

}

TVector3 CalculatePullRivet(Jet lightjet1, const Event & event){
  TVector3 pull(0.0, 0.0, 0.0);
  double PT = lightjet1.pt();
  LorentzVector axis;

  double radiusjet = CalculateRadiusJet(lightjet1);

  for (const PFParticle & candInd : *event.pfparticles){
    if(deltaR(candInd,lightjet1)>radiusjet) continue;
    if(candInd.charge()==0) continue;
    if(candInd.puppiWeight()<0.4) continue;

    // calculate axis                                                                                                                                                                                                                         
    axis += candInd.v4();//*candInd.puppiWeight();                                                                                                                                                                                            
  }

  TVector3 J(axis.Rapidity(), axis.phi(), 0.0);
  // calculate pull                                                                                                                                                                                                                           
  for (const PFParticle & candInd : *event.pfparticles){
    if(deltaR(candInd,lightjet1)>radiusjet) continue;
    if(candInd.charge()==0) continue;
    if(candInd.puppiWeight()<0.4) continue;

    TVector3 ri;
    //ri.SetXYZ(candInd.v4().Rapidity()*candInd.puppiWeight(), candInd.phi()*candInd.puppiWeight(), 0.0);                                                                                                                                     
    ri.SetXYZ(candInd.v4().Rapidity(), candInd.phi(), 0.0);
    ri = ri-J;
    double moduleri = sqrt(ri.x()*ri.x()+ri.y()*ri.y()+ri.z()*ri.z());

    while (ri.y() >  3.1415) ri.SetY(ri.y() - 6.2830);
    while (ri.y() < -3.1415) ri.SetY(ri.y() + 6.2830);
    //pull.SetX(pull.x() + (moduleri * ri.x() * candInd.puppiWeight() * candInd.pt()) / PT);                                                                                                                                                  
    //pull.SetY(pull.y() + (moduleri * ri.y() * candInd.puppiWeight() * candInd.pt()) / PT);                                                                                                                                                  
    pull.SetX(pull.x() + (moduleri * ri.x() * candInd.pt()) / PT);
    pull.SetY(pull.y() + (moduleri * ri.y() * candInd.pt()) / PT);
  }

  return pull;

}

double CalculatePullAngle(Jet jet1, Jet axisjet, const Event & event) {
  TVector3 pull_vector = CalculatePullRivet(jet1, event);
  pull_vector.SetXYZ(1000.*pull_vector.x(), 1000.*pull_vector.y(), 0.);
  double drap = axisjet.v4().Rapidity() - jet1.v4().Rapidity();
  double dphi = DeltaPhi(axisjet.phi(),jet1.phi());
  TVector3 j2_vector(drap, dphi, 0.0);
  //return mapAngleMPiToPi(deltaPhi(pull_vector, j2_vector));                                                                                                                                                                                 
  return DeltaPhi(pull_vector.Phi(), j2_vector.Phi());
}

double DeltaPhiGen(double p1, double p2){
  double pi = 3.1415926535;
  double deltaphi = fabs(p1 - p2);
  if(deltaphi > pi) deltaphi = 2* pi - deltaphi;
  return deltaphi;
}

TVector3 CalculatePullGenRivet(GenJet lightjet1, const Event & event){
  TVector3 pull(0.0, 0.0, 0.0);
  double PT = lightjet1.pt();
  LorentzVector axis;

  for (const auto candInd : lightjet1.genparticles_indices()) {

    if(event.genparticles->at(candInd).charge()==0) continue;

    // calculate axis                                                                                                                                                                                                                         
    axis += event.genparticles->at(candInd).v4();
  }

  TVector3 J(axis.Rapidity(), axis.phi(), 0.0);
  // calculate pull                                                                                                                                                                                                                           
  for (const auto candInd : lightjet1.genparticles_indices()) {
    if(event.genparticles->at(candInd).charge()==0) continue;

    TVector3 ri;
    ri.SetXYZ(event.genparticles->at(candInd).v4().Rapidity(), event.genparticles->at(candInd).phi(), 0.0);
    ri = ri-J;
    double moduleri = sqrt(ri.x()*ri.x()+ri.y()*ri.y()+ri.z()*ri.z());

    while (ri.y() >  3.1415) ri.SetY(ri.y() - 6.2830);
    while (ri.y() < -3.1415) ri.SetY(ri.y() + 6.2830);
    pull.SetX(pull.x() + (moduleri * ri.x() * event.genparticles->at(candInd).pt()) / PT);
    pull.SetY(pull.y() + (moduleri * ri.y() * event.genparticles->at(candInd).pt()) / PT);
  }

  return pull;

}

double CalculatePullGenAngle(GenJet jet1, GenJet axisjet, const Event & event) {
  TVector3 pull_vector = CalculatePullGenRivet(jet1, event);
  pull_vector.SetXYZ(1000.*pull_vector.x(), 1000.*pull_vector.y(), 0.);
  double drap = axisjet.v4().Rapidity() - jet1.v4().Rapidity();
  double dphi = DeltaPhiGen(axisjet.phi(),jet1.phi());
  TVector3 j2_vector(drap, dphi, 0.0);
  //return mapAngleMPiToPi(deltaPhi(pull_vector, j2_vector));                                                                                                                                                                                 
  return DeltaPhiGen(pull_vector.Phi(), j2_vector.Phi());
}

vector<double> FillRapidityPhiPt_Gen(GenParticle candInd, GenJet jet1, GenJet jet2, double weight, double rapidityleadsubjet, double phileadsubjet, double rotationangle, double NormFactor){

  //we set the reference to the center of the leading subjet (among the two we have)                                                                                                                                                          
  double rapiditytransl=(candInd.v4().Rapidity()-rapidityleadsubjet)/NormFactor;
  double phitransl=DeltaPhiGen(candInd.phi(),phileadsubjet)/NormFactor;

  //we apply a rotation such that the two jets are on the horizontal axis                                                                                                                                                                     
  double rapidityrot=rapiditytransl*cos(rotationangle)+phitransl*sin(rotationangle);
  double phirot=phitransl*cos(rotationangle)-rapiditytransl*sin(rotationangle);

  //we normalize the candidate pTs to the pT of the subjet, which they belong to                                                                                                                                                             
  double weightpt=1.;
  if(deltaR(candInd,jet1)<0.4){
    weightpt=candInd.pt()/jet1.pt()*weight;
  }

  if(deltaR(candInd,jet2)<0.4){
    weightpt=candInd.pt()/jet2.pt()*weight;
  }

  vector<double> outputvector;
  outputvector.push_back(rapidityrot);
  outputvector.push_back(phirot);
  outputvector.push_back(weightpt);

  return outputvector;
}
