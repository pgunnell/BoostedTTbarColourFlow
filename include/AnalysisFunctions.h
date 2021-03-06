#pragma once

#include <iostream>
#include <memory>

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"


double CalculateRadiusJet(Jet jet1);
double DeltaPhi(double p1, double p2);
std::vector<double> FillRapidityPhiPt(PFParticle candInd, Jet jet1, Jet jet2, double weight, double rapidityleadsubjet, double phileadsubjet, double rotationangle, double NormFactor);
std::vector<double> GetPull(PFParticle candInd, std::vector<double> Pull, double jetpt, double jetrapidity, double jetphi);
TVector3 CalculatePullRivet(Jet lightjet1, const Event & event);
double CalculatePullAngle(Jet jet1, Jet axisjet, const Event & event);
double CalculatePullGenAngle(GenJet jet1, GenJet axisjet, const Event & event);
double DeltaPhiGen(double p1, double p2);
TVector3 CalculatePullGenRivet(GenJet lightjet1, const Event & event);
std::vector<double> FillRapidityPhiPt_Gen(GenParticle candInd, GenJet jet1, GenJet jet2, double weight, double rapidityleadsubjet, double phileadsubjet, double rotationangle, double NormFactor);

