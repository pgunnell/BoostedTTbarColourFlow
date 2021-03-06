#pragma once

#include <iostream>
#include <memory>

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/JetIds.h"

std::vector<Jet> jets_for_images(const uhh2::Event & event, const TopJet HadronicTopJet);
std::vector<Jet> jets_for_images_withWMass(const uhh2::Event & event, const TopJet HadronicTopJet);

std::vector<double> info_particles(const uhh2::Event & event, const Jet jet_to_cluster);
