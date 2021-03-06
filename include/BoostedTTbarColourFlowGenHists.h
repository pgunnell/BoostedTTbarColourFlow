#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/common/include/JetIds.h"

namespace uhh2examples {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class BoostedTTbarColourFlowGenHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    BoostedTTbarColourFlowGenHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~BoostedTTbarColourFlowGenHists();
protected:
    JetId btag_tight = CSVBTag(CSVBTag::WP_TIGHT);
private:
    uhh2::Event::Handle<std::vector<GenTopJet>>h_genjets_had;
    uhh2::Event::Handle<bool>h_gensel_2;
};

}
