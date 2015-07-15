#ifndef LARLITE_CONTAINEDPI0_CXX
#define LARLITE_CONTAINEDPI0_CXX

#include "containedpi0.h"

namespace larlite {

  bool containedpi0::initialize() {

    return true;
  }
  
  bool containedpi0::analyze(storage_manager* storage) {
        // Bring in the info for the event
        auto mctruth = storage->get_data<event_mctruth>("generator");
            if(!mctruth) {
                print(larlite::msg::kERROR,__FUNCTION__,Form("Did not find specified data product, mctruth!"));
                return false;}// if no mctruth
        auto mcpart = mctruth->at(0).GetParticles();
	// For now this will work on single particle gun... Need to make it work with single gun or bnb_data
        auto mcshower = storage->get_data<event_mcshower>("mcreco");

	
	//If DepContain
	if(_depcontain){
	// For the pi0 to be contained we need to have both photons contained up to a certain percent 
	for(auto const& mcs : *mcshower){
		auto SP = mcs.Start();
		auto ShowerDetProf =  mcs.DetProfile();
		double mccontained = ShowerDetProf.E()/SP.E();
		if(mccontained<_containmentcut) return false;
		}
	}//if Depcontain

	//If BoxContain
	if(_boxcontain){
        // ::geoalgo::AABox  volume(0+_fid,-116.5+_fid,0+_fid,256.35-_fid,116.5-_fid,1036.8-_fid);
	// is this contained? 
	for(auto & p:mcpart) if(p.PdgCode()==111){	
		auto traj = p.Trajectory();
		if(traj[0].X()>256.35-_fid || traj[0].X()<0+_fid) return false;
		if(traj[0].Y()>116.5-_fid || traj[0].Y()<-116.5+_fid) return false;
		if(traj[0].Z()>1036.8-_fid || traj[0].Y()<-0+_fid) return false;}
	}

 
    return true;
  }

  bool containedpi0::finalize() {

    return true;
  }

}
#endif
