/**
 * \file containedpi0.h
 *
 * \ingroup Pi0Ana
 * 
 * \brief Class def header for a class containedpi0
 *
 * @author ryan
 */

/** \addtogroup Pi0Ana

    @{*/

#ifndef LARLITE_CONTAINEDPI0_H
#define LARLITE_CONTAINEDPI0_H

#include "Analysis/ana_base.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mcshower.h"



namespace larlite {
  /**
     \class containedpi0
     User custom analysis class made by SHELL_USER_NAME
   */
  class containedpi0 : public ana_base{
  
  public:

    /// Default constructor
    containedpi0(){ _name="containedpi0"; _fout=0; _containmentcut=0.9; _boxcontain=false;_depcontain = false; _fid=0.0; }

    /// Default destructor
    virtual ~containedpi0(){}

    /** IMPLEMENT in containedpi0.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in containedpi0.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in containedpi0.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();


    void SetContainmentCut(double contain) { _containmentcut = contain ; }
    void SetBoxContain(bool boxcontain) { _boxcontain = boxcontain ; }
    void SetDepContain(bool depcontain) { _depcontain = depcontain ; }
    void SetFidCut(double fid) { _fid = fid ; }
  

  protected:
	double _containmentcut;
	bool _boxcontain;
	bool _depcontain;
	double _fid;
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
