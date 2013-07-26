#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/internal/MinHeap.hh"


#include<cassert>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace opti_details {

  constexpr unsigned short NOWHERE = 62001;

 class OTiledJet {
  public:
   float     eta, phi, kt2, NN_dist;
   unsigned short NN=NOWHERE; 
   unsigned short prev=NOWHERE, next=NOWHERE;
   unsigned short   _jets_index, tile_index=NOWHERE, diJ_posn;
   // routines that are useful in the minheap version of tiled
   // clustering ("misuse" the otherwise unused diJ_posn, so as
   // to indicate whether jets need to have their minheap entries
   // updated).
   inline void label_minheap_update_needed() {diJ_posn = 1;}
   inline void label_minheap_update_done()   {diJ_posn = 0;}
   inline bool minheap_update_needed() const {return diJ_posn==1;}
 };

  /*
  struct OTile {
    /// pointers to neighbouring tiles, including self
    OTile *   begin_tiles[n_tile_neighbours]; 
    /// neighbouring tiles, excluding self
    OTile **  surrounding_tiles; 
    /// half of neighbouring tiles, no self
    OTile **  RH_tiles;  
    /// just beyond end of tiles
    OTile **  end_tiles;
    // index of first and last 
    unsigned short head=NOWHERE, tail=NOWHERE;
    /// sometimes useful to be able to tag a tile
    bool     tagged;    
  };
  */

  template<typename T> T chop(T t) {
     return std::min(T(1.e24),std::max(T(1.e-24),t));
  }

}


//----------------------------------------------------------------------
/// run a tiled clustering, with our minheap for keeping track of the
/// smallest dij
void ClusterSequence::_minheap_optimized_tiled_N2_cluster() {
  using namespace opti_details;

  int n = _jets.size();
  int oriN = n;


  if (n==0) {
    _minheap_faster_tiled_N2_cluster();
    std::cout << "jet done " << _tiles.size() << " " << oriN << " " << _jets.size()
              << " " << _jets[oriN].rap() << "," << _jets[oriN].phi_02pi() << "," << jet_scale_for_algorithm(_jets[oriN])
       	      << " " <<	_jets[_jets.size()-1].rap() << "," << _jets[_jets.size()-1].phi_02pi()  << "," << jet_scale_for_algorithm(_jets[_jets.size()-1])
              << std::endl;
    return;
  }

  _initialise_tiles();
  if (_jets.size()>32000  || _tiles.size() > 62000) return _minheap_faster_tiled_N2_cluster();
  

  OTiledJet briefjets[n];
  
  {
    unsigned int index[n];
    for (int i = 0; i< n; i++) {
      index[i]  = _tile_index(_jets[i].rap(),_jets[i].phi_02pi());
      ++(_tiles[index[i]].nJets);
    }
    

    unsigned int kj[_tiles.size()]; 
    int i=0; int t=0;
    for (auto & tile : _tiles) {
      if(tile.nJets>0){ 
	tile.first=i; i+=tile.nJets;
      }
      kj[t++]=tile.first;
    }
    
    assert(i==n);
    
    
    for (int i = 0; i!=n; ++i) {
      auto k = kj[index[i]];
      auto & j = briefjets[k];
      j.eta = _jets[i].rap();
      j.phi = _jets[i].phi_02pi();
      j.kt2 = chop(jet_scale_for_algorithm(_jets[i]));
      j.NN_dist = _R2;
      j._jets_index=i;
      j.tile_index=index[i];
      if (k!=_tiles[index[i]].first) {
	j.prev = k-1;
	briefjets[k-1].next=k;
      }
      ++kj[index[i]];
    }

    for (unsigned int k=0; k!=_tiles.size(); ++k) assert(kj[k]==_tiles[k].first+_tiles[k].nJets);
    for (int i = 0; i!=n; ++i) assert( briefjets[i].tile_index!=NOWHERE);  
  }


    /*
    auto verify = [&]() {
      for (auto & tile : _tiles) {
	if (0==tile.nJets) assert(tile.first==NOWHERE);
	else {
	  assert(tile.first!=NOWHERE);
	  assert(tile.first < oriN);
	  int k=0;
	  auto p = NOWHERE;
	  for (int iI = tile.first; iI != NOWHERE; iI = briefjets[iI].next) {
	    assert(briefjets[iI].prev==p); p=iI;
	    ++k;
	    if (k<tile.nJets) assert(briefjets[iI].next!=NOWHERE);
	    else assert(briefjets[iI].next==NOWHERE);	    
	  }
	  assert(k==tile.nJets);		    
	}
      }
    };
    */
    // std::cout << "init done " << std::endl;
    // verify();
    // std::cout << " level 1" << std::endl;

  OTiledJet oldB;
  
  // define it locally
  auto bj_diJ = [&](OTiledJet const * const jet)->double  {
    auto kt2 = jet->kt2;
    kt2 = (jet->NN != NOWHERE) ? std::min(kt2,briefjets[jet->NN].kt2) : kt2; 
    return jet->NN_dist * kt2;
  };

  // will be used quite deep inside loops, but declare it here so that
  // memory (de)allocation gets done only once
  vector<int> tile_union(3*n_tile_neighbours);
  
 
  // set up the initial nearest neighbour information
  for (auto const & tile : _tiles) {
    // first do it on this tile
    for (auto jA=tile.first; jA!=tile.first+tile.nJets; ++jA) {
      auto  & jetA = briefjets[jA];
      for (auto jB =tile.first; jB!=jA; ++jB) {
	auto  & jetB = briefjets[jB];
	double dist = _bj_dist(&jetA,&jetB);
	if (dist < jetA.NN_dist) {jetA.NN_dist = dist; jetA.NN = jB;}
	if (dist < jetB.NN_dist) {jetB.NN_dist = dist; jetB.NN = jA;}
      }
    }
    // then do it for RH tiles
    for (Tile ** RTile = tile.RH_tiles; RTile != tile.end_tiles; RTile++) {
      auto  & rtile =  *(*RTile);
      for (int jA=tile.first; jA!=tile.first+tile.nJets; ++jA) {
	auto & jetA = briefjets[jA];
	for (int jB=rtile.first; jB!=rtile.first+rtile.nJets; ++jB) {
	  auto & jetB = briefjets[jB];
	  double dist = _bj_dist(&jetA,&jetB);
	  if (dist < jetA.NN_dist) {jetA.NN_dist = dist; jetA.NN = jB;}
	  if (dist < jetB.NN_dist) {jetB.NN_dist = dist; jetB.NN = jA;}
	}
      }
    }
    // no need to do it for LH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the RH tiles.
  }

  
 

  float diJs[n];
  for (int i = 0; i < n; i++) {
    diJs[i] = bj_diJ(&briefjets[i]);
    briefjets[i].label_minheap_update_done();
  }


  MinHeap<float> minheap(diJs,n);
  // have a stack telling us which jets we'll have to update on the heap
  vector<unsigned short> jets_for_minheap;
  jets_for_minheap.reserve(n); 

  // now run the recombination loop
  int history_location = n-1;

  auto head =  briefjets;

  auto removeFromTile = [&](OTiledJet & jet) {
    --_tiles[jet.tile_index].nJets;
    if(jet.prev==NOWHERE) { // first...
      _tiles[jet.tile_index].first=jet.next;
    } else {
      briefjets[jet.prev].next=jet.next;
    }
    if(jet.next!=NOWHERE) {
      briefjets[jet.next].prev=jet.prev;
    }
    assert(_tiles[jet.tile_index].nJets>=0);
    if (0==_tiles[jet.tile_index].nJets) assert(_tiles[jet.tile_index].first==NOWHERE);
    if (_tiles[jet.tile_index].nJets>0) assert(_tiles[jet.tile_index].first!=NOWHERE);

  };

  while (n > 0) {

    double diJ_min = minheap.minval() *_invR2;
    unsigned short kA = minheap.minloc();
    auto jetA = head + kA;

    assert(jetA->tile_index!=NOWHERE);     

    // do the recombination between A and B
    history_location++;
    auto kB = jetA->NN;
    auto jetB = (kB==NOWHERE) ? nullptr : head + kB;

    if (jetB != nullptr) {

       assert(kB!=NOWHERE);
       assert(kA!=NOWHERE);
	
      // jet-jet recombination
      // If necessary relabel A & B to ensure jetB < jetA, that way if
      // the larger of them == newtail then that ends up being jetA and 
      // the new jet that is added as jetB is inserted in a position that
      // has a future!
      if (kA < kB) {std::swap(jetA,jetB); std::swap(kA,kB);}

      int nn; // new jet index
      _do_ij_recombination_step(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      assert(nn>=oriN);      

      // what was jetB will now become the new jet
      oldB = *jetB;  // take a copy because we will need it...
      assert(oldB.tile_index!=NOWHERE);

     // _bj_remove_from_tiles(jetA);
      //_bj_remove_from_tiles(jetB);
     removeFromTile(*jetA);
     removeFromTile(*jetB);
 

     //      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
                                 // (also registers the jet in the tiling)
     {
       jetB->tile_index=NOWHERE; // just a check
       auto & j = *jetB;
       auto k =  kB;
       j.eta = _jets[nn].rap();
       j.phi = _jets[nn].phi_02pi();
       j.kt2 = chop(jet_scale_for_algorithm(_jets[nn]));
       j.NN_dist = _R2;
       j.NN = NOWHERE;
       j._jets_index=nn;
       j.tile_index=_tile_index(j.eta,j.phi);
       auto & ti = _tiles[j.tile_index];
       j.prev=NOWHERE;
       ++ti.nJets;
       if (ti.nJets>1) {
	 briefjets[ti.first].prev=k;
	 j.next= ti.first;
       } else j.next=NOWHERE;
       ti.first=k;
     }
     assert(jetB->tile_index!=NOWHERE);
     assert(_tiles[jetB->tile_index].first!=NOWHERE);     

    } else {
      // jet-beam recombination
      // std::cout << "jet-beam " << kA << std::endl;
      // get the hist_index
      _do_iB_recombination_step(jetA->_jets_index, diJ_min);
      // _bj_remove_from_tiles(jetA);
      removeFromTile(*jetA);
    }
    
    // remove the minheap entry for jetA
    minheap.remove(kA);
    
    // first establish the set of tiles over which we are going to
    // have to run searches for updated and new nearest-neighbours --
    // basically a combination of vicinity of the tiles of the two old
    // and one new jet.
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union(jetA->tile_index, 
					   tile_union, n_near_tiles);
    if (jetB != nullptr) {
      if (jetB->tile_index != jetA->tile_index) {
	_add_untagged_neighbours_to_tile_union(jetB->tile_index,
					       tile_union,n_near_tiles);
      }
      if (oldB.tile_index != jetA->tile_index && 
	  oldB.tile_index != jetB->tile_index) {

	assert(oldB.tile_index!=NOWHERE);

	// GS: the line below generates a warning that oldB.tile_index
	// may be used uninitialised. However, to reach this point, we
	// ned jetB != NULL (see test a few lines above) and is jetB
	// !=NULL, one would have gone through "oldB = *jetB before
	// (see piece of code ~20 line above), so the index is
	// initialised. We do not do anything to avoid the warning to
	// avoid any potential speed impact.
	_add_untagged_neighbours_to_tile_union(oldB.tile_index,
					       tile_union,n_near_tiles);
      }
      // indicate that we'll have to update jetB in the minheap
      jetB->label_minheap_update_needed();
      jets_for_minheap.push_back(kB);
    }
    
    // verify();

    // Initialise jetB's NN distance as well as updating it for 
    // other particles.
    // Run over all tiles in our union 
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      tile_ptr->tagged = false; // reset tag, since we're done with unions
      // run over all jets in the current tile
      if (tile_ptr->nJets>0)
      for (int iI = tile_ptr->first; iI != NOWHERE; iI = briefjets[iI].next) {
	auto jetI = &briefjets[iI];
	// see if jetI had jetA or jetB as a NN -- if so recalculate the NN
	if (jetI->NN == kA || (jetI->NN == kB && jetB != nullptr)) {
	  jetI->NN_dist = _R2;
	  jetI->NN      = NOWHERE;
	  // label jetI as needing heap action...
	  if (!jetI->minheap_update_needed()) {
	    jetI->label_minheap_update_needed();
	    jets_for_minheap.push_back(iI);}
	  // now go over tiles that are neighbours of I (include own tile)
	  for (Tile ** near_tile  = tile_ptr->begin_tiles; 
	       near_tile != tile_ptr->end_tiles; near_tile++) {
	    // and then over the contents of that tile
	    if ((*near_tile)->nJets>0)
 	    for (auto iJ  = (*near_tile)->first; 
		 iJ != NOWHERE; iJ = briefjets[iJ].next) {
	      auto jetJ = &briefjets[iJ];
	      double dist = _bj_dist(jetI,jetJ);
	      if (dist < jetI->NN_dist && jetJ != jetI) {
		jetI->NN_dist = dist; jetI->NN = iJ;
	      }
	    }
	  }
	}
	// check whether new jetB is closer than jetI's current NN and
	// if jetI is closer than jetB's current (evolving) nearest
	// neighbour. Where relevant update things
	if (jetB != nullptr) {
	  double dist = _bj_dist(jetI,jetB);
	  if (dist < jetI->NN_dist) {
	    if (jetI != jetB) {
	      jetI->NN_dist = dist;
	      jetI->NN = kB;
	      // label jetI as needing heap action...
	      if (!jetI->minheap_update_needed()) {
		jetI->label_minheap_update_needed();
		jets_for_minheap.push_back(iI);}
	    }
	  }
	  if (dist < jetB->NN_dist) {
	    if (jetI != jetB) {
	      jetB->NN_dist = dist;
	      jetB->NN      = iI;}
	  }
	}
      }
    }
    
    // verify();

    // deal with jets whose minheap entry needs updating
    while (jets_for_minheap.size() > 0) {
      auto iI = jets_for_minheap.back(); 
      jets_for_minheap.pop_back();
      minheap.update(iI, bj_diJ(&briefjets[iI]));
      briefjets[iI].label_minheap_update_done();
    }
    n--;
  }
  
  
    std::cout << "jet done " << _tiles.size() << " " << oriN << " " << _jets.size()
              << " " << _jets[oriN].rap() << "," << _jets[oriN].phi_02pi() << "," << jet_scale_for_algorithm(_jets[oriN])
              << " " << _jets[_jets.size()-1].rap() << "," << _jets[_jets.size()-1].phi_02pi() 	<< "," << jet_scale_for_algorithm(_jets[_jets.size()-1])
              << std::endl;
  
}
  
  
FASTJET_END_NAMESPACE
    
