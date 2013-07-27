#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/internal/MinHeap.hh"


#include<cassert>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

#define likely(x) (__builtin_expect(x, true))
#define unlikely(x) (__builtin_expect(x, false))

namespace opti_details {
  
  constexpr unsigned short NOWHERE = 62001;
  
  class OTiledJet {
  public:
    float     eta, phi, kt2=1.e26 /* std::numeric_limits<float>::max()*/, NN_dist=10000.f;
    unsigned short NN=NOWHERE; 
    unsigned short jet_index=NOWHERE, tile_index=NOWHERE;
    bool update=false;
    inline void label_minheap_update_needed() {update=true;}
    inline void label_minheap_update_done()   {update=false;}
    inline bool minheap_update_needed() const {return update;}
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

  unsigned int n = _jets.size();
  unsigned  int oriN = n;


  if (n==0) {
    _minheap_faster_tiled_N2_cluster();
    std::cout << "jet done " << _tiles.size() << " " << oriN << " " << _jets.size()
              << " " << _jets[oriN].rap() << "," << _jets[oriN].phi_02pi() << "," << jet_scale_for_algorithm(_jets[oriN])
       	      << " " <<	_jets[_jets.size()-1].rap() << "," << _jets[_jets.size()-1].phi_02pi()  << "," << jet_scale_for_algorithm(_jets[_jets.size()-1])
              << std::endl;
    return;
  }

  _initialise_tiles();

  //reserve at least one place for tile (maybe two...)
  unsigned int tsize = _tiles.size();
  unsigned int bsize = n+tsize;

  if (n>31000  || bsize > 62000) return _minheap_faster_tiled_N2_cluster();
  

  // will contain max tile size for each tile
  unsigned int mtls[tsize]; 


  OTiledJet briefjets[bsize];
  unsigned short indexNN[2*n]; // redirection of NN s....


  {
    unsigned int index[n];
    for (unsigned int i = 0; i< n; i++) {
      index[i]  = _tile_index(_jets[i].rap(),_jets[i].phi_02pi());
      ++(_tiles[index[i]].nJets);
    }
    
    unsigned int i=0; unsigned int t=0;
    for (auto & tile : _tiles) {
	tile.first=i; i+=(tile.nJets+1); // one more in each tile
	mtls[t++]=tile.first;
    }
    
    assert(i==bsize);
    
    // fil with real jets
    for (unsigned int i = 0; i!=n; ++i) {
      auto k = mtls[index[i]];
      indexNN[i] = k; // point to itlsef...
      auto & j = briefjets[k];
      j.eta = _jets[i].rap();
      j.phi = _jets[i].phi_02pi();
      j.kt2 = chop(jet_scale_for_algorithm(_jets[i]));
      j.NN_dist = _R2;
      j.jet_index=i;
      j.tile_index=index[i];
      ++mtls[index[i]];
    }

    for (unsigned int k=0; k!=_tiles.size(); ++k) assert(mtls[k]==_tiles[k].first+_tiles[k].nJets);
    // for (int i = 0; i!=n; ++i) assert( briefjets[i].tile_index!=NOWHERE);  

    for (unsigned int k=0; k!=tsize; ++k) mtls[k]=_tiles[k].nJets+1;  // this is max size

  }


    
    auto verify = [&]() {
      int ti=0;
      for (auto & tile : _tiles) {
	  assert(tile.first < bsize);
          assert(tile.first+tile.nJets <= bsize);
          assert(tile.nJets <= mtls[ti]);
	  for (int iI = tile.first; iI != tile.first+tile.nJets; ++iI) {
	    assert(briefjets[iI].tile_index==ti);
	  }
         // if (tile.nJets < mtls[ti]) assert(briefjets[tile.first+tile.nJets].tile_index==NOWHERE);
 
          ++ti;
      }
    };
    
    std::cout << "init done " << std::endl;
    verify();
    std::cout << " level 1" << std::endl;

  
  // define it locally
  auto bj_diJ = [&](OTiledJet const * const jet)->float  {
    auto kt2 = jet->kt2;
    kt2 = (jet->NN != NOWHERE) ? std::min(kt2,briefjets[indexNN[jet->NN]].kt2) : kt2; 
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
	auto dist = _bj_dist(&jetA,&jetB);
	if (dist < jetA.NN_dist) {jetA.NN_dist = dist; jetA.NN = jetB.jet_index;}
	if (dist < jetB.NN_dist) {jetB.NN_dist = dist; jetB.NN = jetA.jet_index;}
      }
    }
    // then do it for RH tiles
    for (Tile ** RTile = tile.RH_tiles; RTile != tile.end_tiles; RTile++) {
      auto  & rtile =  *(*RTile);
      for (int jA=tile.first; jA!=tile.first+tile.nJets; ++jA) {
	auto & jetA = briefjets[jA];
	for (int jB=rtile.first; jB!=rtile.first+rtile.nJets; ++jB) {
	  auto & jetB = briefjets[jB];
	  auto dist = _bj_dist(&jetA,&jetB);
	  if (dist < jetA.NN_dist) {jetA.NN_dist = dist; jetA.NN = jetB.jet_index;}
	  if (dist < jetB.NN_dist) {jetB.NN_dist = dist; jetB.NN = jetA.jet_index;}
	}
      }
    }
    // no need to do it for LH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the RH tiles.
  }

  
 

  float diJs[bsize];
  for (unsigned int i = 0; i < bsize; i++) {
    diJs[i] = bj_diJ(&briefjets[i]);
  }


  MinHeap<float> minheap(diJs,bsize);
  // have a stack telling us which jets we'll have to update on the heap
  vector<unsigned short> jets_for_minheap;
  jets_for_minheap.reserve(n); 

  // now run the recombination loop
  int history_location = n-1;

  auto head =  briefjets;

  auto removeFromTile = [&](unsigned short k) {
    auto ti = briefjets[k].tile_index;
    assert(ti<_tiles.size());
    assert(_tiles[ti].nJets>0);
    assert(indexNN[briefjets[k].jet_index] == k);
    // will move last to k...
    --_tiles[ti].nJets;
    auto l = _tiles[ti].first+_tiles[ti].nJets;
    assert(briefjets[l].tile_index==ti);
    assert(k<=l);
    assert(k>=_tiles[ti].first);
    if (l!=k) {
      assert(indexNN[briefjets[l].jet_index] == l);
      briefjets[k]=briefjets[l];
      minheap.update(k,minheap[l]);
      indexNN[briefjets[k].jet_index] = k;
    }
    minheap.remove(l);
    assert(minheap[l]>1.e26);
    briefjets[l].tile_index=NOWHERE;
    briefjets[l].jet_index=NOWHERE;
    if (l!=k) assert(briefjets[k].tile_index==ti);
    assert(_tiles[ti].nJets>=0);
    assert(_tiles[ti].first!=NOWHERE);
  };


  while (n > 0) {

    double diJ_min = minheap.minval() *_invR2;
    unsigned short kA = minheap.minloc();
    auto jetA = head + kA;

    auto tiA = jetA->tile_index;
    auto jiA = jetA->jet_index;

    assert(tiA!=NOWHERE);     
    assert(jiA!=NOWHERE);
     
    assert(_tiles[tiA].nJets>0);
    assert(jiA<_jets.size());


    // do the recombination between A and B
    history_location++;


    auto jiB = jetA->NN;
    auto kB = (jiB!=NOWHERE) ? indexNN[jetA->NN] : NOWHERE;
    auto jetB = (jiB!=NOWHERE) ? head + kB : nullptr;

    unsigned int oldIndex = NOWHERE;

    if likely(jetB != nullptr) {

      assert(jiB==jetB->jet_index);

      assert(kA!=kB);
      assert(kB!=NOWHERE);
      assert(kA!=NOWHERE);
      assert(jetA->jet_index<_jets.size());
      assert(jetB->jet_index<_jets.size());
      assert(_tiles[jetA->tile_index].nJets>0);
      assert(_tiles[jetB->tile_index].nJets>0);


      // jet-jet recombination
      // If necessary relabel A & B to ensure jetB < jetA
      // try to keep active jet close-by

      int nn; // new jet index
      _do_ij_recombination_step(jetA->jet_index, jetB->jet_index, diJ_min, nn);
       assert(nn<int(2*oriN));

      auto eta = _jets[nn].rap();
      auto phi = _jets[nn].phi_02pi();
      auto tin = _tile_index(eta,phi);

      bool inplace=false;
      if (tin==jetB->tile_index) { // in place at kB
	inplace=true;
      } else if (tin==jetA->tile_index) { // in place at kA
	std::swap(jetA,jetB); std::swap(kA,kB); tiA = jetA->tile_index;jiA = jetA->jet_index;
	inplace=true;
      } else {  // in a different tile (if there is space!)
	if (mtls[tin]==_tiles[tin].nJets) {
	  std::cout << "FAILED " << tin << " " << mtls[tin] 
		    << " "  << jetA->tile_index << " "  << jetB->tile_index << std::endl;
	  return _minheap_faster_tiled_N2_cluster();
	}
	if (kA < kB) {
	  std::swap(jetA,jetB); std::swap(kA,kB); tiA = jetA->tile_index;jiA = jetA->jet_index;
	}
      }
      

      // what was jetB will now become the new jet
      oldIndex = jetB->tile_index;  // take a copy because we will need it...
      jiB = jetB->jet_index;
       
    
      assert(oldIndex!=NOWHERE);     
      assert(tiA!=NOWHERE);
      assert(kA!=kB);
      if (tiA==oldIndex) assert(inplace); 

     // _bj_remove_from_tiles(jetA);
      //_bj_remove_from_tiles(jetB);

     assert(tiA == jetA->tile_index);
     assert(kA+head == jetA);
     assert(_tiles[tiA].nJets>0);

     if (tiA==oldIndex)  assert(_tiles[tiA].nJets>1);

     removeFromTile(kA);

     assert(_tiles[oldIndex].nJets>0);
     
     // recompute kb...
     kB = indexNN[jiB];
     jetB = head+kB;

     if (!inplace) removeFromTile(kB);
 

     //      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
                                 // (also registers the jet in the tiling)
     {
       if (!inplace) {
         assert(mtls[tin]>_tiles[tin].nJets);
	 kB = _tiles[tin].first+_tiles[tin].nJets;
	 ++_tiles[tin].nJets;
	 jetB = head+kB;
	 assert(jetB->tile_index==NOWHERE);
	 assert(jetB->jet_index==NOWHERE);
       }

       auto & j = *jetB;
       j.eta = eta;
       j.phi = phi;
       j.kt2 = chop(jet_scale_for_algorithm(_jets[nn]));
       j.NN_dist = _R2;
       j.NN = NOWHERE;
       j.jet_index=nn;
       j.tile_index=tin;
       indexNN[nn]=kB;

     }

      assert(jetB->tile_index==tin);
      assert(_tiles[jetB->tile_index].nJets>0);
  
    } else {
      assert(kA!=NOWHERE);
      assert(jetA->jet_index<_jets.size());


      // jet-beam recombination
      // std::cout << "jet-beam " << kA << std::endl;
      // get the hist_index
      _do_iB_recombination_step(jetA->jet_index, diJ_min);
      // _bj_remove_from_tiles(jetA);
      removeFromTile(kA);
    }
    
    // at this point jetA DOES NOT EXISTS ANYMORE!

    // first establish the set of tiles over which we are going to
    // have to run searches for updated and new nearest-neighbours --
    // basically a combination of vicinity of the tiles of the two old
    // and one new jet.
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union(tiA, 
					   tile_union, n_near_tiles);
    if likely(jetB != nullptr) {
      if (jetB->tile_index != tiA) {
	_add_untagged_neighbours_to_tile_union(jetB->tile_index,
					       tile_union,n_near_tiles);
      }
      if (oldIndex != tiA && 
	  oldIndex != jetB->tile_index) {
	_add_untagged_neighbours_to_tile_union(oldIndex,
					       tile_union,n_near_tiles);
      }
      // indicate that we'll have to update jetB in the minheap
      jetB->label_minheap_update_needed();
      jets_for_minheap.push_back(kB);
    }
    
    verify();

    // Initialise jetB's NN distance as well as updating it for 
    // other particles.
    // Run over all tiles in our union 
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      tile_ptr->tagged = false; // reset tag, since we're done with unions
      // run over all jets in the current tile
      // if likely(tile_ptr->nJets>0)
      for (auto iI = tile_ptr->first; iI !=tile_ptr->first+tile_ptr->nJets; ++iI) {
	auto jetI = &briefjets[iI];
	// see if jetI had jetA or jetB as a NN -- if so recalculate the NN
	if unlikely(jetI->NN == jiA || (jetI->NN == jiB && jetB != nullptr)) {
	  jetI->NN_dist = _R2;
	  jetI->NN      = NOWHERE;
	  // label jetI as needing heap action...
	  if (!jetI->minheap_update_needed()) {
	    jetI->label_minheap_update_needed();
	    jets_for_minheap.push_back(iI);
	  }
	  // now go over tiles that are neighbours of I (include own tile)
	  for (Tile ** near_tile  = tile_ptr->begin_tiles; 
	       near_tile != tile_ptr->end_tiles; near_tile++) {
	    // and then over the contents of that tile
	    // if likely((*near_tile)->nJets>0)
 	    for (auto iJ  = (*near_tile)->first; 
		 iJ != (*near_tile)->first+(*near_tile)->nJets; ++iJ) {
	      auto jetJ = &briefjets[iJ];
	      auto dist = _bj_dist(jetI,jetJ);
	      if (dist < jetI->NN_dist && jetJ != jetI) {
		jetI->NN_dist = dist; jetI->NN =jetJ->jet_index;
	      }
	    }
	  }
	} // end JetI NN recomputation
	// check whether new jetB is closer than jetI's current NN and
	// if jetI is closer than jetB's current (evolving) nearest
	// neighbour. Where relevant update things
	if likely(jetB != nullptr) {
	  auto dist = _bj_dist(jetI,jetB);
	  if unlikely(dist < jetI->NN_dist) {
	    if likely(jetI != jetB) {
	      jetI->NN_dist = dist;
	      jetI->NN = jetB->jet_index;
	      // label jetI as needing heap action...
	      if (!jetI->minheap_update_needed()) {
		jetI->label_minheap_update_needed();
		jets_for_minheap.push_back(iI);
	      }
	    }
	  }
	  if (dist < jetB->NN_dist) {
	    if likely(jetI != jetB) {
	      jetB->NN_dist = dist;
	      jetB->NN      = jetI->jet_index;}
	  }
	}
      }
    }
    
    verify();

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
    
