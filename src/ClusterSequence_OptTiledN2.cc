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

  struct OTiles {
    using Itype = unsigned int; 
    std::vector<Itype> first;
    std::vector<Itype> last;
    std::vector<bool> tag;

    float etaMin, etaMax;
    float deta, dphi;  // inverse!
    unsigned int nEta, nPhi;  // real one  eta is row, phi is column...
    unsigned int head, tail;
    unsigned int head2, tailN;  // last of second row, first of last aw
    unsigned int rsize;  // raw size
    unsigned int goff; // distance to the corresponding phi gard

    unsigned int size() const { return first.size();}

    OTiles(float ietaMin, float ietaMax, 
	   unsigned int inEta,  unsigned int inPhi) :
      etaMin(ietaMin), etaMax(ietaMax),
      nEta(inEta), nPhi(inPhi) {
      auto sz = (nEta+2)*(nPhi+2);
      rsize = nEta+2;

      first.resize(sz,NOWHERE);
      last.resize(sz,0);
      tag.resize(sz,false);

      goff = rsize*nPhi;
      head = rsize+1; tail = sz-rsize-1;
      head2 = head+nEta;
      tailN=tail-nEta;
      deta = float(nEta)/(etaMax-etaMin);
      dphi = float(nPhi)/float(twopi);
    }

    unsigned int index(float eta, float phi) const {
      int ieta = 1 + std::max(0,std::min(int(nEta-1),int(floor((eta-etaMin)*deta))));
      // phi assumed between 0 and 2phi...)
      // assert(phi>=0 && phi<=float(twopi));
      int iphi = 1 + std::min(int(nPhi-1),int(phi*dphi));
      return iphi*rsize + ieta;
    }

  };


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


  // OTiles tiles(_tiles_eta_min,_tiles_eta_max, _tiles_ieta_max-_tiles_ieta_min+1, _n_tiles_phi);
  //   apparenlty deta is < _R2 for the above...
  OTiles tiles(_tiles_eta_min,_tiles_eta_max, std::max(2,_tiles_ieta_max-_tiles_ieta_min), _n_tiles_phi);



  //reserve at least one place for tile (maybe two...)
  unsigned int tsize = tiles.size();

 

  // will contain max tile size for each tile
  unsigned int mtls[tsize]; 

  
  unsigned int index[n];
  for (unsigned int i = 0; i< n; i++) {
    index[i]  = tiles.index(_jets[i].rap(),_jets[i].phi_02pi());
    { auto k = index[i];
      assert(k>=tiles.head); assert(k<tiles.tail);
      k = k%tiles.rsize; assert(k!=0); assert(k!=(tiles.rsize-1));
    }
    ++tiles.last[index[i]];  // for the time being is size...
  }
    
  unsigned int i=0; unsigned int t=tiles.head;
  for (unsigned int ip=0; ip!=tiles.nPhi; ++ip) {
    for (unsigned int ie=0; ie!=tiles.nEta; ++ie) {
      tiles.first[t]=i; i+=std::max(4,int(tiles.last[t]+1)); // one more in each tile or at least 4
      mtls[t]=tiles.first[t]; ++t;
    } 
    t+=2;  //skip the two eta gards
  }
  
  assert((t-1)==tiles.size()-tiles.rsize);
  
  unsigned int bsize = i;
  
  // too big for 16bits?
  if (n>31000  || bsize > 62000) return _minheap_faster_tiled_N2_cluster();
  
  
  OTiledJet briefjets[bsize];
  unsigned short indexNN[2*n]; // redirection of NN s....
  
  
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
  
  // for (unsigned int k=0; k!=_tiles.size(); ++k) assert(mtls[k]==_tiles[k].first+_tiles[k].nJets);
  // for (int i = 0; i!=n; ++i) assert( briefjets[i].tile_index!=NOWHERE);  
  
  
  // fix indexing
  for (unsigned int k=0; k!=tsize; ++k) tiles.last[k]+=tiles.first[k];
  for (unsigned int k=0; k!=tsize; ++k) mtls[k]=std::max(tiles.first[k]+4,tiles.last[k]+1);  // this is max size
  // fill phi gards
  for (unsigned int k=0; k!=tiles.nEta; ++k) { 
    tiles.first[1+k] = tiles.first[tiles.tailN+k]; 
    tiles.last[1+k] =  tiles.last[tiles.tailN+k];
  }
  for (unsigned int k=0; k!=tiles.nEta; ++k) { 
    tiles.first[tiles.tailN+tiles.rsize+k] = tiles.first[tiles.head+k]; 
    tiles.last[tiles.tailN+tiles.rsize+k]  = tiles.last[tiles.head+k];
  }
  
  
  
  /*
  auto verify = [&]() {
    for (unsigned int ti=0; ti<tiles.size(); ++ti) {
      assert(!tiles.tag[ti]);
      auto k = ti%tiles.rsize;
      if (k==0 || k==(tiles.rsize-1)) assert(tiles.first[ti]==tiles.last[ti]);
      assert(tiles.first[ti]<=tiles.last[ti]);
      if (k>0 && k<(tiles.rsize-1))  {
	assert(tiles.first[ti]< bsize);
	assert(tiles.last[ti] <= bsize);
	assert(tiles.last[ti] <= mtls[ti]);
      }
      for (int iI = tiles.first[ti]; iI != tiles.last[ti]; ++iI) {
        auto tt = (ti>=tiles.head) ? ti : ti+tiles.goff;
        tt	= (tt<tiles.tail) ? tt : ti-tiles.goff;
	assert(briefjets[iI].tile_index==tt);
        assert(indexNN[briefjets[iI].jet_index]==iI);
      }
    }
  };
  */
  // std::cout << "init done " << std::endl;
  //verify();
  // std::cout << " level 1" << std::endl;
  
  
  // define it locally
  auto bj_diJ = [&](OTiledJet const * const jet)->float  {
    auto kt2 = jet->kt2;
    kt2 = (jet->NN != NOWHERE) ? std::min(kt2,briefjets[indexNN[jet->NN]].kt2) : kt2; 
    return jet->NN_dist * kt2;
  };

 
 
  // set up the initial nearest neighbour information
  for (auto k = tiles.head; k!=tiles.tail; ++k) {
    if (tiles.first[k]==tiles.last[k]) continue;  // empty...
    // first do it on this tile
    for (auto jA=tiles.first[k]; jA!=tiles.last[k]; ++jA) {
      auto  & jetA = briefjets[jA];
      for (auto jB =tiles.first[k]; jB!=jA; ++jB) {
	auto  & jetB = briefjets[jB];
	auto dist = _bj_dist(&jetA,&jetB);
	if (dist < jetA.NN_dist) {jetA.NN_dist = dist; jetA.NN = jetB.jet_index;}
	if (dist < jetB.NN_dist) {jetB.NN_dist = dist; jetB.NN = jetA.jet_index;}
      }
    }
    // then do it for LH tiles (3 on top + one left
    for (auto kk = (k-tiles.rsize)-1; kk!=(k-tiles.rsize)+2; ++kk) {
      for (auto jA=tiles.first[k]; jA!=tiles.last[k]; ++jA) {
	auto & jetA = briefjets[jA];
	for (auto jB=tiles.first[kk]; jB!=tiles.last[kk]; ++jB) {
	  auto & jetB = briefjets[jB];
	  auto dist = _bj_dist(&jetA,&jetB);
	  if (dist < jetA.NN_dist) {jetA.NN_dist = dist; jetA.NN = jetB.jet_index;}
	  if (dist < jetB.NN_dist) {jetB.NN_dist = dist; jetB.NN = jetA.jet_index;}
	}
      }
    }
    auto kk = k-1;
    for (auto jA=tiles.first[k]; jA!=tiles.last[k]; ++jA) {
      auto & jetA = briefjets[jA];
      for (auto jB=tiles.first[kk]; jB!=tiles.last[kk]; ++jB) {
	auto & jetB = briefjets[jB];
	auto dist = _bj_dist(&jetA,&jetB);
	if (dist < jetA.NN_dist) {jetA.NN_dist = dist; jetA.NN = jetB.jet_index;}
	if (dist < jetB.NN_dist) {jetB.NN_dist = dist; jetB.NN = jetA.jet_index;}
      }
    }
    // no need to do it for RH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the LH tiles.
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
    // assert(ti>=tiles.head && ti<tiles.tail);
    // assert(tiles.last[ti]>tiles.first[ti]);
    // will move last to k...
    --tiles.last[ti];
    auto l =tiles.last[ti];
    // fix gards
    if (ti<tiles.head2) --tiles.last[ti+tiles.goff]; 
    if (ti>=tiles.tailN) --tiles.last[ti-tiles.goff]; 
    // assert(briefjets[l].tile_index==ti);
    // assert(k<=l);
    // assert(k>=tiles.first[ti]);
    if (l!=k) {
      //assert(indexNN[briefjets[l].jet_index] == l);
      briefjets[k]=briefjets[l];
      minheap.update(k,minheap[l]);
      indexNN[briefjets[k].jet_index] = k;
    }
    minheap.remove(l);
    //assert(minheap[l]>1.e26);
    briefjets[l].tile_index=NOWHERE;
    briefjets[l].jet_index=NOWHERE;
    // if (l!=k) assert(briefjets[k].tile_index==ti);
  };


  while (n > 0) {

    auto diJ_min = minheap.minval() *_invR2;
    unsigned short kA = minheap.minloc();
    auto jetA = head + kA;

    auto tiA = jetA->tile_index;
    auto jiA = jetA->jet_index;

    // assert(tiA!=NOWHERE);     
    // assert(jiA!=NOWHERE);
     
    // assert(_tiles[tiA].nJets>0);
    // assert(jiA<_jets.size());


    // do the recombination between A and B
    history_location++;


    auto jiB = jetA->NN;
    bool paired = jiB!=NOWHERE;
    if (!paired) ++jiB; // trick so that jiB is not NOWHERE....
    auto kB = paired ? indexNN[jetA->NN] : NOWHERE;
    auto jetB = paired ? head + kB : nullptr;

    unsigned int oldIndex = paired ? jetB->tile_index : NOWHERE;

    // tiles where modified jets lies
    int ct=1;
    unsigned int ttc[3];
    ttc[0]=tiA;
    if (oldIndex =! tiA) ttc[ct++] = oldIndex;



    if likely(paired) {
	

      // jet-jet recombination
      // If necessary relabel A & B depending where the new jet will endup

      int nn; // new jet index
      _do_ij_recombination_step(jiA, jiB, diJ_min, nn);
       assert(nn<int(2*oriN));

      auto eta = _jets[nn].rap();
      auto phi = _jets[nn].phi_02pi();
      auto tin = tiles.index(eta,phi);

      bool inplace=false;
      if (tin==oldIndex) { // in place at kB
	inplace=true;
      } else if (tin==tiA) { // in place at kA
	std::swap(jetA,jetB); std::swap(kA,kB); tiA = jetA->tile_index;jiA = jetA->jet_index;
	inplace=true;
      } else {  // in a different tile (if there is space!)
	ttc[ct++] = tin;  // a new tile!

	if (mtls[tin]==tiles.last[tin]) {
	  std::cout << "FAILED " << tin << " " << mtls[tin] 
		    << " "  << jetA->tile_index << " "  << jetB->tile_index << std::endl;
	  // will need to re-adjust everything
	  // FIXME this is wrong still will may work...  for test
	  inplace=true; // in place at kB
	}

      }
      

      // jetA will be removed and
      // what was jetB will now become the new jet  

      jiB = jetB->jet_index;  // save old jet index
       
     removeFromTile(kA);
     // recompute kb...
     kB = indexNN[jiB];
     jetB = head+kB;

 
     if (!inplace) removeFromTile(kB);
 

     //      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
                                 // (also registers the jet in the tiling)
     {
       if (!inplace) {
         assert(mtls[tin]>tiles.last[tin]);
	 kB = tiles.last[tin];
	 ++tiles.last[tin];
	 // fix gards
	 if (tin<tiles.head2) ++tiles.last[tin+tiles.goff]; 
	 if (tin>=tiles.tailN) ++tiles.last[tin-tiles.goff]; 

	 jetB = head+kB;
	 // assert(jetB->tile_index==NOWHERE);
	 // assert(jetB->jet_index==NOWHERE);
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

     // indicate that we'll have to update jetB in the minheap
     jetB->label_minheap_update_needed();
     jets_for_minheap.push_back(kB);
   
    } else {
 
     // jet-beam recombination
      // std::cout << "jet-beam " << kA << std::endl;
      // get the hist_index
      _do_iB_recombination_step(jiA, diJ_min);
      // _bj_remove_from_tiles(jetA);
      removeFromTile(kA);
    }
    
    // at this point jetA DOES NOT EXISTS ANYMORE!

    // verify();


    // Initialise jetB's NN distance as well as updating it for 
    // other particles.
    // Run over all tiles in our union 
    

    for (int it=0; it!=ct; ++it) {
      auto row=  ttc[it] - tiles.rsize-1;
      for (int r=0;r!=3;++r) {  // rows
	for (auto kk = row; kk!=row+3; ++kk) { // columns
	  if (tiles.tag[kk]) continue;
	  tiles.tag[kk]=true;
	  if (kk<tiles.head2) tiles.tag[kk+tiles.goff]=true; 
	  if (kk>=tiles.tailN) tiles.tag[kk-tiles.goff]=true; 

	  // run over all jets in the current tile
	  // if likely(tile_ptr->nJets>0)
	  for (auto iI = tiles.first[kk]; iI !=tiles.last[kk]; ++iI) {
	    auto jetI = &briefjets[iI];
	    // see if jetI had jetA or jetB as a NN -- if so recalculate the NN  (jiB is eihter ok or NOWHERE+1)
	    if unlikely(jetI->NN == jiA || (jetI->NN == jiB)) {
		jetI->NN_dist = _R2;
		jetI->NN      = NOWHERE;
		// label jetI as needing heap action...
		if (!jetI->minheap_update_needed()) {
		  jetI->label_minheap_update_needed();
		  jets_for_minheap.push_back(iI);
		}
		// now go over tiles that are neighbours of I (include own tile)
                // kk maybe a guard row!
		auto irow = (jetI->tile_index) -tiles.rsize-1;
		for (int ir=0;ir!=3;++ir) {  // rows
		  for (auto ii = irow; ii!=irow+3; ++ii) { // columns
		    for (auto iJ = tiles.first[ii]; iJ !=tiles.last[ii]; ++iJ) {  // vectorization challange: they are all contiguous in a row. not all valid
		      auto jetJ = &briefjets[iJ];
		      auto dist = _bj_dist(jetI,jetJ);
		      if (dist < jetI->NN_dist && jetJ != jetI) {  // FIXME we should find a way to get dist to itself infinite!
			jetI->NN_dist = dist; jetI->NN =jetJ->jet_index;
		      }
		    }
		  }
		  irow+=tiles.rsize;
		}	    
	      } // end JetI NN recomputation
	    
	    // check whether new jetB is closer than jetI's current NN and
	    // if jetI is closer than jetB's current (evolving) nearest
	    // neighbour. Where relevant update things
	    if likely(paired) {
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
		      jetB->NN      = jetI->jet_index;
		    }
		}
	      }  // end jetB update
	    
	  }
	}
	row+=tiles.rsize;
      }
    }
    // clean tag bit
    for (int it=0; it!=ct; ++it) {
      auto row=  ttc[it] - tiles.rsize-1;
      for (int r=0;r!=3;++r) {
	for (auto kk = row; kk!=row+3; ++kk) {
	  tiles.tag[kk]=false;
	  if (kk<tiles.head2) tiles.tag[kk+tiles.goff]=false; 
	  if (kk>=tiles.tailN) tiles.tag[kk-tiles.goff]=false; 

	}
	row+=tiles.rsize;
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
  
    /// backward compatible printout....
    std::cout << "jet done " << _tiles.size() << " " << oriN << " " << _jets.size()
              << " " << _jets[oriN].rap() << "," << _jets[oriN].phi_02pi() << "," << jet_scale_for_algorithm(_jets[oriN])
              << " " << _jets[_jets.size()-1].rap() << "," << _jets[_jets.size()-1].phi_02pi() 	<< "," << jet_scale_for_algorithm(_jets[_jets.size()-1])
              << std::endl;
  
}
  
  
FASTJET_END_NAMESPACE
    
