#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/internal/MinHeap.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh



//----------------------------------------------------------------------
/// run a tiled clustering, with our minheap for keeping track of the
/// smallest dij
void ClusterSequence::_minheap_optimized_tiled_N2_cluster() {


  _initialise_tiles();

  int n = _jets.size();
  int oirN = n;
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB;
  

  // will be used quite deep inside loops, but declare it here so that
  // memory (de)allocation gets done only once
  vector<int> tile_union(3*n_tile_neighbours);
  
  // initialise the basic jet info 
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    //cout << i<<": "<<jetA->tile_index<<"\n";
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * head = briefjets; // a nicer way of naming start

  // set up the initial nearest neighbour information
  vector<Tile>::const_iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    // first do it on this tile
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    // then do it for RH tiles
    for (Tile ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
      for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
	for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
	  double dist = _bj_dist(jetA,jetB);
	  if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	  if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
	}
      }
    }
    // no need to do it for LH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the RH tiles.
  }

  
  //// now create the diJ (where J is i's NN) table -- remember that 
  //// we differ from standard normalisation here by a factor of R2
  //// (corrected for at the end). 
  //struct diJ_plus_link {
  //  double     diJ; // the distance
  //  TiledJet * jet; // the jet (i) for which we've found this distance
  //                  // (whose NN will the J).
  //};
  //diJ_plus_link * diJ = new diJ_plus_link[n];
  //jetA = head;
  //for (int i = 0; i < n; i++) {
  //  diJ[i].diJ = _bj_diJ(jetA); // kt distance * R^2
  //  diJ[i].jet = jetA;  // our compact diJ table will not be in	     
  //  jetA->diJ_posn = i; // one-to-one corresp. with non-compact jets,
  //                      // so set up bi-directional correspondence here.
  //  jetA++; // have jetA follow i 
  //}

  vector<double> diJs(n);
  for (int i = 0; i < n; i++) {
    diJs[i] = _bj_diJ(&briefjets[i]);
    briefjets[i].label_minheap_update_done();
  }
  MinHeap minheap(diJs);
  // have a stack telling us which jets we'll have to update on the heap
  vector<TiledJet *> jets_for_minheap;
  jets_for_minheap.reserve(n); 

  // now run the recombination loop
  int history_location = n-1;
  while (n > 0) {

    double diJ_min = minheap.minval() *_invR2;
    jetA = head + minheap.minloc();

    // do the recombination between A and B
    history_location++;
    jetB = jetA->NN;

    if (jetB != NULL) {
      // jet-jet recombination
      // If necessary relabel A & B to ensure jetB < jetA, that way if
      // the larger of them == newtail then that ends up being jetA and 
      // the new jet that is added as jetB is inserted in a position that
      // has a future!
      if (jetA < jetB) {std::swap(jetA,jetB);}

      int nn; // new jet index
      _do_ij_recombination_step(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      
      // what was jetB will now become the new jet
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
                                 // (also registers the jet in the tiling)
    } else {
      // jet-beam recombination
      // get the hist_index
      _do_iB_recombination_step(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }

    // remove the minheap entry for jetA
    minheap.remove(jetA-head);

    // first establish the set of tiles over which we are going to
    // have to run searches for updated and new nearest-neighbours --
    // basically a combination of vicinity of the tiles of the two old
    // and one new jet.
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union(jetA->tile_index, 
					   tile_union, n_near_tiles);
    if (jetB != NULL) {
      if (jetB->tile_index != jetA->tile_index) {
	_add_untagged_neighbours_to_tile_union(jetB->tile_index,
					       tile_union,n_near_tiles);
      }
      if (oldB.tile_index != jetA->tile_index && 
	  oldB.tile_index != jetB->tile_index) {
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
      jets_for_minheap.push_back(jetB);
    }


    // Initialise jetB's NN distance as well as updating it for 
    // other particles.
    // Run over all tiles in our union 
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      tile_ptr->tagged = false; // reset tag, since we're done with unions
      // run over all jets in the current tile
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
	// see if jetI had jetA or jetB as a NN -- if so recalculate the NN
	if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
	  jetI->NN_dist = _R2;
	  jetI->NN      = NULL;
	  // label jetI as needing heap action...
	  if (!jetI->minheap_update_needed()) {
	    jetI->label_minheap_update_needed();
	    jets_for_minheap.push_back(jetI);}
	  // now go over tiles that are neighbours of I (include own tile)
	  for (Tile ** near_tile  = tile_ptr->begin_tiles; 
	               near_tile != tile_ptr->end_tiles; near_tile++) {
	    // and then over the contents of that tile
	    for (TiledJet * jetJ  = (*near_tile)->head; 
                            jetJ != NULL; jetJ = jetJ->next) {
	      double dist = _bj_dist(jetI,jetJ);
	      if (dist < jetI->NN_dist && jetJ != jetI) {
		jetI->NN_dist = dist; jetI->NN = jetJ;
	      }
	    }
	  }
	}
	// check whether new jetB is closer than jetI's current NN and
	// if jetI is closer than jetB's current (evolving) nearest
	// neighbour. Where relevant update things
	if (jetB != NULL) {
	  double dist = _bj_dist(jetI,jetB);
	  if (dist < jetI->NN_dist) {
	    if (jetI != jetB) {
	      jetI->NN_dist = dist;
	      jetI->NN = jetB;
	      // label jetI as needing heap action...
	      if (!jetI->minheap_update_needed()) {
		jetI->label_minheap_update_needed();
		jets_for_minheap.push_back(jetI);}
	    }
	  }
	  if (dist < jetB->NN_dist) {
	    if (jetI != jetB) {
	      jetB->NN_dist = dist;
	      jetB->NN      = jetI;}
	  }
	}
      }
    }

    // deal with jets whose minheap entry needs updating
    while (jets_for_minheap.size() > 0) {
      TiledJet * jetI = jets_for_minheap.back(); 
      jets_for_minheap.pop_back();
      minheap.update(jetI-head, _bj_diJ(jetI));
      jetI->label_minheap_update_done();
    }
    n--;
  }

  // final cleaning up;
  delete[] briefjets;

  std::cout << "jet done " << oirN << " " << _jets.size() << std::endl;


}


FASTJET_END_NAMESPACE

