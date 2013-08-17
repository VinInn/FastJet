#include<iostream>
#include<vector>
#include<cmath>
#include<algorithm>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/internal/MinHeap.hh"
#include<cstdlib>
#include<cstring>

#include<cassert>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

#define likely(x) (__builtin_expect(x, true))
#define unlikely(x) (__builtin_expect(x, false))

namespace opti_details {
  
  using FixPoint = signed short;
  constexpr float fp2fix_f = 4096.f;
  constexpr float fix2fp_f = 1./4096.;
  constexpr double fp2fix_d = 4096.;
  constexpr double fix2fp_d = 1./4096.;

  constexpr FixPoint inffix = 32700;
  constexpr FixPoint maxfix = 32000;
  constexpr double maxeta_d = 7.2;
  constexpr float maxeta_f = 7.2;


  // cannot be constexpr!
  inline FixPoint fp2fix(float x) {  return std::rint(x*fp2fix_f);}
  inline FixPoint fp2fix(double x) { return std::rint(x*fp2fix_d);}

  inline FixPoint fp2fixS(float x) {  return fp2fix(std::copysign(std::min(maxeta_f,std::abs(x)),x));}
  inline FixPoint fp2fixS(double x) { return fp2fix(std::copysign(std::min(maxeta_d,std::abs(x)),x));}


  constexpr float fix2f(FixPoint x) { return x*fix2fp_f;}

  constexpr unsigned short NOWHERE = 62001;
  constexpr float pif=M_PI;
  constexpr float twopif=2*pif;

  // constexpr FixPoint pifix= fp2fix(M_PI);
  // constexpr FixPoint twopifix = fp2fix(2*pif);
  constexpr FixPoint pifix=  M_PI*fp2fix_d + 0.5;
  constexpr FixPoint twopifix = 2*M_PI*fp2fix_d + 0.5;
  

  
  class ProtoJet {
  public:
    explicit ProtoJet(int noloc) : NN(noloc), jet_index(noloc) {}
    float kt2=1.e26; 
    FixPoint  eta=maxfix, phi=maxfix; 
    FixPoint  NN_dist=inffix;
    unsigned short NN=62005; 
    unsigned short jet_index=62005, tile_index=NOWHERE;
  };

  inline
  ProtoJet* minit(unsigned int sz,  unsigned int noloc) {
      ProtoJet* jets= ((ProtoJet*)(::malloc(sz*sizeof(ProtoJet))) );
      unsigned char * cj = (unsigned char *)(jets);
      auto last = cj+sz*sizeof(ProtoJet);
      auto half = cj + (last-cj)/2;
      jets[0]=ProtoJet(noloc);
      auto start = cj+sizeof(ProtoJet);
      auto bs = start-cj;
      while (start<half) {
        ::memcpy(start,cj,bs);
        start+=bs;
        bs = start-cj;
      };
      assert(last-start<=start-cj);
      ::memcpy(start,cj,last-start);
      return jets;
  }
  
  // AOS
  class PrJetsAOS {
  public:
    using Itype = unsigned short;
    
    PrJetsAOS(unsigned int sz,  unsigned int noloc) : 
      jets(minit(sz,noloc)), s(sz) {
    }
    
    ~PrJetsAOS() {
      free(jets);
    }

    unsigned int size() const { return s;}
    float kt2(int i) const { return jets[i].kt2;}
    float & kt2(int i) { return jets[i].kt2;}
    FixPoint eta(int i) const { return jets[i].eta;}
    FixPoint & eta(int i) { return jets[i].eta;}
    FixPoint phi(int i) const { return jets[i].phi;}
    FixPoint & phi(int i) { return jets[i].phi;}
    FixPoint NN_dist(int i) const { return jets[i].NN_dist;}
    FixPoint & NN_dist(int i) { return jets[i].NN_dist;}
    Itype NN(int i) const { return jets[i].NN;}
    Itype & NN(int i) { return jets[i].NN;}
    Itype jet_index(int i) const { return jets[i].jet_index;}
    Itype & jet_index(int i) { return jets[i].jet_index;}
    Itype tile_index(int i) const { return jets[i].tile_index;}
    Itype & tile_index(int i) { return jets[i].tile_index;}
    
    
    ProtoJet const & back() const { return jets[s-1];}

    
    FixPoint dist(unsigned int i, unsigned int j) const {
      auto dphi = std::abs(phi(i) - phi(j));
      auto deta = eta(i) - eta(j);
      dphi =  (dphi > pifix) ? twopifix - dphi : dphi;
      // return dphi*dphi + deta*deta;
      return (i==j) ? inffix : dphi*dphi + deta*deta;
    }
    
  
    // valid if we are sure dphi < pi
    FixPoint safeDist(unsigned int i, unsigned int j) const {
      auto dphi = phi(i) - phi(j);
      auto deta = eta(i) - eta(j);
      return (i==j) ? inffix : dphi*dphi + deta*deta;
    }
    
    
    // valid if we are sure dphi < pi and i!=j
    FixPoint safeDist1(unsigned int i, unsigned int j) const {
      auto dphi = phi(i) - phi(j);
      auto deta = eta(i) - eta(j);
      return dphi*dphi + deta*deta;
    }
    
    
    // valid if we are sure dphi > pi
    FixPoint safeDist2(unsigned int i, unsigned int j) const {
      auto dphi = twopifix - std::abs(phi(i) - phi(j));
      auto deta = eta(i) - eta(j);
      // can never be i==j
      return dphi*dphi + deta*deta;
    }
    
    
    void move(unsigned int j, unsigned int i) {
       ::memcpy(&jets[i],&jets[j],sizeof(ProtoJet));
       ::memcpy(&jets[j],&back(),sizeof(ProtoJet));
      // jets[i]=jets[j]; jets[j]=back();
    }
    
    void reset (unsigned int i) {
      ::memcpy(&jets[i],&back(),sizeof(ProtoJet));
      // jets[i]=back();
    }
    
    
    void copy (PrJetsAOS const & o, unsigned int j, unsigned int i) {
      ::memcpy(&jets[i],&o.jets[j],sizeof(ProtoJet));
      // jets[i] = o.jets[j];
    }
    
    void copy (unsigned int j, unsigned int i) {
      copy(*this, j, i);
    }
    
    void swap (PrJetsAOS & o) {
      std::swap(jets,o.jets);
      std::swap(s,o.s);
    }
  
  private:
    
    ProtoJet * jets;
    unsigned int s;
  };
  
  
  // SOA
  class PrJetsSOA {
  public:
    using Itype = unsigned short;
    
    PrJetsSOA(unsigned int sz,  unsigned int noloc) : 
      veta(sz,100.f), vphi(sz,100.f), vkt2(sz,1.e26), vNN_dist(sz,10000.f),
      vNN(sz,noloc), vjet_index(sz,noloc), vtile_index(sz,NOWHERE),
       s(sz) {}
    
    unsigned int size() const { return s;}
    float eta(int i) const { return veta[i];}
    float & eta(int i) { return veta[i];}
    float phi(int i) const { return vphi[i];}
    float & phi(int i) { return vphi[i];}
    float kt2(int i) const { return vkt2[i];}
    float & kt2(int i) { return vkt2[i];}
    float NN_dist(int i) const { return vNN_dist[i];}
    float & NN_dist(int i) { return vNN_dist[i];}
    Itype NN(int i) const { return vNN[i];}
    Itype & NN(int i) { return vNN[i];}
    Itype jet_index(int i) const { return vjet_index[i];}
    Itype & jet_index(int i) { return vjet_index[i];}
    Itype tile_index(int i) const { return vtile_index[i];}
    Itype & tile_index(int i) { return vtile_index[i];}
    
    
    
    float dist(unsigned int i, unsigned int j) const {
      auto dphi = std::abs(phi(i) - phi(j));
      auto deta = eta(i) - eta(j);
      dphi =  (dphi > pif) ? twopif - dphi : dphi;
      // return dphi*dphi + deta*deta;
      return (i==j) ? 10000.f : dphi*dphi + deta*deta;
    }
    
    
    // valid if we are sure dphi < pi
    float safeDist(unsigned int i, unsigned int j) const {
      auto dphi = phi(i) - phi(j);
      auto deta = eta(i) - eta(j);
      return (i==j) ? 10000.f : dphi*dphi + deta*deta;
    }
    
    
    // valid if we are sure dphi < pi and i!=j
    float safeDist1(unsigned int i, unsigned int j) const {
      auto dphi = phi(i) - phi(j);
      auto deta = eta(i) - eta(j);
      return dphi*dphi + deta*deta;
    }
    
    
    // valid if we are sure dphi > pi
    float safeDist2(unsigned int i, unsigned int j) const {
      auto dphi = twopif - std::abs(phi(i) - phi(j));
      auto deta = eta(i) - eta(j);
      // can never be i==j
      return dphi*dphi + deta*deta;
    }
    
    
    void move(unsigned int j, unsigned int i) {
      eta(i)=eta(j); eta(j)=100.f;
      phi(i)=phi(j); phi(j)=100.f;
      kt2(i)=kt2(j); kt2(j)=1.e26;
      NN_dist(i)=NN_dist(j);
      NN(i)=NN(j); NN(j)=vNN.back();
      jet_index(i)= jet_index(j); jet_index(j)=vjet_index.back();
      tile_index(i)= tile_index(j); tile_index(j)=vtile_index.back();
    }
    
    void reset (unsigned int j) {
      eta(j)=100.f;
      phi(j)=100.f;
      kt2(j)=1.e26;
      NN_dist(j)=10000.f;
      NN(j)=vNN.back();
      jet_index(j)=vjet_index.back();
      tile_index(j)=vtile_index.back();
    }

    
    void copy (PrJetsSOA const & o, unsigned int j, unsigned int i) {
      eta(i)=o.eta(j);
      phi(i)=o.phi(j);
      kt2(i)=o.kt2(j); 
      NN_dist(i)=o.NN_dist(j);
      NN(i)=o.NN(j);
      jet_index(i)= o.jet_index(j);
      tile_index(i)= o.tile_index(j); 
    }
  
    void copy (unsigned int j, unsigned int i) {
      copy(*this, j, i);
    }
    
    void swap (PrJetsSOA & o) {
      veta.swap(o.veta);
      vphi.swap(o.vphi);
      vkt2.swap(o.vkt2); 
      vNN_dist.swap(o.vNN_dist);
      vNN.swap(o.vNN);
      vjet_index.swap(o.vjet_index);
      vtile_index.swap(o.vtile_index); 
    }
    
private:
    
    std::vector<float> veta, vphi, vkt2, vNN_dist;
    std::vector<unsigned short> vNN; 
    std::vector<unsigned short> vjet_index,  vtile_index;
    unsigned int s;
  };
  
  using  PrJets = PrJetsAOS;
  
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

      first.resize(sz,0);
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
    /*
    auto half = oriN + (_jets.size()-oriN)/2;
    std::cout << "jet done " << _tiles.size() << " " << oriN << " " << _jets.size()
              << " " << _jets[oriN].rap() << "," << _jets[oriN].phi_02pi() << "," << jet_scale_for_algorithm(_jets[oriN])
       	      << " " <<	_jets[half].rap() << "," << _jets[half].phi_02pi()  << "," << jet_scale_for_algorithm(_jets[half])
       	      << " " <<	_jets[_jets.size()-1].rap() << "," << _jets[_jets.size()-1].phi_02pi()  << "," << jet_scale_for_algorithm(_jets[_jets.size()-1])
              << std::endl;
    */
    return;
  }

  _initialise_tiles();

  // too few tiles to optize deltaphi...
  if (_n_tiles_phi<5)   return _minheap_faster_tiled_N2_cluster();

  auto R2fix = fp2fix(_R2);


  // OTiles tiles(_tiles_eta_min,_tiles_eta_max, _tiles_ieta_max-_tiles_ieta_min+1, _n_tiles_phi);
  //   apparenlty deta is < _R2 for the above...
  OTiles tiles(_tiles_eta_min,_tiles_eta_max, std::max(2,_tiles_ieta_max-_tiles_ieta_min), _n_tiles_phi);

  // std::cout << "tiles " << _tiles_eta_min << "," << _tiles_eta_max << " " << tiles.size() << " " << tiles.rsize << std::endl;


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
      tiles.first[t]=i; i+=std::max(2,int(tiles.last[t]+1)); // one more in each tile or at least 4
      mtls[t]=tiles.first[t]; ++t;
    } 
    //skip the two eta gards
    tiles.first[t]=i; t++;
    tiles.first[t]=i; t++;
  }
  
  assert((t-1)==tiles.size()-tiles.rsize);
  
  unsigned int bsize = i;
  
  // too big for 16bits?
  if (n>31000  || bsize > 62000) return _minheap_faster_tiled_N2_cluster();
  
  
  int NOLOC=2*n;
  PrJets protojets(bsize+1,NOLOC);
  unsigned short indexNN[2*n+1]; // redirection of NN s....
  indexNN[NOLOC]=bsize;  // dummy place for vectorization
  
  // fil with real jets
  auto & j = protojets;
  for (unsigned int i = 0; i!=n; ++i) {
    auto k = mtls[index[i]];
    indexNN[i] = k; // point to itlsef...
    j.eta(k) = fp2fixS(_jets[i].rap());
    j.phi(k) = fp2fix(_jets[i].phi_02pi());
    j.kt2(k) = chop(jet_scale_for_algorithm(_jets[i]));
    j.NN_dist(k) = R2fix;
    j.jet_index(k)=i;
    j.tile_index(k)=index[i];
    ++mtls[index[i]];
  }
  
  // for (unsigned int k=0; k!=_tiles.size(); ++k) assert(mtls[k]==_tiles[k].first+_tiles[k].nJets);
  // for (int i = 0; i!=n; ++i) assert( protojets[i].tile_index!=NOWHERE);  
  
  
  // fix indexing
  for (unsigned int k=0; k!=tsize; ++k) tiles.last[k]+=tiles.first[k];
  // fill phi gards
  for (unsigned int k=0; k!=tiles.rsize; ++k) { 
    tiles.first[k] = tiles.first[tiles.tailN+k-1]; 
    tiles.last[k] =  tiles.last[tiles.tailN+k-1];
  }
  for (unsigned int k=0; k!=tiles.rsize; ++k) { 
    tiles.first[tiles.tailN+tiles.rsize+k-1] = tiles.first[tiles.head+k-1]; 
    tiles.last[tiles.tailN+tiles.rsize+k-1]  = tiles.last[tiles.head+k-1];
  }
  
  for (unsigned int k=0; k!=tsize; ++k) mtls[k]=std::max(tiles.first[k]+2,tiles.last[k]+1);  // this is max size

  
  /*
  auto verify = [&]() {
    assert(protojets.kt2[indexNN[NOLOC]]>_R2);
    assert(protojets.jet_index.back()==NOLOC);
    assert(protojets.NN.back()==NOLOC);

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
	assert(protojets.tile_index[iI]==tt);
        assert(indexNN[protojets.jet_index[iI]]==iI);
        assert(protojets.NN_dist[iI]< 10000.f);
      }
    }
  };
  
  // std::cout << "init done " << std::endl;
  verify();
  std::cout << " level 1" << std::endl;
  */
  
  // define it locally
  auto bj_diJ = [&](unsigned int jet)->float  {
    // NOLOC exits and point to a very large kt2! 
    return fix2f(protojets.NN_dist(jet))*std::min(protojets.kt2(jet),protojets.kt2(indexNN[protojets.NN(jet)])); 
    // auto kt2 = jet->kt2;
    // kt2 = (jet->NN != NOLOC) ? std::min(kt2,protojets[indexNN[jet->NN]].kt2) : kt2; 
    // return jet->NN_dist * kt2;
  };

 
 
  // set up the initial nearest neighbour information
  for (auto k = tiles.head; k!=tiles.tail; ++k) {
    if (tiles.first[k]==tiles.last[k]) continue;  // empty...
    // first do it on this tile
    for (auto jA=tiles.first[k]; jA!=tiles.last[k]; ++jA) {
      for (auto jB =tiles.first[k]; jB!=jA; ++jB) {
	auto dist = protojets.safeDist1(jA,jB);
	if (dist < protojets.NN_dist(jA)) {protojets.NN_dist(jA) = dist; protojets.NN(jA) = protojets.jet_index(jB);}
	if (dist < protojets.NN_dist(jB)) {protojets.NN_dist(jB) = dist; protojets.NN(jB) = protojets.jet_index(jA);}
      }
    }
    // then do it for LH tiles (3 on top + one left
    if (k<tiles.head2) // dphi>pi
      for (auto kk = (k-tiles.rsize)-1; kk!=(k-tiles.rsize)+2; ++kk) {
	for (auto jA=tiles.first[k]; jA!=tiles.last[k]; ++jA) {
	  for (auto jB=tiles.first[kk]; jB!=tiles.last[kk]; ++jB) {
	    auto dist = protojets.safeDist2(jA,jB);
	    if (dist < protojets.NN_dist(jA)) {protojets.NN_dist(jA) = dist; protojets.NN(jA) = protojets.jet_index(jB);}
	    if (dist < protojets.NN_dist(jB)) {protojets.NN_dist(jB) = dist; protojets.NN(jB) = protojets.jet_index(jA);}
	  }
	}
      }
    else
      for (auto kk = (k-tiles.rsize)-1; kk!=(k-tiles.rsize)+2; ++kk) {
	for (auto jA=tiles.first[k]; jA!=tiles.last[k]; ++jA) {
	  for (auto jB=tiles.first[kk]; jB!=tiles.last[kk]; ++jB) {
	    auto dist = protojets.safeDist1(jA,jB);
	    if (dist < protojets.NN_dist(jA)) {protojets.NN_dist(jA) = dist; protojets.NN(jA) = protojets.jet_index(jB);}
	    if (dist < protojets.NN_dist(jB)) {protojets.NN_dist(jB) = dist; protojets.NN(jB) = protojets.jet_index(jA);}
	  }
	}
      }

    auto kk = k-1;
    for (auto jA=tiles.first[k]; jA!=tiles.last[k]; ++jA) {
      for (auto jB=tiles.first[kk]; jB!=tiles.last[kk]; ++jB) {
	auto dist = protojets.safeDist1(jA,jB);
	if (dist < protojets.NN_dist(jA)) {protojets.NN_dist(jA) = dist; protojets.NN(jA) = protojets.jet_index(jB);}
	if (dist < protojets.NN_dist(jB)) {protojets.NN_dist(jB) = dist; protojets.NN(jB) = protojets.jet_index(jA);}
      }
    }
    // no need to do it for RH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the LH tiles.
  }

  
 

  float diJs[bsize];
  for (unsigned int i = 0; i != bsize; i++) {
    diJs[i] = bj_diJ(i);
  }


  MinHeap<float> minheap(diJs,bsize);
  // have a stack telling us which jets we'll have to update on the heap
  vector<unsigned short> jets_for_minheap;
  jets_for_minheap.reserve(n); 

  // now run the recombination loop
  int history_location = n-1;

  auto removeFromTile = [&](unsigned short k) {
    auto ti = protojets.tile_index(k);
    // assert(ti>=tiles.head && ti<tiles.tail);
    // assert(tiles.last[ti]>tiles.first[ti]);
    // will move last to k...
    --tiles.last[ti];
    auto l =tiles.last[ti];
    // fix gards
    if (ti<tiles.head2) --tiles.last[ti+tiles.goff]; 
    if (ti>=tiles.tailN) --tiles.last[ti-tiles.goff]; 
    // assert(protojets[l].tile_index==ti);
    // assert(k<=l);
    // assert(k>=tiles.first[ti]);
    if (l!=k) {
      // assert(indexNN[protojets.jet_index[l]] == l);
      protojets.copy(l,k);
      minheap.update(k,minheap[l]);
      indexNN[protojets.jet_index(k)] = k;
    } 
    protojets.reset(l);
    minheap.remove(l);
    //assert(minheap[l]>1.e26);
    // if (l!=k) assert(protojets.tile_index[k]==ti);
  };


  auto compactify  = [&]() {
    // compactify the array
    PrJets newProtojets(bsize+1,NOLOC);  // can be made smaller...

    unsigned int i=0; unsigned int t=tiles.head;
    for (unsigned int ip=0; ip!=tiles.nPhi; ++ip) {
      for (unsigned int ie=0; ie!=tiles.nEta; ++ie) {
	auto sz = tiles.last[t]-tiles.first[t];
	auto fo =tiles.first[t];
	tiles.first[t]=i; i+=std::max(2,int(sz+1)); // one more in each tile or at least 4
	tiles.last[t] = tiles.first[t]+sz; 
	mtls[t]=i;
	// copy
	auto ki=tiles.first[t];
	for (auto k=fo; k!=fo+sz; ++k) {
	  newProtojets.copy(protojets,k,ki);
	  indexNN[newProtojets.jet_index(ki)]=ki;
	  ki++;
	}
	// zero
	// for (auto k=tiles.last[t]; k!=i; ++k) protojets[k]=ProtoJet();
	++t;
      } 
      //skip the two eta gards
      tiles.last[t]=tiles.first[t]=i; t++;
      tiles.last[t]=tiles.first[t]=i; t++;
    }
    assert(i<=bsize);
    assert((t-1)==tiles.size()-tiles.rsize);
    newProtojets.swap(protojets);
    bsize=i;
    indexNN[NOLOC]=bsize;

    // fill phi gards
    for (unsigned int k=0; k!=tiles.rsize; ++k) { 
      tiles.first[k] = tiles.first[tiles.tailN+k-1]; 
      tiles.last[k] =  tiles.last[tiles.tailN+k-1];
    }
    for (unsigned int k=0; k!=tiles.rsize; ++k) { 
      tiles.first[tiles.tailN+tiles.rsize+k-1] = tiles.first[tiles.head+k-1]; 
      tiles.last[tiles.tailN+tiles.rsize+k-1]  = tiles.last[tiles.head+k-1];
    }
    
    // rebuild heap
    float diJ[i];
    for (unsigned int k = 0; k != i; ++k) {
      diJ[k] = bj_diJ(k);
    }
    minheap = MinHeap<float>(diJ,i);
    // done ???
  };
  
  
  auto nOld = n;
  constexpr int nMin=256; // tsize???

  // verify();
  // std::cout << " level 2" << std::endl;
 


  while (n > 0) {

    auto diJ_min = minheap.minval() *_invR2;
    unsigned short kA = minheap.minloc();

    auto tiA = protojets.tile_index(kA);
    auto jiA = protojets.jet_index(kA);


    //assert(kA!=NOWHERE);
    //assert(tiA!=NOWHERE);     
    //assert(jiA!=NOWHERE);
    //assert(jiA<_jets.size());


    // do the recombination between A and B
    history_location++;


    auto jiB = protojets.NN(kA);
    bool paired = jiB!=NOLOC;
    auto kB = paired ? indexNN[jiB] : NOWHERE;
    if (!paired) ++jiB; // trick so that jiB is not NOLOC....
  
    unsigned int otiB = paired ? protojets.tile_index(kB) : NOWHERE;

    // tiles where modified jets lies
    int ct=1;
    unsigned int ttc[3];
    ttc[0]=tiA;
    if (paired && otiB != tiA) ttc[ct++] = otiB;



    if likely(paired) {
	
      //assert(kB!=NOWHERE);
      //assert(otiB!=NOWHERE);
      //assert(jiB!=NOWHERE);
      //assert(jiB!=(NOWHERE+1));
      //assert(jiB<_jets.size());

      

      // jet-jet recombination
      // If necessary relabel A & B depending where the new jet will endup

      int nn; // new jet index
      _do_ij_recombination_step(jiA, jiB, diJ_min, nn);
       assert(nn<int(2*oriN));

      auto eta = _jets[nn].rap();
      auto phi = _jets[nn].phi_02pi();
      auto tin = tiles.index(eta,phi);

      bool inplace=false;
      if (tin==otiB) { // in place at kB
	inplace=true;
      } else if (tin==tiA) { // in place at kA
	std::swap(kA,kB); tiA = protojets.tile_index(kA); jiA = protojets.jet_index(kA);  jiB = protojets.jet_index(kB); 
	inplace=true;
      } else {  // in a different tile (if there is space!)
	ttc[ct++] = tin;  // a new tile!
	
	if (mtls[tin]==tiles.last[tin]) {
	  std::cout << "FAILED " << tin << " " << mtls[tin] 
		    << " "  << protojets.tile_index(kA) << " "  << protojets.tile_index(kB) << std::endl;
	  // will need to re-adjust everything
	  compactify();
	  kA = indexNN[jiA];
	}
      }
      
      //assert(ct<3);
      //assert(kA!=kB);
      //assert(kA!=NOWHERE);
      //assert(tiA!=NOWHERE);
      //assert(jiA!=NOWHERE);
      //assert(jiA!=(NOWHERE+1));
      //assert(jiA<_jets.size());


      // jetA will be removed and
      // what was jetB will now become the new jet  

       
     removeFromTile(kA);
     // recompute kb...
     kB = indexNN[jiB];

 
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

      }

       auto & j = protojets;
       j.eta(kB) = fp2fixS(eta);
       j.phi(kB) = fp2fix(phi);
       j.kt2(kB) = chop(jet_scale_for_algorithm(_jets[nn]));
       j.NN_dist(kB) = R2fix;
       j.NN(kB) = NOLOC;
       j.jet_index(kB)=nn;
       j.tile_index(kB)=tin;
       indexNN[nn]=kB;

     }

     // indicate that we'll have to update jetB in the minheap
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
	    bool mod=false;
            bool isB = iI==kB;
	    //  if (protojets.tile_index[iI]>tsize) std::cout << "??? k " << kk << " " << iI << " " << tiles.first[kk] << " " << tiles.last[kk] << " " << protojets.jet_index[iI] << " " << protojets.NN[iI] << std::endl;
            // if (protojets.jet_index[iI]>_jets.size()) std::cout << "??? j " << kk << " " << iI << " " << tiles.first[kk] << " " << tiles.last[kk] << " " << protojets.jet_index[iI] << " " << protojets.NN[iI] << std::endl;
	    // see if jetI had jetA or jetB as a NN -- if so recalculate the NN  (jiB is eihter ok or NOLOC+1)
	    if unlikely(protojets.NN(iI) == jiA || (protojets.NN(iI) == jiB)) {
		// now go over tiles that are neighbours of I (include own tile)
		auto ndist=R2fix;
		auto nind = NOLOC;
                // kk maybe a guard row!
		auto irow = protojets.tile_index(iI) -tiles.rsize-1;
		for (int ir=0;ir!=3;++ir) {  // rows
		  if (irow<tiles.rsize || irow>tiles.tail) // deltaphi>pi 
		    for (auto ii = irow; ii!=irow+3; ++ii) { // columns
		      for (auto iJ = tiles.first[ii]; iJ !=tiles.last[ii]; ++iJ) { 
			auto dist = protojets.safeDist2(iI,iJ);
			nind =  (dist<ndist) ?  protojets.jet_index(iJ) : nind;
			ndist = (dist<ndist) ? dist : ndist;
		      }
		    }
		  else if(ir!=1)
		    for (auto ii = irow; ii!=irow+3; ++ii) { // columns
		      for (auto iJ = tiles.first[ii]; iJ !=tiles.last[ii]; ++iJ) { 
			auto dist = protojets.safeDist1(iI,iJ);
			nind =  (dist<ndist) ?  protojets.jet_index(iJ) : nind;
			ndist = (dist<ndist) ? dist : ndist;
		      }
		    }
		  else  // central row
		    for (auto ii = irow; ii!=irow+3; ++ii) { // columns
		      for (auto iJ = tiles.first[ii]; iJ !=tiles.last[ii]; ++iJ) { 
			auto dist = protojets.safeDist(iI,iJ);  // iJ can be == iI ...
			nind =  (dist<ndist) ?  protojets.jet_index(iJ) : nind;
			ndist = (dist<ndist) ? dist : ndist;
		      }
		    }
		  irow+=tiles.rsize;
		}  

		  // if (tiles.last[irow+2]-tiles.first[irow]>bsize) std::cout << "irow " << irow << " " << tiles.first[irow] << " " << tiles.last[irow+2] << std::endl;
		  // if (tiles.last[irow+2]<tiles.first[irow]) std::cout << "irow " << irow << " " << tiles.first[irow] << " " << tiles.last[irow+2] << std::endl;
		  
		  //for (auto iJ = tiles.first[ii]; iJ !=tiles.last[ii]; ++iJ) {  // vectorization challange (when minloc pattern will vectorize!): they are all contiguous in a row. not all valid
		  //for (auto iJ = tiles.first[irow]; iJ !=tiles.last[irow+2]; ++iJ) {
		  // auto dist = protojets.dist(iI,iJ);
		  //	nind =  (dist<ndist) ?  protojets.jet_index(iJ) : nind;
		  //	ndist = (dist<ndist) ? dist : ndist;
		  // if (dist < protojets.NN_dist[iI] && iJ != iI) {  // FIXME we should find a way to get dist to itself infinite!
		  // if (protojets.jet_index[iI]>_jets.size()) std::cout << "??? d " << dist << " " << protojets.NN_dist[iI] << " " << protojets.jet_index[iJ] << std::endl;
		  // protojets.NN_dist[iI] = dist; protojets.NN[iI] = protojets.jet_index[iJ];
		  //}


		protojets.NN_dist(iI) = ndist;
		protojets.NN(iI)   = nind;
		// label jetI as needing heap action...
                mod=true;		
	      } // end JetI NN recomputation
	    
	    // check whether new jetB is closer than jetI's current NN and
	    // if jetI is closer than jetB's current (evolving) nearest
	    // neighbour. Where relevant update things
	    if likely(paired &!isB) {
		auto dist = protojets.dist(iI,kB);
		if unlikely(dist < protojets.NN_dist(iI)) {
		    protojets.NN_dist(iI) = dist;
		    protojets.NN(iI) = protojets.jet_index(kB);
		    // label jetI as needing heap action...
		    mod = true;
		  }
		if (dist < protojets.NN_dist(kB)) {
		  protojets.NN_dist(kB)=dist;;
		  protojets.NN(kB) = protojets.jet_index(iI);
		}
	      }  // end jetB update

	    if (mod) jets_for_minheap.push_back(iI);	    

	  } // end iI loop
	} // end kk loop
	row+=tiles.rsize;
      } // end "row" loop
    }   // end tiles loop
	  
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
      minheap.update(iI, bj_diJ(iI));
    }
    n--;
    

    if (n>nMin && n < nOld/2) {
      // std::cout << "compact " << n << std::endl;
      nOld = n;
      compactify();
      // verify();
    }
  }
  /*
  /// backward compatible printout....
  auto half = oriN + (_jets.size()-oriN)/2;
  std::cout << "jet done " << _tiles.size() << " " << oriN << " " << _jets.size()
	    << " " << _jets[oriN].rap() << "," << _jets[oriN].phi_02pi() << "," << jet_scale_for_algorithm(_jets[oriN])
	    << " " <<	_jets[half].rap() << "," << _jets[half].phi_02pi()  << "," << jet_scale_for_algorithm(_jets[half])
 	    << " " << _jets[_jets.size()-1].rap() << "," << _jets[_jets.size()-1].phi_02pi() 	<< "," << jet_scale_for_algorithm(_jets[_jets.size()-1])
	    << std::endl;
  */
}
  
  
FASTJET_END_NAMESPACE
    
