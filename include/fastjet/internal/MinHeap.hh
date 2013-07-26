//STARTHEADER
// $Id: MinHeap.hh 2577 2011-09-13 15:11:38Z salam $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#ifndef __FASTJET_MINHEAP__HH__
#define __FASTJET_MINHEAP__HH__

#include<vector>
#include<cassert>
#include<memory>
#include<limits>
#include "fastjet/internal/base.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//======================================================================
/// \if internal_doc
/// @ingroup internal
/// \class MinHeap
/// A class which provides a "heap"-like structure that allows
/// access to a the minimal value of a dynamically changing set of numbers
/// \endif
template<typename T>
class MinHeap {
public:
  /// construct a MinHeap from the vector of values, allowing future
  /// expansion to a maximum size max_size;
  MinHeap (const std::vector<T> & values, unsigned int max_size) :
    _heap(max_size) {_initialise(&values.front(), values.size());}

  /// constructor in which the the maximum size is the size of the values array
  MinHeap (const std::vector<T> & values) :
    _heap(values.size()) {_initialise(&values.front(), values.size());}
  

  MinHeap(T const * values, unsigned int size) : _heap(size) {
    _initialise(values,size);
  }

  /// return the location of the minimal value on the heap
  inline unsigned int minloc() const { return _heap[0].minloc;}
  
  /// return the minimal value on the heap
    inline T       minval() const {return _heap[minloc()].value;};

  inline T operator[](int i) const {return _heap[i].value;};

  /// remove the value at the specified location (i.e. replace it with
  /// the largest possible value).
  void remove(unsigned int loc) {
    update(loc,std::numeric_limits<T>::max());};

  /// update the value at the specified location
  void update(unsigned int, T);

private:

  struct ValueLoc{
    T value;
    unsigned int minloc;
  };
      
  std::vector<ValueLoc> _heap;

  void _initialise(T const * values, unsigned int size);

};


//----------------------------------------------------------------------
/// construct the MinHeap; structure will be as follows:
///   . _heap[0].minloc points to globally smallest entry
///     _heap[1].minloc points to smallest entry in one half of heap
///     _heap[2].minloc points to smallest entry in other half of heap
///
///   . for _heap[i], the "parent" is to be found at (i-1)/2
template<typename T>
void MinHeap<T>::_initialise(T const * values, unsigned int size){
  
  // fill the high-range of the heap with the largest possible value
  // (minloc of each entry is itself)
  for (unsigned i = size; i < _heap.size(); i++) {
    _heap[i].value = std::numeric_limits<T>::max();
    _heap[i].minloc = i;
  }

  // fill the rest of the heap with the actual values
  // (minloc of each entry is itself)
  for (unsigned i = 0; i < size; i++) {
    _heap[i].value = values[i];
    _heap[i].minloc = i;
  }
  
  // now adjust the minlocs so that everything is OK...
  for (unsigned i = _heap.size()-1; i > 0; i--) {
    ValueLoc & parent = _heap[(i-1)/2];
    ValueLoc const & here   = _heap[i];
    if (_heap[here.minloc].value < _heap[parent.minloc].value) {
      parent.minloc = here.minloc;
    }
  }
  //cout << minloc() << " "<<sqrt(minval())<<endl;
  //cout << sqrt(_heap[47].value)<<endl;
  //cout << sqrt(_heap[48].value)<<endl;
  //cout << sqrt(_heap[25].value)<<endl;
}


//----------------------------------------------------------------------
template<typename T>
void MinHeap<T>::update(unsigned int loc, T new_value) {
  

  assert(loc < _heap.size());
  auto start = loc;
  auto & vstart = _heap[start];

  // if the minloc is somewhere below us and our value is no smaller
  // than the previous value, we can't possibly change the minloc
  if (vstart.minloc != start && !(new_value < _heap[vstart.minloc].value)) {
    vstart.value = new_value;
    //std::cout << "                     had easy exit\n";
    return;
  }

  // update the value and put a temporary location
  vstart.value = new_value;
  vstart.minloc = start;
  // warn that a change has been made at this place
  bool change_made = true;
  // to make sure we don't go off edge...
  auto heap_end = _heap.size();

  // now work our way up the heap
  while(change_made) {
    auto here = loc;
    auto & vhere = _heap[here];

    // if we were pointing to start, then we must re-initialise things
    if (vhere.minloc == start) {
      vhere.minloc = here; change_made = true;
    }

    // now compare current location to children (at 2*loc+1, 2*loc+2)
    auto child = 2*loc+1;
    {
      auto const & vchild = _heap[child];
      if (child < heap_end && _heap[vchild.minloc].value < _heap[vhere.minloc].value ) {
	vhere.minloc = vchild.minloc;
	change_made = true;
      }
    }
    {
      child++;
      auto const & vchild = _heap[child];
      if (child < heap_end && _heap[vchild.minloc].value < _heap[vhere.minloc].value ) {
	vhere.minloc = vchild.minloc;
	change_made = true;
      }
    }
    // then we move up (loc ->(loc-1)/2) if there's anywhere to go 
    if (loc == 0) {break;}
    loc = (loc-1)/2;
  }

}


FASTJET_END_NAMESPACE

#endif // __FASTJET_MINHEAP__HH__
