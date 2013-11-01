#ifndef __FASTJET__VORONOI_H__
#define __FASTJET__VORONOI_H__

//STARTHEADER
// $Id: Voronoi.hh 2686 2011-11-14 09:28:22Z soyez $
//
// Copyright (c) 1994 by AT&T Bell Laboratories (see below)
//
//
//----------------------------------------------------------------------
// This file is included as part of FastJet but was mostly written by
// S. Fortune in C, put into C++ with memory management by S
// O'Sullivan, and with further interface and memeory management
// modifications by Gregory Soyez.
//
// Permission to use, copy, modify, and distribute this software for
// any purpose without fee is hereby granted, provided that this
// entire notice is included in all copies of any software which is or
// includes a copy or modification of this software and in all copies
// of the supporting documentation for such software. THIS SOFTWARE IS
// BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED WARRANTY.
// IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY REPRESENTATION
// OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY OF THIS
// SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
//
//----------------------------------------------------------------------
//ENDHEADER


/*
* The author of this software is Steven Fortune.  
* Copyright (c) 1994 by AT&T Bell Laboratories.
* Permission to use, copy, modify, and distribute this software for any
* purpose without fee is hereby granted, provided that this entire notice
* is included in all copies of any software which is or includes a copy
* or modification of this software and in all copies of the supporting
* documentation for such software.
* THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
* WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
* REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
* OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
*/

/* 
* This code was originally written by Stephan Fortune in C code.  I,
* Shane O'Sullivan, have since modified it, encapsulating it in a C++
* class and, fixing memory leaks and adding accessors to the Voronoi
* Edges.  Permission to use, copy, modify, and distribute this
* software for any purpose without fee is hereby granted, provided
* that this entire notice is included in all copies of any software
* which is or includes a copy or modification of this software and in
* all copies of the supporting documentation for such software.  THIS
* SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
* WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
* REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE
* MERCHANTABILITY OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR
* PURPOSE.
*/

/*
 * This code, included in the FastJet distribution, was originally
 * written by Stephan Fortune in C and adapted to C++ by Shane
 * O'Sullivan under the terms repported above.
 *
 * Below are the list of changes implemented by the FastJet authors:
 *
 * 2011-11-14  Gregory Soyez  <soyez@fastjet.fr>
 * 
 *      * removed 'plot' and 'triangulate' (were always 0)
 *      * removed unused plot functions (openpl, circle, range, 
 *        out_bisector, out_ep, out_vertex, out_site, out_triple)
 *      * removed unused 'VPoint p' in 'intersect'
 * 
 * 
 * 2011-07-22  Gregory Soyez  <soyez@fastjet.fr>
 * 
 *      * replaced Point by VPoint (to avoid any potential conflict
 *        with an already existing class Point in FastJet
 * 
 * 
 * 2008-04-01  Gregory Soyez  <soyez@fastjet.fr>
 * 
 *      * declared ystar volatile in HalfEdge (apparently fixes a bug
 *        related to VD computations with points on a grid)
 * 
 * 
 * 2007-05-07  Gregory Soyez  <soyez@fastjet.fr>
 * 
 *      * put the code in the fastjet namespace
 * 
 *      * replaced float by double
 * 
 *      * generateVoronoi() takes a vector of Point instead of 2
 *        pointers
 * 
 *      * added info about the parent sites to GraphEdge
 * 
 *      * removed condition on minimal distance between sites
 * 
 */

#include "fastjet/LimitedWarning.hh"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DELETED -2
#define le 0
#define re 1

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/**
 * \if internal_doc
 * @ingroup internal
 * \class VPoint
 * class to handle a 2d point
 * \endif
 */
class VPoint{
public:
  /// defailt ctor
  VPoint() : x(0.0), y(0.0) {}

  /// ctor with initialisation
  VPoint(double _x, double _y) : x(_x), y(_y) {}

  /// addition
  inline VPoint operator + (const VPoint &p) const{
    return VPoint(x+p.x, y+p.y);
  }

  /// subtraction
  inline VPoint operator - (const VPoint &p) const{
    return VPoint(x-p.x, y-p.y);
  }

  /// scalar multiplication
  inline VPoint operator * (const double t) const{
    return VPoint(x*t, y*t);
  }

  /// vector coordinates
  double x,y;
};


/// norm of a vector
inline double norm(const VPoint p){
  return p.x*p.x+p.y*p.y;
}


/// 2D vector product
inline double vector_product(const VPoint &p1, const VPoint &p2){
  return p1.x*p2.y-p1.y*p2.x;
}


/// scalar product
inline double scalar_product(const VPoint &p1, const VPoint &p2){
  return p1.x*p2.x+p1.y*p2.y;
}


/**
 * \if internal_doc
 * @ingroup internal
 * \class GraphEdge
 * handle an edge of the Voronoi Diagram.
 * \endif
 */
class GraphEdge{
public:
  /// coordinates of the extreme points
  double x1,y1,x2,y2;

  /// indices of the parent sites that define the edge
  int point1, point2;

  /// pointer to the next edge
  GraphEdge* next;
};


/**
 * \if internal_doc
 * @ingroup internal
 * \class Site
 * structure used both for particle sites and for vertices.
 * \endif
 */
class Site{
 public:
  VPoint	coord;
  int sitenbr;
  int refcnt;
};



class Freenode{
public:
  Freenode *nextfree;
};


class FreeNodeArrayList{
public:
  Freenode* memory;
  FreeNodeArrayList* next;
};


class Freelist{
public:
  Freenode *head;
  int nodesize;
};

class Edge{
public:
  double a,b,c;
  Site *ep[2];
  Site *reg[2];
  int edgenbr;
};


class Halfedge{
public:
  Halfedge *ELleft, *ELright;
  Edge *ELedge;
  int ELrefcnt;
  char ELpm;
  Site *vertex;
  volatile double ystar;
  Halfedge *PQnext;
};

/**
 * \if internal_doc
 * @ingroup internal
 * \class VoronoiDiagramGenerator
 * Shane O'Sullivan C++ version of Stephan Fortune Voronoi diagram
 * generator
 * \endif
 */
class VoronoiDiagramGenerator{
public:
  inline VoronoiDiagramGenerator();
  inline ~VoronoiDiagramGenerator();

  inline bool generateVoronoi(std::vector<VPoint> *_parent_sites,
		       double minX, double maxX, double minY, double maxY, 
		       double minDist=0);

  inline void resetIterator(){
    iteratorEdges = allEdges;
  }

  bool getNext(GraphEdge **e){
    if(iteratorEdges == 0)
      return false;
    
    *e = iteratorEdges;
    iteratorEdges = iteratorEdges->next;
    return true;
  }
  
  std::vector<VPoint> *parent_sites;
  int n_parent_sites;

private:
  inline void cleanup();
  inline void cleanupEdges();
  inline char *getfree(Freelist *fl);	
  inline Halfedge *PQfind();
  inline int PQempty();
	
  Halfedge **ELhash;

  inline Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
  inline Halfedge *HEcreate(Edge *e,int pm);
  
  inline VPoint PQ_min();
  inline Halfedge *PQextractmin();	
  inline void freeinit(Freelist *fl,int size);
  inline void makefree(Freenode *curr,Freelist *fl);
  inline void geominit();
  inline void plotinit();

  // GS: removed the unused (always ==0) argument
  bool voronoi(/*int triangulate*/);

  inline void ref(Site *v);
  inline void deref(Site *v);
  inline void endpoint(Edge *e,int lr,Site * s);

  inline void ELdelete(Halfedge *he);
  inline Halfedge *ELleftbnd(VPoint *p);
  inline Halfedge *ELright(Halfedge *he);
  inline void makevertex(Site *v);
  
  inline void PQinsert(Halfedge *he,Site * v, double offset);
  inline void PQdelete(Halfedge *he);
  inline bool ELinitialize();
  inline void ELinsert(Halfedge *lb, Halfedge *newHe);
  inline Halfedge * ELgethash(int b);
  inline Halfedge *ELleft(Halfedge *he);
  inline Site *leftreg(Halfedge *he);
  inline bool PQinitialize();
  inline int PQbucket(Halfedge *he);
  inline void clip_line(Edge *e);
  inline char *myalloc(unsigned n);
  inline int right_of(Halfedge *el,VPoint *p);

  inline Site *rightreg(Halfedge *he);
  inline Edge *bisect(Site *s1, Site *s2);
  inline double dist(Site *s,Site *t);

  // GS: 'p' is unused and always ==0 (see also comment by
  //     S. O'Sullivan in the source file), so we remove it
  inline Site *intersect(Halfedge *el1, Halfedge *el2 /*, VPoint *p=0*/);

  inline Site *nextone();

  inline void pushGraphEdge(double x1, double y1, double x2, double y2, 
		     Site *s1, Site *s2);

  // Gregory Soyez: unused plotting methods
  // void openpl();
  // void circle(double x, double y, double radius);
  // void range(double minX, double minY, double maxX, double maxY);
  // 
  // void out_bisector(Edge *e);
  // void out_ep(Edge *e);
  // void out_vertex(Site *v);
  // void out_site(Site *s);
  // 
  // void out_triple(Site *s1, Site *s2,Site * s3);

  Freelist hfl;
  Halfedge *ELleftend, *ELrightend;
  int ELhashsize;
  
  int sorted, debug;
  double xmin, xmax, ymin, ymax, deltax, deltay;
  
  Site *sites;
  int nsites;
  int siteidx;
  int sqrt_nsites;
  int nvertices;
  Freelist sfl;
  Site *bottomsite;
  
  int nedges;
  Freelist efl;
  int PQhashsize;
  Halfedge *PQhash;
  int PQcount;
  int PQmin;
  
  int ntry, totalsearch;
  double pxmin, pxmax, pymin, pymax, cradius;
  int total_alloc;
  
  double borderMinX, borderMaxX, borderMinY, borderMaxY;
  
  FreeNodeArrayList* allMemoryList;
  FreeNodeArrayList* currentMemoryBlock;
  
  GraphEdge* allEdges;
  GraphEdge* iteratorEdges;
  
  double minDistanceBetweenSites;

  static LimitedWarning _warning_degeneracy;
};

inline
int scomp(const void *p1,const void *p2)
{
  VPoint *s1 = (VPoint*)p1, *s2=(VPoint*)p2;
  if(s1->y < s2->y) return(-1);
  if(s1->y > s2->y) return(1);
  if(s1->x < s2->x) return(-1);
  if(s1->x > s2->x) return(1);
  return(0);
}

/* return a single in-storage site */
Site * VoronoiDiagramGenerator::nextone()
{
  Site *s;
  if(siteidx < nsites)
    {	
      s = &sites[siteidx];
      siteidx += 1;
      return(s);
    }
  else	
    return( (Site *)NULL);
}



VoronoiDiagramGenerator::VoronoiDiagramGenerator(){
  siteidx = 0;
  sites = NULL;
  parent_sites = NULL;

  allMemoryList = new FreeNodeArrayList;
  allMemoryList->memory = NULL;
  allMemoryList->next = NULL;
  currentMemoryBlock = allMemoryList;
  allEdges = NULL;
  iteratorEdges = 0;
  minDistanceBetweenSites = 0;

  ELhash = NULL;
  PQhash = NULL;
}

VoronoiDiagramGenerator::~VoronoiDiagramGenerator(){
  cleanup();
  cleanupEdges();

  if (allMemoryList != NULL)
    delete allMemoryList;
}



bool VoronoiDiagramGenerator::generateVoronoi(std::vector<VPoint> *_parent_sites,
					      double minX, double maxX, 
					      double minY, double maxY, 
					      double minDist){
  cleanup();
  cleanupEdges();
  int i;
  double x, y;

  minDistanceBetweenSites = minDist;

  parent_sites = _parent_sites;

  nsites = n_parent_sites = parent_sites->size();
  debug = 1;
  sorted = 0; 
  freeinit(&sfl, sizeof (Site));
		
  //sites = (Site *) myalloc(nsites*sizeof( *sites));
  sites = (Site *) myalloc(nsites*sizeof( Site));
  //parent_sites = (Site *) myalloc(nsites*sizeof( Site));
 
  if (sites == 0)
    return false;

  xmax = xmin = (*parent_sites)[0].x;
  ymax = ymin = (*parent_sites)[0].y;
  
  for(i=0;i<nsites;i++){
    x = (*parent_sites)[i].x;
    y = (*parent_sites)[i].y;
    sites[i].coord.x = x;
    sites[i].coord.y = y;
    sites[i].sitenbr = i;
    sites[i].refcnt  = 0;
    
    if (x<xmin)
      xmin=x;
    else if (x>xmax)
      xmax=x;
    
    if (y<ymin)
      ymin=y;
    else if (y>ymax)
      ymax=y;
  }
	
  qsort(sites, nsites, sizeof (*sites), scomp);

  // Gregory Soyez
  // 
  // check if some of the particles are degenerate to avoid a crash.
  //
  // At the moment, we work under the assumption that they will be
  // "clustered" later on, so we just keep the 1st one and discard the
  // others
  unsigned int offset=0;
  for (int is=1;is<nsites;is++){
    if (sites[is].coord.y==sites[is-1].coord.y && sites[is].coord.x==sites[is-1].coord.x){
      offset++;
    } else if (offset>0){
      sites[is-offset] = sites[is];
    }
  }

  if (offset>0){
    nsites-=offset;
    _warning_degeneracy.warn("VoronoiDiagramGenerator: two (or more) particles are degenerate in rapidity and azimuth, Voronoi cell assigned to the first of each set of degenerate particles.");
  }

  siteidx = 0;
  geominit();
  double temp = 0;
  if(minX > maxX){
    temp = minX;
    minX = maxX;
    maxX = temp;
  }
  if(minY > maxY){
    temp = minY;
    minY = maxY;
    maxY = temp;
  }
  borderMinX = minX;
  borderMinY = minY;
  borderMaxX = maxX;
  borderMaxY = maxY;
	
  siteidx = 0;
  voronoi();

  return true;
}

bool VoronoiDiagramGenerator::ELinitialize(){
  int i;
  freeinit(&hfl, sizeof(Halfedge));
  ELhashsize = 2 * sqrt_nsites;
  //ELhash = (Halfedge **) myalloc ( sizeof *ELhash * ELhashsize);
  ELhash = (Halfedge **) myalloc ( sizeof(Halfedge*)*ELhashsize);
  
  if(ELhash == 0)
    return false;
  
  for(i=0; i<ELhashsize; i +=1) ELhash[i] = (Halfedge *)NULL;
  ELleftend = HEcreate((Edge*) NULL, 0);
  ELrightend = HEcreate((Edge*) NULL, 0);
  ELleftend->ELleft = (Halfedge*) NULL;
  ELleftend->ELright = ELrightend;
  ELrightend->ELleft = ELleftend;
  ELrightend->ELright = (Halfedge*) NULL;
  ELhash[0] = ELleftend;
  ELhash[ELhashsize-1] = ELrightend;

  return true;
}


Halfedge* VoronoiDiagramGenerator::HEcreate(Edge *e,int pm){
  Halfedge *answer;
  answer = (Halfedge*) getfree(&hfl);
  answer->ELedge = e;
  answer->ELpm = pm;
  answer->PQnext = (Halfedge*) NULL;
  answer->vertex = (Site*) NULL;
  answer->ELrefcnt = 0;
  return answer;
}


void VoronoiDiagramGenerator::ELinsert(Halfedge *lb, Halfedge *newHe)
{
  newHe->ELleft = lb;
  newHe->ELright = lb->ELright;
  (lb->ELright)->ELleft = newHe;
  lb->ELright = newHe;
}


/* Get entry from hash table, pruning any deleted nodes */
Halfedge* VoronoiDiagramGenerator::ELgethash(int b){
  Halfedge *he;
	
  if ((b<0) || (b>=ELhashsize)) 
    return (Halfedge *) NULL;

  he = ELhash[b]; 
  if ((he==(Halfedge*) NULL) || (he->ELedge!=(Edge*) DELETED)) 
    return he;
	
  /* Hash table points to deleted half edge.  Patch as necessary. */
  ELhash[b] = (Halfedge*) NULL;
  if ((he->ELrefcnt -= 1) == 0) 
    makefree((Freenode*)he, &hfl);
  return (Halfedge*) NULL;
}	

Halfedge * VoronoiDiagramGenerator::ELleftbnd(VPoint *p){
  int i, bucket;
  Halfedge *he;
	
  /* Use hash table to get close to desired halfedge */
  // use the hash function to find the place in the hash map that this
  // HalfEdge should be
  // Gregory Soyez: the original code was 
  //
  //   bucket = (int)((p->x - xmin)/deltax * ELhashsize);	
  //   // make sure that the bucket position in within the range of the
  //   //hash array
  //   if (bucket<0) bucket =0;
  //   if (bucket>=ELhashsize) bucket = ELhashsize - 1;
  //
  // but this runs the risk of having a overflow which would 
  // cause bucket to be truncated at 0 instead of ELhashsize - 1
  // (or vice-versa)
  // We fix this by performing the test immediately on the double
  // We put in a extra bit of margin to be sure the conversion does
  // not play dirty tricks on us

  //const double &px = p->x;
  if (p->x < xmin){ bucket=0;}
  else if (p->x >= xmax){ bucket = ELhashsize - 1;}
  else{
    bucket= (int)((p->x - xmin)/deltax * ELhashsize);
    if (bucket>=ELhashsize) bucket = ELhashsize - 1;  // the lower cut should be robust
  }

  he = ELgethash(bucket);

  // if the HE isn't found, search backwards and forwards in the hash
  // map for the first non-null entry
  if (he == (Halfedge*) NULL){   
    for(i=1;1;i++){	
      if ((he=ELgethash(bucket-i)) != (Halfedge*) NULL) 
	break;
      if ((he=ELgethash(bucket+i)) != (Halfedge*) NULL) 
	break;
    };
    totalsearch += i;
  };
  ntry += 1;
  /* Now search linear list of halfedges for the correct one */
  if ((he==ELleftend)  || (he != ELrightend && right_of(he,p))){
    do{
      he = he->ELright;
    } while (he!=ELrightend && right_of(he,p));
      // keep going right on the list until either the end is reached,
      // or you find the 1st edge which the point
      he = he->ELleft;	//isn't to the right of
  } else 
    // if the point is to the left of the HalfEdge, then search left
    // for the HE just to the left of the point
    do{
      he = he->ELleft;
    } while ((he!=ELleftend )&& (!right_of(he,p)));
		
  /* Update hash table and reference counts */
  if ((bucket > 0) && (bucket <ELhashsize-1)){	
    if(ELhash[bucket] != (Halfedge *) NULL){
      ELhash[bucket]->ELrefcnt -= 1;
    }
    ELhash[bucket] = he;
    ELhash[bucket]->ELrefcnt += 1;
  };

  return he;
}


/* This delete routine can't reclaim node, since pointers from hash
   table may be present.   */
void VoronoiDiagramGenerator::ELdelete(Halfedge *he){
  (he->ELleft)->ELright = he->ELright;
  (he->ELright)->ELleft = he->ELleft;
  he->ELedge = (Edge*) DELETED;
}


Halfedge* VoronoiDiagramGenerator::ELright(Halfedge *he){
  return he->ELright;
}


Halfedge* VoronoiDiagramGenerator::ELleft(Halfedge *he){
  return he->ELleft;
}


Site * VoronoiDiagramGenerator::leftreg(Halfedge *he){
  if (he->ELedge == (Edge*) NULL) 
    return(bottomsite);
  return (he->ELpm == le) 
    ? he->ELedge->reg[le] 
    : he->ELedge->reg[re];
}

Site * VoronoiDiagramGenerator::rightreg(Halfedge *he){
  if (he->ELedge == (Edge*) NULL)
    // if this halfedge has no edge, return the bottom site (whatever
    // that is)
    return bottomsite;
  
  // if the ELpm field is zero, return the site 0 that this edge
  // bisects, otherwise return site number 1
  return (he->ELpm == le)
    ? he->ELedge->reg[re] 
    : he->ELedge->reg[le];
}

void VoronoiDiagramGenerator::geominit(){
  double sn;
  
  freeinit(&efl, sizeof(Edge));
  nvertices = 0;
  nedges = 0;
  sn = (double)nsites+4;
  sqrt_nsites = (int)sqrt(sn);
  deltay = ymax - ymin;
  deltax = xmax - xmin;
}


Edge * VoronoiDiagramGenerator::bisect(Site *s1, Site *s2){
  double dx,dy,adx,ady;
  Edge *newedge;	

  newedge = (Edge*) getfree(&efl);
	
  newedge->reg[0] = s1; //store the sites that this edge is bisecting
  newedge->reg[1] = s2;
  ref(s1); 
  ref(s2);

  // to begin with, there are no endpoints on the bisector - it goes
  // to infinity
  newedge->ep[0] = (Site*) NULL;
  newedge->ep[1] = (Site*) NULL;
  
  // get the difference in x dist between the sites
  dx = s2->coord.x - s1->coord.x;
  dy = s2->coord.y - s1->coord.y;

  // make sure that the difference is positive
  adx = dx>0 ? dx : -dx;
  ady = dy>0 ? dy : -dy;

  // get the slope of the line
  newedge->c = (double)(s1->coord.x * dx + s1->coord.y * dy
			+ (dx*dx + dy*dy)*0.5);

  if (adx>ady){	
    //set formula of line, with x fixed to 1
    newedge->a = 1.0; newedge->b = dy/dx; newedge->c /= dx;
  } else {	
    //set formula of line, with y fixed to 1
    newedge->b = 1.0; newedge->a = dx/dy; newedge->c /= dy;
  }
	
  newedge->edgenbr = nedges;
  nedges++;

  return newedge;
}


// create a new site where the HalfEdges el1 and el2 intersect - note
// that the VPoint in the argument list is not used, don't know why
// it's there
//
// Gregory Soyez: removed the uinused point p
Site* VoronoiDiagramGenerator::intersect(Halfedge *el1, Halfedge *el2
					 /*, VPoint *p*/){
  Edge *e1,*e2, *e;
  Halfedge *el;
  double d, xint, yint;
  int right_of_site;
  Site *v;
	
  e1 = el1->ELedge;
  e2 = el2->ELedge;
  if ((e1==(Edge*)NULL) || (e2==(Edge*)NULL)) 
    return ((Site*) NULL);

  // if the two edges bisect the same parent, return null
  if (e1->reg[1] == e2->reg[1]) 
    return (Site*) NULL;

  // Gregory Soyez:
  //	
  // if the 2 parents are too close, the intersection is going to be
  // computed from the "long edges" of the triangle which could causes
  // large rounding errors. In this case, use the bisector of the 2
  // parents to find the interaction point
  // 
  // The following replaces 
  //   d = e1->a * e2->b - e1->b * e2->a;
  //   if (-1.0e-10<d && d<1.0e-10) 
  //     return (Site*) NULL;
  //   	
  //   xint = (e1->c*e2->b - e2->c*e1->b)/d;
  //   yint = (e2->c*e1->a - e1->c*e2->a)/d;

  double dx = e2->reg[1]->coord.x - e1->reg[1]->coord.x;
  double dy = e2->reg[1]->coord.y - e1->reg[1]->coord.y;
  double dxref = e1->reg[1]->coord.x - e1->reg[0]->coord.x;
  double dyref = e1->reg[1]->coord.y - e1->reg[0]->coord.y;

  if (dx*dx + dy*dy < 1e-14*(dxref*dxref+dyref*dyref)){
    // make sure that the difference is positive
    double adx = dx>0 ? dx : -dx;
    double ady = dy>0 ? dy : -dy;
    
    // get the slope of the line
    double a,b;
    double c = (double)(e1->reg[1]->coord.x * dx + e1->reg[1]->coord.y * dy
			+ (dx*dx + dy*dy)*0.5);
    
    if (adx>ady){
      a = 1.0; b = dy/dx; c /= dx;
    } else {
      b = 1.0; a = dx/dy; c /= dy;
    }

    d = e1->a * b - e1->b * a;
    if (-1.0e-10<d && d<1.0e-10) {
      return (Site*) NULL;
    }
	
    xint = (e1->c*b - c*e1->b)/d;
    yint = (c*e1->a - e1->c*a)/d;
    
  } else {	
    d = e1->a * e2->b - e1->b * e2->a;
    if (-1.0e-10<d && d<1.0e-10) {
      return (Site*) NULL;
    }
	
    xint = (e1->c*e2->b - e2->c*e1->b)/d;
    yint = (e2->c*e1->a - e1->c*e2->a)/d;
  }
  // end of Gregory Soyez's modifications

  volatile double local_y1 = e1->reg[1]->coord.y;
  volatile double local_y2 = e2->reg[1]->coord.y;
	
  if( (local_y1 < local_y2) ||
      ((local_y1 == local_y2) &&
       (e1->reg[1]->coord.x <  e2->reg[1]->coord.x)) ){	
    el = el1; 
    e = e1;
  } else {	
    el = el2; 
    e = e2;
  }
	
  right_of_site = xint >= e->reg[1]->coord.x;
  if ((right_of_site && el->ELpm == le) || (!right_of_site && el->ELpm == re)) 
    return (Site*) NULL;
	
  // create a new site at the point of intersection - this is a new
  // vector event waiting to happen
  v = (Site*) getfree(&sfl);
  v->refcnt = 0;
  v->coord.x = xint;
  v->coord.y = yint;
  return v;
}

//HERE

/* returns 1 if p is to right of halfedge e */
int VoronoiDiagramGenerator::right_of(Halfedge *el,VPoint *p)
{
  Edge *e;
  Site *topsite;
  int right_of_site, above, fast;
  double dxp, dyp, dxs, t1, t2, t3, yl;
	
  e = el->ELedge;
  topsite = e->reg[1];
  right_of_site = p->x > topsite->coord.x;
  if(right_of_site && el->ELpm == le) return(1);
  if(!right_of_site && el->ELpm == re) return (0);
	
  if (e->a == 1.0)
    {	dyp = p->y - topsite->coord.y;
    dxp = p->x - topsite->coord.x;
    fast = 0;
    if ((!right_of_site & (e->b<0.0)) | (right_of_site & (e->b>=0.0)) )
      {	above = dyp>= e->b*dxp;	
      fast = above;
      }
    else 
      {	above = p->x + p->y*e->b > e-> c;
      if(e->b<0.0) above = !above;
      if (!above) fast = 1;
      };
    if (!fast)
      {	dxs = topsite->coord.x - (e->reg[0])->coord.x;
      above = e->b * (dxp*dxp - dyp*dyp) <
	dxs*dyp*(1.0+2.0*dxp/dxs + e->b*e->b);
      if(e->b<0.0) above = !above;
      };
    }
  else  /*e->b==1.0 */
    {	yl = e->c - e->a*p->x;
    t1 = p->y - yl;
    t2 = p->x - topsite->coord.x;
    t3 = yl - topsite->coord.y;
    above = t1*t1 > t2*t2 + t3*t3;
    };
  return (el->ELpm==le ? above : !above);
}


void VoronoiDiagramGenerator::endpoint(Edge *e,int lr,Site * s)
{
  e->ep[lr] = s;
  ref(s);
  if(e->ep[re-lr]== (Site *) NULL) 
    return;

  clip_line(e);

  deref(e->reg[le]);
  deref(e->reg[re]);
  makefree((Freenode*)e, &efl);
}


double VoronoiDiagramGenerator::dist(Site *s,Site *t)
{
  double dx,dy;
  dx = s->coord.x - t->coord.x;
  dy = s->coord.y - t->coord.y;
  return (double)(sqrt(dx*dx + dy*dy));
}


void VoronoiDiagramGenerator::makevertex(Site *v)
{
  v->sitenbr = nvertices;
  nvertices += 1;
  //GS unused plot: out_vertex(v);
}


void VoronoiDiagramGenerator::deref(Site *v)
{
  v->refcnt -= 1;
  if (v->refcnt == 0 ) 
    makefree((Freenode*)v, &sfl);
}

void VoronoiDiagramGenerator::ref(Site *v)
{
  v->refcnt += 1;
}

//push the HalfEdge into the ordered linked list of vertices
void VoronoiDiagramGenerator::PQinsert(Halfedge *he,Site * v, double offset)
{
  Halfedge *last, *next;
	
  he->vertex = v;
  ref(v);
  he->ystar = (double)(v->coord.y + offset);
  last = &PQhash[PQbucket(he)];
  while ((next = last->PQnext) != (Halfedge *) NULL &&
	 (he->ystar  > next->ystar  ||
	  (he->ystar == next->ystar && v->coord.x > next->vertex->coord.x)))
    {	
      last = next;
    };
  he->PQnext = last->PQnext; 
  last->PQnext = he;
  PQcount += 1;
}

//remove the HalfEdge from the list of vertices 
void VoronoiDiagramGenerator::PQdelete(Halfedge *he)
{
  Halfedge *last;
	
  if(he->vertex != (Site *) NULL)
    {	
      last = &PQhash[PQbucket(he)];
      while (last->PQnext != he) 
	last = last->PQnext;

      last->PQnext = he->PQnext;
      PQcount -= 1;
      deref(he->vertex);
      he->vertex = (Site *) NULL;
    };
}

int VoronoiDiagramGenerator::PQbucket(Halfedge *he)
{
  // Gregory Soyez: the original code was 
  //
  //   bucket = (int)((he->ystar - ymin)/deltay * PQhashsize);
  //   if (bucket<0) bucket = 0;
  //   if (bucket>=PQhashsize) bucket = PQhashsize-1 ;
  //   if (bucket < PQmin) PQmin = bucket;
  //   return(bucket);
  //
  // but this runs the risk of having a overflow which would 
  // cause bucket to be truncated at 0 instead of PQhashsize-1
  // (or vice-versa)
  // We fix this by performing the test immediately on the double
  // We put in a extra bit of margin to be sure the conversion does
  // not play dirty tricks on us

  int bucket;
	
  double hey = he->ystar;
  if (hey < ymin){ bucket = 0;}
  else if (hey >= ymax){ bucket = PQhashsize-1;}
  else {
    bucket = (int)((hey - ymin)/deltay * PQhashsize);
    if (bucket>=PQhashsize) bucket = PQhashsize-1 ;
  }

  if (bucket < PQmin) PQmin = bucket;
  return(bucket);
}



int VoronoiDiagramGenerator::PQempty()
{
  return(PQcount==0);
}


VPoint VoronoiDiagramGenerator::PQ_min()
{
  VPoint answer;
	
  while(PQhash[PQmin].PQnext == (Halfedge *)NULL) {PQmin += 1;};
  answer.x = PQhash[PQmin].PQnext->vertex->coord.x;
  answer.y = PQhash[PQmin].PQnext->ystar;
  return (answer);
}

Halfedge * VoronoiDiagramGenerator::PQextractmin()
{
  Halfedge *curr;
	
  curr = PQhash[PQmin].PQnext;
  PQhash[PQmin].PQnext = curr->PQnext;
  PQcount -= 1;
  return(curr);
}


bool VoronoiDiagramGenerator::PQinitialize()
{
  int i; 
	
  PQcount = 0;
  PQmin = 0;
  PQhashsize = 4 * sqrt_nsites;
  //PQhash = (Halfedge *) myalloc(PQhashsize * sizeof *PQhash);
  PQhash = (Halfedge *) myalloc(PQhashsize * sizeof(Halfedge));

  if(PQhash == 0)
    return false;

  for(i=0; i<PQhashsize; i+=1) PQhash[i].PQnext = (Halfedge *)NULL;

  return true;
}


void VoronoiDiagramGenerator::freeinit(Freelist *fl,int size)
{
  fl->head = (Freenode *) NULL;
  fl->nodesize = size;
}

char * VoronoiDiagramGenerator::getfree(Freelist *fl)
{
  int i; 
  Freenode *t;

  if(fl->head == (Freenode *) NULL)
    {	
      t =  (Freenode *) myalloc(sqrt_nsites * fl->nodesize);

      if(t == 0)
	return 0;
		
      currentMemoryBlock->next = new FreeNodeArrayList;
      currentMemoryBlock = currentMemoryBlock->next;
      currentMemoryBlock->memory = t;
      currentMemoryBlock->next = 0;

      for(i=0; i<sqrt_nsites; i+=1) 	
	makefree((Freenode *)((char *)t+i*fl->nodesize), fl);		
    };
  t = fl->head;
  fl->head = (fl->head)->nextfree;
  return((char *)t);
}



void VoronoiDiagramGenerator::makefree(Freenode *curr,Freelist *fl)
{
  curr->nextfree = fl->head;
  fl->head = curr;
}

void VoronoiDiagramGenerator::cleanup()
{
  if(sites != NULL){
    free(sites);
    sites = 0;
  }

  FreeNodeArrayList* current=NULL, *prev=NULL;

  current = prev = allMemoryList;

  while(current->next!=NULL){
    prev = current;
    current = current->next;
    free(prev->memory);
    delete prev;
    prev = NULL;
  }

  if(current!=NULL){
    if (current->memory!=NULL){
      free(current->memory);
    }
    delete current;
  }

  allMemoryList = new FreeNodeArrayList;
  allMemoryList->next = NULL;
  allMemoryList->memory = NULL;
  currentMemoryBlock = allMemoryList;

  if (ELhash!=NULL)
    free(ELhash);

  if (PQhash!=NULL)
    free(PQhash);
}

void VoronoiDiagramGenerator::cleanupEdges()
{
  GraphEdge* geCurrent = NULL, *gePrev = NULL;
  geCurrent = gePrev = allEdges;

  while(geCurrent!=NULL){
    gePrev = geCurrent;
    geCurrent = geCurrent->next;
    delete gePrev;
  }

  allEdges = 0;

}

void VoronoiDiagramGenerator::pushGraphEdge(double x1, double y1, double x2, double y2,
					    Site *s1, Site *s2)
{
  GraphEdge* newEdge = new GraphEdge;
  newEdge->next = allEdges;
  allEdges = newEdge;
  newEdge->x1 = x1;
  newEdge->y1 = y1;
  newEdge->x2 = x2;
  newEdge->y2 = y2;

  newEdge->point1 = s1->sitenbr;
  newEdge->point2 = s2->sitenbr;
}


char * VoronoiDiagramGenerator::myalloc(unsigned n)
{
  char *t=0;	
  t=(char*)malloc(n);
  total_alloc += n;
  return(t);
}


// unused plot functions
//
// /* for those who don't have Cherry's plot */
// /* #include <plot.h> */
// void VoronoiDiagramGenerator::openpl(){}
// void VoronoiDiagramGenerator::circle(double x, double y, double radius){}
// void VoronoiDiagramGenerator::range(double minX, double minY, double maxX, double maxY){}
// 
// 
// 
// void VoronoiDiagramGenerator::out_bisector(Edge *e)
// {
// 	
// 
// }
// 
// 
// void VoronoiDiagramGenerator::out_ep(Edge *e)
// {
// 	
// 	
// }
// 
// void VoronoiDiagramGenerator::out_vertex(Site *v)
// {
// 	
// }
// 
// 
// void VoronoiDiagramGenerator::out_site(Site *s)
// {
//   // Gregory Soyez: 
//   //   plot was always 0 so the expression below was always false
//   //   and even if it was not, 'circle' does nothing!
//   //
//   // if(!triangulate & plot & !debug)
//   //   circle (s->coord.x, s->coord.y, cradius);
// 	
// }
// 
// 
// void VoronoiDiagramGenerator::out_triple(Site *s1, Site *s2,Site * s3)
// {
// 	
// }



void VoronoiDiagramGenerator::plotinit()
{
  double dx,dy,d;
	
  dy = ymax - ymin;
  dx = xmax - xmin;
  d = (double)(( dx > dy ? dx : dy) * 1.1);
  pxmin = (double)(xmin - (d-dx)/2.0);
  pxmax = (double)(xmax + (d-dx)/2.0);
  pymin = (double)(ymin - (d-dy)/2.0);
  pymax = (double)(ymax + (d-dy)/2.0);
  cradius = (double)((pxmax - pxmin)/350.0);
  //GS unused: openpl();
  //GS unused: range(pxmin, pymin, pxmax, pymax);
}


void VoronoiDiagramGenerator::clip_line(Edge *e)
{
  Site *s1, *s2;
  double x1=0,x2=0,y1=0,y2=0; //, temp = 0;

  x1 = e->reg[0]->coord.x;
  x2 = e->reg[1]->coord.x;
  y1 = e->reg[0]->coord.y;
  y2 = e->reg[1]->coord.y;

  //if the distance between the two points this line was created from is less than 
  //the square root of 2, then ignore it
  //TODO improve/remove
  //if(sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) < minDistanceBetweenSites)
  //  {
  //    return;
  //  }
  pxmin = borderMinX;
  pxmax = borderMaxX;
  pymin = borderMinY;
  pymax = borderMaxY;

  if(e->a == 1.0 && e ->b >= 0.0)
    {	
      s1 = e->ep[1];
      s2 = e->ep[0];
    }
  else 
    {
      s1 = e->ep[0];
      s2 = e->ep[1];
    };
	
  if(e->a == 1.0)
    {
      y1 = pymin;
      if (s1!=(Site *)NULL && s1->coord.y > pymin)
	{
	  y1 = s1->coord.y;
	}
      if(y1>pymax) 
	{
	  //	printf("\nClipped (1) y1 = %f to %f",y1,pymax);
	  y1 = pymax;
	  //return;
	}
      x1 = e->c - e->b * y1;
      y2 = pymax;
      if (s2!=(Site *)NULL && s2->coord.y < pymax) 
	y2 = s2->coord.y;

      if(y2<pymin) 
	{
	  //printf("\nClipped (2) y2 = %f to %f",y2,pymin);
	  y2 = pymin;
	  //return;
	}
      x2 = (e->c) - (e->b) * y2;
      if (((x1> pxmax) & (x2>pxmax)) | ((x1<pxmin)&(x2<pxmin))) 
	{
	  //printf("\nClipLine jumping out(3), x1 = %f, pxmin = %f, pxmax = %f",x1,pxmin,pxmax);
	  return;
	}
      if(x1> pxmax)
	{	x1 = pxmax; y1 = (e->c - x1)/e->b;};
      if(x1<pxmin)
	{	x1 = pxmin; y1 = (e->c - x1)/e->b;};
      if(x2>pxmax)
	{	x2 = pxmax; y2 = (e->c - x2)/e->b;};
      if(x2<pxmin)
	{	x2 = pxmin; y2 = (e->c - x2)/e->b;};
    }
  else
    {
      x1 = pxmin;
      if (s1!=(Site *)NULL && s1->coord.x > pxmin) 
	x1 = s1->coord.x;
      if(x1>pxmax) 
	{
	  //printf("\nClipped (3) x1 = %f to %f",x1,pxmin);
	  //return;
	  x1 = pxmax;
	}
      y1 = e->c - e->a * x1;
      x2 = pxmax;
      if (s2!=(Site *)NULL && s2->coord.x < pxmax) 
	x2 = s2->coord.x;
      if(x2<pxmin) 
	{
	  //printf("\nClipped (4) x2 = %f to %f",x2,pxmin);
	  //return;
	  x2 = pxmin;
	}
      y2 = e->c - e->a * x2;
      if (((y1> pymax) & (y2>pymax)) | ((y1<pymin)&(y2<pymin))) 
	{
	  //printf("\nClipLine jumping out(6), y1 = %f, pymin = %f, pymax = %f",y2,pymin,pymax);
	  return;
	}
      if(y1> pymax)
	{	y1 = pymax; x1 = (e->c - y1)/e->a;};
      if(y1<pymin)
	{	y1 = pymin; x1 = (e->c - y1)/e->a;};
      if(y2>pymax)
	{	y2 = pymax; x2 = (e->c - y2)/e->a;};
      if(y2<pymin)
	{	y2 = pymin; x2 = (e->c - y2)/e->a;};
    };
	
  //printf("\nPushing line (%f,%f,%f,%f)",x1,y1,x2,y2);
  //fprintf(stdout, "Line with vertices (%f,%f) and (%f,%f)\n", 
  //	e->reg[0]->coord.x, e->reg[1]->coord.x, e->reg[0]->coord.y, e->reg[1]->coord.y);
  pushGraphEdge(x1,y1,x2,y2,e->reg[0],e->reg[1]);
}



FASTJET_END_NAMESPACE

#endif
