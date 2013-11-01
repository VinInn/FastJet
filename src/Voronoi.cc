//STARTHEADER
// $Id: Voronoi.cc 2767 2011-11-24 21:43:09Z salam $
//
// Copyright (c) 1994 by AT&T Bell Laboratories (see below)
//
//
//----------------------------------------------------------------------
// This file is included as part of FastJet but was mostly written by
// S. Fortune in C, put into C++ with memory management by S
// O'Sullivan, and with further interface and memory management
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
 * O'Sullivan under the terms reported above.
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
 * 2011-06-28  Gregory Soyez  <soyez@fastjet.fr>
 * 
 *      * added support for situations with degenerate particles (we just
 *        discard the particles degenerate wiht an already existing
 *        one which is perfectly sufficient for our needs)
 *      * in 'VoronoiDiagramGenerator::intersect', improved the numerical
 *        precision in cases where the 2 parents are nearly degenerate
 * 
 * 
 * 2011-06-14  Gregory Soyez  <soyez@fastjet.fr>
 * 
 *      * fixed a potential overflow bug in VoronoiDiagramGenerator::PQbucket
 * 
 * 
 * 2007-05-07  Gregory Soyez  <soyez@fastjet.fr>
 * 
 *      * fied a few memory leaks
 *
 *      * put the code in the fastjet namespace
 * 
 *      * replaced float by double
 * 
 *      * generateVoronoi() takes a vector of Point instead of 2
 *        pointers
 * 
 *      * added info about the parent sites to GraphEdge (and clip_line)
 * 
 *      * removed condition on minimal distance between sites
 */

#include <stdio.h>
#include "fastjet/internal/Voronoi.hh"

using namespace std;

FASTJET_BEGIN_NAMESPACE

LimitedWarning VoronoiDiagramGenerator::_warning_degeneracy;


/* implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax,
   deltax, deltay (can all be estimates).
   Performance suffers if they are wrong; better to make nsites,
   deltax, and deltay too big than too small.  (?) */

bool VoronoiDiagramGenerator::voronoi()
{
  Site *newsite, *bot, *top, *temp, *p;
  Site *v;
  VPoint newintstar;
  int pm;
  Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
  Edge *e;
	
  PQinitialize();
  bottomsite = nextone();
  //GS unused plot: out_site(bottomsite);
  bool retval = ELinitialize();

  if(!retval)
    return false;
	
  newsite = nextone();
  while(1)
    {

      if(!PQempty()) 
	newintstar = PQ_min();
		
      //if the lowest site has a smaller y value than the lowest vector intersection, process the site
      //otherwise process the vector intersection		
      if (newsite != (Site *)NULL  && (PQempty() || newsite->coord.y < newintstar.y
				       || (newsite->coord.y == newintstar.y && newsite->coord.x < newintstar.x)))
	{/* new site is smallest - this is a site event*/
	  //GS unused plot: out_site(newsite);						//output the site
	  lbnd = ELleftbnd(&(newsite->coord));				//get the first HalfEdge to the LEFT of the new site
	  rbnd = ELright(lbnd);						//get the first HalfEdge to the RIGHT of the new site
	  bot = rightreg(lbnd);						//if this halfedge has no edge, , bot = bottom site (whatever that is)
	  e = bisect(bot, newsite);					//create a new edge that bisects 
	  bisector = HEcreate(e, le);					//create a new HalfEdge, setting its ELpm field to 0			
	  ELinsert(lbnd, bisector);					//insert this new bisector edge between the left and right vectors in a linked list	
	    
	  if ((p = intersect(lbnd, bisector)) != (Site *) NULL) 	//if the new bisector intersects with the left edge, remove the left edge's vertex, and put in the new one
	    {	
	      PQdelete(lbnd);
	      PQinsert(lbnd, p, dist(p,newsite));
	    };
	  lbnd = bisector;						
	  bisector = HEcreate(e, re);					//create a new HalfEdge, setting its ELpm field to 1
	  ELinsert(lbnd, bisector);					//insert the new HE to the right of the original bisector earlier in the IF stmt
	    
	  if ((p = intersect(bisector, rbnd)) != (Site *) NULL)	//if this new bisector intersects with the
	    {	
	      PQinsert(bisector, p, dist(p,newsite));			//push the HE into the ordered linked list of vertices
	    };
	  newsite = nextone();	
	}
      else if (!PQempty()) /* intersection is smallest - this is a vector event */			
	{	
	  lbnd = PQextractmin();						//pop the HalfEdge with the lowest vector off the ordered list of vectors				
	  llbnd = ELleft(lbnd);						//get the HalfEdge to the left of the above HE
	  rbnd = ELright(lbnd);						//get the HalfEdge to the right of the above HE
	  rrbnd = ELright(rbnd);						//get the HalfEdge to the right of the HE to the right of the lowest HE 
	  bot = leftreg(lbnd);						//get the Site to the left of the left HE which it bisects
	  top = rightreg(rbnd);						//get the Site to the right of the right HE which it bisects
	    
	  //GS unused plot: out_triple(bot, top, rightreg(lbnd));		//output the triple of sites, stating that a circle goes through them
	    
	  v = lbnd->vertex;						//get the vertex that caused this event
	  makevertex(v);							//set the vertex number - couldn't do this earlier since we didn't know when it would be processed
	  endpoint(lbnd->ELedge,lbnd->ELpm,v);	//set the endpoint of the left HalfEdge to be this vector
	  endpoint(rbnd->ELedge,rbnd->ELpm,v);	//set the endpoint of the right HalfEdge to be this vector
	  ELdelete(lbnd);							//mark the lowest HE for deletion - can't delete yet because there might be pointers to it in Hash Map	
	  PQdelete(rbnd);							//remove all vertex events to do with the  right HE
	  ELdelete(rbnd);							//mark the right HE for deletion - can't delete yet because there might be pointers to it in Hash Map	
	  pm = le;								//set the pm variable to zero
	    
	  if (bot->coord.y > top->coord.y)		//if the site to the left of the event is higher than the Site
	    {										//to the right of it, then swap them and set the 'pm' variable to 1
	      temp = bot; 
	      bot = top; 
	      top = temp; 
	      pm = re;
	    }
	  e = bisect(bot, top);					//create an Edge (or line) that is between the two Sites. This creates
	  //the formula of the line, and assigns a line number to it
	  bisector = HEcreate(e, pm);				//create a HE from the Edge 'e', and make it point to that edge with its ELedge field
	  ELinsert(llbnd, bisector);				//insert the new bisector to the right of the left HE
	  endpoint(e, re-pm, v);					//set one endpoint to the new edge to be the vector point 'v'.
	  //If the site to the left of this bisector is higher than the right
	  //Site, then this endpoint is put in position 0; otherwise in pos 1
	  deref(v);								//delete the vector 'v'
	    
	  //if left HE and the new bisector don't intersect, then delete the left HE, and reinsert it 
	  if((p = intersect(llbnd, bisector)) != (Site *) NULL)
	    {	
	      PQdelete(llbnd);
	      PQinsert(llbnd, p, dist(p,bot));
	    };
	    
	  //if right HE and the new bisector don't intersect, then reinsert it 
	  if ((p = intersect(bisector, rrbnd)) != (Site *) NULL)
	    {	
	      PQinsert(bisector, p, dist(p,bot));
	    };
	}
      else break;
    };

	


  for(lbnd=ELright(ELleftend); lbnd != ELrightend; lbnd=ELright(lbnd))
    {	
      e = lbnd->ELedge;

      clip_line(e);
    };

  //cleanup();

  return true;
	
}




FASTJET_END_NAMESPACE
