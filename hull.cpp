// hull.cpp
//
// constructing the convex hull of a SimplicialComplex
// -------------------------------------------------------------------
// This file is part of Geode
// Copyright (C) 2003, 2006 Jesse H. Willett
// email: jhw@lmi.net, jhwjhw@gmail.com
// web:   http://users.lmi.net/~jhw
//
// Geode is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// The file LICENSE should be included with this distribution, containing
// the full text of the GNU General Public License.
//
// There may be restrictions to this software's release under the GPL: please
// see also the file README.

#include "hull.hpp"
#include "simplicial.hpp"
#include <string>
#include <iostream>
using std::cerr;
using std::string;
using std::vector;
using __gnu_cxx::hash_set;
using __gnu_cxx::hash_map;

const static bool debug = false;

ConvexHull::ConvexHull (SimplicialComplex &sci, double volumeSignTolerancei) :
  sc(sci),
  volumeSignTolerance(volumeSignTolerancei)  
{
  base = NULL;
}

bool ConvexHull::isProcessed (const Vertex &v) const
{
  return processed.find(const_cast<Vertex *>(&v)) != processed.end();
}

bool ConvexHull::isOnHull (const Vertex &v) const
{
  return onHull.find(const_cast<Vertex *>(&v)) != onHull.end();
}

bool ConvexHull::isVisible (const Face &f) const
{
  return visible.find(const_cast<Face *>(&f)) != visible.end();
}


/*---------------------------------------------------------------------
	MakeFace creates a new face structure from three vertices (in ccw
	order).  It returns a pointer to the face.
	---------------------------------------------------------------------*/
Face &SimplicialComplex::MakeFace (Material &material,
																	 Vertex *v0,
																	 Vertex *v1,
																	 Vertex *v2)
{
	/* Create edges of the initial triangle. */
	Edge &e0 = MakeEdge(material,*v0,*v1);
	Edge &e1 = MakeEdge(material,*v1,*v2);
	Edge &e2 = MakeEdge(material,*v2,*v0);
	
	/* Create face for triangle. */
	Face &f = MakeNullFace();
	f.edge[0]   = &e0;  f.edge[1]   = &e1; f.edge[2]   = &e2;
	f.vertex[0] =  v0;  f.vertex[1] =  v1; f.vertex[2] =  v2;
	
	/* Link edges to face. */
	e0.face[0] = e1.face[0] = e2.face[0] = &f;
	
	return f;
}
Face &SimplicialComplex::MakeFace (Face &fold,
																	 Vertex *v0,
																	 Vertex *v1,
																	 Vertex *v2)
{
	/* Copy from fold, in reverse order. */
	Edge &e0 = *fold.edge[2];
	Edge &e1 = *fold.edge[1];
	Edge &e2 = *fold.edge[0];

	e0.vertex[0] = v0;              e0.vertex[1] = v1;
	e1.vertex[0] = v1;              e1.vertex[1] = v2;
	e2.vertex[0] = v2;              e2.vertex[1] = v0;
	
	/* Create face for triangle. */
	Face &f = MakeNullFace();
	f.edge[0]   = &e0;  f.edge[1]   = &e1; f.edge[2]   = &e2;
	f.vertex[0] =  v0;  f.vertex[1] =  v1; f.vertex[2] =  v2;
	
	/* Link edges to face. */
	e0.face[0] = e1.face[0] = e2.face[0] = &f;
	
	return f;
}

bool ConvexHull::isAdvancingHull ()
{
  return hullIt != sc.vertices.end();
}

void ConvexHull::StartHull ()
{
  // clear any preexisting edges or faces
  sc.edges.clear();
  sc.faces.clear();
  sc.materials.clear();
  base = &sc.MakeMaterial(  0.0, 0.20, true,  0.1);

  // experimental shuffle with vertices:
  if (false) {
    vector<Vertex> arr;
    for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++) {
      Vertex &v = it->second;
      arr.push_back(v);
    }

    random_shuffle(arr.begin(),arr.end());

    sc.vertices.clear();
    for (unsigned int i = 0; i < arr.size(); i++) {
      arr[i].objID = sc.currentObjID++;
      sc.vertices[arr[i].objID] = arr[i];
    }
  }

  // clear any preexisting state (in case this object is being reused)
  processed.clear();
  duplicate.clear();
  for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++)
    onHull.insert(&it->second);

  // set up initial conditions:
  Vertex *nextAddOne = DoubleTriangle();
  if (nextAddOne == NULL) {
    hullIt = sc.vertices.end();
    return;
  }
  else if (isProcessed(*nextAddOne))
    throw "nextAddOne is processed";
  else {
    if (debug) cerr << "nextAddOne id: " << nextAddOne->objID << "\n";
    processed.insert(nextAddOne);
    if (debug) cerr << "calling AddOne() on nextAddOne\n";
    MarkVisiblity(*nextAddOne);
    if (!isOnHull(*nextAddOne))
      throw "unexpected nextAddOne not on hull";
    AddOne(*nextAddOne);
    if (debug) cerr << "calling CleanUp() on nextAddOne\n";
    // nextAddOne is specified to be noncoplanar:
    // (NOT REALLY: this one is known to be on hull right now
    if (CleanUp(nextAddOne->objID)) throw "nextAddOne needs erasure";
    if (debug) cerr << "clear of CleanUp()on nextAddOne\n";
  }

  hullIt = sc.vertices.begin();
}

void ConvexHull::AdvanceHull ()
{
  if (hullIt == sc.vertices.end()) return;

  if (debug) cerr << "next vertex: " << hullIt->first << "\n";

  while (hullIt != sc.vertices.end() && isProcessed(hullIt->second))
    hullIt++;

  if (hullIt == sc.vertices.end()) return;

  Vertex &v = hullIt->second;

  processed.insert(&v);

  MarkVisiblity(v);

  try {
    if (isOnHull(v))
      AddOne(v);
  } 
  catch (const char *str) {
    cerr << "exception in AddOne(): " << str << "\n";
    throw str;
  }

  // need to be careful about deleting v if CleanUp() wants to
  if (CleanUp(v.objID)) {
    if (debug) cerr << "erasing v at top level: " << v.objID << "\n";
    sc.vertices.erase(hullIt++);
  }
  else
    hullIt++;
}

/*---------------------------------------------------------------------
	AddOne is passed a vertex.  It first determines all faces visible from 
	that point.  If none are visible then the point is marked as not 
	onhull.  Next is a loop over sc.edges.  If both faces adjacent to an edge
	are visible, then the edge is marked for deletion.  If just one of the
	adjacent faces is visible then a new face is constructed.
	---------------------------------------------------------------------*/

void ConvexHull::MarkVisiblity (Vertex &p)
{
  /* Mark faces visible from p. */
  bool vis = false;
  for (FaceMapIt f = sc.faces.begin(); f != sc.faces.end(); f++) {
    double vol = VolumeSign(f->second,p);
    if (vol < 0.0) {
      visible.insert(&f->second);
      vis = true;
    }
  }
  if (!vis)
    onHull.erase(onHull.find(&p));
}

bool ConvexHull::AddOne (Vertex &p)
{
  // Mark edges in interior of visible region for deletion.
  // Erect a newface based on each border edge.
  // NOTE: This only needs to be considered over pre-existing geometry.
  //       The original code also scanned any edges added by MakeConeFace()
  //       -jhw

  for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++) {
    Edge &e = it->second;
    if (e.face[0] == NULL) throw "bogus early NULL face[0]";
    if (e.face[1] == NULL) throw "bogus early NULL face[1]";
  }

  hash_set<Vertex *> borderOfVisibleRegionV;
  hash_set<Edge *> borderOfVisibleRegionE;

  const int upperBound = sc.currentObjID;

  hash_map<int,Edge *> sort;
  for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++)
    sort[it->first] = &it->second;

  // OK, the problem is with mutation of the underlying sc.edges during
  // the loop.  This problem probably exists elsewhere, such as in
  // SubdivideFaces, etc.
  //    - yep, confirmed.  icosa | hull | subdivide 3 fails similarly.
  // This is bad... 
  //    - it might mean lots of code subtly depends on ordering
  //    - it might mean hash_map iterators are unstable under insertion
  // In either case, an automated test suite is called for!
  
  for (hash_map<int,Edge *>::iterator it = sort.begin(); 
       it != sort.end(); it++) {
    Edge &e = *it->second;
    if (e.objID >= upperBound) continue;
    if (e.face[0] == NULL) throw "bogus forward NULL face[0]";
    if (e.face[1] == NULL) throw "bogus forward NULL face[1]";
    // e interior: mark for deletion.
    // e border: make a new face. 
    if (isVisible(*e.face[0]) && isVisible(*e.face[1])) {
      deadEdges.insert(&e);
    }
    else if (isVisible(*e.face[0]) || isVisible(*e.face[1])) {
      borderOfVisibleRegionE.insert(&e);
      borderOfVisibleRegionV.insert(e.vertex[0]);
      borderOfVisibleRegionV.insert(e.vertex[1]);
      newFaces[&e] = &MakeConeFace(e, p);
    }
  }

  if (debug || 
      borderOfVisibleRegionV.size() != 
      borderOfVisibleRegionE.size()) {
    cerr << "\n";
    cerr << "vertices: " << 
      borderOfVisibleRegionV.size() << "\n";
    cerr << "   edges: " << 
      borderOfVisibleRegionE.size() << "\n";
    cerr << "   diff: " << 
      (long(borderOfVisibleRegionE.size()) -
       long(borderOfVisibleRegionV.size())) 
				 << "\n";
  }

  for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++) {
    Edge &e = it->second;
    if (e.face[0] == NULL) throw "bogus late NULL face[0]";
    if (e.face[1] == NULL) throw "bogus late NULL face[1]";
  }

  if (debug) cerr << "done\n";

  return true;
}

/*---------------------------------------------------------------------
	MakeConeFace makes a new face and two new edges between the 
	edge and the point that are passed to it. It returns a pointer to
	the new face.
	---------------------------------------------------------------------*/
Face &ConvexHull::MakeConeFace (Edge &e, Vertex &p)
{
  // make two new edges (if don't already exist):
  for (int i = 0; i < Edge::numVertices; ++i) {
    if (duplicate.find(e.vertex[i]) != duplicate.end()) continue;
    duplicate[e.vertex[i]] = &sc.MakeEdge(*base,*e.vertex[i],p);
  }

  // make the new face
  Face &new_face = sc.MakeNullFace();   
  new_face.edge[0] = &e;
  new_face.edge[1] = duplicate[e.vertex[0]];
  new_face.edge[2] = duplicate[e.vertex[1]];
  MakeCcw(new_face, e, p); 
        
  // set the adjacent face pointers, if unset
  for (int i = 0; i < Edge::numVertices; ++i) {
    for (int j = 0; j < Edge::numFaces; ++j) {
      Edge *dup = duplicate[e.vertex[i]];
      if (dup->face[j] != NULL) continue;
      dup->face[j] = &new_face;
      // only one NULL link should be set to new_face
      break;
    }
  }

  return new_face;
}

/*---------------------------------------------------------------------
	MakeCcw puts the vertices in the face structure in counterclock wise 
	order.  We want to store the vertices in the same 
	order as in the visible face.  The third vertex is always p.

	Although no specific ordering of the edges of a face are used
	by the code, the following condition is maintained for each face f:
	one of the two endpoints of f->edge[i] matches f->vertex[i]. 
	But note that this does not imply that f->edge[i] is between
	f->vertex[i] and f->vertex[(i+1)%3].  (Thanks to Bob Williamson.)
	---------------------------------------------------------------------*/
void ConvexHull::MakeCcw (Face &f, Edge &e, Vertex &p)
{
	int    i;    /* Index of e->endpoint[0] in fv. */
      
	Face *fv;   /* The visible face adjacent to e */
	if  (isVisible(*e.face[0]))
		fv = e.face[0];
	else fv = e.face[1];
       
	/* Set vertex[0] & [1] of f to have the same orientation
		 as do the corresponding vertices of fv. */ 
	for ( i=0; fv->vertex[i] != e.vertex[0]; ++i )
		;
	/* Orient f the same as fv. */
	if ( fv->vertex[ (i+1) % 3 ] != e.vertex[1] ) {
		f.vertex[0] = e.vertex[1];  
		f.vertex[1] = e.vertex[0];    
	}
	else {                               
		f.vertex[0] = e.vertex[0];   
		f.vertex[1] = e.vertex[1];   
		// swap
		Edge *s = f.edge[1];
		f.edge[1] = f.edge[2];
		f.edge[2] = s;
	}
	/* This swap is tricky. e is edge[0]. edge[1] is based on endpt[0],
		 edge[2] on endpt[1].  So if e is oriented "forwards," we
		 need to move edge[1] to follow [0], because it precedes. */
   
	f.vertex[2] = &p;
}

double ConvexHull::VolumeSign (const Face &f, const Vertex &p)
{
  double vol = f.volume(p.pos);
  /* The volume should be an integer. */
  if      ( vol >  volumeSignTolerance)  return  1.0;
  else if ( vol < -volumeSignTolerance)  return -1.0;
  else                    return  0.0;
}

/*---------------------------------------------------------------------
	CleanUp goes through each data structure list and clears all
	flags and NULLs out some pointers.  The order of processing
	(edges, faces, vertices) is important.
	---------------------------------------------------------------------*/
bool ConvexHull::CleanUp (int currentVertexID)
{
  CleanEdges();
  CleanFaces();
  return CleanVertices(currentVertexID);
}

/*---------------------------------------------------------------------
	CleanEdges runs through the edge list and cleans up the structure.
	If there is a newface then it will put that face in place of the 
	visible face and NULL out newface. It also deletes so marked sc.edges.
	---------------------------------------------------------------------*/
void ConvexHull::CleanEdges ()
{
  if (false && debug) cerr << "CleanEdges starting...";

  /* Integrate the newface's into the data structure. */
  /* Check every edge. */
  for (EdgeMapIt e = sc.edges.begin(); e != sc.edges.end(); e++) {
    if (newFaces.find(&e->second) == newFaces.end()) continue;
    if (isVisible(*e->second.face[0]))
      e->second.face[0] = newFaces[&e->second];
    else
      e->second.face[1] = newFaces[&e->second];

  }
  newFaces.clear();

  /* Delete any edges marked for deletion. */
  for (hash_set<Edge *>::iterator it = deadEdges.begin(); 
       it != deadEdges.end(); it++)
    sc.edges.erase(sc.edges.find((*it)->objID));
  deadEdges.clear();

  if (false && debug) cerr << "done\n";
}

/*---------------------------------------------------------------------
	CleanFaces runs through the face list and deletes any face marked visible.
	---------------------------------------------------------------------*/
void ConvexHull::CleanFaces ()
{
  if (false && debug) cerr << "CleanFaces starting...";
  for (FaceMapIt f = sc.faces.begin(); f != sc.faces.end(); ) {
    if (!isVisible(f->second)) {
      f++;
      continue;
    }
    if (false && debug) cerr << "CleanFaces erasing: " << f->first << "\n";
    sc.faces.erase(f++);
  }
  visible.clear();
  if (false && debug) cerr << "done\n";
}

/*---------------------------------------------------------------------
	CleanVertices runs through the vertex list and deletes the 
	vertices that are marked as processed but are not incident to any 
	undeleted sc.edges. 
	---------------------------------------------------------------------*/
bool ConvexHull::CleanVertices (int currentVertexID)
{
  bool need_to_delete_current = false;

  // mark all vertices incident to some undeleted edge as on the hull
  for (EdgeMapIt e = sc.edges.begin(); e != sc.edges.end(); e++) {
    onHull.insert(e->second.vertex[0]);
    onHull.insert(e->second.vertex[1]);
  }
	
  // delete all vertices that have been processed but are not on the hull
  // reset flags for the rest
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); ) {
    if (false && debug)
      cerr << "CleanVertices considering: " << v->first << "\n";
    if (!isProcessed(v->second) || isOnHull(v->second)) {
      v++;
      if (false && debug) cerr << "   onhull or unprocessed\n";
    }
    else if (currentVertexID != v->first) {
      if (debug) cerr << "   CleanVertices erasing: " << v->first << "\n";
      sc.vertices.erase(v++);
    }
    else {
      if (debug) cerr << "   CleanVertices need to delete current: " 
											<< v->first << "\n";
      need_to_delete_current = true;
      v++;
    }
  }

  duplicate.clear();
  for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++)
    onHull.insert(&it->second);

  return need_to_delete_current;
}

/*---------------------------------------------------------------------
	Consistency runs through the edge list and checks that all
	adjacent faces have their endpoints in opposite order.  This verifies
	that the vertices are in counterclockwise order.
	---------------------------------------------------------------------*/
void SimplicialComplex::Consistency ()
{
  EdgeMapIt it = edges.begin();
  for (; it != edges.end(); it++) {
    Edge &e = it->second;

    if (e.face[0] == NULL) continue;
    if (e.face[1] == NULL) continue;

    int i, j;

    /* find index of endpoint[0] in adjacent face[0] */
    for (i = 0; e.face[0]->vertex[i] != e.vertex[0]; ++i )
      ;

    /* find index of endpoint[0] in adjacent face[1] */
    for (j = 0; e.face[1]->vertex[j] != e.vertex[0]; ++j )
      ;

    /* check if the endpoints occur in opposite order */
    if ( !( e.face[0]->vertex[ (i+1) % 3 ] ==
						e.face[1]->vertex[ (j+2) % 3 ] ||
						e.face[0]->vertex[ (i+2) % 3 ] ==
						e.face[1]->vertex[ (j+1) % 3 ] )  )
      break;
  }
  if (it != edges.end())
    cerr << "Checks: edges are NOT consistent.\n";
  else
    cerr << "Checks: edges consistent.\n";
}

// checks that every vertex is "behind" all the faces:
//    - this doesn't guarantee full convexity
//    - things like one big face with a bunch of vertices
//      all on the proper side pass this test
//    - you also need to guarantee sphere-ness
//    - you maybe also need to guarantee connectedness
bool ConvexHull::isConvex (const SimplicialComplex &sc,
													 const double volumeSignTolerance)
{
  for (const_FaceMapIt f = sc.faces.begin(); f != sc.faces.end(); f++) {
    for (const_VertexMapIt v = sc.vertices.begin();
				 v != sc.vertices.end(); v++) {
      double vol = f->second.volume(v->second.pos);
      if (vol < -volumeSignTolerance)
				return false;
    }
  }
  return sc.faces.size() > 0;
}

/*---------------------------------------------------------------------
	CheckEuler checks Euler's relation, as well as its implications when
	all faces are known to be triangles.  Only prints positive information
	when debug is true, but always prints negative information.
	---------------------------------------------------------------------*/
void SimplicialComplex::CheckEuler (int V, int E, int F)
{
  cerr << "Checks: V, E, F = " << V << " " << E << " " << F << ":\t";

  if ( (V - E + F) != 2 )
    cerr << "Checks: V-E+F != 2\n";
  else
    cerr << "V-E+F = 2\t";

  if ( F != (2 * V - 4) )
    cerr << "Checks: F=" << F << " != 2V-4=" <<(2*V-4) << "; V=" << V << "\n";
  else
    cerr << "F = 2V-4\t";
   
  if ( (2 * E) != (3 * F) )
    cerr << "Checks: 2E=" << (2*E) << " != 3F=" << (3*F) 
				 << "; E=" << E << ", F=" << F<< "\n";
  else
    cerr << "2E = 3F\n";
}

/*-------------------------------------------------------------------
	Checks that, for each face, for each i={0,1,2}, the [i]th vertex of
	that face is either the [0]th or [1]st endpoint of the [ith] edge of
	the face.
	-------------------------------------------------------------------*/
void SimplicialComplex::CheckEndpts ()
{
	bool error = false;
	if (Face::numVertices != Face::numEdges)
		throw "bogus Face::numEdges/numVertices";
	const int lim = Face::numVertices;
	for (FaceMapIt f = faces.begin(); f != faces.end(); f++) {
		for(int i = 0; i < lim; ++i) {
			Vertex *v = f->second.vertex[i];
			Edge *e = f->second.edge[i];
			if (v != e->vertex[0] && v != e->vertex[1] ) {
				error = true;
				cerr << "CheckEndpts: Error!\n";
				cerr << "edges : (" 
						 << e->vertex[0]->objID << ", " 
						 << e->vertex[1]->objID << ")\n";
			}
		}
	}

	if (error)
		cerr << "Checks: ERROR found and reported above.\n";
	else
		cerr << "Checks: All vertex of all edges of all faces check.\n";
}

/*------------------------------------------------------------------
  EdgeOrderOnFaces: puts e0 between v0 and v1, e1 between v1 and v2,
  e2 between v2 and v0 on each face.  This should be unnecessary, alas.
  Not used in code, but useful for other purposes.
	------------------------------------------------------------------*/
void SimplicialComplex::EdgeOrderOnFaces ()
{
  for (FaceMapIt f = faces.begin(); f != faces.end(); f++) {
    for (int i = 0; i < Face::numEdges; i++) {
      if (!(((f->second.edge[i]->vertex[0] == f->second.vertex[i]) &&
             (f->second.edge[i]->vertex[1] == f->second.vertex[(i+1)%3])) ||
            ((f->second.edge[i]->vertex[1] == f->second.vertex[i]) &&
             (f->second.edge[i]->vertex[0] == f->second.vertex[(i+1)%3])))) {
				cerr << "jive edge order\n";
        /* Change the order of the edges on the face: */
        for (int j = 0; j < Face::numEdges; j ++) {
          /* find the edge that should be there */
          if (((f->second.edge[j]->vertex[0] == f->second.vertex[i]) &&
               (f->second.edge[j]->vertex[1] == f->second.vertex[(i+1)%3])) ||
              ((f->second.edge[j]->vertex[1] == f->second.vertex[i]) &&
               (f->second.edge[j]->vertex[0] == f->second.vertex[(i+1)%3]))) {
            /* Swap it with the one erroneously put into its place: */
            Edge *new_edge = f->second.edge[i];
            f->second.edge[i] = f->second.edge[j];
            f->second.edge[j] = new_edge;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------
	DoubleTriangle builds the initial double triangle.  It first finds 3 
	noncollinear points and makes two faces out of them, in opposite order.
	It then finds a fourth point that is not coplanar with that face.  The  
	vertices are stored in the face structure in counterclockwise order so 
	that the volume between the face and the point is negative. Lastly, the
	3 new faces to the fourth point are constructed and the data structures
	are cleaned up. 
	---------------------------------------------------------------------*/
Vertex *ConvexHull::DoubleTriangle ()
{
  if (debug) cerr << "DoubleTriangle() called\n";

  if (debug) cerr << "DoubleTriangle() part 1\n";

  // find 3 noncollinear points:
  vector<int> arr;
  for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++)
    arr.push_back(it->first);

  Vertex   *v0, *v1, *v2;
  unsigned int i = 0;
  for (; i+2 < sc.vertices.size(); i++) {
    v0 = &sc.vertices[arr[i+0]];
    v1 = &sc.vertices[arr[i+1]];
    v2 = &sc.vertices[arr[i+2]];
    if (!collinear(v0->pos,v1->pos,v2->pos)) break;
  }
  if (i+2 >= sc.vertices.size()) return NULL;

  if (debug) cerr << "DoubleTriangle() part 2\n";

  /* Mark the vertices as processed. */
  processed.insert(v0);
  processed.insert(v1);
  processed.insert(v2);
   
  if (debug) cerr << "DoubleTriangle() part 3\n";

  /* Create the two "twin" faces */
  Face &f0 = sc.MakeFace(*base, v0, v1, v2);
  Face &f1 = sc.MakeFace(f0,    v2, v1, v0);

  if (debug) cerr << "DoubleTriangle() part 4\n";

  /* Link adjacent face fields. */
  f0.edge[0]->face[1] = &f1;
  f0.edge[1]->face[1] = &f1;
  f0.edge[2]->face[1] = &f1;
  f1.edge[0]->face[1] = &f0;
  f1.edge[1]->face[1] = &f0;
  f1.edge[2]->face[1] = &f0;
	
  if (debug) cerr << "DoubleTriangle() part 5\n";

  /* Find a fourth, noncoplanar point to form tetrahedron. */
  /* Put it at front of vertices to ensure it will be the first added. */
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    if (debug) cerr << "DoubleTriangle() part 5a\n";
    if (isProcessed(v->second)) continue;
    if (debug) cerr << "DoubleTriangle() part 5b\n";
    double vol = VolumeSign(f0, v->second);
    if (debug) cerr << "DoubleTriangle() part 5c\n";
    if (vol == 0.0) continue;
    if (debug) cerr << "DoubleTriangle() done\n";
    return &v->second;
  }

  if (debug) cerr << "DoubleTriangle:  All points are coplanar!\n";
  return NULL;
}

void ConvexHull::FinishHull ()
{
  sc.EdgeOrderOnFaces();

  // clean out any vertices not on edges:
  hash_set<Vertex *> alive;
  for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++) {
    alive.insert(it->second.vertex[0]);
    alive.insert(it->second.vertex[1]);
  }
  hash_set<Vertex *> dead;
  for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++) {
    if (alive.find(&it->second) != alive.end()) continue;
    dead.insert(&it->second);
  }
  for (hash_set<Vertex *>::iterator it = dead.begin(); 
       it != dead.end(); it++)
    sc.vertices.erase(sc.vertices.find((*it)->objID));
}

