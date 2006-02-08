// hull.hpp
// 
// An object embodying the convex hull algorithm, exposing it
// for pretty animation, and sequestering its special-purpose
// data from the main SimplicialComplex geometric data.
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

#ifndef HULL_HPP
#define HULL_HPP

// ugly stuff:
#ifndef NULL
#define NULL 0;
#endif

#include "simplicial.hpp"

class ConvexHull
{
  SimplicialComplex &sc;
  __gnu_cxx::hash_set<Vertex *> processed;     // Vertices already processed
  __gnu_cxx::hash_set<Vertex *> onHull;        // Vertices currently on hull
  __gnu_cxx::hash_map<Vertex *,Edge *> duplicate; // from Vertices to new Edges
  __gnu_cxx::hash_set<Face *> visible;         // Faces visible from current pt
  __gnu_cxx::hash_set<Edge *> deadEdges;       // Edges scheduled for removal
  __gnu_cxx::hash_map<Edge *,Face *> newFaces; // boundary Edges' new Faces
  VertexMapIt hullIt;
  Material *base;
  const double volumeSignTolerance;

public:
  ConvexHull (SimplicialComplex &sci, double volumeSignTolerancei = 0.1);

  // ConstructHull() runs the whole algorithm. 
  // The full method definition is included here to demonstrate 
  // how to unpeel the algorithm for, say, animation or debugging.
  void   ConstructHull ()
  {
    StartHull();
    while (isAdvancingHull()) AdvanceHull();
    FinishHull();
  }

  // StartHull() through FinishHull() expose partial steps in the algorithm
  void StartHull ();
  bool isAdvancingHull ();
  void AdvanceHull ();
  void FinishHull ();

  // exposed for pretty rendering
  bool   isProcessed (const Vertex &v) const;
  bool   isOnHull (const Vertex &v) const;
  bool   isVisible (const Face &f) const;

  // exposed for checks
  static bool isConvex (const SimplicialComplex &sc,
												const double volumeSignTolerance = 0.1);

private:
  void    MarkVisiblity (Vertex &p);
  bool	  AddOne (Vertex &p);

  Vertex  *DoubleTriangle ();  // returns next one to add
  double  VolumeSign(const Face &f, const Vertex &p);
  Face	  &MakeConeFace (Edge &e, Vertex &p);
  void    MakeCcw (Face &f, Edge &e, Vertex &p);
  bool    CleanUp (int currentVertexID);
  void    CleanEdges ();
  void    CleanFaces ();
  bool    CleanVertices (int currentVertexID);
};

#endif // #ifndef HULL_HPP
