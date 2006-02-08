// simplicial.hpp
// 
// basic geometry data structures and some constructive algorithms
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

#ifndef SIMPLICIAL_HPP
#define SIMPLICIAL_HPP

// ugly stuff:
#ifndef NULL
#define NULL 0;
#endif

#include "algebra.hpp"
#include "physics.hpp"

#include <ext/hash_set>
#include <ext/hash_map>
#include <vector>
#include <list>
#include <map>

class Vertex;
class Edge;
class Face;
class Material;
class SimplicialComplex;

class Vertex
{
public:
  // saved fields
  int     objID;        // unattractive tagging for read-write support
  Vector  pos;          // position of the vertex in meters
  Vector  vel;          // motion of the vertex in meters per second
  double  mass;         // scalar on force for delta-v in kilograms
  bool    pinned;       // true iff the vertex cannot move

  // unsaved fields: support for numeric integrators
  Vector  ni_force;      // accumulated force in newtons
  Vector  ni_vel;        // accumulated velocity in meters per seond
  Vector  ni_k1;         // for ImprovedEuler and RungeKutta
  Vector  ni_k2;         // for ImprovedEuler and RungeKutta
  Vector  ni_k3;         // for RungeKutta
  Vector  ni_k4;         // for RungeKutta

  Vertex (int objIDi = -1);
};

class Material
{
public:
  // saved fields
  int     objID;        // unattractive tagging for read-write support
  double  restLength;   // rest length of this material in meters
  double  springCoef;   // spring coef of this material in kg/s^2 = N/m
  double  damper;       // damping coef of this material in kg/s  = Ns/m
  bool    compressible; // true iff this material resists compression

  Material (int objIDi = -1, 
						double restLengthi = 0.0,
						double springCoefi = 0.0,
						bool compressiblei = false, 
						double damperi = 0.0);
};

class Edge
{
public:
  enum { 
    numVertices = 2,  // fundamental to concept of Edge
    numFaces = 2      // soon to become variable
  };

  // saved fields
  int             objID;         // unattractive tagging for read-write support
  Material        *material;
  Vertex          *vertex[numVertices];
  Face            *face[numFaces];
  double          springFactor;  // how much "spring" there is in ???

  // unsaved fields
  int faceID[numFaces]; // unfortunate read() assist for circularity of ref

  bool isCompressed () const;
  double length () const;

  Edge (Material *m = NULL, int objIDi = -1);
};

class Face
{
public:
  enum { 
    numVertices = 3,  // fundamental to concept of Face
    numEdges = 3      // fundamental to concept of Face
  };

  // saved fields
  int     objID;        // unattractive tagging for read-write support
  Vertex  *vertex[numVertices];
  Edge    *edge[numEdges];

  Face (int objIDi = -1);

  double volume (const Vector &v) const;  // volume of tetra Face + v
  double area   () const;                 // area of Face
  Vector areaNormal () const;             // normal to Face, mag area of Face
  Vector unitNormal () const;             // normal to Face, mag 1
  Vector center () const;                 // average over vertex

  // find the vertex on this not on e:
  //   - result undefined in iff e not in edges:
  Vertex &vertexOpposite (const Edge &e);
  const Vertex &vertexOpposite (const Edge &e) const;
};

namespace __gnu_cxx
{
  template<> struct hash<Vertex *>
  {
    size_t operator()(Vertex *x) const { return reinterpret_cast<int>(x); }
  };
  template<> struct hash<Edge *>
  {
    size_t operator()(Edge *x) const { return reinterpret_cast<int>(x); }
  };
  template<> struct hash<Face *>
  {
    size_t operator()(Face *x) const { return reinterpret_cast<int>(x); }
  };
}

typedef __gnu_cxx::hash_map<int,Vertex> VertexMap;
typedef __gnu_cxx::hash_map<int,Material> MaterialMap;
typedef std::map<int,Face> FaceMap;  // NOTE: hash_map breaks Make*() and hull
typedef std::map<int,Edge> EdgeMap;  // NOTE: hash_map breaks Make*() and hull
typedef VertexMap::iterator VertexMapIt;
typedef EdgeMap::iterator EdgeMapIt;
typedef FaceMap::iterator FaceMapIt;
typedef MaterialMap::iterator MaterialMapIt;
typedef VertexMap::const_iterator const_VertexMapIt;
typedef EdgeMap::const_iterator const_EdgeMapIt;
typedef FaceMap::const_iterator const_FaceMapIt;
typedef MaterialMap::const_iterator const_MaterialMapIt;

class SimplicialComplex
{
public:
  // saved: basic geometric data
  int         currentObjID;
  VertexMap   vertices;
  EdgeMap     edges;
  FaceMap     faces;
  MaterialMap materials;
  Physics     physics;

  SimplicialComplex ();

  // various correctness checks:
  void	  EdgeOrderOnFaces  ();
  void    CheckEuler(int V, int E, int F);
  void	  Consistency ();
  void	  CheckEndpts  ();
  void	  myConsistency ();

  // some geometric manipulations:
  void   InvertThroughOrigin ();
  void   CenterOnOrigin ();
  void   BisectEdges ();
  void   SubdivideFaces ();
  void   DiamondSubdivide ();
  void   AddCrossScaffolding (bool copyOriginal);
  void   AddDualScaffolding ();
  void   Sphericize ();
  void   DoublePos ();
  void   HalvePos ();
  bool   Equiangulate ();

  // some geometric/physics manipulators:
  void   Harmonize   ();
  void   Harmonize   (Material &material);
  void   DeHarmonize ();
  void   DeHarmonize (Material &material);

  // some geometric constructions:
  enum   TorusType { left = 0, right = 1, modi = 2, modj = 3, modij = 4 };
  void   MakeTorus (TorusType type, int n, int m);
  void   MakePinnedWeight ();
  void   MakeTwicePinnedWeight ();
  void   MakeBinaryWeight ();
  void   MakeTrinaryWeight ();
  void   MakeTetrahedron ();
  void   MakeOctahedron ();
  void   MakeIcosahedron ();
  void   MakeSpiral (int num);
  void   MakeSphere (int num);
  void   MakeCube (int num);
  void   MakeSuspensionBridge (size_t numTowers, size_t sectionsPerSpan);
  void   MakeTensegrityFoo ();
  void   MakePiston ();
  void   MakeUnitCube (double scale);

public:
  // stuff public for ConvexHull
  Vertex   &MakeVertex (const Vector &v);
  Face     &MakeNullFace ();
  Material &MakeMaterial (double restLengthi,
													double springCoefi,
													bool compressiblei,
													double damperi);
  Edge     &MakeEdge (Material &material, Vertex &a, Vertex &b);
  Face     &MakeFace (Material &material, Vertex *v0, Vertex *v1, Vertex *v2);
  Face     &MakeFace (Face &fold,         Vertex *v0, Vertex *v1, Vertex *v2);
  
  // stuff public for ConvexHull & IO
  void    clear ();

private:
  typedef __gnu_cxx::hash_map<Edge *,Vertex *>                  BisectorMap;
  typedef std::list<Edge *>                                     ChildList;
  typedef __gnu_cxx::hash_set<Edge *>                           ChildSet;
  typedef __gnu_cxx::hash_map<Edge *,std::pair<Edge *,Edge *> > ChildMap;
  void    bisect (Edge &e,
									BisectorMap &bisectors,
									ChildMap    &childMap,
									ChildSet    &childSet);
  void    nsect (Edge         &e,
								 size_t       num,
								 ChildList    &childList);
  void    DivisionCleanup (double oldKE, 
													 int numParents,
													 int oldNumV,
													 BisectorMap &bisectors,
													 ChildSet    &childSet);
};

#endif // #ifndef SIMPLICIAL_HPP
