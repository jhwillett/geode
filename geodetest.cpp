// geodetest.cpp
//
// automated test suite for geode project
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

#include "simplicial.hpp"
#include "hull.hpp"
#include "io.hpp"
#include "integrator.hpp"
#include <math.h>
#include <time.h>
#include <sstream>
#include <string>
#include <vector>
using std::cerr;
using std::cin;
using std::cout;
using std::ostringstream;
using std::istringstream;
using std::string;
using std::vector;

void testIOStability (const SimplicialComplex &sc);
void testIOStability (const Material &m);
void testConvexHullness (const SimplicialComplex &sc);

bool operator== (const SimplicialComplex &sc0, const SimplicialComplex &sc1);
bool operator== (const Physics &p0, const Physics &p1);
bool operator== (const Face &f0, const Face &f1);
bool operator== (const Edge &e0, const Edge &e1);
bool operator== (const Material &m0, const Material &m1);
bool operator== (const Vertex &v0, const Vertex &v1);
bool operator== (const Vector &v0, const Vector &v1);

bool operator!= (const SimplicialComplex &sc0, const SimplicialComplex &sc1)
{
  return !(sc0==sc1);
}
bool operator!= (const Physics &p0, const Physics &p1)
{
  return !(p0==p1);
}
bool operator!= (const Face &f0, const Face &f1)
{
  return !(f0==f1);
}
bool operator!= (const Edge &e0, const Edge &e1)
{
  return !(e0==e1);
}
bool operator!= (const Material &m0, const Material &m1)
{
  return !(m0==m1);
}
bool operator!= (const Vertex &v0, const Vertex &v1)
{
  return !(v0==v1);
}
bool operator!= (const Vector &v0, const Vector &v1)
{
  return !(v0==v1);
}

void testConvexHullness (const SimplicialComplex &sc)
{
  // NOTE: passing Euler's formula does not guarantee spherical genus
  int V = sc.vertices.size();
  int E = sc.edges.size();
  int F = sc.faces.size();
  int euler = V-E+F;
  if (euler != 2) {
    cerr << "V:     " << V << "\n";
    cerr << "E:     " << E << "\n";
    cerr << "F:     " << F << "\n";
    cerr << "V-E+F: " << euler << "\n";
    throw "not spherical genus";
  }

  // NOTE: passing isConvex() doesn't guarantee convexity
  if (!ConvexHull::isConvex(sc))
    throw "not convex";

  // TODO: some kind of connectedness test?  isManifold()?
}

void geodetest ()
{
  // some commands take random seed: use time
  //srandom((unsigned int)time(NULL));
  srandom(0);

  testVectors();

  {
    cerr << "testing raw Material io\n";
    Material m0;
    testIOStability(m0);
    Material m1 (37,100.0,101.0,false,0.13254);
    testIOStability(m1);
  }

  {
    cerr << "testing empty construction and destruction\n";
    SimplicialComplex sc;
  }

  {
    cerr << "testing empty w/ io stability\n";
    SimplicialComplex sc;
    testIOStability(sc);
  }

  {
    cerr << "testing unit cube area, etc\n";
    SimplicialComplex sc;
    double size = 1.0;
    double tolerance = 0.0000001;
    for (int i = 0; i < 7; i++) {
      sc.MakeUnitCube(size);
      ConvexHull(sc).ConstructHull();
      sc.physics.CalculateEpiphenomena();
      double expLength = 12*size + 6*size*sqrt(2);
      double errLength = sc.physics.totalLength - expLength;
      double expArea   = 6*size*size;
      double errArea   = sc.physics.totalArea - expArea;
      double expVolume = size*size*size;
      double errVolume = sc.physics.totalVolume - expVolume;
      if (false) {
				cerr << "   size: " << size << "\n";
				cerr << "      length: " << sc.physics.totalLength << "\n";
				cerr << "         exp: " << expLength << "\n";
				cerr << "         err: " << errLength << "\n";
				cerr << "      area:   " << sc.physics.totalArea << "\n";
				cerr << "         exp: " << expArea << "\n";
				cerr << "         err: " << errArea << "\n";
				cerr << "      volume: " << sc.physics.totalVolume << "\n";
				cerr << "         exp: " << expVolume << "\n";
				cerr << "         err: " << errVolume << "\n";
      }
      if (fabs(errLength) > tolerance * expLength) throw "bunk length";
      if (fabs(errArea)   > tolerance * expArea)   throw "bunk area";
      if (fabs(errVolume) > tolerance * expVolume) throw "bunk volume";
      size *= 2.21345; // random-ish number for noise in check
    }
  }

  {
    cerr << "testing tetra io stability\n";
    SimplicialComplex sc;
    sc.MakeTetrahedron();
    testIOStability(sc);
  }

  {
    cerr << "testing pin2\n";
    SimplicialComplex sc;
    sc.MakeTwicePinnedWeight();
    testIOStability(sc);
  }

  {
    cerr << "testing pin2 | bisect 5\n";
    SimplicialComplex sc;
    sc.MakeTwicePinnedWeight();
    sc.BisectEdges();
    sc.BisectEdges();
    sc.BisectEdges();
    sc.BisectEdges();
    sc.BisectEdges();
    testIOStability(sc);
  }

  {
    cerr << "testing icosa\n";
    SimplicialComplex sc;
    sc.MakeIcosahedron();
    testIOStability(sc);
  }

  {
    cerr << "testing tetra hull\n";
    SimplicialComplex sc;
    sc.MakeTetrahedron();
    testIOStability(sc);
    ConvexHull(sc).ConstructHull();
    testConvexHullness(sc);
    testIOStability(sc);
  }

  {
    cerr << "testing octa hull\n";
    SimplicialComplex sc;
    sc.MakeOctahedron();
    ConvexHull(sc).ConstructHull();
    testConvexHullness(sc);
    testIOStability(sc);
  }

  {
    cerr << "testing icosa hull\n";
    SimplicialComplex sc;
    sc.MakeIcosahedron();
    ConvexHull(sc).ConstructHull();
    testConvexHullness(sc);
    testIOStability(sc);
  }

  if (true) {
    cerr << "testing torus left 10 10\n";
    SimplicialComplex sc;
    sc.MakeTorus(SimplicialComplex::left,10,10);
    testIOStability(sc);
    sc.AddCrossScaffolding(true);
    testIOStability(sc);
  }

  if (true) {
    cerr << "testing subdivide 4 hull\n";
    SimplicialComplex sc;
    sc.MakeIcosahedron();
    ConvexHull(sc).ConstructHull();
    testConvexHullness(sc);
    sc.SubdivideFaces();
    sc.SubdivideFaces();
    sc.AddCrossScaffolding(true);
    sc.SubdivideFaces();
    sc.SubdivideFaces();
    testIOStability(sc);
  }

  vector<int> sizes;
  sizes.push_back(50);
  //sizes.push_back(300);
  //sizes.push_back(500);
  for (size_t i = 0; i < sizes.size(); i++) {
    if (true) {
      cerr << "testing cube " << sizes[i] << " hull\n";
      SimplicialComplex sc;
      sc.MakeCube(sizes[i]);
      ConvexHull(sc).ConstructHull();
      testConvexHullness(sc);
      testIOStability(sc);
    }
    if (true) {
      cerr << "testing sphere " << sizes[i] << " hull\n";
      SimplicialComplex sc;
      sc.MakeSphere(sizes[i]);
      ConvexHull(sc).ConstructHull();
      testConvexHullness(sc);
      testIOStability(sc);
    }
    if (true) {
      cerr << "testing spiral " << sizes[i] << " hull\n";
      SimplicialComplex sc;
      sc.MakeSpiral(sizes[i]);
      ConvexHull(sc).ConstructHull();
      testConvexHullness(sc);
      testIOStability(sc);
    }
  }

  if (true) {
    cerr << "testing big nasty spaceship SC\n";
    SimplicialComplex sc;
    sc.MakeIcosahedron();
    sc.DiamondSubdivide();
    ConvexHull(sc).ConstructHull();
    sc.AddCrossScaffolding(true);
    sc.AddDualScaffolding();
    sc.SubdivideFaces();
    sc.SubdivideFaces();
    sc.SubdivideFaces();
    testIOStability(sc);
    cerr << "testing big nasty spaceship Forward Euler io stability\n";
    int numSteps = 3;
    double dt = dt;
    for (int i = 0; i < numSteps; i++)
      ForwardEuler(sc).step(dt);
    testIOStability(sc);
    cerr << "testing big nasty spaceship Improved Euler io stability\n";
    for (int i = 0; i < numSteps; i++)
      ImprovedEuler(sc).step(dt);
    testIOStability(sc);
    cerr << "testing big nasty spaceship Runge-Kutta euler io stability\n";
    for (int i = 0; i < numSteps; i++)
      RungeKutta(sc).step(dt);
    testIOStability(sc);
  }

  cerr << "done\n";
}

void testIOStability (const Material &m)
{
  string str1;
  {
    ostringstream o;
    o << m;
    str1 = o.str();
  }
  Material m1;
  string str2;
  {
    istringstream i(str1);
    ostringstream o;
    i >> m1;
    o << m1;
    str2 = o.str();
  }
  if (m != m) {
    cerr << "equality instability: m != m\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "equality instability: m != m";
  }
  if (m1 != m1) {
    cerr << "equality instability: m1 != m1\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "equality instability: m1 != m1";
  }
  if (m != m1) {
    cerr << "io instability: m != m1\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "io instability: m != m1";
  }
  if (m1 != m) {
    cerr << "io instability: m1 != m\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "io instability: m1 != m";
  }
  if (true) {
    if (str1 != str2) {
      if (false) {
				cerr << str1;
				cerr << str2;
      }
      cerr << "io instability initial iteration (hash_map?)\n";
      //throw "io instability initial iteration";
    }
  }
}

void testIOStability (const SimplicialComplex &sc)
{
  string str1;
  {
    ostringstream o;
    o << sc;
    str1 = o.str();
  }
  SimplicialComplex sc1;
  string str2;
  {
    istringstream i(str1);
    ostringstream o;
    i >> sc1;
    o << sc1;
    str2 = o.str();
  }
  if (sc != sc) {
    cerr << "equality instability: sc != sc\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "equality instability: sc != sc";
  }
  if (sc1 != sc1) {
    cerr << "equality instability: sc1 != sc1\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "equality instability: sc1 != sc1";
  }
  if (sc != sc1) {
    cerr << "io instability: sc != sc1\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "io instability: sc != sc1";
  }
  if (sc1 != sc) {
    cerr << "io instability: sc1 != sc\n";
    if (true) {
      cerr << str1 << "\n";
      cerr << str2 << "\n";
    }
    throw "io instability: sc1 != sc";
  }
  if (true) {
    if (str1 != str2) {
      if (false) {
				cerr << str1;
				cerr << str2;
      }
      cerr << "io instability initial iteration (hash_map?)\n";
      //throw "io instability initial iteration";
    }
  }
}

bool operator== (const Physics &p0, const Physics &p1)
{
  if (p0.viscosityActive != p1.viscosityActive) return false;
  if (p0.viscosity != p1.viscosity) return false;
  if (p0.centeringActive != p1.centeringActive) return false;
  if (p0.centering != p1.centering) return false;
  if (p0.springsActive != p1.springsActive) return false;
  if (p0.springDampingActive != p1.springDampingActive) return false;
  if (p0.gravityActive != p1.gravityActive) return false;
  if (p0.gravity != p1.gravity) return false;
  if (p0.internalPressureActive != p1.internalPressureActive) return false;
  if (p0.moles != p1.moles) return false;
  if (p0.internalTemperature != p1.internalTemperature) return false;
  if (p0.externalPressureActive != p1.externalPressureActive) return false;
  if (p0.externalPressure != p1.externalPressure) return false;
  if (p0.externalTemperature != p1.externalTemperature) return false;
  if (p0.vertexLimitActive != p1.vertexLimitActive) return false;
  if (p0.vertexLimit != p1.vertexLimit) return false;
  if (p0.currentGamma != p1.currentGamma) return false;
  if (p0.possibleGammas.size() != p1.possibleGammas.size()) return false;
  for (size_t i = 0; i < p0.possibleGammas.size(); i++) {
    if (p0.possibleGammas[i] != p1.possibleGammas[i]) return false;
  }

  return true;
}

bool operator== (const SimplicialComplex &sc0, const SimplicialComplex &sc1)
{
  if (sc0.currentObjID != sc1.currentObjID) return false;
  if (sc0.physics != sc1.physics) return false;

  if (true) {
    for (const_VertexMapIt it = sc0.vertices.begin(); 
				 it != sc0.vertices.end(); it++) {
      if (sc1.vertices.find(it->first) == sc1.vertices.end()) return false;
      if (sc1.vertices.find(it->first)->second != it->second) return false;
    }
    for (const_VertexMapIt it = sc1.vertices.begin(); 
				 it != sc1.vertices.end(); it++) {
      if (sc0.vertices.find(it->first) == sc0.vertices.end()) return false;
      if (sc0.vertices.find(it->first)->second != it->second) return false;
    }
  }

  if (true) {
    for (const_MaterialMapIt it = sc0.materials.begin(); 
				 it != sc0.materials.end(); it++) {
      if (sc1.materials.find(it->first) == sc1.materials.end()) {
				cerr << "failed to find Material id: " << it->first << "\n";
				return false;
      }
      if (sc1.materials.find(it->first)->second != it->second) return false;
    }
    for (const_MaterialMapIt it = sc1.materials.begin(); 
				 it != sc1.materials.end(); it++) {
      if (sc0.materials.find(it->first) == sc0.materials.end()) {
				cerr << "failed to find Material id: " << it->first << "\n";
				return false;
      }
      if (sc0.materials.find(it->first)->second != it->second) return false;
    }
  }

  if (true) {
    for (const_EdgeMapIt it = sc0.edges.begin(); 
				 it != sc0.edges.end(); it++) {
      if (sc1.edges.find(it->first) == sc1.edges.end()) return false;
      if (sc1.edges.find(it->first)->second != it->second) return false;
    }
    for (const_EdgeMapIt it = sc1.edges.begin(); 
				 it != sc1.edges.end(); it++) {
      if (sc0.edges.find(it->first) == sc0.edges.end()) return false;
      if (sc0.edges.find(it->first)->second != it->second) return false;
    }
  }

  if (true) {
    for (const_FaceMapIt it = sc0.faces.begin(); 
				 it != sc0.faces.end(); it++) {
      if (sc1.faces.find(it->first) == sc1.faces.end()) return false;
      if (sc1.faces.find(it->first)->second != it->second) return false;
    }
    for (const_FaceMapIt it = sc1.faces.begin(); 
				 it != sc1.faces.end(); it++) {
      if (sc0.faces.find(it->first) == sc0.faces.end()) return false;
      if (sc0.faces.find(it->first)->second != it->second) return false;
    }
  }

  return true;
}

// handy helper for checking NULL-or-objID match
template<class A> bool checkID (const A *a0, const A *a1)
{
  if ((a0 == NULL) != (a1 == NULL)) return false;
  if ((a0 == NULL) || (a1 == NULL)) return true;
  return (a0->objID == a1->objID);
}

bool operator== (const Face &f0, const Face &f1)
{
  if (f0.objID != f1.objID) return false;

  for (int i = 0; i < Face::numVertices; i++)
    if (!checkID(f0.vertex[i],f1.vertex[i])) return false;
  for (int i = 0; i < Face::numEdges; i++)
    if (!checkID(f0.edge[i],f1.edge[i])) return false;

  return true;
}

bool operator== (const Edge &e0, const Edge &e1)
{
  if (e0.objID != e1.objID) return false;
  if (e0.springFactor != e1.springFactor) return false;

  if (!checkID(e0.material,e1.material)) return false;
  for (int i = 0; i < Edge::numVertices; i++)
    if (!checkID(e0.vertex[i],e1.vertex[i])) return false;
  for (int i = 0; i < Edge::numFaces; i++)
    if (!checkID(e0.face[i],e1.face[i])) return false;
  return true;
}

bool operator== (const Material &m0, const Material &m1)
{
  if (m0.objID != m1.objID) return false;
  if (m0.restLength != m1.restLength) return false;
  if (m0.springCoef != m1.springCoef) return false;
  if (m0.damper != m1.damper) return false;
  if (m0.compressible != m1.compressible) return false;
  return true;
}

bool operator== (const Vertex &v0, const Vertex &v1)
{
  if (v0.objID != v1.objID) return false;
  if (v0.pos != v1.pos) return false;
  if (v0.vel != v1.vel) return false;
  if (v0.pinned != v1.pinned) return false;
  if (v0.mass != v1.mass) return false;
  return true;
}
bool operator== (const Vector &v0, const Vector &v1)
{
  return (v0.x == v1.x && v0.y == v1.y && v0.z == v1.z);
}

