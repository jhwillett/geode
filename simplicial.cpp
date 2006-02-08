// simplicial.cpp
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

#include "simplicial.hpp"
#include <set>
#include <list>
#include <iostream>
#include <ext/hash_set>

using std::set;
using std::list;
using __gnu_cxx::hash_set;
using __gnu_cxx::hash_map;
using std::pair;
using std::make_pair;
using std::vector;
using std::cerr;

const static bool debug = false;
const double PI = 3.14159265358979;

SimplicialComplex::SimplicialComplex () :
  currentObjID(0),
  physics(*this)
{
}

Vertex &SimplicialComplex::MakeVertex (const Vector &v)
{
  int id = currentObjID++;
  Vertex &vert = (vertices[id] = Vertex(id));
  vert.pos = v;
  return vert;
}

Edge &SimplicialComplex::MakeEdge (Material &material, Vertex &a, Vertex &b)
{
  int id = currentObjID++;
  Edge &e = (edges[id] = Edge(&material,id));
  e.vertex[0] = &a;
  e.vertex[1] = &b;
  return e;
}

Face &SimplicialComplex::MakeNullFace ()
{
  int id = currentObjID++;
  faces[id] = Face(id);
  return faces[id];
}

Material &SimplicialComplex::MakeMaterial (double restLengthi,
																					 double springCoefi,
																					 bool compressiblei,
																					 double damperi)
{
  int id = currentObjID++;
  Material &m = (materials[id] = Material(id,
																					restLengthi,
																					springCoefi,
																					compressiblei,
																					damperi));
  return m;
}

Vertex::Vertex (int objIDi) : 
  pos(0,0,0),
  vel(0,0,0),
  ni_force(0,0,0),
  ni_vel(0,0,0),
  ni_k1(0,0,0),
  ni_k2(0,0,0),
  ni_k3(0,0,0),
  ni_k4(0,0,0)
{
  objID = objIDi;
  mass = 1.0;
  pinned = false;
}

Material::Material (int objIDi,
										double restLengthi,
										double springCoefi,
										bool compressiblei,
										double damperi)
{
  objID = objIDi;
  restLength = restLengthi;
  springCoef = springCoefi;
  compressible = compressiblei;
  damper = damperi;
}

Edge::Edge (Material *m, int objIDi)
{
  for (int i = 0; i < Edge::numFaces; i++) {
    face[i] = NULL;
    faceID[i] = -1;
  }
  for (int i = 0; i < Edge::numVertices; i++) {
    vertex[i] = NULL;
  }
  material = m;
  objID = objIDi;
  springFactor = 1.0;
}

bool Edge::isCompressed () const
{
  const Vertex &v0 = *vertex[0];
  const Vertex &v1 = *vertex[1];
  const double r = material->restLength * springFactor;
  const Vector L = v0.pos - v1.pos;
  const double l = L.magnitude();
  return l < r;
}

double Edge::length () const
{
  return (vertex[1]->pos - vertex[0]->pos).magnitude();
}

Face::Face (int objIDi)
{
  for (int i = 0; i < Face::numEdges; i++)
    edge[i] = NULL;
  for (int i = 0; i < Face::numVertices; i++)
    vertex[i] = NULL;
  objID = objIDi;
}

double Face::volume (const Vector &v) const
{
  // for conveniance, work as if v were the origin:
  Vector a = vertex[0]->pos - v;
  Vector b = vertex[1]->pos - v;
  Vector c = vertex[2]->pos - v;

  // vector scalar product yields volume of parallelpiped
  // formed by the three vectors (a, b, c):
  double parallelpipedVolume = a * (b ^ c);

  // tetrahedron is 1/6 of that parallelpiped:
  return parallelpipedVolume / 6.0;
}

double Face::area () const
{
  // vector normal yields:
  //   - direction normal to (v1-v0, v2-v0)
  //   - magnitude equal to the area of the parallelogram 
  //     formed by the two vectors (v1-v0, v2-v0):
  const Vector norm = normal(vertex[0]->pos,vertex[1]->pos,vertex[2]->pos);

  // triangle is 1/2 of that parallelogram:
  return norm.magnitude() / 2.0;
}

Vector Face::areaNormal () const
{
  // vector normal yields:
  //   - direction normal to (v1-v0, v2-v0)
  //   - magnitude equal to the area of the parallelogram 
  //     formed by the two vectors (v1-v0, v2-v0):
  const Vector norm = normal(vertex[0]->pos,vertex[1]->pos,vertex[2]->pos);

  // triangle is 1/2 of that parallelogram:
  return norm / 2.0;
}

Vector Face::unitNormal () const
{
  // vector normal yields returns the area of the parallelogram
  // formed by the two vectors (v1-v0, v2-v0):
  const Vector norm = normal(vertex[0]->pos,vertex[1]->pos,vertex[2]->pos);

  // unit normal has unit length:
  return norm / norm.magnitude();
}

Vector Face::center () const
{
  // returns average of vertex[0..2]:
  return (vertex[0]->pos + vertex[1]->pos + vertex[2]->pos) / 3.0;
}

Vertex &Face::vertexOpposite (const Edge &e)
{
  for (int i = 0; i < Face::numVertices; i++) {
    if (vertex[i] == e.vertex[0]) continue;
    if (vertex[i] == e.vertex[1]) continue;
    return *vertex[i];
  }
  throw "really unexpected result";
}

const Vertex &Face::vertexOpposite (const Edge &e) const
{
  for (int i = 0; i < Face::numVertices; i++) {
    if (vertex[i] == e.vertex[0]) continue;
    if (vertex[i] == e.vertex[1]) continue;
    return *vertex[i];
  }
  throw "really unexpected result";
}

void SimplicialComplex::AddCrossScaffolding (bool copyOriginal)
{
  if (debug) cerr << "in addcs\n";

  // make materials
  Material &cross = MakeMaterial(100.0, 0.20, true,  0.1);
  Material &scaff = MakeMaterial(  0.0, 0.20, true,  0.1);

  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    Edge &e = it->second;
    Face *f0 = e.face[0];
    Face *f1 = e.face[1];
    if (f0 == NULL || f1 == NULL) continue;

    // find the vertex for each adjacent face which is *not* on this edge
    Vertex &v0 = f0->vertexOpposite(e);
    Vertex &v1 = f1->vertexOpposite(e);

    // link those two vertices:
    MakeEdge(cross,v0,v1);

    // copy the original edge in scaff, if requested
    if (copyOriginal)
      MakeEdge(scaff,*e.vertex[0],*e.vertex[1]);
  }
  if (debug) cerr << "done addcs\n";
}

void SimplicialComplex::AddDualScaffolding ()
{
  if (debug) cerr << "in addds\n";

  Material &strut = MakeMaterial(0.0, 0.20, true,  0.1);
  Material &scaff = MakeMaterial(0.0, 0.20, true,  0.1);

  hash_map<Face *,Vertex *> sputniks;

  for (FaceMapIt it = faces.begin(); it != faces.end(); it++) {
    Face &f = it->second;
    const Vector center = f.center();
    const Vector unitNormal = f.unitNormal();

    Vertex &newp = MakeVertex(center - unitNormal*20);

    MakeEdge(strut,*f.vertex[0],newp);
    MakeEdge(strut,*f.vertex[1],newp);
    MakeEdge(strut,*f.vertex[2],newp);

    sputniks[&f] = &newp;
  }

  if (debug) cerr << "built faces\n";

  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    Edge &e = it->second;
    if (e.face[0] == NULL) continue;
    if (e.face[1] == NULL) continue;
    if (sputniks.find(e.face[0]) == sputniks.end()) throw "sputnik 0";
    if (sputniks.find(e.face[1]) == sputniks.end()) throw "sputnik 1";
    MakeEdge(scaff,*sputniks[e.face[0]],*sputniks[e.face[1]]);
  }

  if (debug) cerr << "done addds\n";
}

void SimplicialComplex::bisect (Edge &e,
																BisectorMap &bisectors,
																ChildMap &childMap,
																ChildSet &childSet)
{
  // sort of a sanity check: algorithms using bisect() depend on this
  // rewrite them before removing this check!
  if (bisectors.find(&e) != bisectors.end() || 
      childMap.find(&e) != childMap.end())
    throw "tried to bisect bisected edge!";

  // create bisecting vertex:
  Vertex &bisector = MakeVertex((e.vertex[0]->pos + e.vertex[1]->pos) / 2);
  bisector.vel = (e.vertex[0]->vel + e.vertex[1]->vel) / 2;
  bisector.mass = (e.vertex[0]->mass + e.vertex[1]->mass) / 2;
  bisectors[&e] = &bisector;

  // create the 2 new edges:
  Edge &childA = MakeEdge(*e.material,*e.vertex[0],bisector);
  childA.springFactor = e.springFactor;

  Edge &childB = MakeEdge(*e.material,bisector,*e.vertex[1]);
  childB.springFactor = e.springFactor;

  childMap[&e] = make_pair(&childA,&childB);
  childSet.insert(&childA);
  childSet.insert(&childB);
}

void SimplicialComplex::BisectEdges ()
{
  const int upperBound = currentObjID;
  int oldNumV = vertices.size();
  double oldKE = physics.kineticEnergy;

  // first subdivide each edge
  int numParents = 0;
  ChildMap childMap;
  ChildSet childSet;
  BisectorMap bisectors;
  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    if (it->first >= upperBound) continue;
    Edge &e = it->second;
    bisect(e,bisectors,childMap,childSet);
    numParents += 1;
  }

  DivisionCleanup(oldKE,numParents,oldNumV,bisectors,childSet);
}

void solveSpringEnergyQuadratic (Edge &e,
																 int numParents,
																 int numChildren)
{
}

void SimplicialComplex::DivisionCleanup (double oldKE,
																				 int numParents,
																				 int oldNumV,
																				 BisectorMap &bisectors,
																				 ChildSet &childSet)
{
  physics.CalculateEpiphenomena();

  // scale the masses properly so as to conserve mass:
  int newNumV = vertices.size();
  double massFactor = (1.0 * oldNumV) / newNumV;
  for (VertexMapIt v = vertices.begin(); v != vertices.end(); v++)
    v->second.mass *= massFactor;

  // scale the velocities properly so as to conserve energy:
  double curKE = physics.kineticEnergy;
  if (curKE != 0.0) {
    double keFactor = sqrt(oldKE/curKE);
    for (VertexMapIt v = vertices.begin(); v != vertices.end(); v++)
      v->second.vel *= keFactor;
  }

  // scale the material factors properly so as to conserve energy:
  //    - note, the parents are going away soon
  //    - the magic X represents the assumption that all edges bisected
  //    - it would be 3 following edge trisection
  const double X = 2;
  //    - for nonzero restLength, energy conservation depends on everything!
  //    - we need pre-stressed materials
  //       - this code still doesn't work if the edges aren't all the same
  //         length, etc.  Such as when subdividing a spinning structure.
  //         - the code is wrong-ish
  //         - as right as it is, the code shouldn't resolve for each edge
  for (ChildSet::iterator it = childSet.begin(); it != childSet.end(); it++) {
    Edge &e = **it;

    const double a = e.springFactor;
    const double r = e.material->restLength;
    const double l = (e.vertex[1]->pos - e.vertex[0]->pos).magnitude();
    const double p = numParents;
    const double c = childSet.size();

    const double A = r * r;
    const double B = - ( 2*l*r + (p/(c*a))*(X*l-r*a)*(X*l-r*a));
    const double C = l * l;

    double sF;
    if (A == 0) {
      sF = -C/B;
    }
    else {
      double det = B*B - 4*A*C;
      if (det < 0.0)
				throw "no real roots!";
      double sF0 = (-B + sqrt(det)) / (2*A);
      double sF1 = (-B - sqrt(det)) / (2*A);
      
      if (false) cerr << "a:   "  << a << "\n";
      if (false) cerr << "sF0: " << sF0 << "\n";
      if (false) cerr << "sF1: " << sF1 << "\n";

      // The two roots represent the energy balance in the two available
      // modes, compression and tension.  In a material which supports 
      // compression (r > 0), you could have the same potential energy in
      // either mode.
      //   - one root is "same energy, same mode"
      //   - other root is "same energy, complimentary mode"
      // The following falsed-out code would alternate between modes
      // each division:
      if (false) {
				double test0 = (l*X-r*a)*(l-r*sF0);
				if (false) cerr << "test0: " << test0 << "\n";
				double test1 = (l*X-r*a)*(l-r*sF1);
				if (false) cerr << "test1: " << test1 << "\n";
				if (test0 < 0.0)
					sF = sF0;
				else
					sF = sF1;
      }
      else {
				// we want the mode-preserving root:
				double test0 = (l*X-r*a)*(l-r*sF0);
				if (false) cerr << "test0: " << test0 << "\n";
				if (test0 > 0.0)
					sF = sF0;
				else
					sF = sF1;
      }
    }

    if (debug) cerr << "sF:  " << sF << "\n";
      
    if (sF < 0.0)
      throw "root out of range!";

    e.springFactor = sF;
  }

  // remove all original edges: those with bisectors
  for (EdgeMapIt it = edges.begin(); it != edges.end(); ) {
    Edge &e = it->second;
    if (bisectors.find(&e) == bisectors.end())
      it++;
    else
      edges.erase(it++);
  }
}

void SimplicialComplex::SubdivideFaces ()
{
  const int upperBound = currentObjID;
  int oldNumV = vertices.size();
  double oldKE = physics.kineticEnergy;

  // first subdivide each edge which is on a face,
  // keeping track of each bisector:
  int numParents = 0;
  BisectorMap bisectors;
  ChildMap childMap;
  ChildSet childSet;
  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    if (it->first >= upperBound) continue;
    Edge &e = it->second;
    if (e.face[0] == NULL && e.face[1] == NULL) continue;
    bisect(e,bisectors,childMap,childSet);
    numParents += 1;
  }

  // create the 3 new edges corresponding to the 3 bisectors corresponding
  // to each face's edges
  // also create the 4 new faces thereby implied
  const int upperBound2 = currentObjID;
  for (FaceMapIt it = faces.begin(); it != faces.end(); it++) {
    if (it->first >= upperBound2) continue;
    Face &f = it->second;
    Vertex *va = f.vertex[0];
    Vertex *vb = f.vertex[1];
    Vertex *vc = f.vertex[2];
    if (debug) cerr << "yo1\n";
    if (va == NULL) throw "dang";
    if (vb == NULL) throw "dang";
    if (vc == NULL) throw "dang";
    Edge *ea = f.edge[0];
    Edge *eb = f.edge[1];
    Edge *ec = f.edge[2];
    if (ea == NULL) throw "ea == NULL";
    if (eb == NULL) throw "eb == NULL";
    if (ec == NULL) throw "ec == NULL";
    if (ea->vertex[0] != va && ea->vertex[1] != va) throw "yo3a";
    if (ea->vertex[0] != vb && ea->vertex[1] != vb) throw "yo3b";
    if (eb->vertex[0] != vb && eb->vertex[1] != vb) throw "yo3c";
    if (eb->vertex[0] != vc && eb->vertex[1] != vc) throw "yo3d";
    if (ec->vertex[0] != vc && ec->vertex[1] != vc) throw "yo3f";
    if (ec->vertex[0] != va && ec->vertex[1] != va) throw "yo3g";

    if (bisectors.find(f.edge[0]) == bisectors.end()) throw "yo4aa";
    if (bisectors.find(f.edge[1]) == bisectors.end()) throw "yo4ba";
    if (bisectors.find(f.edge[2]) == bisectors.end()) throw "yo4ca";

    Vertex *v0 = bisectors[f.edge[0]];
    Vertex *v1 = bisectors[f.edge[1]];
    Vertex *v2 = bisectors[f.edge[2]];
    if (v0 == NULL) throw "yo4a";
    if (v1 == NULL) throw "yo4b";
    if (v2 == NULL) throw "yo4c";

    Edge &e0 = MakeEdge(*ea->material,*v0,*v1);
    Edge &e1 = MakeEdge(*eb->material,*v1,*v2);
    Edge &e2 = MakeEdge(*ec->material,*v2,*v0);
    e0.springFactor = ea->springFactor;
    e1.springFactor = eb->springFactor;
    e2.springFactor = ec->springFactor;
    childSet.insert(&e0);
    childSet.insert(&e1);
    childSet.insert(&e2);

    Edge *ea0, *eb0;
    if (va == childMap[ea].first->vertex[0] || 
				va == childMap[ea].first->vertex[0]) {
      ea0 = childMap[ea].first;
      eb0 = childMap[ea].second;
    }
    else {
      ea0 = childMap[ea].second;
      eb0 = childMap[ea].first;
    }

    Edge *eb1, *ec1;
    if (vb == childMap[eb].first->vertex[0] || 
				vb == childMap[eb].first->vertex[0]) {
      eb1 = childMap[eb].first;
      ec1 = childMap[eb].second;
    }
    else {
      eb1 = childMap[eb].second;
      ec1 = childMap[eb].first;
    }

    Edge *ec2, *ea2;
    if (vc == childMap[ec].first->vertex[0] || 
				vc == childMap[ec].first->vertex[0]) {
      ec2 = childMap[ec].first;
      ea2 = childMap[ec].second;
    }
    else {
      ec2 = childMap[ec].second;
      ea2 = childMap[ec].first;
    }

    // I know this is really convoluted...it's just getting all the
    // face/edge/vertex relationships linked up right.

    if (true)
			{
				Face &fx = MakeNullFace();
				fx.edge[0]   = &e0;
				fx.edge[1]   = &e1;
				fx.edge[2]   = &e2;
				fx.vertex[0] =  v0;
				fx.vertex[1] =  v1;
				fx.vertex[2] =  v2;
			}

    if (true)
			{
				Face &fa = MakeNullFace();
				fa.edge[0]   = ea0;
				fa.edge[1]   = &e2;
				fa.edge[2]   = ea2;
				fa.vertex[0] =  va;
				fa.vertex[1] =  v0;
				fa.vertex[2] =  v2;
				if (debug) cerr << "edge ids: " 
												<< fa.edge[0]->objID << " " 
												<< fa.edge[1]->objID << " " 
												<< fa.edge[2]->objID;
			}

    if (true)
			{
				Face &fb = MakeNullFace();
				fb.edge[0]   = eb1;
				fb.edge[1]   = &e0;
				fb.edge[2]   = eb0;
				fb.vertex[0] =  vb;
				fb.vertex[1] =  v1;
				fb.vertex[2] =  v0;
				if (debug) cerr << "edge ids: " 
												<< fb.edge[0]->objID << " " 
												<< fb.edge[1]->objID << " " 
												<< fb.edge[2]->objID;
			}

    if (true)
			{
				Face &fc = MakeNullFace();
				fc.edge[0]   = ec2;
				fc.edge[1]   = &e1;
				fc.edge[2]   = ec1;
				fc.vertex[0] =  vc;
				fc.vertex[1] =  v2;
				fc.vertex[2] =  v1;
				if (debug) cerr << "edge ids: " 
												<< fc.edge[0]->objID << " " 
												<< fc.edge[1]->objID << " " 
												<< fc.edge[2]->objID;
			}
   
    if (debug) cerr << "yoX: face " << f.objID << "  ok\n";
  }

  // remove all original faces: those with 3 edges with bisectors
  if (debug) cerr << "yo remove face\n";
  for (FaceMapIt it = faces.begin(); it != faces.end(); ) {
    Face &f = it->second;
    if (bisectors.find(f.edge[0]) == bisectors.end() ||
				bisectors.find(f.edge[1]) == bisectors.end() ||
				bisectors.find(f.edge[2]) == bisectors.end()) {
      it++;
      continue;
    }
    if (debug) cerr << "erasing face: " << it->first << "\n";
    faces.erase(it++);
  }
  if (debug) cerr << "yo remove edge\n";

  DivisionCleanup(oldKE,numParents,oldNumV,bisectors,childSet);

  // relink all edge's face links 'cause they're messed up:
  if (debug) cerr << "yo relink face links\n";
  for (EdgeMapIt e = edges.begin(); e != edges.end(); e++) {
    e->second.face[0] = NULL;
    e->second.face[1] = NULL;
  }
  if (debug) cerr << "yo relink some more\n";
  for (FaceMapIt f = faces.begin(); f != faces.end(); f++) {
    if (debug) cerr << "   on face: " << f->first << "\n";
    for (int i = 0; i < Face::numEdges; i++) {
      if (debug) cerr << "      ed: " << i << "\n";
      if (f->second.edge[i]->face[0] == NULL) {
				if (f->second.edge[i]->face[1] != NULL) {
					throw "complaint 1";
				}
				f->second.edge[i]->face[0] = &f->second;
				continue;
      }
      else {
				if (f->second.edge[i]->face[1] != NULL) {
					throw "complaint 2";
				}
				f->second.edge[i]->face[1] = &f->second;
      }
    }
  }

  // sanity check
  if (debug) cerr << "yo sanity 1\n";
  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    Edge &e = it->second;
    if (bisectors.find(&e) != bisectors.end()) throw "yo sanity1a";
    if (childMap.find(&e) != childMap.end()) throw "yo sanity1a";
  }
  if (debug) cerr << "yo sanity 2\n";
  for (FaceMapIt it = faces.begin(); it != faces.end(); it++) {
    Face &f = it->second;
    if (bisectors.find(f.edge[0]) != bisectors.end()) throw "yo sanity2a";
    if (bisectors.find(f.edge[1]) != bisectors.end()) throw "yo sanity2a";
    if (bisectors.find(f.edge[2]) != bisectors.end()) throw "yo sanity2a";
  }
  if (debug) cerr << "yo clear\n";
}

void SimplicialComplex::clear ()
{
  // clear all preexisting geometry
  vertices.clear();
  materials.clear();
  edges.clear();
  faces.clear();
}

void SimplicialComplex::MakeTwicePinnedWeight ()
{
  clear(); 

  Material &base = MakeMaterial(0.0, 0.20, true,  0.1);

  // define locations on half of a cube's face's centers
  Vertex &v0 = MakeVertex(Vector(   0,   0,   0));
  Vertex &v1 = MakeVertex(Vector( 100,   0,   0));
  Vertex &v2 = MakeVertex(Vector(-100,   0,   0));

  v1.pinned = true;
  v2.pinned = true;

  MakeEdge(base,v0,v1);
  MakeEdge(base,v0,v2);

  CenterOnOrigin();
}


void SimplicialComplex::MakePinnedWeight ()
{
  clear();
  
  Material &base = MakeMaterial(0.0, 0.20, true,  0.1);
  
  // define locations on half of a cube's face's centers
  Vertex &v0 = MakeVertex(Vector(100,   0,   0));
  Vertex &v1 = MakeVertex(Vector(  0,   0,   0));
  v0.vel = Vector( 0,  50 * sqrt(2.0),   0);
  v1.pinned = true;

  MakeEdge(base,v0,v1);

  //CenterOnOrigin();
}

void SimplicialComplex::MakeBinaryWeight ()
{
  clear();

  // make a material
  Material &base = MakeMaterial(  0.0, 0.20, true,  0.1);
  
  // define locations on half of a cube's face's centers
  Vertex &v0 = MakeVertex(Vector(  50,   0,   0));
  Vertex &v1 = MakeVertex(Vector( -50,   0,   0));
  v0.vel = Vector( 0,  50,   0);
  v1.vel = Vector( 0, -50,   0);

  MakeEdge(base,v0,v1);

  CenterOnOrigin();
}

void SimplicialComplex::MakeTrinaryWeight ()
{
  clear();

  const int n = 3;
  
  Material &base = MakeMaterial(0.0, 0.20, true,  0.1);

  // define locations on half of a cube's face's centers
  const double scale   = 30;
  const double angle   = 2 * PI / n;

  Vertex *v[n];

  for (int i = 0; i < n; i++) {
		v[i] = &MakeVertex(Vector(cos(i*angle),sin(i*angle),0)*scale);
		v[i]->vel = Vector(-sin(i*angle),cos(i*angle),0)*scale;  
  }

  for (int i = 0; i < n; i++)
    MakeEdge(base,*v[i],*v[(i+1)%n]);

  CenterOnOrigin();
}

void SimplicialComplex::CenterOnOrigin ()
{
  // center of mass is moved to origin
  Vector ave(0.0,0.0,0.0);
  for (VertexMapIt v = vertices.begin(); 
       v != vertices.end(); v++) {
    ave += v->second.pos;
  }
  ave /= vertices.size();
  for (VertexMapIt v = vertices.begin(); 
       v != vertices.end(); v++) {
    v->second.pos -= ave;
  }
}

void SimplicialComplex::InvertThroughOrigin ()
{
  for (VertexMapIt v = vertices.begin(); v != vertices.end(); v++) {
    v->second.pos *= -1.0;
    v->second.vel *= -1.0;
  }
}

void SimplicialComplex::MakeTorus (TorusType type, int ni, int mi)
{
  const int n = max(ni,2);
  const int m = max(mi,2);

  clear();

  // make materials
  Material &base  = MakeMaterial(  0.0, 0.20, true,  0.1);

  // How nice that we can exploit the torus's modulo-modulo
  // symmetry with two nested loops!  Easier than sphere.

  const double scale = 10;
  const double angle = 2 * PI / m;

  vector<vector<Vertex *> > v;
  
  v.resize(n);
  for (int i = 0; i < n; i++) {
    v[i].resize(m);
    for (int j = 0; j < m; j++)
      v[i][j] = &MakeVertex(Vector(cos(j*angle),sin(j*angle),i)*scale);
  }

  // this makes a bunch of redundant edges:
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {

      bool t = false;
      if (type == left) t = true;
      else if (type == right) t = false;
      else if (type == modi) t = (0 == i%2);
      else if (type == modj) t = (0 == j%2);
      else if (type == modij) t = (0 == (i+j)%2);
      else throw "unrecognized TorusType";

      if (t) {
				MakeFace(base,
								 v[(i+1)%n][(j+1)%m],
								 v[(i+1)%n][(j+0)%m],
								 v[(i+0)%n][(j+0)%m]);
				MakeFace(base,
								 v[(i+1)%n][(j+1)%m],
								 v[(i+0)%n][(j+0)%m],
								 v[(i+0)%n][(j+1)%m]);
      }
      else {
				MakeFace(base,
								 v[(i+1)%n][(j+0+m)%m],
								 v[(i+0)%n][(j+0+m)%m],
								 v[(i+0)%n][(j+1+m)%m]);
				MakeFace(base,
								 v[(i+1)%n][(j+0+m)%m],
								 v[(i+0)%n][(j+1+m)%m],
								 v[(i+1)%n][(j+1+m)%m]);
      }
    }
  }

  cerr << "num edges: " << edges.size() << "\n";

  set<pair<int,int> > redundantEdges;
  for (EdgeMapIt i0 = edges.begin(); i0 != edges.end(); i0++) {
    Edge &e0 = i0->second;
    for (EdgeMapIt i1 = i0; i1 != edges.end(); i1++) {
      if (i0 == i1) continue;
      Edge &e1 = i1->second;
      if (e0.vertex[0] != e1.vertex[0] && e0.vertex[0] != e1.vertex[1])
				continue;
      if (e0.vertex[1] != e1.vertex[0] && e0.vertex[1] != e1.vertex[1])
				continue;
      redundantEdges.insert(make_pair(i0->first,i1->first));
    }
  }
      
  cerr << "redundant: " << redundantEdges.size() << "\n";

  for (set<pair<int,int> >::iterator i = redundantEdges.begin();
       i != redundantEdges.end();
       i++) {
    Edge &e0 = edges[i->first];
    Edge &e1 = edges[i->second];
    if (e0.face[0] == NULL) throw "unexpected null face";
    if (e1.face[0] == NULL) throw "unexpected null face";
    if (e0.face[1] != NULL) throw "unexpected non-null face";
    if (e1.face[1] != NULL) throw "unexpected non-null face";
    Face *f = e1.face[0];
    e0.face[1] = f;
    if (f->edge[0] == &e1)
      f->edge[0] = &e0;
    else if (f->edge[1] == &e1)
      f->edge[1] = &e0;
    else if (f->edge[2] == &e1)
      f->edge[2] = &e0;
    else
      throw "unexpected bad face linkage";
    edges.erase(edges.find(e1.objID));
  }

  cerr << "num edges: " << edges.size() << "\n";

  CenterOnOrigin();
  myConsistency();
}

void SimplicialComplex::MakeTetrahedron ()
{
  clear();
  
  Material &cross = MakeMaterial(100.0, 0.20, true,  0.1);
  
  // define locations on half of a cube's face's centers
  Vertex &v0 = MakeVertex(Vector( 0,  0,  0));
  Vertex &v1 = MakeVertex(Vector(50, 50,  0));
  Vertex &v2 = MakeVertex(Vector(50,  0, 50));
  Vertex &v3 = MakeVertex(Vector( 0, 50, 50));

  // link all combos:
  MakeEdge(cross,v0,v1);
  MakeEdge(cross,v0,v2);
  MakeEdge(cross,v0,v3);
  MakeEdge(cross,v1,v2);
  MakeEdge(cross,v2,v3);
  MakeEdge(cross,v3,v1);

  CenterOnOrigin();
}

void SimplicialComplex::MakeOctahedron ()
{
  clear();

  Material &cross = MakeMaterial(100.0, 0.20, true,  0.1);
  
  // define locations centered around origin
  Vertex *v[3][2];
  v[0][0] = &MakeVertex(Vector( -25,   0,   0));
  v[0][1] = &MakeVertex(Vector(  25,   0,   0));
  v[1][0] = &MakeVertex(Vector(   0, -25,   0));
  v[1][1] = &MakeVertex(Vector(   0,  25,   0));
  v[2][0] = &MakeVertex(Vector(   0,   0, -25));
  v[2][1] = &MakeVertex(Vector(   0,   0,  25));

  // edges from "top" vertex to "equatorial" vertices:
  MakeEdge(cross,*v[0][0],*v[1][0]);
  MakeEdge(cross,*v[0][0],*v[1][1]);
  MakeEdge(cross,*v[0][0],*v[2][0]);
  MakeEdge(cross,*v[0][0],*v[2][1]);

  // edges from "bottom" vertex to "equatorial" vertices:
  MakeEdge(cross,*v[0][1],*v[1][0]);
  MakeEdge(cross,*v[0][1],*v[1][1]);
  MakeEdge(cross,*v[0][1],*v[2][0]);
  MakeEdge(cross,*v[0][1],*v[2][1]);

  // edges along "equator":
  MakeEdge(cross,*v[1][0],*v[2][0]);
  MakeEdge(cross,*v[1][0],*v[2][1]);
  MakeEdge(cross,*v[1][1],*v[2][0]);
  MakeEdge(cross,*v[1][1],*v[2][1]);
}

void SimplicialComplex::MakeIcosahedron ()
{
  clear();

  Material &cross = MakeMaterial(100.0, 0.20, true,  0.1);
  
  // magic numbers from Pythagoreans:
  static const int numVertices = 12;
  static const int numEdges = 30;
  //static const int numFaces = 20;

  // heavy use of magic numbers pulled from:
  //    OpenGL Programming Guide, 2nd ed.
  //    by Woo, Neider, Davis, Shreiner
  static const double X = 0.525731112119133606;
  static const double Z = 0.850650808352039932;
  static const double vdata[numVertices][3] = {
    // coordinates of vertices, centered on origin
    {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
    {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
    {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
  };
  // index of vertex triples (ref into vdata) to make faces: unused
  /*
		static const int findices[numFaces][3] = {
    {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
    {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
    {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
    {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
    };
  */

  // more magic numbers:
  static const int eindices[numEdges][2] = {
    // index of vertex pairs (ref into vdata) to make edges
    { 0,  1}, { 0,  4}, { 0,  6}, { 0,  9}, { 0, 11},
    { 1,  4}, { 1,  8}, { 1,  6}, { 1, 10}, 
    { 2,  3}, { 2,  5}, { 2,  7}, { 2,  9}, { 2, 11},
    { 3,  5}, { 3,  7}, { 3,  8}, { 3, 10},
    { 4,  5}, { 4,  8}, { 4,  9},
    { 5,  9}, { 5,  8},
    { 6,  7}, { 6, 10}, { 6, 11},
    { 7, 10}, { 7, 11},
    { 8, 10},
    { 9, 11},
  };
  static const double scale = 25.0; // make it bigger

  // make vertices according to magic numbers:
  int vertexIDtable[numVertices];
  for (int i = 0; i < numVertices; i++) {
    vertexIDtable[i] = MakeVertex(Vector(scale*vdata[i][0],
																				 scale*vdata[i][1],
																				 scale*vdata[i][2]))
      .objID;
  }

  // link up appropriate edges:
  int edgeIDtable[numEdges];
  for (int i = 0; i < numEdges; i++) {
    edgeIDtable[i] = 
      MakeEdge(cross,
							 vertices[vertexIDtable[eindices[i][0]]],
							 vertices[vertexIDtable[eindices[i][1]]])
      .objID;
  }
}

void SimplicialComplex::DiamondSubdivide ()
{
  // NOTE: This method has unnecessary quadraticisms and cubism.
  // NOTE: It's really slow, but rarely used and (almost) deprecated.

  const int upperBound = currentObjID;
  int oldNumV = vertices.size();
  double oldKE = physics.kineticEnergy;

  // remove any preexisting faces 'cause they'll get screwed up:
  for (EdgeMapIt e = edges.begin(); e != edges.end(); e++) {
    e->second.face[0] = NULL;
    e->second.face[1] = NULL;
  }
  faces.clear();

  // bisect each Edge (creating bisector, child[0], child[1])
  int numParents = 0;
  BisectorMap bisectors;
  ChildMap childMap;
  ChildSet childSet;
  for (EdgeMapIt e = edges.begin(); e != edges.end(); e++) {
    if (e->first >= upperBound) continue;
    bisect(e->second,bisectors,childMap,childSet);
    numParents += 1;
  }

  // for each original vertex...
  for (VertexMapIt v = vertices.begin(); v != vertices.end(); v++) {
    if (v->first >= upperBound) continue;

    if (debug) cerr << "vertex: " << v->first << "\n";

    // for each original edge intersecting that vertex...
    for (EdgeMapIt e0 = edges.begin(); e0 != edges.end(); e0++) {
      if (e0->first >= upperBound) continue;
      if (bisectors.find(&e0->second) == bisectors.end()) throw "dang";
      if (e0->second.vertex[0] != &v->second &&
					e0->second.vertex[1] != &v->second)
				continue;

      // for each distinct original edge intersecting that vertex...
      for (EdgeMapIt e1 = edges.begin(); e1 != edges.end(); e1++) {
				if (e1->first >= upperBound) continue;
				if (bisectors.find(&e1->second) == bisectors.end()) throw "dang";
				if (e1->first == e0->first) continue;
				if (e1->second.vertex[0] != &v->second &&
						e1->second.vertex[1] != &v->second)
					continue;
				if (e1->first < e0->first) continue; // avoid dupes
				if (e1->second.vertex[0] != e0->second.vertex[0] &&
						e1->second.vertex[0] != e0->second.vertex[1] &&
						e1->second.vertex[1] != e0->second.vertex[0] &&
						e1->second.vertex[1] != e0->second.vertex[1])
					continue;
	
				// if there is a third distinct original edge forming a triangle...
				for (EdgeMapIt e2 = edges.begin(); e2 != edges.end(); e2++) {
					if (e2->first >= upperBound) continue;
					if (bisectors.find(&e2->second) == bisectors.end()) throw "dang";
					if (e2->first == e0->first) continue;
					if (e2->first == e1->first) continue;
					if (e2->second.vertex[0] == &v->second ||
							e2->second.vertex[1] == &v->second)
						continue;
					if (e2->second.vertex[0] != e0->second.vertex[0] &&
							e2->second.vertex[0] != e0->second.vertex[1] &&
							e2->second.vertex[1] != e0->second.vertex[0] &&
							e2->second.vertex[1] != e0->second.vertex[1])
						continue;
					if (e2->second.vertex[0] != e1->second.vertex[0] &&
							e2->second.vertex[0] != e1->second.vertex[1] &&
							e2->second.vertex[1] != e1->second.vertex[0] &&
							e2->second.vertex[1] != e1->second.vertex[1])
						continue;

					// create an edge from the bisectors of the first two edges:
					Edge &e = MakeEdge(*e0->second.material,
														 *bisectors.find(&e0->second)->second,
														 *bisectors.find(&e1->second)->second);
					e.springFactor = e0->second.springFactor;
					childSet.insert(&e);
				}
      }
    }
  }

  DivisionCleanup(oldKE,numParents,oldNumV,bisectors,childSet);
}

void SimplicialComplex::Sphericize ()
{
  // calculate center
  Vector center(0,0,0);
  for (VertexMapIt v = vertices.begin(); v != vertices.end(); v++) {
    center += v->second.pos;
  }
  center /= vertices.size();

  // calculate average radius
  double radius = 0.0;
  for (VertexMapIt v = vertices.begin(); v != vertices.end(); v++) {
    radius += (v->second.pos-center).magnitude();
  }
  radius /= vertices.size();

  // project every vertex to be length radius:
  for (VertexMapIt v = vertices.begin(); v != vertices.end(); v++) {
    if (v->second.pinned) continue;
    v->second.pos *= radius / v->second.pos.magnitude();
  }
}

void SimplicialComplex::Harmonize ()
{
  for (MaterialMapIt it = materials.begin(); it != materials.end(); it++)
    Harmonize(it->second);
}

void SimplicialComplex::Harmonize (Material &material)
{
  // find the average length and springFactor for all edges in category:
  double aveLength = 0.0;
  double aveSpringFactor = 0.0;
  int num = 0;
  for (EdgeMapIt e = edges.begin(); e != edges.end(); e++) {
    if (e->second.material != &material) continue;
    Vector v = e->second.vertex[1]->pos - e->second.vertex[0]->pos;
    aveLength += v.magnitude();
    aveSpringFactor += e->second.springFactor;
    num++;
  }
  if (num == 0) return;
  aveLength /= num;
  aveSpringFactor /= num;

  // scale each edge's springFactor to the material's relative length:
  for (EdgeMapIt e = edges.begin(); e != edges.end(); e++) {
    if (e->second.material != &material) continue;
    Vector v = e->second.vertex[1]->pos - e->second.vertex[0]->pos;
    e->second.springFactor = aveSpringFactor * (v.magnitude() / aveLength);
  }
}

void SimplicialComplex::DeHarmonize ()
{
  for (MaterialMapIt it = materials.begin(); it != materials.end(); it++)
    DeHarmonize(it->second);
}

void SimplicialComplex::DeHarmonize (Material &material)
{
  // find the average springFactor for all edges in category:
  double aveSpringFactor = 0.0;
  int num = 0;
  for (EdgeMapIt e = edges.begin(); e != edges.end(); e++) {
    if (e->second.material != &material) continue;
    aveSpringFactor += e->second.springFactor;
    num++;
  }
  if (num == 0) return;
  aveSpringFactor /= num;

  cerr << "deharmonize: " << aveSpringFactor << "\n";

  // set each edge's springFactor to the average
  for (EdgeMapIt e = edges.begin(); e != edges.end(); e++) {
    if (e->second.material != &material) continue;
    e->second.springFactor = aveSpringFactor;
  }
}

void SimplicialComplex::myConsistency ()
{
  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    Edge &e = it->second;
    
    // check objID in map:
    if (e.objID != it->first)
      throw "Edge objID bogus";

    // check cross-links:
    if (e.face[0] != NULL) {
      Face *f = e.face[0];
      if (f->edge[0] != &e && f->edge[1] != &e && f->edge[2] != &e)
				throw "Edge face[0] not relinked";
    }
    if (e.face[1] != NULL) {
      Face *f = e.face[1];
      if (f->edge[0] != &e && f->edge[1] != &e && f->edge[2] != &e)
				throw "Edge face[1] not relinked";
    }
  }

  for (FaceMapIt it = faces.begin(); it != faces.end(); it++) {
    Face &f = it->second;
    
    // check objID in map:
    if (f.objID != it->first)
      throw "Face objID bogus";

    // check cross-links:
    if (f.edge[0] != NULL) {
      Edge *e = f.edge[0];
      if (e->face[0] != &f && e->face[1] != &f)
				throw "Face edge[0] not relinked";
      if (e->vertex[0] != f.vertex[0] &&
					e->vertex[0] != f.vertex[1] &&
					e->vertex[0] != f.vertex[2] &&
					e->vertex[1] != f.vertex[0] &&
					e->vertex[1] != f.vertex[1] &&
					e->vertex[1] != f.vertex[2])
				throw "Face edge[0] bogus vertices";
    }
    if (f.edge[1] != NULL) {
      Edge *e = f.edge[1];
      if (e->face[0] != &f && e->face[1] != &f)
				throw "Face edge[1] not relinked";
      if (e->vertex[0] != f.vertex[0] &&
					e->vertex[0] != f.vertex[1] &&
					e->vertex[0] != f.vertex[2] &&
					e->vertex[1] != f.vertex[0] &&
					e->vertex[1] != f.vertex[1] &&
					e->vertex[1] != f.vertex[2])
				throw "Face edge[0] bogus vertices";
    }
    if (f.edge[2] != NULL) {
      Edge *e = f.edge[2];
      if (e->face[0] != &f && e->face[1] != &f)
				throw "Face edge[2] not relinked";
      if (e->vertex[0] != f.vertex[0] &&
					e->vertex[0] != f.vertex[1] &&
					e->vertex[0] != f.vertex[2] &&
					e->vertex[1] != f.vertex[0] &&
					e->vertex[1] != f.vertex[1] &&
					e->vertex[1] != f.vertex[2])
				throw "Face edge[0] bogus vertices";
    }
  }
}


void SimplicialComplex::MakeSpiral (int num)
{
  clear();

  /*
    -------------------------------------------------------------------------
    Taken from spiral.cpp:
    -------------------------------------------------------------------------
    Some after market mods by jhw: castrated, for instance.
    -------------------------------------------------------------------------
    This program will generate a given number of spiral points uniformly 
    distributed on the surface of a sphere. The number of points is given on 
    the command line as the first parameter.  Thus `spiral 100' will generate
    100 points on the surface of a sphere, and output them to stdout.
    A number of different command-line flags are provided to set the
    radius of the sphere, control the output format, or generate points on
    an ellipsoid.  The definition of the flags is printed if the program is
    run without arg0uments: `spiral'.
    The idea behind the algorithm is that one can cut the globe with 
    N horizontal planes spaced 2/(N-1) units apart, forming N circles of 
    latitude on the sphere, each latitude containing one spiral point.  To
    obtain the kth spiral point, one proceeds upward from the (k-1)st point
    (theta(k-1), phi(k-1)) along a great circle to the next latitude and 
    travels counterclockwise along ti for a fixed distance to arrive at the 
    kth point (theta(k), phi(k)).
    The default output is integers, rounded from the floating point
    computation.  The rounding implies that some points will fall outside
    the sphere, and some inside.  If all are required to be inside, then
    the calls to irint() should be removed.
    The flags -a, -b, -c are used to set ellipsoid axis lengths.
    Note that the points are not uniformly distributed on the ellipsoid: they
    are uniformly distributed on the sphere and that is scaled to an
    ellipsoid. random() is used to generate random numbers, seeded with
    time().
    How to compile:
    gcc -o spiral spiral.c -lm

    Reference: E.B. Saff and A.B.J. Kuijlaars,
    Distributing Many Points on a Sphere,
    The Mathematical Intelligencer, 19(1), Winter (1997);

    Written by Joseph O'Rourke and Min Xu, June 1997.
    Used in the textbook, "Computational Geometry in C."
    Questions to orourke@cs.smith.edu.
    --------------------------------------------------------------------
    This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
    redistributed in its entirety provided that this copyright notice is
    not removed.
    --------------------------------------------------------------------
  */
    
  int k;                /* index */
  double phi1, phi, theta, h, x, y, z;
  double R = 100.0;	/* default radius */
  int r;		/* true radius */
  int r1, r2, r3;	/* ellipsoid axis lengths */
  
  r = (int)R;
  r1 = r2 = r3 = r;

  phi1 = 0.0;

  MakeVertex(Vector(0,0,-1*r3));

  for ( k = 2; k <= num - 1; k ++ ) {
    /* Generate a random point on a sphere of radius 1. */
    h = -1 + 2 * ( k - 1 ) / ( double )( num - 1 );
    theta = acos ( h );

    if ( theta < 0 || theta > M_PI )
      throw "theta out of range error in spiral";

    phi = phi1 + 3.6 / ( sqrt ( ( double )num * ( 1 - h * h ) ) ); 
    phi = fmod ( phi, 2 * M_PI );
    phi1 = phi;

    x = cos ( phi ) * sin ( theta );
    y = sin ( phi ) * sin ( theta );
    /* z = cos ( theta ); But z==h, so: */
    z = h;
    MakeVertex(Vector(r1*x,r2*y,r3*z));
  }

  MakeVertex(Vector(0,0,1*r3));
}

void SimplicialComplex::MakeSphere (int num)
{
  clear();

  /*
    -------------------------------------------------------------------------
    taken from sphere.cpp
    -------------------------------------------------------------------------
    Some after market mods by jhw: castrated for instance
    -------------------------------------------------------------------------
    This program will generate a given number of points uniformly distributed
    on the surface of a sphere. The number of points is given on the command
    line as the first parameter.  Thus `sphere 100' will generate 100 points 
    on the surface of a sphere, and output them to stdout.
    A number of different command-line flags are provided to set the 
    radius of the sphere, control the output format, or generate points on 
    an ellipsoid.  The definition of the flags is printed if the program is 
    run without arg0uments: `sphere'.
    The idea behind the algorithm is that for a sphere of radius r, the 
    area of a zone of width h is always 2*pi*r*h, regardless of where the
    sphere is sliced.  The implication is that the z-coordinates of random
    points on a sphere are uniformly distributed, so that x and y can always
    be generated by a given z and a given angle.
    The default output is integers, rounded from the floating point 
    computation.  The rounding implies that some points will fall outside
    the sphere, and some inside.  If all are required to be inside, then
    the calls to irint() should be removed.  
    The flags -a, -b, -c are used to set ellipsoid axis lengths.  
    Note that the points are not uniformly distributed on the ellipsoid:
    they are uniformly distributed on the sphere and that is scaled to an
    ellipsoid.
    random() is used to generate random numbers, seeded with time().
    How to compile:
    gcc -o sphere sphere.c -lm

    Reference: J. O'Rourke, Computational Geometry Column 31,
    Internat. J. Comput. Geom. Appl. 7 379--382 (1997);
    Also in SIGACT News, 28(2):20--23 (1997), Issue 103.

    Written by Joseph O'Rourke and Min Xu, June 1997.
    Used in the textbook, "Computational Geometry in C."
    Questions to orourke@cs.smith.edu.
    --------------------------------------------------------------------
    This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
    redistributed in its entirety provided that this copyright notice is
    not removed.
    --------------------------------------------------------------------
  */
  double x, y, z, w, t;
  double R = 100.0;	/* default radius */
  int r;		/* true radius */
  int r1, r2, r3;	/* ellipsoid axis lengths */

  r = (int)R;
  r1 = r2 = r3 = r;

  const double resolution = 2147483647.0; // related to maximum int
 
  while (num--) {
    /* Generate a random point on a sphere of radius 1. */
    /* the sphere is sliced at z, and a random point at angle t
       generated on the circle of intersection. */
    z = 2.0 * random() / resolution - 1.0;
    t = 2.0 * M_PI * random() / resolution;
    w = sqrt( 1 - z*z );
    x = w * cos( t );
    y = w * sin( t );
    
    MakeVertex(Vector(r1*x,r2*y,r3*z));
  }
}
 
void SimplicialComplex::MakeCube (int num)
{
  clear();

  const double resolution = 2147483647.0; // related to maximum int
  while (num--) {
    /* Generate a random point on a unit cube. */
    double x = 2.0 * random() / resolution - 1.0;
    double y = 2.0 * random() / resolution - 1.0;
    double z = 2.0 * random() / resolution - 1.0;
    MakeVertex(Vector(x,y,z) * 50.0);
  }
}

void SimplicialComplex::MakeSuspensionBridge (size_t numTowers,
																							size_t sectionsPerSpan)
{
  clear();

  const double stay2TowerLength  = 100.0;
  const double tower2TowerLength = 100.0;
  const double towerHeight = 40.0;
  const double towerWidth  = 40.0;

  const double totalLength = 
    tower2TowerLength * (double(numTowers)-1.0) + stay2TowerLength * 2;

  double hpos = -totalLength/2;

  // bridge is build from some custom:
  Material &post    = MakeMaterial(125.0, 0.90, true,  0.7);
  Material &cable   = MakeMaterial(  0.0, 0.90, true,  0.7);
  Material &strand  = MakeMaterial(  0.0, 0.10, true,  0.7);
  Material &road    = MakeMaterial(tower2TowerLength, 0.90, true,  0.7);

  vector<Vertex *> cableTies;
  vector<Vertex *> roadTies;

  // lay cable stay at start of bridge:
  Vertex &stay0 = MakeVertex(Vector(hpos,towerHeight,0));
  stay0.pinned = true;
  cableTies.push_back(&stay0);
  roadTies.push_back(&stay0);

  // lay towers evenly between start and end of bridge:
  hpos += stay2TowerLength;
  for (size_t i = 0; i < numTowers; i++) {
    Vertex *t  = &MakeVertex(Vector( hpos, towerHeight,             0));
    Vertex *b0 = &MakeVertex(Vector( hpos,           0,  towerWidth/2));
    Vertex *b1 = &MakeVertex(Vector( hpos,           0, -towerWidth/2));
    b0->pinned = true;
    b1->pinned = true;
    MakeEdge(post,*b0,*t);
    MakeEdge(post,*b1,*t);
    cableTies.push_back(t);
    Vertex *r  = &MakeVertex(Vector( hpos, towerHeight*0.6, 0));
    roadTies.push_back(r);
    hpos += tower2TowerLength;
  }
  hpos -= tower2TowerLength;
  hpos += stay2TowerLength;

  // lay cable stay at end of bridge:
  Vertex &stay1 = MakeVertex(Vector( hpos, towerHeight, 0));
  stay1.pinned = true;
  cableTies.push_back(&stay1);
  roadTies.push_back(&stay1);

  // run cable between cable stays:
  vector<Edge *> cables;
  for (size_t i = 1; i < cableTies.size(); i++)
    cables.push_back(&MakeEdge(cable,*cableTies[i-1],*cableTies[i]));

  // run road between cable stays:
  vector<Edge *> roads;
  for (size_t i = 1; i < roadTies.size(); i++)
    roads.push_back(&MakeEdge(road,*roadTies[i-1],*roadTies[i]));

  // note num vertices for mass rescaling, below:
  const int oldNumV = vertices.size();

  // refine the cables and roads:
  ChildList newCables;
  for (size_t i = 0; i < cables.size(); i++) {
    nsect(*cables[i],sectionsPerSpan, newCables);
    edges.erase(edges.find(cables[i]->objID));
  }
  ChildList newRoads;
  for (size_t i = 0; i < roads.size(); i++) {
    nsect(*roads[i],sectionsPerSpan, newRoads);
    edges.erase(edges.find(roads[i]->objID));
  }

  if (debug) cerr << cables.size() << " " << roads.size() << "\n";
  if (debug) cerr << newCables.size() << " " << newRoads.size() << "\n";
  
  // link strands from the cable to the road, but not at ends:
  if (newCables.size() != newRoads.size())
    throw "cable/road size mismatch";
  ChildList::iterator cit = newCables.begin();
  ChildList::iterator rit = newRoads.begin();
  for (;
       cit != newCables.end() && rit != newRoads.end();
       cit++, rit++) {
    MakeEdge(strand,*(*cit)->vertex[0],*(*rit)->vertex[0]);
  }

  // remass everything:
  const int newNumV = vertices.size();
  for (VertexMapIt it = vertices.begin(); it != vertices.end(); it++)
    it->second.mass *= (1.0*oldNumV)/(1.0*newNumV);

  // sortof prestress everything:
  if (true) {
    Harmonize(post);
    Harmonize(cable);
    Harmonize(road);
    Harmonize(strand);
  }
  CenterOnOrigin();
}

// nsect() doesn't work like bisect()/DivisionCleanup()
//   - handles multiple divisions in one step
//   - tracks children for calling function
void SimplicialComplex::nsect (Edge &e, size_t num, ChildList &childList)
{
  if (num < 1)
    throw "bogus nsect() arg";

  Vector along = e.vertex[1]->pos - e.vertex[0]->pos;

  if (debug) cerr << "num: " << num << "\n";
  if (debug) cerr << "len: " << along.magnitude() << "\n";

  vector<Vertex *> vertex;
  vertex.push_back(e.vertex[0]);
  for (size_t i = 0; i+1 < num; i++) {
    Vector v = e.vertex[0]->pos;
    v += along * (i+1.0) / num;
    vertex.push_back(&MakeVertex(v));
  }
  vertex.push_back(e.vertex[1]);

  if (vertex.size() != num+1)
    throw "nsect(): unexpected number of vertex";

  for (size_t i = 0; i < vertex.size()-1; i++) {
    Edge &newE = MakeEdge(*e.material,*vertex[i],*vertex[i+1]);
    newE.springFactor = e.springFactor / num;
    childList.push_back(&newE);
  }
}

void SimplicialComplex::MakePiston ()
{
  MakeTetrahedron();
  if (vertices.size() < 4)
    throw "too few vertices";
  VertexMapIt it = vertices.begin();
  it++;
  for (; it != vertices.end(); it++)
    it->second.pinned = true;
  CenterOnOrigin();
}

void SimplicialComplex::MakeUnitCube (double scale)
{
  if (scale <= 0.0) 
    throw "inappropriate scale for MakeUnitCube()";
  clear();
  MakeVertex(Vector(    0,    0,    0));
  MakeVertex(Vector(scale,    0,    0));
  MakeVertex(Vector(    0,scale,    0));
  MakeVertex(Vector(    0,    0,scale));
  MakeVertex(Vector(    0,scale,scale));
  MakeVertex(Vector(scale,    0,scale));
  MakeVertex(Vector(scale,scale,    0));
  MakeVertex(Vector(scale,scale,scale));
  CenterOnOrigin();
}

void SimplicialComplex::DoublePos ()
{
  for (VertexMapIt it = vertices.begin(); it != vertices.end(); it++)
    it->second.pos *= 2.0;
  physics.CalculateEpiphenomena();
}

void SimplicialComplex::HalvePos ()
{
  for (VertexMapIt it = vertices.begin(); it != vertices.end(); it++)
    it->second.pos /= 2.0;
  physics.CalculateEpiphenomena();
}

bool SimplicialComplex::Equiangulate ()
{
  // this idea is borrowed from SurfaceEvolver:
  //    - for each edge whith has two neighboring faces (all of them for now)
  //         - if switching the diagonal makes the triangles more equiangular
  //             - do it
  // hypothesis: can do this based on edge lengths...

  cerr << "equiangulate:\n";
  
  typedef hash_set<int> hsi;
  hsi edgesToDel;
  hsi facesToDel;

  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    Edge &e = it->second;
    if (e.face[0] == NULL || e.face[1] == NULL) continue;
    Face &f0 = *e.face[0];
    Face &f1 = *e.face[1];

    Vertex *a = e.vertex[0];
    Vertex *b = e.vertex[1];

    // find the vertex for each adjacent face which is *not* on this edge
    Vertex *v0 = &f0.vertexOpposite(e);
    Vertex *v1 = &f1.vertexOpposite(e);

    if ((v1->pos - v0->pos).magnitude() >= e.length()) continue;

    cerr << "   relinking across edge: " << e.objID << "\n";
    facesToDel.insert(f0.objID);
    facesToDel.insert(f1.objID);
    for (int i = 0; i < Face::numEdges; i++) {
      // yep, this includes e.objID twice, but it's a set!
      if (f0.edge[i] != NULL)
				edgesToDel.insert(f0.edge[i]->objID);
      if (f1.edge[i] != NULL)
				edgesToDel.insert(f1.edge[i]->objID);
    }

    // make new faces and edges:
    Face &newF0 = MakeFace(*e.material,v0,v1, b);
    Face &newF1 = MakeFace(*e.material, a,v1,v0);

    // tweak out redundant edge:
    Edge *e0 = NULL;
    Edge *e1 = NULL;
    for (int i = 0; i < Face::numEdges; i++) {
      if (newF0.edge[i]->vertex[0] != a && newF0.edge[i]->vertex[1] != a)
				e0 = newF0.edge[i];
      if (newF1.edge[i]->vertex[0] != a && newF1.edge[i]->vertex[1] != a)
				e1 = newF1.edge[i];
    }
    if (e0 == NULL) throw "dang e0";
    if (e1 == NULL) throw "dang e1";
    if (e0 == e1) throw "dang e0 and e1";

    if (false) {
      bool foundIt = false;
      for (int i = 0; i < Face::numEdges; i++) {
				if (newF1.edge[i] == e1) {
					if (foundIt) throw "found it twice";
					foundIt = true;
					newF1.edge[i] = e0;
				}
      }
      if (!foundIt) throw "didn't find it";
    
      e1->face[0] = e1->face[1] = NULL;
      edgesToDel.insert(e1->objID);
    }
    
    break;  // this method can only handle one at a time so far
  }

  for (EdgeMapIt it = edges.begin(); it != edges.end(); it++) {
    Edge &e = it->second;
    for (int i = 0; i < Edge::numFaces; i++) {
      if (e.face[i] != NULL && 
					facesToDel.find(e.face[i]->objID) != facesToDel.end())
				e.face[i] = NULL;
    }
  }
  for (FaceMapIt it = faces.begin(); it != faces.end(); it++) {
    Face &f = it->second;
    for (int i = 0; i < Face::numEdges; i++) {
      if (f.edge[i] != NULL && 
					edgesToDel.find(f.edge[i]->objID) != edgesToDel.end())
				f.edge[i] = NULL;
    }
  }

  for (hsi::iterator it = edgesToDel.begin(); it != edgesToDel.end(); it++) {
    const int objID = *it;
    cerr << "   deleting edge: " << objID << "\n";
    edges.erase(edges.find(objID));
  }
  for (hsi::iterator it = facesToDel.begin(); it != facesToDel.end(); it++) {
    const int objID = *it;
    cerr << "   deleting face: " << objID << "\n";
    faces.erase(faces.find(objID));
  }

  cerr << "equiangulate done: calling CalculateEpiphenomena()\n";
  physics.CalculateEpiphenomena();

  return (edgesToDel.size() != 0);
}


void SimplicialComplex::MakeTensegrityFoo ()
{
  clear();

  // working on Kenneth Snelson's "Needle Tower"

  // make materials
  Material &base  = MakeMaterial(  0.0, 0.20, true,  0.1);
  Material &cross = MakeMaterial(100.0, 0.20, true,  0.1);

  const double scale = 50;
  const int num = 6;
  const double angle = 2 * PI / num;

  vector<Vector> ref(num,Vector(0,0,0));
  for (int n = 0; n < num; n++)
    ref[n] = Vector( cos(n*angle), sin(n*angle), 0) * scale;

  const int levels = 6;

  vector<vector<Vertex *> > v (num);
  for (int l = 0; l < levels; l++) {
    v[l].resize(num);
    for (int n = 0; n < num; n++) {
      v[l][n] = &MakeVertex(ref[n] + Vector(0,0,l*scale));
      if (l == 0 || l == levels-1)
				v[l][n]->pinned = true;
    }
  }

  // compressive elements:
  for (int l = 1; l < levels; l++) {
    for (int n = 0; n < num; n++) {
      if (n%2 == 0) continue;
      if (l%2 == 0)
				MakeEdge(cross,*v[l-1][n],*v[l][(n+num+1)%num]);
      else
				MakeEdge(cross,*v[l-1][n],*v[l][(n+num-1)%num]);
    }
  }

  // latitudinal tensive elements:
  for (int l = 1; l < levels-1; l++) {
    for (int n = 0; n < num; n++) {
      MakeEdge(base,*v[l][n],*v[l][(n+1)%num]);
    }
  }

  // longtudinal tensive elements:
  for (int l = 1; l < levels; l++) {
    for (int n = 0; n < num; n++) {
      //MakeEdge(base,*v[l-1][n],*v[l][n]);
    }
  }

  if (false) {
    // primitive prism
    Vertex &v0 = MakeVertex(Vector( 0, 0, 0));
    Vertex &v1 = MakeVertex(Vector(10, 0, 0));
    Vertex &v2 = MakeVertex(Vector(10,10, 0));

    Vertex &u0 = MakeVertex(Vector( 0, 0,10));
    Vertex &u1 = MakeVertex(Vector(10, 0,10));
    Vertex &u2 = MakeVertex(Vector(10,10,10));

    MakeEdge(cross,v0,u1);
    MakeEdge(cross,v1,u2);
    MakeEdge(cross,v2,u0);

    MakeEdge(base,v0,u0);
    MakeEdge(base,v1,u1);
    MakeEdge(base,v2,u2);

    MakeEdge(base,v0,v1);
    MakeEdge(base,v1,v2);
    MakeEdge(base,v2,v0);

    MakeEdge(base,u0,u1);
    MakeEdge(base,u1,u2);
    MakeEdge(base,u2,u0);
  }

  CenterOnOrigin();
}

