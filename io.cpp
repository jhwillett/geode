// io.cpp
//
// input/output functionality
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

#include "io.hpp"
#include "simplicial.hpp"
#include "profile.hpp"
#include "physics.hpp"
#include <sstream>
using std::cerr;
using std::istream;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::ios_base;
using std::map;
using __gnu_cxx::hash_map;
using std::vector;

const static bool debug = false;

istream &operator>> (istream &s, Edge &e);
ostream &operator<< (ostream &s, const Edge &e);
istream &operator>> (istream &s, Face &f);
ostream &operator<< (ostream &s, const Face &f);

// stuff which is unfortunate, due to circularity of reference
void read (istream &s, Edge &e, SimplicialComplex &sc);
void fixFaces (Edge &e, FaceMap &f);
void read (istream &s, Face &f, VertexMap &vertices, EdgeMap &edges);

// overriding standard double io: we want every bit preserved:
template<class T> class LosslessIO
{
public:
  // TODO: address big-endian vs. little-endian in LosslessIO<T>

  // reinterpret target as unsigned char array for io:
  unsigned char d[sizeof(T)/sizeof(unsigned char)];

  // constructing these generally is unnecessary: we just do it
  // below as a hook for double-checking some size equivalences.
  LosslessIO ()
  {
    if (sizeof(T) == sizeof(*this)) return;
    cerr << "unexpected LosslessIO size mismatch:\n";
    cerr << "  LosslessIO<T>: " << sizeof(*this) << "\n";
    cerr << "              T: " << sizeof(T)     << "\n";
    throw "unexpected LosslessIO size mismatch";
  }
};
namespace {
  LosslessIO<double> testDouble;
  LosslessIO<float>  testFloat;
}
template<class T> istream &operator>> (istream &s, LosslessIO<T> &l)
{
  char a, b;
  s >> a;
  if (a != '-') throw "unexpeced not -";
  for (size_t i = 0; i < sizeof(T)/sizeof(unsigned char); i++) {
    s >> a >> b;
    l.d[i] = ((a-'A')<<4) | (b-'A'); 
  }
  return s;
}
template<class T> ostream &operator<< (ostream &s, const LosslessIO<T> &l)
{
  char a, b;
  s << "-";
  for (size_t i = 0; i < sizeof(T)/sizeof(unsigned char); i++) {
    a = (l.d[i]>>4) + 'A';
    b = (l.d[i]&0x0F) + 'A';
    s << a << b;
  }
  return s;
}
template<class T> inline LosslessIO<T> &lossless (T &d)
{
  return reinterpret_cast<LosslessIO<T> &>(d);
}
template<class T> inline const LosslessIO<T> &lossless (const T &d)
{
  return reinterpret_cast<const LosslessIO<T> &>(d);
}

istream &eat (istream &is)
{
  char c;
  while (is.get(c)) {
    if (isdigit(c) || c == '-') {
      is.putback(c);
      break;
    }
    //cerr << "eating: " << c << "\n";
  }
  return is;
}

// for reading map<int,T> where T.objID is present:
template<class T> istream &operator>> (istream &s, map<int,T> &m)
{
  eat(s);
  int num;
  s >> num;
  for (int i = 0; i < num; i++) {
    T t;
    s >> t;
    m[t.objID] = t;
  }
  return s;
}
template<class T> ostream &operator<< (ostream &s, const map<int,T> &m)
{
  s << "(map<int,T> " << m.size();
  for (typename map<int,T>::const_iterator it = m.begin(); it != m.end(); it++)
    s << "   " << it->second;
  s << ")";
  return s;
}

// for io of hash_map<int,T> where T.objID is present:
template<class T> istream &operator>> (istream &s, hash_map<int,T> &m)
{
  eat(s);
  int num;
  s >> num;
  for (int i = 0; i < num; i++) {
    T t;
    s >> t;
    m[t.objID] = t;
  }
  return s;
}
template<class T> ostream &operator<< (ostream &s, const hash_map<int,T> &m)
{
  s << "(hash_map<int,T> " << m.size();
  for (typename hash_map<int,T>::const_iterator it = m.begin();
       it != m.end(); it++)
    s << "   " << it->second;
  s << ")";
  return s;
}

// for io of vector<double> (not template: need lossless() here):
istream &operator>> (istream &s, vector<double> &v)
{
  eat(s);
  int num;
  s >> num;
  v.resize(num);
  for (int i = 0; i < num; i++) {
    s >> lossless(v[i]);
  }
  return s;
}
ostream &operator<< (ostream &s, const vector<double> &v)
{
  s << "(vector<double> " << v.size();
  for (size_t i = 0; i < v.size(); i++)
    s << "   " << lossless(v[i]);
  s << ")";
  return s;
}

istream &operator>> (istream &s, Physics &p)
{
  eat(s);
  s >> p.viscosityActive;
  s >> lossless(p.viscosity);
  s >> p.centeringActive;
  s >> lossless(p.centering);
  s >> p.springsActive;
  s >> p.springDampingActive;
  s >> p.internalPressureActive;
  s >> lossless(p.moles);
  s >> lossless(p.idealGasConstant);
  s >> lossless(p.internalTemperature);
  s >> lossless(p.gasK);
  s >> p.externalPressureActive;
  s >> lossless(p.externalPressure);
  s >> lossless(p.externalTemperature);
  s >> p.vertexLimitActive;
  s >> lossless(p.vertexLimit);
  s >> p.gravityActive;
  s >> p.gravity;
  s >> p.currentGamma;
  p.possibleGammas.clear();
  s >> p.possibleGammas;
  return s;
}

istream &operator>> (istream &s, SimplicialComplex &sc)
{
  Profile prof ("SimplicialComplex >>");

  sc.clear();

  eat(s);

  s >> sc.currentObjID;
  s >> sc.physics;
  s >> sc.vertices;
  s >> sc.materials;

  eat(s);
  int numEdges;
  s >> numEdges;
  cerr << "num edges: " << numEdges << "\n";
  for (int i = 0; i < numEdges; i++) {
    Edge e;
    read(s,e,sc);
    sc.edges[e.objID] = e;
  }

  eat(s);
  int numFaces;
  s >> numFaces;
  cerr << "num faces: " << numFaces << "\n";
  for (int i = 0; i < numFaces; i++) {
    Face f;
    read(s,f,sc.vertices,sc.edges);
    sc.faces[f.objID] = f;
  }

  for (EdgeMapIt e = sc.edges.begin(); e != sc.edges.end(); e++)
    fixFaces(e->second,sc.faces);

  // update other fields:
  sc.physics.CalculateEpiphenomena();

  return s;
}

ostream &operator<< (ostream &s, const Physics &p)
{
  s << "(physics ";
  s << " ";
  s << p.viscosityActive;
  s << " ";
  s << lossless(p.viscosity);
  s << " ";
  s << p.centeringActive;
  s << " ";
  s << lossless(p.centering);
  s << " ";
  s << p.springsActive;
  s << " ";
  s << p.springDampingActive;
  s << " ";
  s << p.internalPressureActive;
  s << " ";
  s << lossless(p.moles);
  s << " ";
  s << lossless(p.idealGasConstant);
  s << " ";
  s << lossless(p.internalTemperature);
  s << " ";
  s << lossless(p.gasK);
  s << " ";
  s << p.externalPressureActive;
  s << " ";
  s << lossless(p.externalPressure);
  s << " ";
  s << lossless(p.externalTemperature);
  s << " ";
  s << p.vertexLimitActive;
  s << " ";
  s << lossless(p.vertexLimit);
  s << " ";
  s << p.gravityActive;
  s << " ";
  s << p.gravity;
  s << " ";
  s << p.currentGamma;
  s << " ";
  s << p.possibleGammas;
  s << ")";
  return s;
}

ostream &operator<< (ostream &s, const SimplicialComplex &sc)
{
  Profile prof ("SimplicialComplex <<");
  s << "(simplicialcomplex \n";
  s << sc.currentObjID;
  s << "\n";
  s << sc.physics;
  s << "\n";
  s << sc.vertices;
  s << "\n";
  s << sc.materials;
  s << "\n";
  s << sc.edges;
  s << "\n";
  s << sc.faces;
  s << ")\n";
  return s;
}

ostream &operator<< (ostream &s, const Vector &v)
{
  if (false) {
    cerr << "v: " 
				 << v.x << " " 
				 << v.y << " " 
				 << v.z;
    cerr << "\n";
    cerr << "l: " 
				 << lossless(v.x) << " " 
				 << lossless(v.y) << " "
				 << lossless(v.z);
    cerr << "\n";
  }
  s << "(";
  s << lossless(v.x) << " " << lossless(v.y) << " " << lossless(v.z);
  s << ")";
  return s;
}

istream &operator>> (istream &s, Vector &v) 
{
  eat(s);
  s >> lossless(v.x) >> lossless(v.y) >> lossless(v.z);
  eat(s);
  return s;
}

istream &operator>> (istream &s, Vertex &v) 
{
  eat(s);
  s >> v.objID;
  s >> v.pos;
  s >> v.vel;
  s >> lossless(v.mass);
  s >> v.pinned;
  return s;
}

ostream &operator<< (ostream &s, const Vertex &v)
{
  s << "(";
  s << v.objID;
  s << " ";
  s << v.pos;
  s << " ";
  s << v.vel;
  s << " ";
  s << lossless(v.mass);
  s << " ";
  s << v.pinned;
  s << ")";
  return s;
}

istream &operator>> (istream &s, Material &m)
{
  eat(s);
  s >> m.objID;
  s >> lossless(m.restLength);
  s >> lossless(m.springCoef);
  s >> lossless(m.damper);
  s >> m.compressible;
  return s;
}

ostream &operator<< (ostream &s, const Material &m)
{
  s << "(";
  s << m.objID;
  s << " ";
  s << lossless(m.restLength);
  s << " ";
  s << lossless(m.springCoef);
  s << " ";
  s << lossless(m.damper);
  s << " ";
  s << m.compressible;
  s << ")";
  return s;
}

void read (istream &s, Edge &e, SimplicialComplex &sc)
{
  eat(s);

  s >> e.objID;

  s >> lossless(e.springFactor);

  int temp;
  s >> temp;
  e.material = &sc.materials[temp];

  for (int i = 0; i < Edge::numVertices; i++) {
    s >> temp;
    e.vertex[i] = &sc.vertices[temp];
  }

  for (int i = 0; i < Edge::numFaces; i++) {
    s >> e.faceID[i];
  }
}

void fixFaces (Edge &e, FaceMap &faces)
{
  for (int i = 0; i < Edge::numFaces; i++) {
    if (e.faceID[i] != -1)
      e.face[i] = &faces[e.faceID[i]];
    else 
      e.face[i] = NULL;
  }
}

ostream &operator<< (ostream &s, const Edge &e)
{
  s << "(";
  s << e.objID;
  s << " ";
  s << lossless(e.springFactor);
  s << " ";
  if (e.material == NULL)
    s << -1;
  else
    s << e.material->objID;

  for (int i = 0; i < Edge::numVertices; i++) {
    s << " ";
    s << e.vertex[i]->objID;
  }

  for (int i = 0; i < Edge::numFaces; i++) {
    s << " ";
    if (e.face[i] == NULL)
      s << "-1";
    else
      s << e.face[i]->objID;
  }

  s << ")";
  return s;
}

void read (istream &s, Face &f, VertexMap &vertices, EdgeMap &edges)
{
  eat(s);
  s >> f.objID;
  int temp;
  for (int i = 0; i < Face::numVertices; i++) {
    s >> temp;
    if (temp != -1)
      f.vertex[i] = &vertices[temp];
    else
      f.vertex[i] = NULL;
  }
  for (int i = 0; i < Face::numEdges; i++) {
    s >> temp;
    if (temp != -1)
      f.edge[i] = &edges[temp];
    else
      f.edge[i] = NULL;
  }
}

ostream &operator<< (ostream &s, const Face &f)
{
  s << "(";
  s << f.objID;
  for (int i = 0; i < Face::numVertices; i++) {
    s << " ";
    if (f.vertex[i] == NULL)
      s << -1;
    else
      s << f.vertex[i]->objID;
  }
  for (int i = 0; i < Face::numEdges; i++) {
    s << " ";
    if (f.edge[i] == NULL)
      s << -1;
    else
      s << f.edge[i]->objID;
  }
  s << ")";
  return s;
}



