// geodefunc.cpp
//
// Exposes a bunch of functionality:
//    for reducing compilation, link, and especially Makefile complexity.
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

#include "profile.hpp"
#include "render.hpp"
#include "hull.hpp"
#include "simplicial.hpp"
#include "integrator.hpp"
#include "io.hpp"
#include <math.h>
#include <time.h>
#include <iostream>
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

// it's just one function.... lame excuse for not having a header file:
void geodetest ();

// forward declaration of some functions in this file:
int main (int argc, char *argv[]);
int testmain (int argc, char *argv[]);
void timetest ();

int testmain (int argc, char *argv[])
{
  for (int i = 0; i < argc; i++) {
    cerr << "argv[" << i << "]: " << argv[i] << "\n";
  }
  // geodefunc is meant to be called from a script of the form:
  // 
  //    #/usr/bin/tcsh
  //    geodefunc $0 $1 $2 $3 $4 $5 $6 $7 $8 $9
  // 
  // which means:
  //    - argv[0] is name of this program (usually "geodefunc")
  //    - argv[1] is the full path name of the calling script
  if (argc < 2) {
    cerr << "Must call with at least one argument\n";
    return 1;
  }
  string command = string(argv[1]);
  command = command.substr(1+command.rfind('/'));
  //    - the remaining argv are optional arguments, and we accept many types:
  const int numParameters = 9;
  vector<int>    argInt    (numParameters);
  vector<double> argDouble (numParameters);
  const int numArgs = argc-2;
  for (int i = 0; i < numParameters; i++) {
    argInt   [i] = 0;
    argDouble[i] = 0.0;
    if (i < numArgs) {
      istringstream inInt    (argv[i+2]);
      istringstream inDouble (argv[i+2]);
      inInt    >> argInt   [i];
      inDouble >> argDouble[i];
    }
  }

  cerr << "doing " << command;
  for (int i = 0; i < numParameters; i++)
    cerr << " " << argInt[i] << "(" << argDouble[i] << ")";
  cerr << "\n";

  // some commands take random seed: use time
  //srandom((unsigned int)time(NULL));
  srandom(0);

  if (command == "geode") {
    SimplicialComplex sc;
    cin >> sc;
    Render render (sc, argInt[0]!=1, argInt[1]!=1, argInt[2]!=1, argInt[3]!=0);
    render.run(argc, argv);
  }
  else if (command == "null") {
    SimplicialComplex sc;
    cout << sc;
  }
  else if (command == "testio") {
    SimplicialComplex sc;
    cin >> sc;
    cout << sc;
    Profile::summarize(cerr);
  }
  else if (command == "hull") {
    SimplicialComplex sc;
    cin >> sc;
    ConvexHull(sc).ConstructHull();
    cout << sc;
  }
  else if (command == "physics") {
    SimplicialComplex sc;
    const Time t0;
    cin >> sc;
    const Time t1;
    cerr << "load:     " << (t1-t0) << " sec\n";
    Integrator *integrator;
    switch (argInt[0]) {
    case 0: integrator = new ForwardEuler(sc);  break;
    case 1: integrator = new ImprovedEuler(sc); break;
    case 2: integrator = new RungeKutta(sc);    break;
    default: throw "bogus integrator";
    }
    double step = 0.01;
    try {
      for (int i = 0; i < argInt[1]; i++)
	integrator->step(step);
    }
    catch (...) {
      delete integrator;
      throw;
    }
    const Time t2;
    cerr << "physics:  " << (t2-t1) << " sec\n";
    //Profile::summarize(cerr);
  }
  else if (command == "timetest") {
    timetest();
  }
  else if (command == "injectgas") {
    SimplicialComplex sc;
    cin >> sc;
    sc.physics.InjectGas(argDouble[0],argDouble[1]);
    cout << sc;
  }
  else if (command == "embedgas") {
    SimplicialComplex sc;
    cin >> sc;
    sc.physics.EmbedGas(argDouble[0],argDouble[1]);
    cout << sc;
  }
  else if (command == "equiangulate") {
    SimplicialComplex sc;
    sc.MakeCube(argInt[0]);
    ConvexHull(sc).ConstructHull();
    while (sc.Equiangulate());
    cout << sc;
  }
  else if (command == "energy") {
    SimplicialComplex sc;
    cin >> sc;
    cerr << "kin: " << sc.physics.kineticEnergy << "\n";
    cerr << "grv: " << sc.physics.gravitationalEnergy << "\n";
    cerr << "cnt: " << sc.physics.centeringEnergy << "\n";
    cerr << "spr: " << sc.physics.springEnergy << "\n";
    cerr << "int: " << sc.physics.internalPressureEnergy << "\n";
    cerr << "ext: " << sc.physics.externalPressureEnergy << "\n";
    cerr << "tot: " << sc.physics.totalEnergy << "\n";
    cout << sc;
  }
  else if (command == "mass") {
    SimplicialComplex sc;
    cin >> sc;
    double mass = 0.0;
    Vector velocity(0,0,0);
    Vector momentum(0,0,0);
    double velocity_magnitude = 0.0;
    double momentum_magnitude = 0.0;
    for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++) {
      Vertex &v = it->second;
      mass += v.mass;
      velocity += v.vel;
      momentum += v.vel * v.mass;
      velocity_magnitude += v.vel.magnitude();
      momentum_magnitude += v.vel.magnitude() * v.mass;
    }
    cerr << "mass:               " << mass << "\n";
    cerr << "velocity:           " << velocity << "\n";
    cerr << "momentum:           " << momentum << "\n";
    cerr << "velocity magnitude: " << velocity_magnitude << "\n";
    cerr << "momentum magnitude: " << momentum_magnitude << "\n";
    cout << sc;
  }
  else if (command == "harmonize") {
    SimplicialComplex sc;
    cin >> sc;
    sc.Harmonize();
    cout << sc;
  }
  else if (command == "deharmonize") {
    SimplicialComplex sc;
    cin >> sc;
    sc.DeHarmonize();
    cout << sc;
  }
  else if (command == "bisect") {
    SimplicialComplex sc;
    cin >> sc;
    for (int i = 0; i < argInt[0]; i++)
      sc.BisectEdges();
    cout << sc;
  }
  else if (command == "subdivide") {
    SimplicialComplex sc;
    cin >> sc;
    for (int i = 0; i < argInt[0]; i++)
      sc.SubdivideFaces();
    cout << sc;
  }
  else if (command == "sphericize") {
    SimplicialComplex sc;
    cin >> sc;
    sc.Sphericize();
    cout << sc;
  } 
 else if (command == "subsphericize") {
    SimplicialComplex sc;
    cin >> sc;
    for (int i = 0; i < argInt[0]; i++) {
      sc.SubdivideFaces();
      sc.Sphericize();
      sc.Harmonize();
    }
    cout << sc;
  }
  else if (command == "diamond") {
    SimplicialComplex sc;
    cin >> sc;
    for (int i = 0; i < argInt[0]; i++)
      sc.DiamondSubdivide();
    cout << sc;
  }
  else if (command == "cross") {
    SimplicialComplex sc;
    cin >> sc;
    sc.AddCrossScaffolding(argInt[0]!=0);
    cout << sc;
  }
  else if (command == "dual") {
    SimplicialComplex sc;
    cin >> sc;
    sc.AddDualScaffolding();
    cout << sc;
  }
  else if (command == "sphere") {
    SimplicialComplex sc;
    sc.MakeSphere(argInt[0]);
    cout << sc;
  }
  else if (command == "spiral") {
    SimplicialComplex sc;
    sc.MakeSpiral(argInt[0]);
    cout << sc;
  }
  else if (command == "cube") {
    SimplicialComplex sc;
    sc.MakeCube(argInt[0]);
    cout << sc;
  }
  else if (command == "bridge") {
    SimplicialComplex sc;
    sc.MakeSuspensionBridge(argInt[0],argInt[1]);
    cout << sc;
  }
  else if (command == "tensegrity") {
    SimplicialComplex sc;
    sc.MakeTensegrityFoo();
    cout << sc;
  }
  else if (command == "piston") {
    SimplicialComplex sc;
    sc.MakePiston();
    ConvexHull(sc).ConstructHull();
    cout << sc;
  }
  else if (command == "unitcube") {
    SimplicialComplex sc;
    sc.MakeUnitCube(argDouble[0]);
    ConvexHull(sc).ConstructHull();
    cout << sc;
  }
  else if (command == "tetra") {
    SimplicialComplex sc;
    sc.MakeTetrahedron();
    cout << sc;
  }
  else if (command == "torus") {
    SimplicialComplex sc;
    sc.MakeTorus(static_cast<SimplicialComplex::TorusType>(argInt[0]),
		 argInt[1],
		 argInt[2]);
    cout << sc;
  }
  else if (command == "pin") {
    SimplicialComplex sc;
    sc.MakePinnedWeight();
    cout << sc;
  }
  else if (command == "binary") {
    SimplicialComplex sc;
    sc.MakeBinaryWeight();
    cout << sc;
  }
  else if (command == "trinary") {
    SimplicialComplex sc;
    sc.MakeTrinaryWeight();
    cout << sc;
  }
  else if (command == "pin2") {
    SimplicialComplex sc;
    sc.MakeTwicePinnedWeight();
    cout << sc;
  }
  else if (command == "octa") {
    SimplicialComplex sc;
    sc.MakeOctahedron();
    cout << sc;
  }
  else if (command == "icosa") {
    SimplicialComplex sc;
    sc.MakeIcosahedron();
    cout << sc;
  }
  else if (command == "testo") {
    SimplicialComplex sc;
    while (true) {
      Vector v(0,0,0);
      eat(cin);
      if (!cin) break;
      cin >> v.x;
      eat(cin);
      if (!cin) throw "data not in triples: y missing";
      cin >> v.y;
      eat(cin);
      if (!cin) throw "data not in triples: z missing";
      cin >> v.z;
      sc.MakeVertex(v);
    }    
    cout << sc;
  }
  else if (command == "oldschool") {
    SimplicialComplex sc;
    cin >> sc;
    for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++) {
      Vector &p = it->second.pos;
      cout << int(p.x+0.5) << " " 
	   << int(p.y+0.5) << " " 
	   << int(p.z+0.5) << "\n";
    }
  }
  else if (command == "checks") {
    SimplicialComplex sc;
    cin >> sc;
    sc.CheckEuler(sc.vertices.size(),sc.edges.size(),sc.faces.size());
    sc.EdgeOrderOnFaces();
    sc.CheckEndpts();
    if (!ConvexHull::isConvex(sc))
      cerr << "not convex\n";
    sc.Consistency();
    sc.myConsistency();
  }
  else if (command == "geodetest") {
    geodetest();
  }
  else {
    cerr << "unrecognized command: " << command << "\n";
  }

  return 0;
}



int main (int argc, char *argv[])
{
  try {
    return testmain(argc,argv);
  }
  catch (std::string str) {
    cerr << "exception: " << str << "\n";
  }
  catch (const char *str) {
    cerr << "exception: " << str << "\n";
  }
  catch (...) {
    cerr << "unrecognized exception\n";
  }
}

void wait (double sec);
void fuu ();
void bar ();
void baz ();

void timetest ()
{
  Time start;
  for (int i = 0; i < 3; i++)
  {
    Profile prof ("timetest");
    fuu();
    bar();
    bar();
  }
  cerr << "total time: " << (Time()-start) << "\n";
  Profile::summarize(cerr);
}

void fuu ()
{
  Profile prof ("fuu");
  wait(0.03);
}

void bar ()
{
  Profile prof ("bar");
  wait(0.02);
  baz();
}

void baz ()
{
  Profile prof ("baz");
  wait(0.07);
}

void wait (double sec)
{
  Time start;
  while ((Time()-start).secs() < sec);
}
