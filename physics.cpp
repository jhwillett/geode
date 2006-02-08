// physics.cpp
//
// physics on SimplicialComplexes
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

// notes to me:
// From http://www.citycollegiate.com/viscosity.htm:
//
// Two definitions of surface tension....
//
// "Perpendicular force acting on the uint length of a liquid is called S.T/"
// Surface Tension = force/length
// g = F/L
//
// "Energy per unit area on the surface of a liquid is called S.T."
// s = energy/area
//
// Units of surface tension : N/m or J/m^2
//
// Units of surface energy: J/m
//
// also from: http://hyperphysics.phy-astr.gsu.edu/hbase/surten.html
// 
// Surface energy may have two parts: stretching energy and bending energy.

#include "physics.hpp"
#include "simplicial.hpp"
#include "profile.hpp"
#include <iostream>
using std::cerr;

const double Physics::stpPressure = 101325.0;   // Pascals (N/m^2): wow!
const double Physics::stpTemperature = 273.15;  // Kelvin, eq. to 0 deg.C

static const bool debug = false;

Physics::Physics (SimplicialComplex &sci) : 
  sc(sci),
  gravityActive(false),
  gravity(0,-10,0),
  centerOfMass(0,0,0)
{
  viscosityActive = false;
  viscosity = 0.03;

  centeringActive = true;
  centering = 0.2;

  springsActive = true;
  springDampingActive = true;
  
  surfaceTensionActive = false;
  surfaceTension = 0.01;

  internalPressureActive = true;
  moles = 0.0;
  internalTemperature = 0.0;
  gasK = 0.0;
  idealGasConstant = 8.319;

  possibleGammas.push_back(5.0/3.0); // monoatomic
  possibleGammas.push_back(7.0/5.0); // diatomic, mid temp domain (no vib)
  possibleGammas.push_back(9.0/7.0); // triatomic, mid temp domain (no vib)
  currentGamma = 0;

  externalPressureActive = true;
  externalPressure = 0.00; 
  externalTemperature = 273.15;   // in Kelvin == 0 deg. Centigrade

  vertexLimitActive = true;
  vertexLimit = 1000;

  totalMass = 0;
  pressure = 0.0;
  totalVolume = 0.0;
  totalArea = 0;
  totalLength = 0;

  kineticEnergy = 0;
  gravitationalEnergy = 0;
  centeringEnergy = 0;
  springEnergy = 0;
  surfaceEnergy = 0;
  internalPressureEnergy = 0;
  externalPressureEnergy = 0;
  totalEnergy = 0;

  CalculateEpiphenomena();
}

Physics::~Physics ()
{
}

void Physics::InjectGas (double molesi, double temperaturei)
{
  // inject moles at kelvin
  CalculateEpiphenomena();
  molesi = max(molesi,0.0);
  temperaturei = max(temperaturei,0.0);
  double current_NT = moles * internalTemperature;
  double inject_nt = molesi * temperaturei;
  double total_NT = current_NT + inject_nt;
  moles = moles + molesi;
  internalTemperature = total_NT / moles;
  gasK = 
    moles * idealGasConstant * internalTemperature 
    * pow( totalVolume, gamma()-1.0);
  CalculateEpiphenomena();
}

void Physics::ReleaseGas (double molesi)
{
  CalculateEpiphenomena();
  molesi = max(molesi, 0.0);
  moles -= molesi;
  if (moles <= 0) {
    moles = 0;
    internalTemperature = 0;
  }
  gasK =
    moles * idealGasConstant * internalTemperature 
    * pow( totalVolume, gamma()-1.0);
  CalculateEpiphenomena();
}

void Physics::EmbedGas (double pressurei, double temperaturei)
{
  externalPressure = pressurei;
  externalTemperature = temperaturei;
  CalculateEpiphenomena();
}

void Physics::OpenValve ()
{
  CalculateEpiphenomena();
  internalTemperature = externalTemperature;
  moles = 
    (externalPressure * totalVolume) / 
    (idealGasConstant * externalTemperature);
  gasK =
    moles * idealGasConstant * internalTemperature 
    * pow( totalVolume, gamma()-1.0);
  CalculateEpiphenomena();
}

void Physics::NextGamma ()
{
  CalculateEpiphenomena();
  currentGamma += 1;
  currentGamma %= possibleGammas.size();
  gasK = 
    moles * idealGasConstant * internalTemperature 
    * pow( totalVolume, gamma()-1.0);
  CalculateEpiphenomena();
}

void Physics::AddSpin (double factor)
{
  // add spin velocity around y-axis proportionate to |x,y|
  // NOTE: this is not a force or torque or even angular velocity
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    if (v->second.pinned) continue;
    double scale = factor/100.0;
    v->second.vel.x += scale * (- v->second.pos.y);
    v->second.vel.y += scale * (  v->second.pos.x);
  }
}

void Physics::AddVel (const Vector &vel)
{
  // add velocity 
  // NOTE: this is not a force, just free velocity
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    if (v->second.pinned) continue;
    v->second.vel += vel;
  }
}

void Physics::AddImpact ()
{
  Vector ref (1,0,0);
  for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++) {
    Vertex &v = it->second;
    if (v.pinned) continue;
    while (v.pos * ref >= 0.8 * v.pos.magnitude())
      v.pos += Vector(-10,0,0);
  }
}

void Physics::Jiggle (double range)
{
  // add random vectors of magnitude [0..range] to each vertex
  //    - NOTE: these vectors are skewed toward corners of cube
  //            a spherical gaussian might have better properties
  //            in terms of energy, entropy, and frequency distribution
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    if (v->second.pinned) continue;
    Vector offset (rand()-RAND_MAX/2,
									 rand()-RAND_MAX/2,
									 rand()-RAND_MAX/2);
    offset /= RAND_MAX/2;
    offset *= range;
    v->second.pos += offset;
  }
}

void Physics::Cool (double factor)
{
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    v->second.vel *= factor;
  }
}

void Physics::CalculateForces ()
{
  Profile prof ("CalculateForces");

  // NOTE: most code here has corresponding code in CalculateEpiphenomena()
  // NOTE: (the energy tabulations) which must be kept in sync.

  // some set-up data:
  const Vector zero(0,0,0);

  // per-Vertex effects:
  const Vector centering_accel = centerOfMass * (-centering);
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vertex = v->second;

    // clear all previously accumulated forces:
    vertex.ni_force = zero;

    // centering force to keep the subject on camera:
    //   - center of mass attracted to origin
    //   - force is proportionate to mass (like gravity, so all same delta-v)
    //   - spring-like mechanics (not really gravity)
    if (centeringActive)
      vertex.ni_force += centering_accel * vertex.mass;

    // Real (near-surface of earth) gravity:
    //   - center of mass attracted along gravity
    //   - force is proportionate to mass
    if (gravityActive)
      vertex.ni_force += gravity * vertex.mass;

    // Medium is viscous: slow each vertex some percent
    //    - friction force against velocity
    if (viscosityActive)
      vertex.ni_force -= vertex.ni_vel * viscosity;
  }

  // per-Edge effects:
  for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++) {
    const Edge &e = it->second;
    const Material &material = *e.material;
    Vertex &v0 = *e.vertex[0];
    Vertex &v1 = *e.vertex[1];

    // Edges are springs:
    // each edge imposes a force on each of its vertices:
    //    - force is along the edge
    //    - force is proportionate to edge's magnitude
    //    - dampening is proportionate to dot with velocity
    if (springsActive) {
      double ks = material.springCoef / e.springFactor;
      double r = material.restLength * e.springFactor;
      Vector L = v0.pos - v1.pos;
      double l = L.magnitude();

      double factorS = -ks*(l-r);
      if (!material.compressible && l < r)
				factorS = 0;

      double factorD = 0.0;
      if (springDampingActive) {
				double kd = (material.damper) * e.springFactor;
				factorD = -kd*((v0.ni_vel-v1.ni_vel)*L)/l;
      }

      Vector F = L * ((factorS + factorD) / l);

      v0.ni_force += F;
      v1.ni_force -= F;
    }
  }

  // per-Face effects:
  for (FaceMapIt it = sc.faces.begin(); it != sc.faces.end(); it++) {
    const Face &f = it->second;
    const Vector areaNormal = f.areaNormal();
    
    // There is surface tension:
    if (surfaceTensionActive) {
      const Vector center = f.center();
      for (int i = 0; i < Face::numVertices; i++) {
				const Vector off = center - f.vertex[i]->pos;
				if (false) {
					f.vertex[i]->ni_force += off * surfaceTension;// * 2.0 / 3.0;
				}
				else {
					const double factor = surfaceTension * areaNormal.magnitude();
					f.vertex[i]->ni_force += factor * off/off.magnitude();
				}
      }      
    }

    // The interior and exterior of each face is pressurized:
    //    - force is normal to the face 
    //    - force is proportionate to face area
    //    - internal component inversely prop. to total volume
    //    - external component independent of volume
    //    - both components zeroed when volume is nonpositive for safety
    //    - DANGER: negative interior pressure corresponds to infinite energy
    //    - DANGER: negative exterior pressure corresponds to infinite energy
    //    - DANGER: w/o inside-out correction, positive exterior pressure on
    //              inside-out structures also corresponds to infinite energy

    const double internal_factor = 
      (internalPressureActive && totalVolume > 0.0) 
      ? pressure 
      : 0;
    const double external_factor = 
      (externalPressureActive && totalVolume > 0.0) 
      ? -externalPressure
      : 0;

    Vector force = areaNormal * (internal_factor + external_factor);

    // pressure is on the face: distribute it to the vertices:
    // i.e. the "mass" of this face is the masses at the corners.
    force /= Face::numVertices;
    for (int i = 0; i < Face::numVertices; i++)
      f.vertex[i]->ni_force += force;
  }
}

void Physics::CalculateLimitEffects ()
{
  Profile prof ("CalculateLimitEffects");

  // Edge-of-universe effects:
  //    - cap position values for safety
  //    - position values can kill X Server!
  //    - also null velocity values
  //        - we don't want invisible momentum
  //        - the edge-of-the-universe is *in*elastic!
  if (vertexLimitActive) {
    for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
      double mag = v->second.pos.magnitude();
      if (mag > vertexLimit) {
				v->second.pos *= vertexLimit/mag;
				v->second.vel *= 0;
      }
    }
  }

  // Floor effects:
  //    - cheezy floor collision detection and response:
  //    - this kind of sucks
  if (false & vertexLimitActive) {
    for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
      if (v->second.pos.y >= -100.0) continue;
      v->second.pos.y = -100.0;
      v->second.vel.y = max(0.0,v->second.vel.y);
    }
  }
}

void Physics::CalculateEpiphenomena ()
{
  Profile prof ("CalculateEpiphenomena");

  // sets up caches of expensive values:
  const Vector zero(0,0,0);

  // center of mass stuff:
  centerOfMass = zero;
  totalMass = 0;
  for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++) {
    Vertex &v = it->second;
    if (v.pinned) continue;
    totalMass += v.mass;
    centerOfMass += v.pos * v.mass;
  }
  centerOfMass /= totalMass;

  // totalLength based on all edge's position:
  totalLength = 0;
  for (EdgeMapIt e = sc.edges.begin(); e != sc.edges.end(); e++) {
    const Edge &edge = e->second;
    totalLength += edge.length();
  }
  if (totalLength < 0.0) {
    cerr << "unexpected negative totalLength\n";
    throw "unexpected negative totalLength";
  }

  // totalArea and totalVolume based on all face's position:
  totalVolume = 0;
  totalArea = 0;
  for (FaceMapIt it = sc.faces.begin(); it != sc.faces.end(); it++) {
    const Face &f = it->second;
      
    // Volume of a polyhedron is sum of the (signed) volumes of
    // tetrahedra formed by an arbitrary point and each triangular
    // face.
    totalVolume += f.volume(zero);

    // Area of is sum of the (signed) areas of faces.
    totalArea += f.area();
  }
  if (totalArea < 0.0) {
    cerr << "unexpected negative totalArea\n";
    throw "unexpected negative totalArea";
  }
  
  // calculate current thermo values:
  if (moles <= 0.0 || totalVolume <= 0.0) {
    pressure = 0.0;
    internalTemperature = 0.0;
  }
  else {
    pressure =  gasK * pow(totalVolume, -gamma());
    internalTemperature = 
      (gasK * pow(totalVolume, 1.0-gamma())) / (moles * idealGasConstant);
  }    

  // kinetic energy, the easy one:
  //    - KE = (1/2)mv^2
  //    - factor 1/2 from integrating linear (momentum) to quadratic (energy)
  kineticEnergy = 0.0;
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vertex = v->second;
    kineticEnergy += vertex.mass * (vertex.vel * vertex.vel) / 2.0;
  }

  // gravitational force energy
  //    - arbitrary referance point: origin? vertex limit? floor?
  //    - this term can go negative!
  //    - PE = mgh
  //    - no magic factor: integrating constant (force) to linear (energy)
  gravitationalEnergy = 0.0;
  if (gravityActive) {
    Vector ref = gravity * vertexLimit / gravity.magnitude();
    gravitationalEnergy += 
      (ref-centerOfMass) * gravity * totalMass;
  }

  // centering force energy
  //    - springlike
  //    - PE = (1/2)kmx^2
  //    - factor 1/2 from integrating linear (position) to quadratic (energy)
  centeringEnergy = 0.0;
  if (centeringActive)
    centeringEnergy += 
      centering * totalMass * (centerOfMass * centerOfMass) / 2.0;

  // spring force energy
  //    - PE = (1/2)kx^2
  //    - factor 1/2 from integrating linear (position) to quadratic (energy)
  springEnergy = 0.0;
  for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++) {
    const Edge &e = it->second;
    const Material &material = *e.material;
    Vertex &v0 = *e.vertex[0];
    Vertex &v1 = *e.vertex[1];

    if (springsActive) {
      double ks = material.springCoef / e.springFactor;
      double r = material.restLength * e.springFactor;
      Vector L = v0.pos - v1.pos;
      double l = L.magnitude();

      if (!material.compressible && l < r)
				springEnergy += 0;
      else
				springEnergy += ks * (l-r) * (l-r) / 2.0;
    }
  }

  // surface energy:
  surfaceEnergy = 0.0;
  if (surfaceTensionActive)
    surfaceEnergy = totalArea * surfaceTension;

  // pressurization energy:
  internalPressureEnergy = 0.0;
  if (internalPressureActive && totalVolume > 0.0) {
    internalPressureEnergy = 
      (gasK * pow( totalVolume, 1.0-gamma())) / (gamma()-1.0);
  }

  if (true && internalPressureActive) {
    // check potential work versus internal thermal energy:
    //   - note internal energy cheaper to calculate!
    double eint = 
      moles * idealGasConstant * internalTemperature / (gamma()-1.0);
    if (fabs(eint - internalPressureEnergy) > 0.01)
      cerr << "eint and internalPressureEnergy disagree!!!!!!!\n"
					 << eint << " " << internalPressureEnergy;
  }
  
  externalPressureEnergy = 0.0;
  if (externalPressureActive && totalVolume > 0.0) {
    // NOTE: This looks like magic, but...
    //       - considering F = dE/dr, this is correct
    //       - this results in proper C.o.E with the other forces
    externalPressureEnergy = externalPressure * totalVolume;
  }

  // PRESSURE ASSUMPTIONS:
  //   - gas is ideal
  //      - molecules don't collide with each other
  //      - gamma is now parameterized
  //
  // NECESSARY: 
  //   - expansion is adiabatic 
  //      - heat has nowhere else to go or come from
  //   - can't be isothermic
  //      - that implies arbitrary heat input, not conservative of energy
  //   - can't be isovolumetric
  //      - the container is elastic and gas does no work without delta V
  //   - can't be isobaric
  //      - without some weird Temp <=> Volume thing, prob. not C.o.E.

  totalEnergy = 0.0;
  totalEnergy += kineticEnergy;
  totalEnergy += gravitationalEnergy;
  totalEnergy += centeringEnergy;
  totalEnergy += springEnergy;
  totalEnergy += surfaceEnergy;
  totalEnergy += internalPressureEnergy;
  totalEnergy += externalPressureEnergy;
}


