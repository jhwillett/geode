// physics.hpp
//
// physics module: force representations, etc... 
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

#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "algebra.hpp"
#include <vector>

class SimplicialComplex;

class Physics
{
  SimplicialComplex &sc;

public:
  static const double stpPressure;     // Pascals (N/m^2): wow!
  static const double stpTemperature;  // Kelvin, eq. to 0 deg.C

  // saved fields:
  bool   viscosityActive;
  double viscosity;             // Newtons per (meters per second)
  bool   centeringActive;
  double centering;             // Newtons per (meter kilogram)
  bool   springsActive;
  bool   springDampingActive;
  bool   internalPressureActive;        
  double moles;                 // quantity of gas, in moles
  double idealGasConstant;      // magic number, ideal gas constant
  double internalTemperature;   // in degrees Kelvin
  double gasK;                  // == P*V^gamma, const for adiabatic delta-vol
  bool   externalPressureActive;
  double externalPressure;      // in Pascals (Newtons per (meter^2))
  double externalTemperature;   // in degrees Kelvin
  bool   vertexLimitActive;
  double vertexLimit;           // in meters
  bool   gravityActive;
  Vector gravity;               // constant in pos, vel: N/kg

  bool   surfaceTensionActive;
  double surfaceTension;        // N/m or J/m^2

  std::vector<double> possibleGammas;
  int                 currentGamma;

  // pseudofield: gamma of gas (5/3 if monoatomic, etc)
  inline double gamma () const { return possibleGammas[currentGamma]; }

  // unsaved fields, updated by CalculateEpiphenomena
  double pressure;        // pressure of internal gas, in Pa == N/m^2
  double totalVolume;     // volume of faces, in m^3, may be negative
  double totalArea;       // area of faces, in m^2, always nonnegative
  double totalLength;     // length of edges, in m, always nonnegative
  Vector centerOfMass;    // center of Vertex masses, in m
  double totalMass;       // sum of Vertex masses, in kg

  double kineticEnergy;
  double gravitationalEnergy;
  double centeringEnergy;
  double springEnergy;
  double surfaceEnergy;
  double internalPressureEnergy;
  double externalPressureEnergy;
  double totalEnergy;

public:
  Physics (SimplicialComplex &sci);
  ~Physics ();

  // particular support for integrators:
  void CalculateEpiphenomena ();
  void CalculateForces ();
  void CalculateLimitEffects ();

  void InjectGas (double molesi, double temperaturei);
  void ReleaseGas (double molesi);
  void EmbedGas (double pressurei, double temperaturei);
  void OpenValve ();
  void AddSpin (double factor);
  void AddVel (const Vector &vel);
  void AddImpact ();
  void Jiggle (double range);
  void Cool (double factor);

  void NextGamma ();
};

#endif // #ifndef PHYSICS_HPP
