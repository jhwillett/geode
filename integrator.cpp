// integrator.cpp
//
// collection of numeric integration techniques
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

#include "integrator.hpp"
#include "simplicial.hpp"
#include "profile.hpp"
#include <iostream>
using std::cerr;

static const bool debug = false;

void ForwardEuler::step (double dt)
{
  Profile prof ("ForwardEuler.step()");
  if (debug) cerr << "ForwardEuler.step(" << dt << ") called\n";
  sc.physics.CalculateEpiphenomena();
  // set up initial ni_vel:
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vertex = v->second;
    vertex.ni_vel = vertex.vel;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;

    // apply final calculation:
    vert.vel += vert.ni_force * (dt / vert.mass);
    vert.pos += vert.vel * dt;
  }
  sc.physics.CalculateLimitEffects();
  sc.physics.CalculateEpiphenomena();
}

void ImprovedEuler::step (double dt)
{
  Profile prof ("ImprovedEuler.step()");
  if (debug) cerr << "ImprovedEuler.step(" << dt << ") called\n";
  sc.physics.CalculateEpiphenomena();
  // set up initial ni_vel:
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vertex = v->second;
    vertex.ni_vel = vertex.vel;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k1 = vert.ni_force * (dt / vert.mass);
    vert.ni_vel = vert.vel + vert.ni_k1;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k2 = vert.ni_force * (dt / vert.mass);

    // apply final calculation:
    vert.vel += (vert.ni_k1 + vert.ni_k2) / 2;
    vert.pos += vert.vel * dt;
  }
  sc.physics.CalculateLimitEffects();
  sc.physics.CalculateEpiphenomena();
}

void MidpointEuler::step (double dt)
{
  Profile prof ("MidpointEuler.step()");
  if (debug) cerr << "MidpointEuler.step(" << dt << ") called\n";
  sc.physics.CalculateEpiphenomena();
  // set up initial ni_vel:
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vertex = v->second;
    vertex.ni_vel = vertex.vel;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k1 = vert.ni_force * (dt / vert.mass);
    vert.ni_vel = vert.vel + vert.ni_k1 / 2.0;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k2 = vert.ni_force * (dt / vert.mass);

    // apply final calculation:
    vert.vel += vert.ni_k2;
    vert.pos += vert.vel * dt;
  }
  sc.physics.CalculateLimitEffects();
  sc.physics.CalculateEpiphenomena();
}

void RungeKutta::step (double dt)
{
  Profile prof ("RungeKutta.step()");
  if (debug) cerr << "RungeKutta.step(" << dt << ") called\n";
  sc.physics.CalculateEpiphenomena();
  // set up initial ni_vel:
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vertex = v->second;
    vertex.ni_vel = vertex.vel;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k1 = vert.ni_force * (dt / vert.mass);
    vert.ni_vel = vert.vel + vert.ni_k1 / 2.0;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k2 = vert.ni_force * (dt / vert.mass);
    vert.ni_vel = vert.vel + vert.ni_k2 / 2.0;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k3 = vert.ni_force * (dt / vert.mass);
    vert.ni_vel = vert.vel + vert.ni_k3;
  }
  sc.physics.CalculateForces();
  for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
    Vertex &vert = v->second;
    if (vert.pinned) continue;
    vert.ni_k4 = vert.ni_force * (dt / vert.mass);

    // apply final calculation:
    vert.vel += (vert.ni_k1 + (vert.ni_k2 + vert.ni_k3)*2 + vert.ni_k4) / 6;
    vert.pos += vert.vel * dt;
  }
  sc.physics.CalculateLimitEffects();
  sc.physics.CalculateEpiphenomena();
}
