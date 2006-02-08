// interator.hpp
// 
// collection of numeric integratration techniques
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

#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include <string>
#include "algebra.hpp"

class SimplicialComplex;

class Integrator // abstract class
{
protected:
  SimplicialComplex &sc;
public:
  Integrator (SimplicialComplex &sci) : sc(sci) {}
  virtual std::string name () = 0;   // pure virtual
  virtual void step (double dt) = 0; // pure virtual
};

class ForwardEuler : public Integrator
{
public:
  ForwardEuler (SimplicialComplex &sci) : Integrator(sci) {}
  virtual std::string name () { return "ForwardEuler"; }
  virtual void step (double dt);
};

class ImprovedEuler : public Integrator
{
public:
  ImprovedEuler (SimplicialComplex &sci) : Integrator(sci) {}
  virtual std::string name () { return "ImprovedEuler"; }
  virtual void step (double dt);
};

class MidpointEuler : public Integrator
{
public:
  MidpointEuler (SimplicialComplex &sci) : Integrator(sci) {}
  virtual std::string name () { return "MidpointEuler"; }
  virtual void step (double dt);
};

class RungeKutta : public Integrator
{
public:
  RungeKutta (SimplicialComplex &sci) : Integrator(sci) {}
  virtual std::string name () { return "RungeKutta"; }
  virtual void step (double dt);
};

#endif // #ifndef RENDER_HPP
