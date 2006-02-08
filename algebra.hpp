// algebra.hpp
//
// linear algebra classes for 3D stuff
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

#ifndef ALGEBRA_HPP
#define ALGEBRA_HPP

#include <cmath>

template<class T> const T &min (const T& a, const T&b) 
{
  return (a < b) ? a : b;
}

template<class T> const T &max (const T& a, const T&b) 
{
  return (a > b) ? a : b;
}

// stuff found in algebra.cpp:
void testVectors();

// stuff found in this file (so this file has forward reference)
class Vector;

// vector addition and subtraction:
Vector operator+ (const Vector &a, const Vector &b);
Vector operator- (const Vector &a, const Vector &b);

// vector unary negation:
Vector operator- (const Vector &a);

// vector scalar products:
Vector operator* (const Vector &a, const double d);
Vector operator* (const double d, const Vector &a);
Vector operator/ (const Vector &a, const double d);
Vector operator/ (const double d, const Vector &a);

// vector cross product:
Vector operator^(const Vector &a, const Vector &b);

// vector dot product:
double operator*(const Vector &a, const Vector &b);

// vector normal:
Vector normal (const Vector &v0, const Vector &v1, const Vector &v2);

// test that vectors are collinear:
bool collinear (const Vector &v0,const Vector &v1,const Vector &v2);

class Vector
{
  // mostly from "Physics for Game Developers", David M. Bourg
public:
  // NOTE: this object layout is designed to be compatible with
  // calls like:
  //    Vector v (1,2,3);
  //    glVertex3d((GLdouble *)&v);
  // so changing the ordering or size of these 3 members is HAZARDOUS!
  // Also, changing sizeof(Vector) might cause problems if we have
  // begun passing whole arrays of Vectors to OpenGL as geometry lists,
  // without careful attention.

  double  x;
  double  y;
  double  z;
  
  // We don't want:
  //    Vector(double xi = 0.0, double yi = 0.0, double zi = 0.0) {}
  // because it permits thinks like:
  //    Vector v = 1234;
  // which initializes v to (0,0,0) but plain shouldn't compile.
  // DONE: figure this out, clearly there's something I'm missing.
  // GUESS: each 1-argument constructor implies assignment/conversion op
  // SPOOR: "explicit" keyword, C++PL, Stroustrup, 2nd ed, 11.7.1
  // With "explicit":
  //    Vector v = 1234;
  // won't compile, but
  //    Vector v (1234);
  // yields (1234,0,0).  Bingo.
  //    explicit Vector (double xi = 0.0, double yi = 0.0, double zi = 0.0) {}
  // Still, I think I prefer the no-default-args form anyhow...

  Vector (const double xi, const double yi, const double zi)
  {
    x = xi;
    y = yi;
    z = zi;
  }

  Vector (const Vector &vi)
  {
    x = vi.x;
    y = vi.y;
    z = vi.z;
  }

  double magnitude () const
  {
    return sqrt((*this)*(*this));
  }

  Vector &operator/= (const double d)
  {
    x /= d;
    y /= d;
    z /= d;
    return *this;
  }
  
  Vector &operator*= (const double d)
  {
    x *= d;
    y *= d;
    z *= d;
    return *this;
  }
  
  Vector &operator+= (const Vector &a)
  {
    x += a.x;
    y += a.y;
    z += a.z;
    return *this;
  }
  
  Vector &operator-= (const Vector &a)
  {
    x -= a.x;
    y -= a.y;
    z -= a.z;
    return *this;
  }
};

// vector addition and subtraction:
inline Vector operator+ (const Vector &a, const Vector &b)
{
  return Vector(a.x+b.x,a.y+b.y,a.z+b.z);
}
inline Vector operator- (const Vector &a, const Vector &b)
{
  return Vector(a.x-b.x,a.y-b.y,a.z-b.z);
}

// vector unary negation:
inline Vector operator- (const Vector &a)
{
  return Vector(-a.x,-a.y,-a.z);
}

// vector scalar products:
inline Vector operator* (const Vector &a, const double d)
{
  return Vector(a.x*d,a.y*d,a.z*d);
}
inline Vector operator* (const double d, const Vector &a)
{
  return Vector(a.x*d,a.y*d,a.z*d);
}
inline Vector operator/ (const Vector &a, const double d)
{
  return Vector(a.x/d,a.y/d,a.z/d);
}
inline Vector operator/ (const double d, const Vector &a)
{
  return Vector(a.x/d,a.y/d,a.z/d);
}

// vector cross product:
inline Vector operator^(const Vector &a, const Vector &b)
{
  return Vector(a.y * b.z - a.z * b.y,
								a.z * b.x - a.x * b.z,
								a.x * b.y - a.y * b.x);
}

// vector dot product:
inline double operator*(const Vector &a, const Vector &b)
{
  return  a.x*b.x + a.y*b.y + a.z*b.z;
}

// vector normal:
inline Vector normal (const Vector &v0, const Vector &v1, const Vector &v2)
{
  // taking v0 as origin
  Vector a = v1-v0;
  Vector b = v2-v0;
  return a^b;
}

// test that vectors are collinear:
inline bool collinear (const Vector &v0,const Vector &v1,const Vector &v2)
{
  return 
    ( v2.z - v0.z ) * ( v1.y - v0.y ) -
    ( v1.z - v0.z ) * ( v2.y - v0.y ) == 0
    && ( v1.z - v0.z ) * ( v2.x - v0.x ) -
    ( v1.x - v0.x ) * ( v2.z - v0.z ) == 0
    && ( v1.x - v0.x ) * ( v2.y - v0.y ) -
    ( v1.y - v0.y ) * ( v2.x - v0.x ) == 0  ;
}

#endif // #ifndef ALGEBRA_HPP
