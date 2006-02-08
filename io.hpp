// io.hpp
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

#ifndef IO_HPP
#define IO_HPP

#include <iostream>
class Vector;
class SimplicialComplex;
class Physics;
class Material;
class Vertex;

std::ostream &operator<< (std::ostream &s, const Vector &v); 
std::istream &operator>> (std::istream &s, Vector &v); 

std::ostream &operator<< (std::ostream &s, const SimplicialComplex &sc);
std::istream &operator>> (std::istream &s, SimplicialComplex &sc);

std::ostream &operator<< (std::ostream &s, const Physics &);
std::istream &operator>> (std::istream &s, Physics &p);

std::istream &operator>> (std::istream &s, Vertex &v);
std::ostream &operator<< (std::ostream &s, const Vertex &v);

std::istream &operator>> (std::istream &s, Material &m);
std::ostream &operator<< (std::ostream &s, const Material &m);

// chews up all whitespace and non-digits on input stream:
//    - most primitive concievable non-trivial tokenizer?
std::istream &eat (std::istream &is);

#endif // #ifndev IO_HPP
