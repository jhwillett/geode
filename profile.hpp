// profile.hpp
//
// various home-brew tools for timing code
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

#ifndef PROFILING_HPP
#define PROFILING_HPP

#include <ostream>
#include <list>

class Time
{
  // For now, this is a really cheezy wrapper for <sys/time.h>'s
  // gettimeofday() and struct timeval.  So cheezy as to be an
  // almost unjustified code bloat, though at present the use
  // of constructor semantics does provide syntactic cleanup
  // for mostplaces this is used.
  long sec;
  long usec;
  
public:
  Time (bool poll = true);         // if poll, current time, else 0
  Time (const Time &t);            // copy       (orthodox canonical form)
  Time &operator= (const Time &t); // assignment (orthodox canonical form)
  void repoll ();
  double secs () const;
  friend Time operator-(const Time &a, const Time &b);
  friend std::ostream &operator<< (std::ostream &s, const Time &t);
};

// self-instrumented profiling
//    - not thread-safe
//    - not accurate
//    - but it's mine!
class Profile
{
public:
  class Record;

private:
  // implementation details hidden so I can avoid universal recompiles, etc.
  Record *record;

  // these objects parallel the runtime stack:
  static Record *root;
  static Record *stackHead;
  static Record *historyHead;

  // these are meant to go on your runtime stack:
  void *operator new (size_t);
public:
  Profile (const char *name);
  ~Profile ();
  static void summarize (std::ostream &s);
};

#endif // #ifndef PROFILING_HPP
