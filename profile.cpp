// profile.cpp
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

#include "profile.hpp"
#include <sys/time.h>
#include <iostream>
using std::cerr;
using std::ostream;
using std::ios_base;

// ugly stuff:
#ifndef NULL
#define NULL 0;
#endif

// clock() from <time.h> doesn't cut it for FPS calculations.
// It only measures processor time consumed by the current process,
// which in this case ignores the time spend by, say, the OpenGL
// server process!
// time() from <time.h> only has 1 second precision.
// We kind of want Java's System.currentTimeMillis() or Win's timeGetTime()
// The solution is to use POSIX's gettimeofday() from <std/time.h>

// static checking that gettimeofday() is supported
namespace {
  class CheckGetTimeOfDay
  {
  public:
    CheckGetTimeOfDay () {
      static timeval tv;
      int err = gettimeofday(&tv,NULL); // ignoring timezone stuff
      if (err != 0) {
				cerr << "error on gettimeofday(): " << err << "\n";
#ifdef ENOSYS
				if (err == ENOSYS) {
					cerr << "  (means gettimeofday() is unsupported)\n";
				}
#endif
				throw "error with gettimeofday()";
      }
      if (false) {
				cerr << "just checking static init exception behavior\n";
				throw "just checking static init exception behavior";
      }
    }
  };
  static CheckGetTimeOfDay check;
}

// usec is microseconds is millionths of a second
static const long usec_scale = 1000000;

Time::Time (bool poll)
{
  if (poll) {
    repoll();
  }
  else {
    sec = 0;
    usec = 0;
  }
}
Time::Time (const Time &t)
{
  sec = t.sec;
  usec = t.usec;
}
Time &Time::operator= (const Time &t)
{
  sec = t.sec;
  usec = t.usec;
  return *this;
}

void Time::repoll ()
{
  timeval tv;
  gettimeofday(&tv,0); // ignoring timezones, errors
  sec  = tv.tv_sec;
  usec = tv.tv_usec;
}

ostream &operator<< (ostream &s, const Time &t)
{
  // cheezy way of helping keep columns steady:
  const ios_base::fmtflags old_options = s.flags(ios_base::fixed);
  const double d = t.secs();
  if (d > 0.0)
    s << " " << d;
  else
    s << d;
  s.flags(old_options);
  return s;
}

Time operator-(const Time &a, const Time &b)
{
  Time t;
  t.sec = a.sec - b.sec;
  t.usec = a.usec - b.usec;
  return t;
}

double Time::secs () const
{
  return (usec_scale * sec + usec)/(1.0*usec_scale);
}


class Profile::Record
{
public:
  Time start;
  Time stop;
  const char *name;
  Record *stackNext;
  Record *historyNext;
  Record (const char *namei) :
    start(true),
    stop(false)
  {
    name = namei;
  }
  Time time () const { return stop-start; }
};

Profile::Record *Profile::root = new Profile::Record("");
Profile::Record *Profile::stackHead = Profile::root;
Profile::Record *Profile::historyHead = NULL;

Profile::Profile (const char *name)
{
  if (true) {
    record = NULL;
    return;
  }
  record = new Record(name);
  record->stackNext = stackHead;
  record->historyNext = historyHead;
  stackHead = record;
  historyHead = record;
}

Profile::~Profile ()
{
  if (record == NULL) return;
  record->stop.repoll();
  stackHead = stackHead->stackNext;
}

#include <set>
#include <map>
#include <string>
#include <iostream>
using std::ostream;
using std::string;
using std::multimap;
using std::map;
using std::pair;
using std::make_pair;

typedef multimap<const Profile::Record *, const Profile::Record *> CallMap;
typedef CallMap::const_iterator CallMapIt_const;

void summa (ostream &s,
						string indent,
						const Profile::Record *caller,
						CallMap &calls)
{
  s << indent << caller->name << ": " << caller->time() << "\n";
  indent += "   ";
  pair<CallMapIt_const, CallMapIt_const> range = calls.equal_range(caller);
  for (CallMapIt_const it = range.first; it != range.second; it++) {
    const Profile::Record *callee = it->second;
    summa(s,indent,callee,calls);
  }
}

class Name;
typedef map<string,Name *> NameMap;
class Name
{
public:
  string name;
  double totalTime;
  int    numCalls;
  NameMap subNames;
  Name (const char *namei)
  {
    name = namei;
    totalTime = 0;
    numCalls = 0;
  }
};


void Profile::summarize (ostream &s)
{
  if (true)
    return;

  root->stop.repoll();

  if (false) {
    CallMap calls;
    for (Record *r = historyHead; r != NULL; r = r->historyNext)
      calls.insert(make_pair(r->stackNext,r));
    summa(s,string(""),root,calls);
  }    
  else {
    NameMap names;
    names[""] = new Name("");
    for (Record *r = historyHead; r != NULL; r = r->historyNext) {
      if (names.find(r->name) == names.end())
				names[r->name] = new Name(r->name);
      names[r->name]->totalTime += r->time().secs();
      names[r->name]->numCalls += 1;
    }
    for (Record *r = historyHead; r != NULL; r = r->historyNext) {
      names[r->stackNext->name]->subNames[r->name] = names[r->name];
    }
    for (NameMap::iterator it = names.begin(); it != names.end(); it++) {
      Name *name = it->second;
      cerr << name->name << ": " 
					 << name->totalTime << " " 
					 << name->numCalls << "\n";
      for (NameMap::iterator subIt = name->subNames.begin(); 
					 subIt != name->subNames.end(); subIt++) {
				Name *subName = subIt->second;
				cerr << "   " 
						 << subName->name << ": " 
						 << subName->totalTime << " " 
						 << subName->numCalls << "\n";
      }
    }
    for (NameMap::iterator it = names.begin(); it != names.end(); it++)
      delete it->second;
  }

  // clean up Records:
  while (historyHead != NULL) {
    Record *r = historyHead;
    historyHead = r->historyNext;
    delete r;
  }
  root = stackHead = new Profile::Record("");
  historyHead = NULL;
}



