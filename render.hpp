// render.hpp
// 
// main OpenGL rendering and interface engine
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

#ifndef RENDER_HPP
#define RENDER_HPP

#include <string>
#include <vector>
#include <set>
#include "algebra.hpp"
#include "profile.hpp"

class SimplicialComplex;
class ConvexHull;
class Integrator;

class Render
{
public:
  Render (SimplicialComplex &sci,
					bool renderVertices,
					bool renderEdges,
					bool renderFaces,
					bool renderLights);
  void run (int argc, char *argv[]);
  ~Render ();

private:
  int windowID;     // window ID from GLUT
  int windowWidth;  // window size from GLUT, screen pixels
  int windowHeight; // window size from GLUT, screen pixels

  double xRot;      // camera x position (radians)
  double yRot;      // camera y position (radians)
  double zOff;      // camera z position (openGL units/meters)

  double dtPerStep;
  int    stepsPerFrame;

  SimplicialComplex &sc;
  ConvexHull *hull;

  bool renderingVertices;
  bool renderingEdges;
  bool renderingFaces;
  bool renderingLights;

  bool renderingNormals;
  bool renderingForce;
  bool renderingAxes;

  int integrator;
  bool systemPaused;
  bool truncateText;
  bool numRenderSceneCalls;

  // FPS statistics state
  const int frameRateSamples;
  Time fpsLast;
  int frameCount;
  double frameRate;

  std::vector<Integrator *> integrators;

  // glut callback handlers (not actual callbacks):
  void RenderScene (const bool select = false,
										const int x = 0,
										const int y = 0);
  void KeyPressed (unsigned char key, int x, int y);
  void SpecialKeyPressed (int key, int x, int y);
  void ResizeScene (int width, int height);
  void Menu (int value);
  void Mouse (int button, int state, int x, int y);

  // forwarding statics, suitable for glut callbacks:
  static std::set<Render *> theRenders;
  static void cbRenderScene ();
  static void cbKeyPressed (unsigned char key, int x, int y);
  static void cbSpecialKeyPressed (int key, int x, int y);
  static void cbResizeScene (int width, int height);
  static void cbMenu (int value);
  static void cbMouse (int button, int state, int x, int y);

  int menuID;

  // helpers for RenderScene()
  void updateFPS ();
  std::string generateDisplayText ();
  void renderSC (const bool select, const int x, const int y);
  void renderLittleLine (const Vector &a, const Vector &b);

  // nontrivial UI actions:
  void quit ();
  void takeSnapshot ();

  std::set<int> selectedObjIDs;
};

#endif // #ifndef RENDER_HPP
