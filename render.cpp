// render.cpp
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

#include "render.hpp"
#include "simplicial.hpp"
#include "hull.hpp"
#include "profile.hpp"
#include "integrator.hpp"
#include "io.hpp"
#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

using std::list;
using std::cerr;
using std::ostringstream;
using std::ifstream;
using std::ofstream;
using std::string;

const string WINDOW_TITLE = "Geode by jhw; <esc> to exit";

// TODO: get proper license notice controls in executable
const string LICENSE_NOTICE = 
"Geode version 0.1.1, Copyright (C) 2003 Jesse H. Willett\n"
"Geode comes with ABSOLUTELY NO WARRANTY; for details ???\n"
"This is free software, and you are welcome to redistribute it\n"
"under certain conditions; type ??? for details.\n";

std::set<Render *> Render::theRenders;

void Render::cbRenderScene ()
{
  for (std::set<Render *>::iterator it = theRenders.begin();
       it != theRenders.end(); it++)
    (*it)->RenderScene();
}
void Render::cbKeyPressed (unsigned char key, int x, int y)
{
  for (std::set<Render *>::iterator it = theRenders.begin();
       it != theRenders.end(); it++)
    (*it)->KeyPressed(key,x,y);
}
void Render::cbSpecialKeyPressed (int key, int x, int y)
{
  for (std::set<Render *>::iterator it = theRenders.begin();
       it != theRenders.end(); it++)
    (*it)->SpecialKeyPressed(key,x,y);
}
void Render::cbResizeScene (int width, int height)
{
  for (std::set<Render *>::iterator it = theRenders.begin();
       it != theRenders.end(); it++)
    (*it)->ResizeScene(width,height);
}
void Render::cbMenu (int value)
{
  for (std::set<Render *>::iterator it = theRenders.begin();
       it != theRenders.end(); it++)
    (*it)->Menu(value);
}
void Render::cbMouse (int button, int state, int x, int y)
{
  for (std::set<Render *>::iterator it = theRenders.begin();
       it != theRenders.end(); it++)
    (*it)->Mouse(button,state,x,y);
}

Render::Render (SimplicialComplex &sci,
								bool renderVertices,
								bool renderEdges,
								bool renderFaces,
								bool renderLights) :
  sc(sci),
  renderingVertices(renderVertices),
  renderingEdges(renderEdges),
  renderingFaces(renderFaces),
  renderingLights(renderLights),
  renderingNormals(false),
  renderingForce(false),
  renderingAxes(false),
  frameRateSamples(50)
{
  windowID = -1;
  windowWidth = 500;
  windowHeight = 400;

  xRot   = 0.0f;
  yRot   = 0.0;
  zOff   =-300.0;

  dtPerStep = 0.03;
  stepsPerFrame = 1;

  systemPaused = true;
  truncateText = false;

  hull = NULL;

  frameCount = 0;
  frameRate  = 0;

  numRenderSceneCalls = 0;

  integrators.push_back(new ImprovedEuler(sc));
  integrators.push_back(new MidpointEuler(sc));
  integrators.push_back(new RungeKutta(sc));
  //integrators.push_back(new ForwardEuler(sc));
  integrator = 0;

  menuID = 0;

  theRenders.insert(this);
}

Render::~Render ()
{
  for (size_t i = 0; i < integrators.size(); i++)
    delete integrators[i];
  theRenders.erase(theRenders.find(this));
}

class Color
{
public:
  double r; // red
  double g; // green
  double b; // blue
  double a; // alpha
  Color (double ri, double gi, double bi, double ai)
  {
    r = ri;
    g = gi;
    b = bi;
    a = ai;
  }
  Color (const Color &c)
  {
    r = c.r;
    g = c.g;
    b = c.b;
    a = c.a;
  }
  void activate () const 
  {
    glColor4d(r,g,b,a);
  }
};

const Color BLACK       (0.0, 0.0, 0.0, 1.0);
const Color RED         (1.0, 0.0, 0.0, 1.0);
const Color GREEN       (0.0, 1.0, 0.0, 1.0);
const Color DARK_GREEN  (0.2, 0.7, 0.2, 1.0);
const Color YELLOW      (1.0, 1.0, 0.0, 1.0);
const Color BLUE        (0.0, 0.0, 1.0, 1.0);
const Color MAGENTA     (1.0, 0.0, 1.0, 1.0);
const Color DARK        (0.2, 0.2, 0.2, 1.0);
const Color MUDDY       (0.0, 0.2, 0.2, 1.0);
const Color CYAN        (0.0, 1.0, 1.0, 1.0);
const Color WHITE       (1.0, 1.0, 1.0, 1.0);

// f==0 => all a, f==1 => all b
Color color_combo(const Color &a, const Color &b, double f) {
  f = min(f,1.0);
  f = max(f,0.0);
  return Color(f*b.r + (1.0-f)*a.r,
							 f*b.g + (1.0-f)*a.g,
							 f*b.b + (1.0-f)*a.b,
							 f*b.a + (1.0-f)*a.a);
}

void Render::updateFPS() 
{
  if (++frameCount >= frameRateSamples) {
    double delta = (Time()-fpsLast).secs();
    frameRate = frameCount / delta;
    frameCount = 0;
    fpsLast.repoll();
  }
}

string Render::generateDisplayText ()
{
  ostringstream o;

  if (true) {
    o << "fps: " << frameRate 
      << " f: " << frameCount; 
    o << "\n";
  }

  if (true) {
    o << integrators[integrator]->name() << "\n";
  }

  if (true) {
    o << "dt/s: " << dtPerStep << "\n";
    o << " s/f: " << stepsPerFrame << "\n";
  }

  if (true) {
    o << "kin: " << sc.physics.kineticEnergy << "\n";
    o << "grv: " << sc.physics.gravitationalEnergy << "\n";
    o << "cnt: " << sc.physics.centeringEnergy << "\n";
    o << "spr: " << sc.physics.springEnergy << "\n";
    o << "srf: " << sc.physics.surfaceEnergy << "\n";
    o << "inp: " << sc.physics.internalPressureEnergy << "\n";
    o << "exp: " << sc.physics.externalPressureEnergy << "\n";
  }
  if (true) {
    o << "tot: " << sc.physics.totalEnergy << "\n";
  }

  if (false) {
    o << "P: " << sc.physics.pressure << "\n";
    o << "V: " << sc.physics.totalVolume << "\n";
    o << "N: " << sc.physics.moles << "\n";
    o << "R: " << sc.physics.idealGasConstant << "\n";
    o << "T: " << sc.physics.internalTemperature << "\n";
    o << "g: " << sc.physics.gamma() << "\n";
    o << "K: " << sc.physics.gasK << "\n";
    o << "ext P: " << sc.physics.externalPressure << "\n";
    o << "ext T: " << sc.physics.externalTemperature << "\n";
  }

  if (true) {
    o << "volume: " << sc.physics.totalVolume << "\n";
    o << "area:   " << sc.physics.totalArea << "\n";
    o << "length: " << sc.physics.totalLength << "\n";
  }

  if (truncateText) return o.str();
  
  if (true)
    o << "v " << sc.vertices.size() << ", "
      << "e " << sc.edges.size() << ", "
      << "f " << sc.faces.size() << "\n";

  if (true) {
    o << "mats: ";
    for (MaterialMapIt it = sc.materials.begin();
				 it != sc.materials.end(); it++) {
      int num = 0;
      for (EdgeMapIt e = sc.edges.begin(); e != sc.edges.end(); e++)
				if (e->second.material == &it->second) num++;
      o << num << " ";
    }
    o << "\n";
  }

  if (false) {
    o << "total mass:    " << sc.physics.totalMass << "\n";
    o << "ave mass:    " << (sc.physics.totalMass/sc.vertices.size()) << "\n";
  }

  if (false) {
    double totalSpringFactor  = 0.0;
    double totalLength  = 0.0;
    int num = sc.edges.size();
    for (EdgeMapIt e = sc.edges.begin(); e != sc.edges.end(); e++) {
      totalSpringFactor += e->second.springFactor;
      totalLength += 
				(e->second.vertex[1]->pos-e->second.vertex[0]->pos).magnitude();
    }
    o << "total spring factor:    " << totalSpringFactor << "\n";
    o << "ave spring factor:    " << (totalSpringFactor/num) << "\n";
    o << "total length:    " << totalLength << "\n";
    o << "ave length:    " << (totalLength/num) << "\n";
  }

  if (sc.physics.viscosityActive) {
    o << "1 viscosity: " << sc.physics.viscosity << "\n";
  }

  if (sc.physics.springDampingActive) {
    double totalDamping = 0.0;
    for (MaterialMapIt it = sc.materials.begin(); 
				 it != sc.materials.end(); it++)
      totalDamping += it->second.damper;
    double aveDamping = totalDamping / sc.materials.size();
    o << "2 spring damping: " << aveDamping << "\n";
  }

  if (sc.physics.centeringActive) {
    o << "3 centering" << "\n";
  }

  if (sc.physics.gravityActive) {
    o << "4 gravity" << "\n";
  }

  if (sc.physics.internalPressureActive) {
    o << "6 N: " << sc.physics.moles << "\n";
  }

  if (sc.physics.externalPressureActive) {
    o << "7 external pressure: " << sc.physics.externalPressure << "\n";
  }

  if (sc.physics.vertexLimitActive) {
    o << "8 vertex limit: " << sc.physics.vertexLimit << "\n";
  }

  return o.str();
}

void Render::RenderScene (const bool select, const int x, const int y)
{
  if (false && numRenderSceneCalls++ > 0)
    return; // something weird in 1st frame

  // Enables, disables or otherwise adjusts as 
  // appropriate for our current settings.
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_NORMALIZE);
  if (renderingLights)
    glEnable(GL_LIGHTING);  // is disabled below for rendering text
   
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND); 
  glEnable(GL_DEPTH_TEST); 

  // Need to manipulate the ModelView matrix to move our model around.
  glMatrixMode(GL_MODELVIEW);

  // Reset to 0,0,0; no rotation, no scaling.
  glLoadIdentity(); 

  // Move the object back from the screen.
  glTranslated(0.0,0.0,zOff);

  // Rotate the calculated amount.
  glRotated(xRot,1.0,0.0,0.0);
  glRotated(yRot,0.0,1.0,0.0);

  // Clear the color and depth buffers.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // render geometry:
  {
    if (!systemPaused && !select) {
      for (int i = 0; i < stepsPerFrame; i++)
				integrators[integrator]->step(dtPerStep);
    }
    else {
      sc.physics.CalculateEpiphenomena();
    }

    if (true) {
      renderSC(select,x,y);
    }

    if (false) {
      // what do the glut platonic geometric built-ins look like?
      YELLOW.activate();
      glutWireTetrahedron();
      GREEN.activate();
      glutWireCube(1.0);
      BLUE.activate();
      glutWireOctahedron();
      RED.activate();
      glutWireIcosahedron();
      MAGENTA.activate();
      glutWireDodecahedron();
    }

    if (false) {
      // what does this look like?
      WHITE.activate();
      glutSolidTorus(50.0, 100.0, 10, 20);
    }

    if (false) {
      // what do the other glut geometric built-ins look like?
      YELLOW.activate();
      glutWireCone( 50.0, 100.0, 10, 20);
      BLUE.activate();
      glutWireSphere(100.0, 10, 20);
    }

    if (false) {
      // what does a nonplanar polygon look like?
      glBegin(GL_POLYGON);
      YELLOW.activate();
      glVertex3d(  0,   0,   0);
      glVertex3d(  0, 100,   0);
      glVertex3d(100, 100,-100);
      glVertex3d(100,   0, 100);
      glVertex3d(100,   0,   0);
      glEnd();
    }

    // render axis:
    if (false) {
      Vector o (   0,   0,   0);
      Vector x ( 100,   0,   0);
      Vector y (   0, 100,   0);
      Vector z (   0,   0, 100);
      WHITE.activate();
      glBegin(GL_LINES); 
      renderLittleLine(o,x);
      renderLittleLine(o,y);
      renderLittleLine(o,z);
      glEnd();
      glBegin(GL_TRIANGLES);
      for (int i = 1; i <= 10; i++) {
				//renderLittleTetra(x*i);
				//renderLittleTetra(y*i);
				//renderLittleTetra(z*i);
      }
      glEnd();
    }
  }

  // render some text:
  const string text = generateDisplayText();
  {
    // Move back to the origin (for the text, below).
    glLoadIdentity();

    // We need to change the projection matrix for the text rendering.  
    glMatrixMode(GL_PROJECTION);

    // But we like our current view too; so we save it here.
    glPushMatrix();

    // Now we set up a new projection for the text.
    glLoadIdentity();
    gluOrtho2D(0,windowWidth,0,windowHeight);

    // Lit or textured text looks awful.
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    // We don't want depth-testing either.
    glDisable(GL_DEPTH_TEST); 

    // Now we want to render the paragraph.
   
    // To ease, simply translate up.  Note we're working in screen
    // pixels in this projection.
    glTranslated(6.0,windowHeight - 14,0.0);

    const int lineHeight = 6; // this is kind of voodoo with const int y below
    const double lineWidth = 
      (windowWidth-12) ;//* log(sc.physics.totalEnergy) / 15;
    glLineWidth(lineHeight);
    glBegin(GL_LINES);
    double kin = sc.physics.kineticEnergy / sc.physics.totalEnergy;
    double grv = sc.physics.gravitationalEnergy / sc.physics.totalEnergy;
    double cnt = sc.physics.centeringEnergy / sc.physics.totalEnergy;
    double spr = sc.physics.springEnergy / sc.physics.totalEnergy;
    double srf = sc.physics.surfaceEnergy / sc.physics.totalEnergy;
    double inp = sc.physics.internalPressureEnergy / sc.physics.totalEnergy;
    double exp = sc.physics.externalPressureEnergy / sc.physics.totalEnergy;
    {
      const int y = 24-lineHeight-windowHeight;
      Vector start (            0, y, 0);
      Vector end   (            0, y, 0);

      start.x = end.x;
      end.x += lineWidth*srf;
      MUDDY.activate();
      renderLittleLine(start,end);

      start.x = end.x;
      end.x += lineWidth*kin;
      RED.activate();
      renderLittleLine(start,end);

      start.x = end.x;
      end.x += lineWidth*grv;
      BLUE.activate();
      renderLittleLine(start,end);

      start.x = end.x;
      end.x += lineWidth*cnt;
      YELLOW.activate();
      renderLittleLine(start,end);

      start.x = end.x;
      end.x += lineWidth*spr;
      GREEN.activate();
      renderLittleLine(start,end);

      start.x = end.x;
      end.x += lineWidth*inp;
      MAGENTA.activate();
      renderLittleLine(start,end);

      start.x = end.x;
      end.x += lineWidth*exp;
      WHITE.activate();
      renderLittleLine(start,end);
    }
    glEnd();

    // draw that text!
    WHITE.activate();
    int y = 0;
    int num = 0;
    int limit = text.size();
    glRasterPos2i(6,y);
    y -= 13;
    for(int i = 0; i < limit; i++) {
      if (text[i] == '\n') {
				glRasterPos2i(6,y);
				y -= 13;
				++num;
				continue;
      }
      // NOTE: Without hardware acceleration, our string rendering, which
      // leverages on GLUT routine, is damn slow!  glutStrokeCharacter()
      // isn't really any better (and takes more work to translate to
      // screen).  With acceleration in place, it's all gold.
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13,text[i]);
      //glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,text[i]);
      //glutStrokeCharacter(GLUT_STROKE_ROMAN,text[i]);
    } 

    // Done with this special projection matrix.  Throw it away.
    glPopMatrix();
  }

  // All done drawing.  Let's show it.
  glutSwapBuffers();

  // And collect our statistics.
  updateFPS();
}

// ------
// Callback function called when a normal key is pressed.
void Render::KeyPressed (unsigned char key, int x, int y)
{
  int mods = glutGetModifiers();
  bool shift = mods & GLUT_ACTIVE_SHIFT;
  bool ctrl  = mods & GLUT_ACTIVE_CTRL;
  bool alt   = mods & GLUT_ACTIVE_ALT;

  // small warning: some key combos (alt-tab, alt-space, etc) get directed
  // to the window manager or elsewhere and don't show up in this method.
  if (false) {
    // other than above, these work as expected
    if (shift) cerr << "shift\n";
    if (ctrl) cerr << "ctrl";
    if (alt) cerr << "alt";
  }
  
  switch (key) {

  case 27: // Escape - We're outta here.
    if (shift)
      Profile::summarize(cerr);
    quit();
    break; // exit doesn't return, but anyway...

  case '1':
    sc.physics.viscosityActive = !sc.physics.viscosityActive;
    break;
  case 'q':
    sc.physics.viscosity += 0.03;
    break;
  case 'Q':
    sc.physics.viscosity = max(0.0,sc.physics.viscosity-0.03);
    break;

  case '2':
    sc.physics.springDampingActive = !sc.physics.springDampingActive;
    break;
  case 'w':
    // TODO: I don't exactly understand the need for a damper limit
    //       - but around 6, icosa | hull | cross blows up ForwardEuler(0.3)
    for (MaterialMapIt it = sc.materials.begin();
				 it != sc.materials.end(); it++)
      it->second.damper = min(it->second.damper+0.1,3.0);
    break;
  case 'W':
    for (MaterialMapIt it = sc.materials.begin();
				 it != sc.materials.end(); it++)
      it->second.damper = max(it->second.damper-0.1,0.0);
    break;
  case 's':
    sc.physics.springsActive = !sc.physics.springsActive;
    break;

  case '3':
    sc.physics.centeringActive = !sc.physics.centeringActive;
    break;

  case '4':
    sc.physics.gravityActive = !sc.physics.gravityActive;
    break;

  case '6':
    sc.physics.internalPressureActive = !sc.physics.internalPressureActive;
    break;
  case 'y':
    sc.physics.InjectGas(100,100);
    break;
  case 'Y':
    sc.physics.ReleaseGas(100);
    break;
  case 'h':
    sc.physics.InjectGas(1000,100);
    break;
  case 'H':
    sc.physics.ReleaseGas(1000);
    break;

  case '7':
    sc.physics.externalPressureActive = !sc.physics.externalPressureActive;
    break;
  case 'u':
    sc.physics.externalPressure += 0.01;
    break;
  case 'U':
    sc.physics.externalPressure -= 0.01;
    break;
  case 'j':
    sc.physics.externalPressure += 0.05;
    break;
  case 'J':
    sc.physics.externalPressure -= 0.05;
    break;

  case '8':
    //sc.physics.vertexLimitActive = !sc.physics.vertexLimitActive;
    break;
  case 'i':
		sc.physics.vertexLimit += 100.0;
    break;
  case 'I':
    sc.physics.vertexLimit = max(0.0,sc.physics.vertexLimit-100.0);
    break;
     
  case 'z':
    sc.physics.AddVel(Vector(  0.0,  0.0, 10.0));
    break;
  case 'x':
    sc.physics.AddVel(Vector( 10.0,  0.0,  0.0));
    break;
  case 'c':
    sc.physics.AddVel(Vector(  0.0, 10.0,  0.0));
    break;
  case 'Z':
    sc.physics.AddVel(Vector(  0.0,  0.0,-10.0));
    break;
  case 'X':
    sc.physics.AddVel(Vector(-10.0,  0.0,  0.0));
    break;
  case 'C':
    sc.physics.AddVel(Vector(  0.0,-10.0,  0.0));
    break;

  case 'v':
    sc.InvertThroughOrigin();
    break;

  case 'b':
    sc.physics.AddSpin(1.0);
    break;
  case 'B':
    sc.physics.AddSpin(-1.0);
    break;
     
  case 'k':
    sc.physics.Cool(0.7);
    break;
  case 'K':
    sc.physics.Cool(0.0);
    break;

  case 'n':
    sc.Equiangulate();
    break;

  case 'm':
    sc.physics.Jiggle(1.0);
    break;
  case 'M':
    sc.physics.Jiggle(10.0);
    break;

  case '`':
    sc.Sphericize();
    break;

  case ',':
    sc.Harmonize();
    break;

  case '<':
    sc.DeHarmonize();
    break;

  case '.':
    sc.SubdivideFaces();
    break;

  case '>':
    sc.DiamondSubdivide();
    break;

  case '/':
    if (hull == NULL)
      ConvexHull(sc).ConstructHull();
    break;

  case '|':
    sc.edges.clear();
    sc.faces.clear();
    if (hull != NULL) {
      delete hull;
      hull = NULL;
    }
    break;

  case '[':
    sc.physics.AddImpact();
    break;

  case '\\':
    try {
      if (hull == NULL && sc.faces.size() == 0) {
				hull = new ConvexHull(sc);
				hull->StartHull();
      }
      else if (hull != NULL && hull->isAdvancingHull()) {
				hull->AdvanceHull();
      }
      else if (hull != NULL) {
				hull->FinishHull();
				delete hull;
				hull = NULL;
      }
    }
    catch (const char *str) {
      cerr << "hull exception: " << str << "\n";
    }
    break;

  case '?':
    sc.BisectEdges();
    break;

  case '=':
    dtPerStep += 0.0001;
    break;

  case '-':
    dtPerStep -= 0.0001;
    dtPerStep = max(0.0,dtPerStep);
    break;

  case '+':
    dtPerStep += 0.01;
    break;

  case '_':
    dtPerStep -= 0.01;
    dtPerStep = max(0.0,dtPerStep);
    break;

  case '~':
    takeSnapshot();
    break;

  default:
    cerr << "KP: No action for " << key << "\n";
    break;
  }
}

void Render::quit ()
{
  glutDestroyWindow(windowID);
  exit(0);
}
void Render::takeSnapshot ()
{
  // find unused filename of form "snapshot.n"
  string str;
  for (int i = 0; true; i++) {
    ostringstream o;
    o << "snapshot." << i;
    str = o.str();
    ifstream s(str.c_str());
    if (!s) break;
  }
  ofstream s(str.c_str());
  s << sc;
  s.close();
}


// ------
// Callback Function called when a special key is pressed.
void Render::SpecialKeyPressed (int key, int x, int y)
{
  int mods = glutGetModifiers();
  bool shift = mods & GLUT_ACTIVE_SHIFT;
  bool ctrl  = mods & GLUT_ACTIVE_CTRL;
  bool alt   = mods & GLUT_ACTIVE_ALT;

  // small warning: some key combos (alt-tab, alt-space, etc) get directed
  // to the window manager or elsewhere and don't show up in this method.
  if (false) {
    // other than above, these work as expected
    if (shift) cerr << "shift\n";
    if (ctrl) cerr << "ctrl";
    if (alt) cerr << "alt";
  }
  
  switch (key) {
  case GLUT_KEY_F1:
    systemPaused = !systemPaused;
    break;
  case GLUT_KEY_F2:
    truncateText = !truncateText;
    break;
  case GLUT_KEY_F3:
    numRenderSceneCalls = 0;
    break;
  case GLUT_KEY_F4:
    sc.physics.NextGamma();
    break;
  case GLUT_KEY_F5:
    sc.DoublePos();
    break;
  case GLUT_KEY_F6:
    sc.HalvePos();
    break;
  case GLUT_KEY_F7:
    sc.physics.OpenValve();
    break;
  case GLUT_KEY_F8:
    // set to STP:
    sc.physics.externalPressure = Physics::stpPressure;
    sc.physics.externalTemperature = Physics::stpTemperature;
    for (VertexMapIt it = sc.vertices.begin(); it != sc.vertices.end(); it++)
      it->second.mass = Physics::stpPressure;
    break;

  case GLUT_KEY_F9:
    if (shift)
      renderingForce = !renderingForce;
    else
      renderingVertices = !renderingVertices;
    break;
  case GLUT_KEY_F10:
    if (shift)
      renderingAxes = !renderingAxes;
    else
      renderingEdges = !renderingEdges;
    break;
  case GLUT_KEY_F11:
    if (shift)
      renderingNormals = !renderingNormals;
    else
      renderingFaces = !renderingFaces;
    break;
  case GLUT_KEY_F12:
    renderingLights = !renderingLights;
    break;

  case GLUT_KEY_PAGE_DOWN:
    zOff -= 30.00;
    break;
  case GLUT_KEY_PAGE_UP:
    zOff += 30.00;
    break;

  case GLUT_KEY_INSERT:
    integrator += 1;
    integrator %= integrators.size();
    break;

  case GLUT_KEY_HOME:
    for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
      v->second.mass += 0.1;
    }
    break;
  case GLUT_KEY_END:
    for (VertexMapIt v = sc.vertices.begin(); v != sc.vertices.end(); v++) {
      v->second.mass -= 0.1;
    }
    break;

  case GLUT_KEY_UP:
    xRot -= 10.0;
    break;
  case GLUT_KEY_DOWN:
    xRot += 10.0;
    break;
  case GLUT_KEY_LEFT:
    yRot -= 10;
    break;
  case GLUT_KEY_RIGHT:
    yRot += 10;
    break;

  default:
    cerr << "SKP: No action for " << key << "\n";
    break;
  }
}

// ------
// Callback routine executed whenever our window is resized.  Lets us
// request the newly appropriate perspective projection matrix for 
// our needs.  Try removing the gluPerspective() call to see what happens.
void Render::ResizeScene (int width, int height)
{
  // Let's not core dump, no matter what.
  if (height == 0)
    height = 1;

  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (true) {
    gluPerspective(45.0,(GLdouble)width/(GLdouble)height,1.0,1000000.0);
  }
  else {
    // NOTE: zooming doesn't work anymore!
    glOrtho(  -0.5*width,     0.5*width,   // left,  right
							-0.5*height,    0.5*height,  // top,   bottom
							0001.0,   1000000.0);        // near,  far
  }

  glMatrixMode(GL_MODELVIEW);

  windowWidth  = width;
  windowHeight = height;
}

// ------
// The "main" function.  Inits OpenGL.  Calls our own init function,
// then passes control onto OpenGL.
void Render::run (int argc, char *argv[])
{
  cerr << LICENSE_NOTICE;

  glutInit(&argc, argv);

  // To see OpenGL drawing, take out the GLUT_DOUBLE request.
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(windowWidth, windowHeight);

  // Open a window 
  windowID = glutCreateWindow(WINDOW_TITLE.c_str());

  // Register the callback function to do the drawing. 
  glutDisplayFunc(&cbRenderScene);

  // If there's nothing to do, draw.
  glutIdleFunc(&cbRenderScene);

  // It's a good idea to know when our window's resized.
  glutReshapeFunc(&cbResizeScene);

  // And let's get some keyboard input.
  glutKeyboardFunc(&cbKeyPressed);
  glutSpecialFunc(&cbSpecialKeyPressed);

  // And let's get some mouse input.
  glutMouseFunc(&cbMouse);

  // OK, OpenGL's ready to go.

  // do everything needed before losing control to the OpenGL event loop:

  // set up lighting
  {
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position0[] = {  1.0, 1.0, 1.0, 0.0 };
    GLfloat light_position1[] = { -1.0, 1.0, 1.0, 0.0 };
    GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat lmodel_ambient[] = { 0.1, 0.1, 0.1, 1.0 };
    //GLfloat emission_light[] = { 0.1, 0.1, 0.1, 1.0 };
    glShadeModel(GL_SMOOTH);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
    //glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission_light);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
    glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, white_light);
    glLightfv(GL_LIGHT1, GL_SPECULAR, white_light);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    glColorMaterial(GL_FRONT,GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    if (renderingLights)
      glEnable(GL_LIGHTING);
  }

  // Depth to clear depth buffer to; type of test.
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClearDepth(1.0);
  glDepthFunc(GL_LESS); 

  // set up the menu
  {
    menuID = glutCreateMenu(cbMenu);
    glutSetMenu(menuID); // should be redundant with glutCreateMenu() above
    glutAddMenuEntry("save snapshot",0);
    glutAddMenuEntry("quit",1);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }

  // Load up the correct perspective matrix; using a callback directly.
  cbResizeScene(windowWidth,windowHeight);

  // Pass off control to OpenGL.
  // Above functions are called as appropriate.
  glutMainLoop();
}

void Render::Menu (int value)
{
  // TODO: better modularization of menu & key commands (esp. IDs)
  cerr << "Menu(" << value << ") called\n";
  switch (value) {
  case 0:
    takeSnapshot();
    break;
  case 1:
    quit();
    break;
  }
}

void Render::Mouse (int button, int state, int x, int y)
{
  const bool debug = false;

  if (button != GLUT_LEFT_BUTTON) return;
  if (state != GLUT_DOWN) return;
  if (debug) 
    cerr << "Mouse(" << button << "," << state<<","<<x<<","<<y<<");\n";

  const GLsizei size = 100;
  GLuint buf[size];
  glSelectBuffer(size,buf);
  glRenderMode(GL_SELECT);
  RenderScene(true,x,y);
  const GLint hits = glRenderMode(GL_RENDER);
  if (debug) 
    cerr << "num hits: " << hits << "\n";
  GLuint *ptr = buf;
  selectedObjIDs.clear();;
  for (int i = 0; i < hits; i++) {
    GLuint numNames = *ptr;
    ptr++;
    float z1 = (float)*ptr/0x7fffffff;
    ptr++;
    float z2 = (float)*ptr/0x7fffffff;
    ptr++;
    if (debug) cerr << "   num of names for hit: " << numNames << "\n";
    if (debug) cerr << "      z1 is:   " << z1 << "\n";
    if (debug) cerr << "      z2 is:   " << z2 << "\n";
    if (debug) cerr << "      name is: ";
    for (GLuint j = 0; j < numNames; j++) {
      GLuint name = *ptr;
      ptr++;
      if (debug) {
				if (sc.vertices.find(name) != sc.vertices.end())
					cerr << "vertex(" << name << ")";
				else if (sc.edges.find(name) != sc.edges.end())
					cerr << "edge(" << name << ")";
				else if (sc.faces.find(name) != sc.faces.end())
					cerr << "face(" << name << ")";
				else
					cerr << "unknown(" << name << ")";
				cerr << " ";
      }
      selectedObjIDs.insert(name);
    }
    if (debug) cerr << "\n";
  }
}

void Render::renderLittleLine (const Vector &a, const Vector &b)
{
  glVertex3dv((GLdouble *)&a);
  glVertex3dv((GLdouble *)&b);
}

void Render::renderSC (const bool select, const int x, const int y)
{
  const bool debug = false;

  if (select) {
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);

    if (debug) {
      for (int i = 0; i < 4; i++)
				cerr << "vp: " << viewport[i] << "\n";
    }

    glInitNames();
    glPushName(0);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPickMatrix((GLdouble)x,(GLdouble)(viewport[3]-y),5.0,5.0,viewport);

    if (debug) {
      cerr << "     x: " << x << "\n";
      cerr << "     y: " << y << "\n";
      cerr << "v[3]-y: " << (viewport[3]-y) << "\n";
    }
    gluPerspective(45.0,
									 (GLdouble)viewport[2]/(GLdouble)viewport[3],
									 1.0,
									 1000000.0);
    glMatrixMode(GL_MODELVIEW);
  }

  if (renderingAxes/* && sc.physics.gravityActive*/) {
    CYAN.activate();
    glLineWidth(1.0);
    glBegin(GL_LINES); 
    renderLittleLine(Vector(0,0,0),sc.physics.gravity);
    glEnd();
  }
  
  if (renderingVertices) {
    for (VertexMapIt it = sc.vertices.begin(); 
				 it != sc.vertices.end(); it++) {
      Vertex &v = it->second;

      // for lighting: normal of vertex (taken from center of mass)
      Vector n = v.pos - sc.physics.centerOfMass;
      glNormal3dv((GLdouble *)&n); // using GL_NORMALIZE

      if (selectedObjIDs.find(v.objID) != selectedObjIDs.end())
				WHITE.activate();
      else if (hull != NULL && hull->isProcessed(v))
				GREEN.activate();
      else
				YELLOW.activate();

      glLoadName(v.objID);
      glPointSize(1.0);
      glBegin(GL_POINTS); 
      glVertex3dv((GLdouble *)&v.pos);
      glEnd();

      if (renderingForce && !v.pinned) {
				YELLOW.activate();
				glLineWidth(1.0);
				glBegin(GL_LINES); 
				glVertex3dv((GLdouble *)&v.pos);
				Vector end = v.pos + v.ni_force * 10;
				glVertex3dv((GLdouble *)&end);
				glEnd();
      }
    }
  }
	
  // this must be before face rendering for last render correctness (if used)
  if (renderingEdges) {
    Vector n (0,0,0);
    for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++) {
      Edge &e = it->second;

      if (selectedObjIDs.find(e.objID) != selectedObjIDs.end())
				WHITE.activate();
      else if (!e.isCompressed())
				RED.activate();
      else
				BLUE.activate();

      glLoadName(e.objID);

      glLineWidth(1.0);
      glBegin(GL_LINES); 
      n =  e.vertex[0]->pos - sc.physics.centerOfMass;
      glNormal3dv((GLdouble *)&n); // using GL_NORMALIZE
      glVertex3dv((GLdouble *)&e.vertex[0]->pos);
      n =  e.vertex[1]->pos - sc.physics.centerOfMass;
      glNormal3dv((GLdouble *)&n); // using GL_NORMALIZE
      glVertex3dv((GLdouble *)&e.vertex[1]->pos);
      glEnd();
    }
  }

  if (hull != NULL) {
    for (EdgeMapIt it = sc.edges.begin(); it != sc.edges.end(); it++) {
      Edge &e = it->second;

      glLineWidth(1.0);
      glBegin(GL_LINES); 
      if (e.face[0] == NULL || e.face[1] == NULL) {
				RED.activate();
				renderLittleLine(e.vertex[0]->pos,e.vertex[1]->pos);
      }
      else if (hull->isVisible(*e.face[0]) != 
							 hull->isVisible(*e.face[1])) {
				YELLOW.activate();
				renderLittleLine(e.vertex[0]->pos,e.vertex[1]->pos);
				GREEN.activate();
				//renderLittleLine(e.vertex[0]->pos,hull->hullIt->second.pos);
				//renderLittleLine(e.vertex[1]->pos,hull->hullIt->second.pos);
      }
      glEnd();
    }
  }

  if (renderingFaces) {
    Vector n (0,0,0);
    for (FaceMapIt it = sc.faces.begin(); it != sc.faces.end(); it++) {
      const Face &f = it->second;
      const Vector areaNormal = f.areaNormal();
      glNormal3dv((GLdouble *)&areaNormal); // using GL_NORMALIZE

      bool gourand = false;
      if (selectedObjIDs.find(f.objID) != selectedObjIDs.end())
				WHITE.activate();
      else if (hull == NULL || hull->isVisible(f))
				gourand = true;
      else
				MUDDY.activate();

      glLoadName(f.objID);

      glBegin(GL_TRIANGLES);
      if (gourand) RED.activate();
      n = f.vertex[0]->pos-sc.physics.centerOfMass;
      glNormal3dv((GLdouble *)&n);
      glVertex3dv((GLdouble *)&f.vertex[0]->pos);
      if (gourand) GREEN.activate();
      n = f.vertex[1]->pos-sc.physics.centerOfMass;
      glNormal3dv((GLdouble *)&n);
      glVertex3dv((GLdouble *)&f.vertex[1]->pos);
      if (gourand) BLUE.activate();
      n = f.vertex[2]->pos-sc.physics.centerOfMass;
      glNormal3dv((GLdouble *)&n);
      glVertex3dv((GLdouble *)&f.vertex[2]->pos);
      glEnd(); 

      if (renderingNormals) {
				const Vector center = f.center();
				glLineWidth(1.0);
				glBegin(GL_LINES); 
				CYAN.activate();
				renderLittleLine(center,center+areaNormal/100.0);
				glEnd();
      }
    }
  }
  if (select) {
    glPopName();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glFlush();
  }
}
