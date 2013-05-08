// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

// Copyright 2002-2009, Daniel Fontijne, University of Amsterdam -- fontijne@science.uva.nl

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gl_includes.h"


#include "glwindow.h"
#include "util.h"
#include "osdep.h"
#include "object.h"
#include "draw.h"
#include "geosphere.h"
#include "uistate.h"
#include "state.h"

int drawVector(const GAIM_FLOAT tail[3], const GAIM_FLOAT dir[3], GAIM_FLOAT scale) {
	TubeDraw &T = gui_state->m_tubeDraw;

	GLboolean l;
	const GAIM_FLOAT rotStep = 2 * M_PI / 32;
	GAIM_FLOAT z;
	e3ga rt, n, rti;
	float rotM[16];

	if (scale == 0.0) return 0;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	if (tail)
		glTranslated(tail[0], tail[1], tail[2]);
	glScaled(scale, scale, scale);

	// draw the stick of the vector
	glGetBooleanv(GL_LIGHTING, &l);
	glDisable(GL_LIGHTING);
	T.begin(GL_LINES);
	T.vertex3d(0.0, 0.0, 0.0);
	T.vertex3d(0.96 * dir[0], 0.96 * dir[1], 0.96 * dir[2]);
	T.end();
/*	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.96 * dir[0], 0.96 * dir[1], 0.96 * dir[2]);
	glEnd();*/
	if (l) glEnable(GL_LIGHTING);

	glPushMatrix();

	// translate to head of vector
	glTranslated(dir[0], dir[1], dir[2]);

	if (gui_state->m_vectorHeadSize >= 0.0) { // temp test for fixed vector head size
		glScaled(1.0 / scale, 1.0 / scale, 1.0 / scale);
		glScaled(gui_state->m_vectorHeadSize, gui_state->m_vectorHeadSize, gui_state->m_vectorHeadSize);
	}
	else {
		glScaled(-gui_state->m_vectorHeadSize, -gui_state->m_vectorHeadSize, -gui_state->m_vectorHeadSize);
		if (scale > 1.0) {
			glScaled(1.0 / scale, 1.0 / scale, 1.0 / scale);
			double s = pow(scale, 0.75);
			glScaled(s, s, s);
		}
	}


	// rotate e3 to vector direction
	e3gaRve3(rt, e3ga(GRADE1, dir[0], dir[1], dir[2]));
	e3gaRotorToOpenGLMatrix(rt, rotM);
	glMultMatrixf(rotM);

//	glDisable(GL_CULL_FACE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
/*
	// ---If normals are not required:	
	glBegin(GL_TRIANGLE_FAN);
	glVertex3d(0.0, 0.0, 0.0);			
	for (z = 0; z <= 2 * M_PI + 0.01; z += 2 * M_PI / 16) {
		// no normals?
		glVertex3d(0.1 * cos(z), 0.1 * sin(z), -0.25);
	}
	glEnd();*/

	// ---Otherwise:
	rt = (-0.5 * rotStep * e3ga::e1 ^ e3ga::e2).exp();
	rti = rt.inverse(); n = e3ga::e1;
	glBegin(GL_QUADS);
	for (z = 0; z < 2 * M_PI; z += rotStep) {
		glNormal3dv(n[GRADE1]);
		glVertex3d(0.1 * cos(z), 0.1 * sin(z), -0.25);
		glVertex3d(0.0, 0.0, 0.00);
		n = (rt * n * rti)(GRADE1);
		glNormal3dv(n[GRADE1]);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3d(0.1 * cos(z + rotStep), 0.1 * sin(z + rotStep), -0.25);
	}
	glEnd();

	n = e3ga::e1;
	glBegin(GL_TRIANGLE_FAN);
	glNormal3dv((-(e3ga&)e3ga::e3)[GRADE1]);
	glVertex3d(0.0f, 0.0f, -0.25);
	for (z = 0; z <= 2 * M_PI + 1e-3; z += rotStep) {
		glVertex3d(0.1 * cos(z), 0.1 * sin(z), -0.25);
	}
	glEnd();

	glDisable(GL_CULL_FACE);

	glPopMatrix();


	// draw the stick of the vector
/*	glGetBooleanv(GL_LIGHTING, &l);
	glDisable(GL_LIGHTING);
//	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.98 * dir[0], 0.98 * dir[1], 0.98 * dir[2]);
	glEnd();
	if (l) glEnable(GL_LIGHTING);*/

	glPopMatrix();

	return 0;
}


int drawBivector(const GAIM_FLOAT base[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], 
				 const GAIM_FLOAT ortho2[3], GAIM_FLOAT scale, int method /*= DRAW_BV_CIRCLE*/, 
				 int flags /* = 0 */, object *o /*= NULL*/) {
	TubeDraw &T = gui_state->m_tubeDraw;
	const GAIM_FLOAT rotStep = 2 * M_PI / 64;
	GAIM_FLOAT x, y;
	e3ga rt;
	float rotM[16];

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	if (base) glTranslated(base[0], base[1], base[2]);

//	printf("method = %d %d\n", method, DRAW_BV_PARALLELOGRAM);
//	printf("scale = %f,\n", scale);
	// scale bivector

	if ((method != DRAW_BV_PARALLELOGRAM) &&
		(method != DRAW_BV_PARALLELOGRAM_NO_VECTORS)) {
		glScaled(scale, scale, scale);

		// rotate e3 to normal direction
		e3gaRve3(rt, e3ga(GRADE1, normal[0], normal[1], normal[2])); 
		e3gaRotorToOpenGLMatrix(rt, rotM);
		glMultMatrixf(rotM);
	}
	else {
		e3ga o1(GRADE1, ortho1);
		e3ga o2(GRADE1, ortho2);
		// scale is based on _circle_, re-scale for square:
		double size = scale * scale * M_PI / sqrt((o1 ^ o2).norm_a());
		double _scale = sqrt(size);
//		printf("Would scale by %f (scale = %f)\n", _scale, scale);
//		printf("o1: %f %f %f\n", ortho1[0], ortho1[1], ortho1[2]);
//		printf("o2: %f %f %f\n", ortho2[0], ortho2[1], ortho2[2]);
		glScaled(_scale, _scale, _scale);
	}

	switch(method) {
	case DRAW_BV_CIRCLE:
		if (o->fgColor(3) > 0.0)  {
			// draw the filled-in circle (back)
			glNormal3d(0.0, 0.0, 1.0);
			glBegin(GL_TRIANGLE_FAN);
			glVertex3d(0.0, 0.0, 0.0);
			for (x = 0; x < M_PI * 2.0 + 0.001; x += rotStep)
				glVertex2d(cos(x), sin(x));
			glEnd();

			// draw the filled-in circle (front)
			glNormal3d(0.0, 0.0, -1.0);
			if (flags & 0x01) glPolygonMode(GL_FRONT,  GL_LINE);
			glBegin(GL_TRIANGLE_FAN);
			glVertex3d(0.0, 0.0, 0.0);
			for (x = 0; x < M_PI * 2.0 + 0.001; x += rotStep)
				glVertex2d(-cos(x), sin(x));
			glEnd();
			if (flags & 0x01) glPolygonMode(GL_FRONT,  GL_FILL);
		}

		// draw the outline
		glDisable(GL_LIGHTING);
		if (o) o->setOlColor();
		T.begin(GL_LINE_LOOP);
		for (x = 0; x < M_PI * 2.0; x += rotStep)
			T.vertex2d(cos(x), sin(x));
		T.end();
		/*	glBegin(GL_LINE_LOOP);
		for (x = 0; x < M_PI * 2.0; x += rotStep)
			glVertex2d(cos(x), sin(x));
		glEnd();*/

		if (flags & 0x01) {
			// draw 6 little 'hooks' along the edge of the circle
			for (x = 0; x < 6; x++) {
				T.begin(GL_LINES);
				T.vertex2d(1.0, 0.0);
				T.vertex2d(1.0, -0.5);
				T.end();
				/*glBegin(GL_LINES);
				glVertex2d(1.0, 0.0);
				glVertex2d(1.0, -0.3);
				glEnd();*/
				glRotated(360 / 6, 0.0, 0.0, 1.0);
			}

			// draw a normal vector
/*			glBegin(GL_LINES);
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3d(0.0, 0.0, 1.0);
			glEnd();*/
		}

		break;
	case DRAW_BV_PARALLELOGRAM:
	case DRAW_BV_PARALLELOGRAM_NO_VECTORS:
		// front
		glNormal3dv(normal);
		glBegin(GL_QUADS);
		glVertex3d(0.0, 0.0, 0.0);
		glVertex3dv(ortho1);
		glVertex3d(ortho1[0] + ortho2[0], ortho1[1] + ortho2[1], ortho1[2] + ortho2[2]);
		glVertex3dv(ortho2);
		glEnd();

		// back
		/*if ((method != DRAW_BV_PARALLELOGRAM) && 
			(flags & 0x01)) glPolygonMode(GL_FRONT,  GL_LINE);*/
		if (flags & 0x01) glPolygonMode(GL_FRONT,  GL_LINE);
		glNormal3d(-normal[0], -normal[1], -normal[2]);
		glBegin(GL_QUADS);
		glVertex3dv(ortho2);
		glVertex3d(ortho1[0] + ortho2[0], ortho1[1] + ortho2[1], ortho1[2] + ortho2[2]);
		glVertex3dv(ortho1);
		glVertex3d(0.0, 0.0, 0.0);
		glEnd();
		if (flags & 0x01) glPolygonMode(GL_FRONT,  GL_FILL);

		// vectors
		if (method != DRAW_BV_PARALLELOGRAM_NO_VECTORS) {
			glPolygonMode(GL_FRONT_AND_BACK, (o->m_drawMode & OD_WIREFRAME) ? GL_LINE : GL_FILL);
			if (o) o->setBgColor();
			drawVector(NULL, ortho1, 1.0);
			if (o->m_drawMode & OD_SHADE) glEnable(GL_LIGHTING); // has been turned off by drawVector()
			drawVector(ortho1, ortho2, 1.0);

			// outline 
			/*
			glDisable(GL_LIGHTING);
			if (o) o->setOlColor();
			//glBegin(GL_LINE_LOOP);
			glVertex3d(0.0, 0.0, 0.0);
			glVertex3dv(ortho1);
			glVertex3d(ortho1[0] + ortho2[0], ortho1[1] + ortho2[1], ortho1[2] + ortho2[2]);
			glVertex3dv(ortho2);
			glEnd();*/

			// draw a normal vector
			if (flags & 0x01) {
				T.begin(GL_LINES);
				T.vertex3d(0.0, 0.0, 0.0);
				T.vertex3dv(normal);
				T.end();
/*				glBegin(GL_LINES);
				glVertex3d(0.0, 0.0, 0.0);
				glVertex3dv(normal);
				glEnd();*/
			}
		}

		break;
	case DRAW_BV_CROSS:
	case DRAW_BV_CURLYTAIL:
		glDisable(GL_LIGHTING);
		T.begin(GL_LINE_STRIP);
		for (y = 0; y <= M_PI *2 + 0.001; y += M_PI * 2 / 64)
			T.vertex2d(-sqrt(y / (2 * M_PI)) * sin(y), sqrt(y / (2 * M_PI)) * cos(y));
		T.end();
		break;
	case DRAW_BV_SWIRL:
		glDisable(GL_LIGHTING);
		for (x = 0; x < 4; x++) {
			T.begin(GL_LINE_STRIP);
			for (y = 0; y <= M_PI / 2 + 0.001; y += M_PI / (2 * 16)) {
				T.vertex2d((1.0 - sin(y)), cos(y));
			}
			T.end();
			glRotated(90, 0.0, 0.0, 1.0);
		}
		break;
	default:
		glPopMatrix();
		return -1;
	}

	glPopMatrix();

	return 0;
}

int drawTriVector(const GAIM_FLOAT base[3], GAIM_FLOAT scale, GAIM_FLOAT vector[][3], int method /*= DRAW_TV_SPHERE*/, int flags /*= 0*/, object *o /*= NULL*/) {
	TubeDraw &T = gui_state->m_tubeDraw;

	GAIM_FLOAT scaleSign = (scale < 0.0) ? -1.0 : 1.0;
	if ((method == DRAW_TV_PARALLELEPIPED) ||
		(method == DRAW_TV_PARALLELEPIPED_NO_VECTORS)) {
		if (vector == NULL) drawTriVector(base, scale, NULL, DRAW_TV_SPHERE, flags, o);

		// adjust scale for cube
		scale = scaleSign * pow(scaleSign * scale, 1.0 / 3.0);
	}
	else {
		// adjust scale for sphere
		scale = scaleSign * pow(scaleSign * scale / ((4.0/3.0) * M_PI), 1.0 / 3.0);
	}

//	printf("Draw scale: %f\n", scale);

	int i, j, vtxIdx;
	GAIM_FLOAT z, s = (scale < 0.0) ? -1.0 : 1.0, f;
	const GAIM_FLOAT zMax = 4 * M_PI;
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	if (base) glTranslated(base[0], base[1], base[2]);
	glScaled(fabs(scale), fabs(scale), fabs(scale));

	switch(method) {
	case DRAW_TV_SPHERE:
		gsDraw(g_sphereSphere, ((flags & 0x01) ? s * 0.1 : 0.0));
		break;
	case DRAW_TV_CROSS:
	case DRAW_TV_CURLYTAIL:
		// mag & ori nog niet helemaal goed....
		if (flags & 0x01) glScaled(s, fabs(s), fabs(s));
		glDisable(GL_LIGHTING);
		T.begin(GL_LINE_STRIP);
		for (z = -zMax; z <= zMax; z += 0.1) {
			f = z / zMax;
			f = sqrt(1 - f * f);
			T.vertex3d(/*s * */1.0 * f * cos(z), 1.0 * f * sin(z), 1.0 * z / zMax);
		}
		T.end();
		break;

	case DRAW_TV_PARALLELEPIPED:
	case DRAW_TV_PARALLELEPIPED_NO_VECTORS:
		if (flags & 0x01) glScaled(s, fabs(s), fabs(s));

		// maybe define vertices of the 'cube'
		{
			GAIM_FLOAT vertex[8][3];
			int vertexVectors[8][3] = 
			{{-1, -1, -1}, // -
			{0, -1, -1},  // 0
			{0, 1, -1}, // 0 + 1
			{1, -1, -1}, // 1
			{2, -1, -1}, // 2
			{0, 2, -1},  // 0 + 2
			{0, 1, 2}, // 0 + 1 + 2
			{1, 2, -1}}; // 1 + 2
			
			int polygon[6][4] = 
			{{0, 1, 5, 4},
			{0, 4, 7, 3},
			{4, 5, 6, 7},
			{1, 2, 6, 5},
			{6, 2, 3, 7},
			{0, 3, 2, 1}};

			for (i = 0; i < 8; i++) {
				vertex[i][0] = vertex[i][1] = vertex[i][2] = 0.0;
				for (j = 0; j < 3; j++) {
					if ( (vtxIdx = vertexVectors[i][j]) < 0) continue;
					vertex[i][0] += vector[vtxIdx][0];
					vertex[i][1] += vector[vtxIdx][1];
					vertex[i][2] += vector[vtxIdx][2];
				}
			}

			glBegin(GL_QUADS);
			for (i = 0; i < 6; i++) {
				e3ga v1(GRADE1, vertex[polygon[i][0]]);
				e3ga v2(GRADE1, vertex[polygon[i][1]]);
				e3ga v3(GRADE1, vertex[polygon[i][3]]);
				e3ga normal(scaleSign * ((v2 - v1) ^ (v3 - v1)).dual());
				glNormal3dv(normal[GRADE1]);
				if (scale >= 0.0)
					for (j = 0; j < 4; j++) glVertex3dv(vertex[polygon[i][j]]);
				else for (j = 3; j >= 0; j--) glVertex3dv(vertex[polygon[i][j]]);
				
			}
			glEnd();

			if (method == DRAW_TV_PARALLELEPIPED) {
				// draw the vectors
				if (o) o->setBgColor();
				drawVector(vertex[0], vector[0], 1.0);
				if (o->m_drawMode & OD_SHADE) glEnable(GL_LIGHTING); // has been turned off by drawVector()
				drawVector(vertex[1], vector[1], 1.0);
				if (o->m_drawMode & OD_SHADE) glEnable(GL_LIGHTING); // has been turned off by drawVector()
				drawVector(vertex[2], vector[2], 1.0);
				if (o->m_drawMode & OD_SHADE) glEnable(GL_LIGHTING); // has been turned off by drawVector()
			}
		}

		break;
	default:
		glPopMatrix();
		return -1;
	}

	glPopMatrix();

	return 0;
}

void glApplyRotor(e3ga rotor)
{
  float matrix[16];
  e3gaRotorToOpenGLMatrix(rotor, matrix);
  glMultMatrixf(matrix);
}

int drawRegulus(e3ga axis, double slant) {
  TubeDraw &T = gui_state->m_tubeDraw;
  e3ga plane = axis.dual();
  e3ga rotor;

  // OpenGL boilerplate
  glMatrixMode(GL_MODELVIEW);
  glDisable(GL_LIGHTING);
  glPushMatrix();

  // green used for 3-blades
  glColor3d(0, 1, 0);

  // rotate e3 to the axis/plane normal
  e3gaRve3(rotor, axis);
  glApplyRotor(rotor);

  // draw a selection of the lines in the regulus
  for (double i = 0; i < 2 * M_PI; i += M_PI / 8) {
    // rotate by i radians around the axis
    rotor = cos(i / 2) - plane * sin(i / 2);
    glApplyRotor(rotor);

    // rotate by the slant angle in the e2 ^ e3 plane
    glPushMatrix();
    rotor = cos(slant / 2) - (e3ga::e2 ^ e3ga::e3) * sin(slant / 2);
    glApplyRotor(rotor);

    T.begin(GL_LINES);
    T.vertex3d(1, 0, -30);
    T.vertex3d(1, 0, 30);
    T.end();
    glPopMatrix();
  }

  glPopMatrix();
  return 0;
}

int drawLine(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], GAIM_FLOAT magnitude, int method /*= DRAW_LINE_HOOKS */, int flags /*= 0*/, object *o /*= NULL*/) {
	TubeDraw &T = gui_state->m_tubeDraw;
  GAIM_FLOAT z, c, stepSize = 0.1, scaleConst = g_state->m_clipDistance * sqrt(2.0);
  e3ga e3gaR;
  float rotM[16];


	glMatrixMode(GL_MODELVIEW);
  glDisable(GL_LIGHTING);
  glPushMatrix();
  if (point) glTranslated(point[0], point[1], point[2]);

  // rotate e3 to line direction
  e3gaRve3(e3gaR, e3ga(GRADE1, normal));
  e3gaRotorToOpenGLMatrix(e3gaR, rotM);
  glMultMatrixf(rotM);

  // draw line
  T.begin(GL_LINE_STRIP);
  for (z = -scaleConst; z <= scaleConst; z += stepSize * scaleConst)
    T.vertex3d(0.0, 0.0, z);
  T.end();
  /*glBegin(GL_LINE_STRIP);
    for (z = -scaleConst; z <= scaleConst; z += stepSize * scaleConst)
    glVertex3d(0.0, 0.0, z);
    glEnd();*/

  // draw 'orientation'
  if (flags & 0x01) { 
    // todo: flip orientation of this stuff
    switch (method) {
      case DRAW_LINE_CURVE:
        // option 1: some kind of curve along the line
        if (o != NULL && o->m_drawMode & OD_MAGNITUDE)
          glScaled(0.5 * fabs(magnitude), 0.5 * fabs(magnitude), 0.5 * fabs(magnitude));
        else glScaled(0.5, 0.5, 0.5);

        T.begin(GL_LINE_STRIP);
        for (z = -scaleConst; z < scaleConst; z += 1.0 / 32)
          T.vertex3d(sin(z) * sin(z * M_PI * 2), sin(z) * cos(z * M_PI * 2), z);
        T.end();
        /*glBegin(GL_LINE_STRIP);
          for (z = -scaleConst; z < scaleConst; z += 1.0 / 32)
          glVertex3d(sin(z) * sin(z * M_PI * 2), sin(z) * cos(z * M_PI * 2), z);
          glEnd();*/
        break;
      case DRAW_LINE_CURLYTAIL:
        // option 2: curly tails
        glTranslated(0.0, 0.0, -scaleConst);
        for (c = 0.0; c <= 1.0; c += stepSize) {
          glPushMatrix();
          // if weight: scale
          if (o != NULL && o->m_drawMode & OD_MAGNITUDE)
            glScaled(0.5 * fabs(magnitude), 0.5 * fabs(magnitude), 0.5 * fabs(magnitude));
          else glScaled(0.5, 0.5, 0.5);

          T.begin(GL_LINE_STRIP);
          for (z = 0.0; z < 1.0; z += 1.0 / 64)
            T.vertex3d(sqrt(z) * sin(z * M_PI * 2), sqrt(z) * cos(z * M_PI * 2), -2.0 * (-0.5 + z) * stepSize * scaleConst);
          T.end();
          /*glBegin(GL_LINE_STRIP);
            for (z = 0.0; z < 1.0; z += 1.0 / 64)
            glVertex3d(sqrt(z) * sin(z * M_PI * 2), sqrt(z) * cos(z * M_PI * 2), -2.0 * (-0.5 + z) * stepSize * scaleConst);
            glEnd();*/

          glPopMatrix();
          glTranslated(0.0, 0.0, 2.0 * stepSize * scaleConst);
        }
        break;
      case DRAW_LINE_HOOKS:
      case DRAW_LINE_HOOKCROSSES:
        // option 3: little hooks 
        glTranslated(0.0, 0.0, -scaleConst);
        for (c = 0.0; c < 1.0; c += stepSize) {
          glPushMatrix();
          // if weight: scale
          if (o != NULL && o->m_drawMode & OD_MAGNITUDE)
            glScaled(0.5 * fabs(magnitude), 0.5 * fabs(magnitude), 0.5 * fabs(magnitude));
          else glScaled(0.5, 0.5, 0.5);

          T.begin(GL_LINE_STRIP);
          T.vertex3d(-0.25, 0.0, -1.0);
          T.vertex3d(0.0, 0.0, 0.0);
          T.vertex3d(0.25, 0.0, -1.0);
          T.end();
          /*glBegin(GL_LINE_STRIP);
            glVertex3d(-0.25, 0.0, -1.0);
            glVertex3d(0.0, 0.0, 0.0);
            glVertex3d(0.25, 0.0, -1.0);
            glEnd();*/

          glPopMatrix();
          glRotated(90, 0.0, 0.0, 1.0);
          glTranslated(0.0, 0.0, 2.0 * stepSize * scaleConst);
        }
        break;
    }
  }

  glPopMatrix();

  return 0;
}


int drawPoint(const GAIM_FLOAT point[3], GAIM_FLOAT weight, int flags, object *o) {
  glDisable(GL_LIGHTING);
  glPolygonMode(GL_FRONT_AND_BACK, (o && o->m_drawMode & OD_WIREFRAME) ? GL_LINE : GL_FILL);
  glPushMatrix();
  glTranslated(point[0], point[1], point[2]);
  glScaled(gui_state->m_pointSize, gui_state->m_pointSize, gui_state->m_pointSize);
  if (o && o->m_drawMode & OD_MAGNITUDE)
    glScaled(fabs(weight), fabs(weight), fabs(weight));

  gsDraw(g_pointSphere, (o && o->m_drawMode & OD_ORI) ? 0.01 * weight : 0.0);
  glPopMatrix();
  return 0;
}


int drawPlane(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], GAIM_FLOAT weight, int method, int flags, object *o) {
  TubeDraw &T = gui_state->m_tubeDraw;
	GAIM_FLOAT x, y;
	int s;

	GAIM_FLOAT stepSize = 0.1;
	GAIM_FLOAT scaleConst = g_state->m_clipDistance * sqrt(2.0);
	GAIM_FLOAT scaleMag = 1.0;
  for (s = 0; s < 2; s++) { // draw both front and back side individually
    if (s == 0) { // front
      glPolygonMode(GL_FRONT, (o && o->m_drawMode & OD_WIREFRAME) ? GL_LINE : GL_FILL);
      glNormal3dv(normal); 
    }
    else { // back
      glPolygonMode(GL_FRONT, ((o && o->m_drawMode & OD_WIREFRAME) || (o && o->m_drawMode & OD_ORI)) ? GL_LINE : GL_FILL);
      glNormal3d(-normal[0], -normal[1], -normal[2]); 
    }
    for (y = -scaleConst; y < scaleConst - stepSize * scaleConst; y += stepSize * scaleConst) {
      glBegin(GL_QUAD_STRIP);
      for (x = -scaleConst; x < scaleConst; x += stepSize * scaleConst) {
        glVertex3d(
            point[0] + x * ortho1[0] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[0],
            point[1] + x * ortho1[1] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[1],
            point[2] + x * ortho1[2] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[2]);
        glVertex3d(
            point[0] + x * ortho1[0] + ((s == 1) ? (y + stepSize * scaleConst) : y) * ortho2[0],
            point[1] + x * ortho1[1] + ((s == 1) ? (y + stepSize * scaleConst) : y) * ortho2[1],
            point[2] + x * ortho1[2] + ((s == 1) ? (y + stepSize * scaleConst) : y) * ortho2[2]);
      }
      glEnd();
    }
  }

  if (o && o->m_drawMode & OD_ORI) { // draw normals
    if (o && o->m_drawMode & OD_MAGNITUDE) scaleMag *= weight;
    glDisable(GL_LIGHTING);
    //T.begin(GL_LINE_STRIP);
    T.begin(GL_LINES);
    for (y = -scaleConst; y <= scaleConst; y += stepSize * scaleConst) {
      for (x = -scaleConst; x <= scaleConst; x += stepSize * scaleConst) {
        T.vertex3d(
            point[0] + x * ortho1[0] + y * ortho2[0],
            point[1] + x * ortho1[1] + y * ortho2[1],
            point[2] + x * ortho1[2] + y * ortho2[2]);
        T.vertex3d(
            point[0] + x * ortho1[0] + y * ortho2[0] + scaleMag *  normal[0],
            point[1] + x * ortho1[1] + y * ortho2[1] + scaleMag *  normal[1],
            point[2] + x * ortho1[2] + y * ortho2[2] + scaleMag *  normal[2]);
      }
    }
    T.end();
  }
  return 0;
}


int drawCircle(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], GAIM_FLOAT radius, GAIM_FLOAT weight, int method /*= DRAW_CIRCLE_HOOKS*/, int flags /*= 0*/, object *o /*= NULL*/) {
	TubeDraw &T = gui_state->m_tubeDraw;
  GAIM_FLOAT x;
  e3ga e3gaR;
  float rotM[16];

  glDisable(GL_LIGHTING);
  glPushMatrix();
  // translate to center, scalar to radius 
  if (point) glTranslated(point[0], point[1], point[2]);
  double scale = radius;
  //glScaled(m_int.m_scalar[0], m_int.m_scalar[0], m_int.m_scalar[0]);

  // rotate e3 to plane normal
  e3gaRve3(e3gaR, e3ga(GRADE1, normal[0], normal[1], normal[2]));
  e3gaRotorToOpenGLMatrix(e3gaR, rotM);
  glMultMatrixf(rotM);

  // draw circle
  {
    T.begin(GL_LINE_LOOP);
    for (x = 0; x < M_PI * 2; x += (M_PI * 2) / 64)
      T.vertex2d(scale * sin(x), scale * cos(x));
    T.end();
  }

  /*		glBegin(GL_LINE_LOOP);
        for (x = 0; x < M_PI * 2; x += (M_PI * 2) / 64)
        glVertex2d(sin(x), cos(x));
        glEnd();

*/

  // draw 6 little 'hooks' along the edge of the circle
  if (o != NULL && o->m_drawMode & OD_ORI) {
    for (x = 0; x < 6; x++) {
      T.begin(GL_LINES);
      T.vertex2d(scale * 1.0, 0.0);
      T.vertex2d(scale * 1.0, scale * -((o->m_drawMode & OD_MAGNITUDE) ? fabs(weight) : 1.0) * 0.3);
      T.end();
      /*glBegin(GL_LINES);
        glVertex2d(1.0, 0.0);
        glVertex2d(1.0, -((m_drawMode & OD_MAGNITUDE) ? fabs(weight) : 1.0) * 0.3);
        glEnd();*/
      glRotated(360 / 6, 0.0, 0.0, 1.0);
    }

    // draw a normal vector (removed)
    /*	glBegin(GL_LINES);
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(0.0, 0.0, (m_drawMode & OD_MAGNITUDE) ? fabs(weight) : 1.0);
        glEnd();*/
  }

  glPopMatrix();

  return 0;
}


int drawIdealLine(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], GAIM_FLOAT weight, int method /*= DRAW_IDEAL_LINE_HOOKS*/, int flags /*= 0*/, object *o /*= NULL*/) {
  switch (method) {
    case DRAW_IDEAL_LINE_HOOKS:
      return drawCircle(point, normal, 1, weight, method, flags, o);
    case DRAW_IDEAL_LINE_RADIUS:
      return drawCircle(point, normal, (o && o->m_drawMode & OD_MAGNITUDE) ? weight : 1, 1, method, flags, o);
    default:
      return drawLine(point, normal, weight, method, (flags | OD_ORI) ? 0x01 : 0, o);
  }
  return drawLine(point, normal, weight, method, (flags | OD_ORI) ? 0x01 : 0, o);
}


int drawLinePencil(const GAIM_FLOAT center[3], GAIM_FLOAT weight, const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], int method /*= DRAW_LINE_PENCIL*/, int flags /*= 0*/, object *o /*= NULL*/) {
	TubeDraw &T = gui_state->m_tubeDraw;
  const GAIM_FLOAT rotStep = 2 * M_PI / 64;
  GAIM_FLOAT x, scaleConst = g_state->m_clipDistance * sqrt(2.0);
  e3ga rt;
  float rotM[16];
  GAIM_FLOAT ori_pt[3], ori_dir[3];
  int orientation = 1;

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glDisable(GL_LIGHTING);
  glTranslated(center[0], center[1], center[2]);

  if (o && o->m_drawMode & OD_MAGNITUDE) {
    //weight = sqrt(fabs(weight)) / M_PI;
    glScaled(weight, weight, weight);
  }
  else {
    //glScaled(scaleConst, scaleConst, scaleConst);
  }

  // rotate e3 to pencil normal
  e3gaRve3(rt, e3ga(GRADE1, normal));
  e3gaRotorToOpenGLMatrix(rt, rotM);
  glMultMatrixf(rotM);

  if (o && o->fgColor(3) > 0.0) {
    if (o && o->m_drawMode & OD_ORI) {
      glEnable(GL_LIGHTING);
      orientation *= (0 < ortho1[0]) - (ortho1[0] <= 0);
      orientation *= (0 < ortho1[1]) - (ortho1[1] <= 0);
      orientation *= (0 < ortho1[2]) - (ortho1[2] <= 0);
      orientation *= (0 < ortho2[0]) - (ortho2[0] <= 0);
      orientation *= (0 < ortho2[1]) - (ortho2[1] <= 0);
      orientation *= (0 < ortho2[2]) - (ortho2[2] <= 0);
      printf("center: %2.2f, %2.2f, %2.2f\n", center[0], center[1], center[2]);
      printf("normal: %2.2f, %2.2f, %2.2f\n", normal[0], normal[1], normal[2]);
      printf("ortho1: %2.2f, %2.2f, %2.2f\n", ortho1[0], ortho1[1], ortho1[2]);
      printf("ortho2: %2.2f, %2.2f, %2.2f\n", ortho2[0], ortho2[1], ortho2[2]);
      printf("ori: %d\n",orientation);
      for (x = 0; x < M_PI * 2 + 0.001; x += 4 *rotStep) {
        ori_dir[0] = -cos(x);
        ori_dir[1] = sin(x);
        ori_pt[0] = (orientation > 0) * cos(x);
        ori_pt[1] = (orientation > 0) * -sin(x);
        drawVector(ori_pt, ori_dir, 1);
      }
    }
    else {
      T.begin(GL_LINE_STRIP);
      for (x = 0; x < M_PI * 2.0 + 0.001; x += rotStep) {
        T.vertex3d(-cos(x), sin(x), 0);
        T.vertex3d(0.0, 0.0, 0.0);
      }
      T.end();
    }
  }

  glPopMatrix();


  return 0;
}


int drawCirclePencil(const GAIM_FLOAT center[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], GAIM_FLOAT weight, int method /*= DRAW_CIRCLE_HOOKS*/, int flags /*= 0*/, object *o /*= NULL*/) {
  const GAIM_FLOAT rotStep = 360 / 32.0;
  GAIM_FLOAT i;
  GAIM_FLOAT negnormal[3];
  negnormal[0] = -normal[0];
  negnormal[1] = -normal[1];
  negnormal[2] = -normal[2];
  int drawMode;

  if (o && o->m_drawMode & OD_ORI) {
    drawVector(negnormal, normal, 2 * weight);
  }
  for (i = 0; i < 180; i += rotStep) { // only need to draw half the circles
    glPushMatrix();
    glRotated(i, normal[0], normal[1], normal[2]);
    drawIdealLine(center, ortho1, weight, method, flags, NULL);
    glPopMatrix();
  }

  return 0;
}


int drawRuledPlane(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], GAIM_FLOAT weight, int method /*= DRAW_RULED_PLANE*/, int flags /*= 0*/, object *o /*= NULL*/) {
  TubeDraw &T = gui_state->m_tubeDraw;
	GAIM_FLOAT x, y;
	int s;

	GAIM_FLOAT stepSize = 0.02;
	GAIM_FLOAT scaleConst = g_state->m_clipDistance * sqrt(2.0);
	GAIM_FLOAT scaleMag = 1.0;
  T.begin(GL_LINES);
  for (y = -scaleConst; y < scaleConst - stepSize * scaleConst; y += stepSize * scaleConst) {
    T.vertex3d(
        point[0] + -scaleConst * ortho1[0] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[0],
        point[1] + -scaleConst * ortho1[1] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[1],
        point[2] + -scaleConst * ortho1[2] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[2]);
    T.vertex3d(
        point[0] + scaleConst * ortho1[0] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[0],
        point[1] + scaleConst * ortho1[1] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[1],
        point[2] + scaleConst * ortho1[2] + ((s == 0) ? (y + stepSize * scaleConst) : y) * ortho2[2]);
  }
  T.end();

  if (o && o->m_drawMode & OD_ORI) { // draw normals
    if (o && o->m_drawMode & OD_MAGNITUDE) scaleMag *= weight;
    glDisable(GL_LIGHTING);
    T.begin(GL_LINES);
    for (y = -scaleConst; y <= scaleConst; y += stepSize * scaleConst) {
      for (x = -scaleConst; x <= scaleConst; x += stepSize * scaleConst) {
        T.vertex3d(
            point[0] + x * ortho1[0] + y * ortho2[0],
            point[1] + x * ortho1[1] + y * ortho2[1],
            point[2] + x * ortho1[2] + y * ortho2[2]);
        T.vertex3d(
            point[0] + x * ortho1[0] + y * ortho2[0] + scaleMag *  normal[0],
            point[1] + x * ortho1[1] + y * ortho2[1] + scaleMag *  normal[1],
            point[2] + x * ortho1[2] + y * ortho2[2] + scaleMag *  normal[2]);
      }
    }
    T.end();
  }
  return 0;
}


int drawScrew(const GAIM_FLOAT point[3], const GAIM_FLOAT direction[3], GAIM_FLOAT weight, GAIM_FLOAT pitch, int method /*= DRAW_SCREW_SPIRAL */, int flags /*= 0*/, object *o /*= NULL*/) {
  TubeDraw &T = gui_state->m_tubeDraw;
  GAIM_FLOAT x, scaleConst = g_state->m_clipDistance * sqrt(2.0);
  double z;
  e3ga e3gaR;
  float rotM[16];
  int stepSize = 64;
  double scale = ((o && o->m_drawMode & OD_MAGNITUDE) ? weight : 1.0) / 2.0;
  double arrowpos[3], arrowdir[3];

  glDisable(GL_LIGHTING);
  glPushMatrix();
  // translate to center, scale to weight 
  if (point) glTranslated(point[0], point[1], point[2]);
  //glScaled(m_int.m_scalar[0], m_int.m_scalar[0], m_int.m_scalar[0]);

  // rotate e3 to plane normal
  e3gaRve3(e3gaR, e3ga(GRADE1, direction[0], direction[1], direction[2]));
  e3gaRotorToOpenGLMatrix(e3gaR, rotM);
  glMultMatrixf(rotM);

  // draw 1 or more spirals
  switch (method) {
    case DRAW_SCREW_SPIRAL:
      T.begin(GL_LINE_STRIP);
      for (x = 0, z = -pitch / 2.0; x < M_PI * 2; x += (M_PI * 2) / stepSize, z += pitch / stepSize) {
        T.vertex3d(scale * sin(x), -scale * cos(x), z);
      }
      T.end();
      if ((flags & 0x01) || (o && o->m_drawMode & OD_ORI)) { // OD_ORI
        glEnable(GL_LIGHTING);
        arrowpos[0] = scale * sin(x - (M_PI * 2) / stepSize);
        arrowpos[1] = -scale * cos(x - (M_PI * 2) / stepSize);
        arrowpos[2] = z;
        arrowdir[0] = scale * (sin(x) - sin(x - (M_PI * 2) / stepSize));
        arrowdir[1] = scale * (cos(x) - cos(x - (M_PI * 2) / stepSize));
        arrowdir[2] = pitch / stepSize;
        z = 10 * sqrt(arrowdir[0] * arrowdir[0] + arrowdir[1] * arrowdir[1] + arrowdir[2] * arrowdir[2]);
        arrowdir[0] /= z;
        arrowdir[1] /= z;
        arrowdir[2] /= z;
        drawVector(arrowpos, arrowdir, 1.0);
      }
      break;
    case DRAW_SCREW_LINE:
      T.begin(GL_LINE_STRIP);
      for (x = 0, z = -scaleConst; z <= scaleConst; x += (M_PI * 2) / stepSize, z += pitch / stepSize) {
        T.vertex3d(scale * sin(x), -scale * cos(x), z);
      }
      T.end();
      if ((flags & 0x01) || (o && o->m_drawMode & OD_ORI)) { // OD_ORI
        glEnable(GL_LIGHTING);
        arrowpos[0] = scale * sin(x - (M_PI * 2) / stepSize);
        arrowpos[1] = -scale * cos(x - (M_PI * 2) / stepSize);
        arrowpos[2] = z;
        arrowdir[0] = scale * (sin(x) - sin(x - (M_PI * 2) / stepSize));
        arrowdir[1] = -scale * (cos(x) - cos(x - (M_PI * 2) / stepSize));
        arrowdir[2] = pitch / stepSize;
        z = 10 * sqrt(arrowdir[0] * arrowdir[0] + arrowdir[1] * arrowdir[1] + arrowdir[2] * arrowdir[2]);
        arrowdir[0] /= z;
        arrowdir[1] /= z;
        arrowdir[2] /= z;
        for (arrowpos[2] = scaleConst; arrowpos[2] >= -scaleConst; arrowpos[2] -= pitch) {
          drawVector(arrowpos, arrowdir, 1.0);
        }
      }
      break;
  }

  glPopMatrix();

  return 0;
}


int drawLinePair(const GAIM_FLOAT point1[3], const GAIM_FLOAT direction1[3], const GAIM_FLOAT point2[3], const GAIM_FLOAT direction2[3], const GAIM_FLOAT weight, int method /*= DRAW_LINE_HOOKS*/, int flags /*= 0*/, object *o /*= NULL*/) {
  GAIM_FLOAT orientation[3];
  int ori = 1;

  ori *= (0 < direction1[0]) - (direction1[0] <= 0);
  ori *= (0 < direction1[1]) - (direction1[1] <= 0);
  ori *= (0 < direction1[2]) - (direction1[2] <= 0);
  ori *= (0 < direction2[0]) - (direction2[0] <= 0);
  ori *= (0 < direction2[1]) - (direction2[1] <= 0);
  ori *= (0 < direction2[2]) - (direction2[2] <= 0);

  if (o && o->m_drawMode & OD_ORI) {
    if (ori < 0) {
      printf("switched\n");
      orientation[0] = point1[0] - point2[0];
      orientation[1] = point1[1] - point2[1];
      orientation[2] = point1[2] - point2[2];
      drawVector(point2, orientation, 1.0);
    }
    else {
      orientation[0] = point2[0] - point1[0];
      orientation[1] = point2[1] - point1[1];
      orientation[2] = point2[2] - point1[2];
      drawVector(point1, orientation, 1.0);
    }
  }
  return drawLine(point1, direction1, weight, DRAW_LINE_HOOKS, flags, o) + drawLine(point2, direction2, weight, DRAW_LINE_HOOKS, flags, o);
}


int drawDualLinePair(const GAIM_FLOAT point1[3], const GAIM_FLOAT direction1[3], const GAIM_FLOAT point2[3], const GAIM_FLOAT direction2[3], const GAIM_FLOAT weight, int method /*= DRAW_DUAL_LINE_PAIR*/, int flags /*= 0*/, object *o /*= NULL*/) {
  TubeDraw &T = gui_state->m_tubeDraw;
  int i, j, resolution = 15;
  GAIM_FLOAT focalpoints[2][1 + 2 * resolution][3],
             tmp[3],
             z,
             c,
             scaleConst = g_state->m_clipDistance * sqrt(2.0),
             focalStepSize = scaleConst / resolution,
             lineStepSize = 0.1;
  e3ga e3gaR;
  float rotM[16];

  GAIM_FLOAT normalize1 = (direction1[0] * direction1[0]) + (direction1[1] * direction1[1]) + (direction1[2] * direction1[2]);
  GAIM_FLOAT normalize2 = (direction2[0] * direction2[0]) + (direction2[1] * direction2[1]) + (direction2[2] * direction2[2]);
  for (i = 0 ; i <= 2 * resolution; ++i) {
    for (j = 0; j < 3; ++j) {
      focalpoints[0][i][j] = (point1[j] + ((i - resolution) * focalStepSize * direction1[j])) / normalize1;
      focalpoints[1][i][j] = (point2[j] + ((i - resolution) * focalStepSize * direction2[j])) / normalize2;
    }
  }

	glMatrixMode(GL_MODELVIEW);

  for (i = 0; i <= 2 * resolution; ++i) {

    if (!(o && o->m_drawMode & OD_ORI)) {
      glPushMatrix();
      glTranslated(focalpoints[0][i][0], focalpoints[0][i][1], focalpoints[0][i][2]);
    }

    for (j = 0; j <= 2 * resolution; ++j) {
      tmp[0] = focalpoints[1][j][0] - focalpoints[0][i][0];
      tmp[1] = focalpoints[1][j][1] - focalpoints[0][i][1];
      tmp[2] = focalpoints[1][j][2] - focalpoints[0][i][2];
      c = sqrt((tmp[0] * tmp[0]) + (tmp[1] * tmp[1]) + (tmp[2] * tmp[2]));
      tmp[0] /= c;
      tmp[1] /= c;
      tmp[2] /= c;

      if (o && o->m_drawMode & OD_ORI) {
        glPolygonMode(GL_FRONT_AND_BACK, (o && o->m_drawMode & OD_WIREFRAME) ? GL_LINE : GL_FILL);
        drawVector(focalpoints[0][i], tmp, c);
      }
      else {
        glPushMatrix();
        e3gaRve3(e3gaR, e3ga(GRADE1, tmp));
        e3gaRotorToOpenGLMatrix(e3gaR, rotM);
        glMultMatrixf(rotM);
        T.begin(GL_LINE_STRIP);
        T.vertex3d(0.0, 0.0, 0.0);
        T.vertex3d(0.0, 0.0, c);
        T.end();
        glPopMatrix();
      }
    }

    if (!(o && o->m_drawMode & OD_ORI)) {
      glPopMatrix();
    }
  }
  return 0;
}


int drawIdealPoint(const GAIM_FLOAT direction[3], const GAIM_FLOAT weight, int method /*= DRAW_IDEAL_POINT*/, int flags /*= 0*/, object *o /*= NULL*/) {
  GAIM_FLOAT dir1[3], dir2[3];

  if (o && o->m_drawMode & OD_MAGNITUDE) {
    dir1[0] = weight * direction[0];
    dir1[1] = weight * direction[1];
    dir1[2] = weight * direction[2];
    dir2[0] = -weight * direction[0];
    dir2[1] = -weight * direction[1];
    dir2[2] = -weight * direction[2];
  }
  else {
    dir1[0] = direction[0];
    dir1[1] = direction[1];
    dir1[2] = direction[2];
    dir2[0] = -direction[0];
    dir2[1] = -direction[1];
    dir2[2] = -direction[2];
  }

  drawPoint(dir1, 1.0, 0, o);
  drawPoint(dir2, 1.0, 0, o);

  if (o && o->m_drawMode & OD_ORI) {
    drawVector(dir2, dir1, 2.0);
  }
  return 0;
}
