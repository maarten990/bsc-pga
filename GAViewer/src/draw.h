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

#ifndef _DRAW_FUNCTIONS_H_
#define _DRAW_FUNCTIONS_H_

#include "gaincludes.h"

class object;



// some drawing functions

/*
Draws a vector (a line with a 3D arrow at its head)

The tail of the vector will be at 'tail'. This argument
can be NULL if you simply want to draw at the origin.

The head of the vector will be at 'tail' + scale * 'dir'.

The arrow will scale with 'scale'
*/
int drawVector(const GAIM_FLOAT tail[3], const GAIM_FLOAT dir[3], GAIM_FLOAT scale);

#define DRAW_BV_CIRCLE 0
#define DRAW_BV_PARALLELOGRAM 1
#define DRAW_BV_PARALLELOGRAM_NO_VECTORS 2
#define DRAW_BV_CROSS 3
#define DRAW_BV_CURLYTAIL 4
#define DRAW_BV_SWIRL 5

/*
Draws a bivector at 'base'. 
'base' can be NULL if you want to draw at the origin.

The bivector is specified by its 'normal', 
and two factors ('factor1' and 'factor2').

The entire bivector will be scaled according to 'scale'

Several drawing methods can be selected using the 'method' argument,
which can have any 'DRAW_BV_XXX' value

If you specify an object 'o' this can be used to set
the appropriate colors (i.e. foreground, outline)

'flags' can be used to specify some extra options.
Currently:
0x01: draw something related to the orientation (if possible)
*/
int drawBivector(const GAIM_FLOAT base[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT factor1[3], const GAIM_FLOAT factor2[3], 
				 GAIM_FLOAT scale, int method = DRAW_BV_CIRCLE, int flags = 0, object *o = NULL);

#define DRAW_TV_SPHERE 0
#define DRAW_TV_CROSS 1
#define DRAW_TV_CURLYTAIL 2
#define DRAW_TV_PARALLELEPIPED 3
#define DRAW_TV_PARALLELEPIPED_NO_VECTORS 4
/*
Draws a trivector at 'base'. 
'base' can be NULL if you want to draw at the origin.

The entire trivector will be scaled according to 'scale'

Several drawing methods can be selected using the 'method' argument,
which can have any 'DRAW_TV_XXX' value

If you specify an object 'o' this can be used to set
the appropriate colors (i.e. foreground, outline)

'flags' can be used to specify some extra options.
Currently:
0x01: draw something related to the orientation (if possible)
*/
//int drawTriVector(const GAIM_FLOAT base[3], GAIM_FLOAT scale, int method = DRAW_TV_SPHERE, int flags = 0, object *o = NULL);
int drawTriVector(const GAIM_FLOAT base[3], GAIM_FLOAT scale, GAIM_FLOAT vector[][3], int method = DRAW_TV_SPHERE, int flags = 0, object *o = NULL);


#define DRAW_LINE_CURVE 0
#define DRAW_LINE_CURLYTAIL 1
#define DRAW_LINE_HOOKS 2
#define DRAW_LINE_HOOKCROSSES 3
/*
Draws a line trough 'point' with direction 'vector', and weight 'magnitude'.

Several drawing methods can be selected using the 'method' argument, which can have any 'DRAW_LINE_XXX' value.

If you specify an object 'o' this can be used to set 
the approriate colors (i.e. foreground, outline)

'flags' can be used to specify some extra options.
Currently:
0x01: draw something related to the orientation (if possible)
 */
int drawLine(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], GAIM_FLOAT magnitude, int method = DRAW_LINE_HOOKS, int flags = 0, object *o = NULL);


int drawPoint(const GAIM_FLOAT point[3], GAIM_FLOAT weight, int flags, object *o);


#define DRAW_PLANE 0
int drawPlane(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], GAIM_FLOAT weight, int method = DRAW_PLANE, int flags = 0, object *o = NULL);


#define DRAW_CIRCLE_HOOKS 0
int drawCircle(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], GAIM_FLOAT radius, GAIM_FLOAT weight, int method = DRAW_CIRCLE_HOOKS, int flags = 0, object *o = NULL);


#define DRAW_IDEAL_LINE_RADIUS 4
#define DRAW_IDEAL_LINE_HOOKS 5
int drawIdealLine(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], GAIM_FLOAT weight, int method = DRAW_IDEAL_LINE_HOOKS, int flags = 0, object *o = NULL);


#define DRAW_PENCIL 0
int drawLinePencil(const GAIM_FLOAT center[3], GAIM_FLOAT weight, const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], int method = DRAW_PENCIL, int flags = 0, object *o = NULL);


int drawCirclePencil(const GAIM_FLOAT center[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], GAIM_FLOAT weight, int method = DRAW_CIRCLE_HOOKS, int flags = 0, object *o = NULL);


#define DRAW_RULED_PLANE 1
int drawRuledPlane(const GAIM_FLOAT point[3], const GAIM_FLOAT normal[3], const GAIM_FLOAT ortho1[3], const GAIM_FLOAT ortho2[3], GAIM_FLOAT weight, int method = DRAW_RULED_PLANE, int flags = 0, object *o = NULL);


#define DRAW_SCREW_SPIRAL 0
#define DRAW_SCREW_LINE 1
int drawScrew(const GAIM_FLOAT point[3], const GAIM_FLOAT direction[3], GAIM_FLOAT weight, GAIM_FLOAT pitch, int method = DRAW_SCREW_SPIRAL, int flags = 0, object *o = NULL);


#define DRAW_LINE_PAIR 0
int drawLinePair(const GAIM_FLOAT point1[3], const GAIM_FLOAT direction1[3], const GAIM_FLOAT point2[3], const GAIM_FLOAT direction2[3], const GAIM_FLOAT weight, int method = DRAW_LINE_PAIR, int flags = 0, object *o = NULL);


#define DRAW_DUAL_LINE_PAIR 0
int drawDualLinePair(const GAIM_FLOAT point1[3], const GAIM_FLOAT direction1[3], const GAIM_FLOAT point2[3], const GAIM_FLOAT direction2[3], const GAIM_FLOAT weight, int method = DRAW_DUAL_LINE_PAIR, int flags = 0, object *o = NULL);


#define DRAW_IDEAL_POINT 0
int drawIdealPoint(const GAIM_FLOAT direction[3], const GAIM_FLOAT weight, int method = DRAW_IDEAL_POINT, int flags = 0, object *o = NULL);

/* Applies a rotor to the currently selected OpenGL matrix*/
void glApplyRotor();

int drawRegulus(e3ga &axis, double slant, double point[],
                double axis1[], double axis2[]);
#endif /* _DRAW_FUNCTIONS_H_ */
