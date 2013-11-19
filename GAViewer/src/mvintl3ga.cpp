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
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <map>
#include <assert.h>

#include "object.h"
#include "mvint.h"
#include "util.h"
#include "state.h"
#include "console/consolevariable.h"

#include <eigen2/Eigen/Core>
#include <eigen2/Eigen/Eigen>
#include <eigen2/Eigen/src/QR/EigenSolver.h>

USING_PART_OF_NAMESPACE_EIGEN

#define p3gaToL3ga(x) (consoleVariable("", x).castToL3ga()->l3())
#define l3gaToP3ga(x) (consoleVariable("", x).castToP3ga()->p())

double factorize_blade(const l3ga &B, int grade, l3ga (&factors)[7]);
int is_parallel(l3ga &a, l3ga &b, double epsilon);
bool only_coordinates_set(const l3ga &a, GAIM_FLOAT epsilon, int grade, int coordinate);
bool only_coordinates_set(const l3ga &a, GAIM_FLOAT epsilon, int grade, std::vector<int> coordinates);


p3ga direction(p3ga x);
p3ga moment(p3ga x);
p3ga cross_product(p3ga x, p3ga y); 
p3ga point_on_line(p3ga x, GAIM_FLOAT t);

MatrixXd versorToMatrix(const l3ga &R);
int regulusParameters(const l3ga &X, VectorXd *mainAxis, VectorXd *axis1,
                      VectorXd *axis2);
l3ga vectorToNullGA(const VectorXd &vec);
VectorXd nullGAToVector(const l3ga &ga);
void l3gaLineDirectionOffset(GAIM_FLOAT *offset, GAIM_FLOAT *direction, l3ga line);

int mvInt::interpret(const l3ga &X, int creationFlags /* = 0*/) {
  const GAIM_FLOAT epsilon = 1e-6; // rather arbitrary limit on fp-noise
  l3ga tmp, factors[7];
  mvInt tmpInt, interpret[7];
  int i;
  double s, x, y, z;
  int null_vectors = 0;
  int ideal_lines = 0;
  int intersection_count = 0;
  m_type = 0;
  m_type = MVI_L3GA;

  //TODO: dualize? (copied from mvint[cp]3ga.cpp; why?) 
  

  // check for zero euclidean norm multivectors
  if (X.norm_a() < epsilon) { // *********** zero blade
    m_type |= MVI_ZERO;
    m_valid = 1;
    return 0;
  }

  // determine the type of the multivectors
  int grade, type = X.mvType(&grade, epsilon); // temp todo test epsilon

  grade = (int)(log((double)grade) / log((double)2) + 0.45);
  if (creationFlags & OC_BLADE) type = GA_BLADE; // force blade interpretation
  else if (creationFlags & OC_VERSOR) type = GA_VERSOR; // force versor interpretation
  // maybe allow force null, dual, and other properties?

  if (type == GA_MULTIVECTOR) {
    m_type |= MVI_UNKNOWN;
    m_valid = 0;

    return 0;
  }

  if (type == GA_BLADE) { // ************************ all blades *********************
    GAIM_FLOAT X2 = (X * X).scalar(), weight2;
    s = factorize_blade(X, grade, factors);
    if (grade > 1) {
      m_scalar[0] = 1;
      m_scalar[1] = 1;
      for (i = 0; i < grade && factors[i].grade(); ++i) {
        if (i == 0) {
          factors[i] *= s;
        }
        interpret[i].interpret(factors[i]);

        if (fabs((factors[i] << factors[i]).scalar()) < epsilon) {
          ++null_vectors;
        }
        if (interpret[i].m_type & MVI_IDEAL) {
          ++ideal_lines;
        }
        m_scalar[0] *= interpret[i].m_scalar[0];
        m_scalar[1] *= interpret[i].m_scalar[1];
      }
      m_scalar[0] = fabs(m_scalar[0]);
    }

    switch (grade) {
      case 0: // ******************** scalar
        m_type |= MVI_SCALAR;
        //printf("scalar\n");
        m_scalar[0] = X.scalar();
        break;
      case 1: // ******************** line, screw, kine
        if (fabs(X2) < epsilon) {
          if (fabs(X[GRADE1][L3GA_E01]) < epsilon &&
              fabs(X[GRADE1][L3GA_E02]) < epsilon &&
              fabs(X[GRADE1][L3GA_E03]) < epsilon) {
            m_type |= MVI_IDEAL_LINE;
            //printf("ideal line\n");
            /*
            scalar 0: weight
            scalar 1: orientation
            vector 0: normal/reciprocal direction
            */
            m_scalar[0] = sqrt(X[GRADE1][L3GA_E23] * X[GRADE1][L3GA_E23] + X[GRADE1][L3GA_E31] * X[GRADE1][L3GA_E31] + X[GRADE1][L3GA_E12] * X[GRADE1][L3GA_E12]);

            m_vector[0][0] = X[GRADE1][L3GA_E23] / m_scalar[0];
            m_vector[0][1] = X[GRADE1][L3GA_E31] / m_scalar[0];
            m_vector[0][2] = X[GRADE1][L3GA_E12] / m_scalar[0];
          }
          else {
            m_type |= MVI_LINE;
            //printf("line\n");
            /*
            scalar 0: weight
            scalar 1: orientation
            point 0: point closest to origin
            vector 0: direction of line
            */


            // 6D-vector = [d, m]
            // direction = d
            // moment = m
            // distance = crossprod(m, d)
            // weight = sqrt(distance^2)
            m_scalar[0] = sqrt(X[GRADE1][L3GA_E01] * X[GRADE1][L3GA_E01] + X[GRADE1][L3GA_E02] * X[GRADE1][L3GA_E02] + X[GRADE1][L3GA_E03] * X[GRADE1][L3GA_E03]);

            m_vector[0][0] = X[GRADE1][L3GA_E01] / m_scalar[0];
            m_vector[0][1] = X[GRADE1][L3GA_E02] / m_scalar[0];
            m_vector[0][2] = X[GRADE1][L3GA_E03] / m_scalar[0];

            m_point[0][0] = ((X[GRADE1][L3GA_E12] * m_vector[0][1]) - (X[GRADE1][L3GA_E31] * m_vector[0][2])) / m_scalar[0];
            m_point[0][1] = ((X[GRADE1][L3GA_E23] * m_vector[0][2]) - (X[GRADE1][L3GA_E12] * m_vector[0][0])) / m_scalar[0];
            m_point[0][2] = ((X[GRADE1][L3GA_E31] * m_vector[0][0]) - (X[GRADE1][L3GA_E23] * m_vector[0][1])) / m_scalar[0];
          }
        }
        else {
          // (ideal) screw
          m_type |= MVI_SCREW;
          //printf("screw\n");

          /*
          scalar 0: weight
          scalar 1: orientation
          scalar 2: pitch (translation distance over 1 rotation).
          point 0: point closest to origin on the screw axis, if not ideal
          vector 0: direction of screw axis

          If pitch > 0, screw is called right-handed, otherwise left-handed.
          */

          m_scalar[0] = sqrt(fabs(X2));

          m_scalar[2] = (X[GRADE1][L3GA_E01] * X[GRADE1][L3GA_E23] + X[GRADE1][L3GA_E02] * X[GRADE1][L3GA_E31] + X[GRADE1][L3GA_E03] * X[GRADE1][L3GA_E12]) / (X[GRADE1][L3GA_E01] * X[GRADE1][L3GA_E01] + X[GRADE1][L3GA_E02] * X[GRADE1][L3GA_E02] + X[GRADE1][L3GA_E03] * X[GRADE1][L3GA_E03]);

          tmp = X - (m_scalar[2] * (X[GRADE1][L3GA_E01] * l3ga::e23 + X[GRADE1][L3GA_E02] * l3ga::e31 + X[GRADE1][L3GA_E03] * l3ga::e12));

          tmpInt.interpret(tmp);

          m_vector[0][0] = tmpInt.m_vector[0][0];
          m_vector[0][1] = tmpInt.m_vector[0][1];
          m_vector[0][2] = tmpInt.m_vector[0][2];

          m_point[0][0] = tmpInt.m_point[0][0];
          m_point[0][1] = tmpInt.m_point[0][1];
          m_point[0][2] = tmpInt.m_point[0][2];
        }

        m_scalar[1] = 1;
        m_scalar[1] *= (0 < m_vector[0][0]) - (m_vector[0][0] <= 0);
        m_scalar[1] *= (0 < m_vector[0][1]) - (m_vector[0][1] <= 0);
        m_scalar[1] *= (0 < m_vector[0][2]) - (m_vector[0][2] <= 0);

        m_valid = 1;
        break;
      case 2: // ******************** line pencil, skew line pair, 'line tangent', 'dual regulus pencil'
        if (null_vectors == 2) { 
          if (fabs(X2) < epsilon) {
            // Detect if lines are parallel.
            // Line L1=(a:b) and L2=(x:y) are parallel when:
            // c = a[i] / x[i]
            if (ideal_lines == null_vectors) {
              m_type |= MVI_IDEAL_LINE_PENCIL;
              printf("pencil of ideal lines\n");
              /*
              scalar 0: weight
              scalar 1: orientation
              vector 0: normal
              vector 1: orthogonal to vector 0 and 2
              vector 2: orthogonal to vector 0 and 1
              */
              //m_scalar[0] = s;

              z = sqrt((factors[0][GRADE1][L3GA_E23] * factors[0][GRADE1][L3GA_E23]) + (factors[0][GRADE1][L3GA_E31] * factors[0][GRADE1][L3GA_E31]) + (factors[0][GRADE1][L3GA_E12] * factors[0][GRADE1][L3GA_E12]));
              m_vector[1][0] = -factors[0][GRADE1][L3GA_E23] / z;
              m_vector[1][1] = -factors[0][GRADE1][L3GA_E31] / z;
              m_vector[1][2] = -factors[0][GRADE1][L3GA_E12] / z;

              z = sqrt((factors[1][GRADE1][L3GA_E23] * factors[1][GRADE1][L3GA_E23]) + (factors[1][GRADE1][L3GA_E31] * factors[1][GRADE1][L3GA_E31]) + (factors[1][GRADE1][L3GA_E12] * factors[1][GRADE1][L3GA_E12]));
              m_vector[2][0] = factors[1][GRADE1][L3GA_E23] / z;
              m_vector[2][1] = factors[1][GRADE1][L3GA_E31] / z;
              m_vector[2][2] = factors[1][GRADE1][L3GA_E12] / z;

              m_vector[0][0] = (m_vector[1][1] * m_vector[2][2]) - (m_vector[1][2] * m_vector[2][1]);
              m_vector[0][1] = (m_vector[1][2] * m_vector[2][0]) - (m_vector[1][0] * m_vector[2][2]);
              m_vector[0][2] = (m_vector[1][0] * m_vector[2][1]) - (m_vector[1][1] * m_vector[2][0]);
            }
            else if (is_parallel(factors[0], factors[1], epsilon)) {
              m_type |= MVI_RULED_PLANE;
              printf("ruled plane\n");
              /*
              scalar 0: weight;
              scalar 1: orientation
              point 0: point on plane closest to chosen origin
              vector 0: normal
              vector 1: line direction
              vector 2: orthogonal to vector 0 and 1
              */
              //m_scalar[0] = s;

              if (interpret[0].m_type & MVI_IDEAL) {
                factors[0] = p3gaToL3ga(l3gaToP3ga(factors[0]) - (p3ga::e0 ^ direction(l3gaToP3ga(factors[1]))));
              }
              if (interpret[1].m_type & MVI_IDEAL) {
                factors[1] = p3gaToL3ga(l3gaToP3ga(factors[1]) - (p3ga::e0 ^ direction(l3gaToP3ga(factors[0]))));
              }
              int t0 = 0, t1 = 0;
              p3ga p0 = point_on_line(l3gaToP3ga(factors[0]), t0),
                   p1 = point_on_line(l3gaToP3ga(factors[1]), t0),
                   p = point_on_line(p0^p1, t1);

              m_point[0][0] = p[GRADE1][P3GA_E1] / p[GRADE1][P3GA_E0];
              m_point[0][1] = p[GRADE1][P3GA_E2] / p[GRADE1][P3GA_E0];
              m_point[0][2] = p[GRADE1][P3GA_E3] / p[GRADE1][P3GA_E0];

              p0 = p0 - p1;
              p0 /= sqrt((p0 * p0).scalar());
              m_vector[2][0] = p0[GRADE1][P3GA_E1];
              m_vector[2][1] = p0[GRADE1][P3GA_E2];
              m_vector[2][2] = p0[GRADE1][P3GA_E3];

              p1.meet(l3gaToP3ga(factors[0]), l3gaToP3ga(factors[1]));
              p1 /= sqrt((p1 * p1).scalar());
              m_vector[1][0] = p1[GRADE1][P3GA_E1];
              m_vector[1][1] = p1[GRADE1][P3GA_E2];
              m_vector[1][2] = p1[GRADE1][P3GA_E3];

              p = cross_product(p0, p1);
              p /= sqrt((p * p).scalar());
              m_vector[0][0] = p[GRADE1][P3GA_E1];
              m_vector[0][1] = p[GRADE1][P3GA_E2];
              m_vector[0][2] = p[GRADE1][P3GA_E3];
            }
            else {
              m_type |= MVI_LINE_PENCIL;
              printf("line pencil\n");
              /*
              scalar 0: weight
              point 0: center
              vector 0: normal
              vector 1: orthogonal to vector 0 and 2
              vector 2: orthogonal to vector 0 and 1
              */
              //m_scalar[0] = (interpret[0].m_scalar[0] * interpret[1].m_scalar[0]);

              // point of intersection is the meet of the two homogeneous lines!
              p3ga a1 = (consoleVariable("", factors[0]).castToP3ga())->p(),
                   a2 = (consoleVariable("", factors[1]).castToP3ga())->p();
              p3ga intersection_point;
              intersection_point.meet(a1, a2);

              m_point[0][0] = intersection_point[GRADE1][P3GA_E1] / intersection_point[GRADE1][P3GA_E0];
              m_point[0][1] = intersection_point[GRADE1][P3GA_E2] / intersection_point[GRADE1][P3GA_E0];
              m_point[0][2] = intersection_point[GRADE1][P3GA_E3] / intersection_point[GRADE1][P3GA_E0];

              // normal to direction is:
              //    e3dual(factors1_d ^ factors2_d) = cross_product(factor1_d, factor2_d) 
              m_vector[0][0] = ((factors[0][GRADE1][L3GA_E02] * factors[1][GRADE1][L3GA_E03]) - (factors[0][GRADE1][L3GA_E03] * factors[1][GRADE1][L3GA_E02]));
              m_vector[0][1] = ((factors[0][GRADE1][L3GA_E03] * factors[1][GRADE1][L3GA_E01]) - (factors[0][GRADE1][L3GA_E01] * factors[1][GRADE1][L3GA_E03]));
              m_vector[0][2] = ((factors[0][GRADE1][L3GA_E01] * factors[1][GRADE1][L3GA_E02]) - (factors[0][GRADE1][L3GA_E02] * factors[1][GRADE1][L3GA_E01]));

              // normalize
              z = sqrt((m_vector[0][0] * m_vector[0][0]) + (m_vector[0][1] * m_vector[0][1]) + (m_vector[0][2] * m_vector[0][2]));
              m_vector[0][0] /= z;
              m_vector[0][1] /= z;
              m_vector[0][2] /= z;

              z = sqrt((factors[0][GRADE1][L3GA_E01] * factors[0][GRADE1][L3GA_E01]) + (factors[0][GRADE1][L3GA_E02] * factors[0][GRADE1][L3GA_E02]) + (factors[0][GRADE1][L3GA_E03] * factors[0][GRADE1][L3GA_E03]));
              m_vector[1][0] = factors[0][GRADE1][L3GA_E01] / z;
              m_vector[1][1] = factors[0][GRADE1][L3GA_E02] / z;
              m_vector[1][2] = factors[0][GRADE1][L3GA_E03] / z;

              m_vector[2][0] = (m_vector[0][1] * m_vector[1][2]) - (m_vector[0][2] * m_vector[1][1]);
              m_vector[2][1] = (m_vector[0][2] * m_vector[1][0]) - (m_vector[0][0] * m_vector[1][2]);
              m_vector[2][2] = (m_vector[0][0] * m_vector[1][1]) - (m_vector[0][1] * m_vector[1][0]);
            }
          }
          else {
            m_type |= MVI_LINE_PAIR;
            printf("line pair\n");
            /*
            scalar 0: weight;
            scalar 1: orientation
            point 0: point closest to origin of line 1
            vector 0: direction of line 1
            point 1: point closest to origin of line 2
            vector 1: direction of line 2
            */
            if (interpret[0].m_type & MVI_IDEAL) {
              m_type |= MVI_IDEAL_LINE_PAIR;
              /*
              scalar 0: weight;
              scalar 1: orientation
              vector 0: normal to the ideal line
              point 1: point closest to origin of the real line
              vector 1: direction of the real line
              */

              // TODO: Better name for this blade.
              // TODO: second, real line has always the same orientation.
            }
            else if (interpret[1].m_type & MVI_IDEAL) {
              m_type |= MVI_IDEAL_LINE_PAIR;
              tmpInt = interpret[0];
              interpret[0] = interpret[1];
              interpret[1] = tmpInt;
            }

            //m_scalar[0] = s;
            m_point[0][0] = interpret[0].m_point[0][0];
            m_point[0][1] = interpret[0].m_point[0][1];
            m_point[0][2] = interpret[0].m_point[0][2];
            m_vector[0][0] = interpret[0].m_vector[0][0];
            m_vector[0][1] = interpret[0].m_vector[0][1];
            m_vector[0][2] = interpret[0].m_vector[0][2];

            m_point[1][0] = interpret[1].m_point[0][0];
            m_point[1][1] = interpret[1].m_point[0][1];
            m_point[1][2] = interpret[1].m_point[0][2];
            m_vector[1][0] = interpret[1].m_vector[0][0];
            m_vector[1][1] = interpret[1].m_vector[0][1];
            m_vector[1][2] = interpret[1].m_vector[0][2];
          }
        }
        else {
          m_type |= MVI_UNKNOWN;
          printf("grade 2 unknown; \n\tX²:\t%2.2f\n\tf0:\t%s\n\tf1:\t%s\n\t#0:\t%d\n", X2, factors[0].string(), factors[1].string(), null_vectors);
          m_valid = 0;

          if (null_vectors == 1) {
            printf("parabolic pencil of linear complexes / line tangent?\n");
          }
          else if (null_vectors == 0) {
            printf("elliptic pencil of linear complexes / dual regulus pencil?\n");
          }
        }
        break;
      case 3: // ******************** line bundle/point, fied of lines/plane, regulus, double wheel pencil, ...
        intersection_count = (fabs(lcont(factors[0], factors[1]).scalar()) < epsilon)
          + (fabs(lcont(factors[1], factors[2]).scalar()) < epsilon)
          + (fabs(lcont(factors[2], factors[0]).scalar()) < epsilon);

        if (only_coordinates_set(X, epsilon, GRADE3, L3GA_E23_E31_E12)) {
          m_type |= MVI_IDEAL_PLANE;
          printf("ideal plane\n");
          /*
          scalar 0: weight
          scalar 1: orientation
          */

          m_valid = 1;
        }
        else if (!((interpret[0].m_type & MVI_DUAL) || (interpret[1].m_type & MVI_DUAL) || (interpret[2].m_type & MVI_DUAL)) &&
                 intersection_count == 3 && (X*X).scalar() == 0) {
          // point of intersection is the meet of the two homogeneous lines!
          p3ga a1 = (consoleVariable("", factors[0]).castToP3ga())->p(),
               a2 = (consoleVariable("", factors[1]).castToP3ga())->p(),
               a3 = (consoleVariable("", factors[2]).castToP3ga())->p();
          p3ga intersection1, intersection2, intersection3;
          intersection1.meet(a1, a2);
          intersection2.meet(a2, a3);
          intersection3.meet(a3, a1);
          int isIdeal = 0;

          // Normalize intersection points
          if (fabs(intersection1[GRADE1][P3GA_E0]) >= epsilon) {
            intersection1 = intersection1 / intersection1[GRADE1][P3GA_E0];
          }
          else {
            isIdeal += 1;
            z = sqrt(intersection1[GRADE1][P3GA_E1] * intersection1[GRADE1][P3GA_E1] + intersection1[GRADE1][P3GA_E2] * intersection1[GRADE1][P3GA_E2] + intersection1[GRADE1][P3GA_E3] * intersection1[GRADE1][P3GA_E3]);
            intersection1 /= z;
          }
            
          if (fabs(intersection2[GRADE1][P3GA_E0]) >= epsilon) {
            intersection2 = intersection2 / intersection2[GRADE1][P3GA_E0];
          }
          else {
            isIdeal += 1;
            z = sqrt(intersection2[GRADE1][P3GA_E1] * intersection2[GRADE1][P3GA_E1] + intersection2[GRADE1][P3GA_E2] * intersection2[GRADE1][P3GA_E2] + intersection2[GRADE1][P3GA_E3] * intersection2[GRADE1][P3GA_E3]);
            intersection2 /= z;
          }
            
          if (fabs(intersection3[GRADE1][P3GA_E0]) >= epsilon) {
            intersection3 = intersection3 / intersection3[GRADE1][P3GA_E0];
          }
          else {
            isIdeal += 1;
            z = sqrt(intersection3[GRADE1][P3GA_E1] * intersection3[GRADE1][P3GA_E1] + intersection3[GRADE1][P3GA_E2] * intersection3[GRADE1][P3GA_E2] + intersection3[GRADE1][P3GA_E3] * intersection3[GRADE1][P3GA_E3]);
            intersection3 /= z;
          }
            

          if (isIdeal == 3 && (X * X).scalar() == 0) {
            m_type |= MVI_IDEAL_POINT;
            //printf("ideal point\n");
            /*
            scalar 0: weight
            scalar 1: orientation
            vector 0: direction
            */
            m_vector[0][0] = m_scalar[1] * intersection2[GRADE1][P3GA_E1];
            m_vector[0][1] = m_scalar[1] * intersection2[GRADE1][P3GA_E2];
            m_vector[0][2] = m_scalar[1] * intersection2[GRADE1][P3GA_E3];
            m_valid = 1;
          }

          else if (fabs((intersection1 - intersection2).norm_a()) < epsilon && 
              fabs((intersection2 - intersection3).norm_a()) < epsilon &&
              fabs((intersection3 - intersection1).norm_a()) < epsilon) {
            // a1, a2 and a3 intersect in one point:
            m_type |= MVI_POINT;
            printf("brush/point\n");
            /*
            scalar 0: weight
            scalar 1: orientation
            point 0: position
            */

            m_point[0][0] = intersection1[GRADE1][P3GA_E1] / intersection1[GRADE1][P3GA_E0];
            m_point[0][1] = intersection1[GRADE1][P3GA_E2] / intersection1[GRADE1][P3GA_E0];
            m_point[0][2] = intersection1[GRADE1][P3GA_E3] / intersection1[GRADE1][P3GA_E0];
            m_valid = 1;
          }
          else {
            m_type |= MVI_PLANE;
            printf("plane\n");

            tmpInt.interpret(intersection1 ^ intersection2 ^ intersection3);

            m_point[0][0] = tmpInt.m_point[0][0];
            m_point[0][1] = tmpInt.m_point[0][1];
            m_point[0][2] = tmpInt.m_point[0][2];

            m_vector[0][0] = tmpInt.m_vector[0][0];
            m_vector[0][1] = tmpInt.m_vector[0][1];
            m_vector[0][2] = tmpInt.m_vector[0][2];

            m_vector[1][0] = tmpInt.m_vector[1][0];
            m_vector[1][1] = tmpInt.m_vector[1][1];
            m_vector[1][2] = tmpInt.m_vector[1][2];

            m_vector[2][0] = tmpInt.m_vector[2][0];
            m_vector[2][1] = tmpInt.m_vector[2][1];
            m_vector[2][2] = tmpInt.m_vector[2][2];

            m_valid = 1;
          }
        }
		else if (intersection_count == 2 &&
			     (X * X).scalar() < epsilon && (X * X).scalar() > -epsilon) {
          m_type |= MVI_LINE_PENCIL_PAIR;
          printf("line pencil pair\n");
          /*
          scalar 0: weight
          scalar 1: orientation
          point 0: center of pencil 1
          point 1: center of pencil 2
          vector 0: normal of pencil 1
          vector 1: normal of pencil 2
          vector 2: common line of pencil 1 and 2, orthogonal to vector 0 and 1
          */
          m_valid = 0;

          printf("f0: %s\nf1: %s\nf2: %s\n", factors[0].string(), factors[1].string(), factors[2].string());
          printf("f0.f1: %s\nf1.f2: %s\nf2.f0: %s\n", (factors[0] << factors[1]).string(), (factors[1] << factors[2]).string(), (factors[2] << factors[0]).string());
          if (fabs(lcont(factors[0], factors[1]).scalar()) >= epsilon) {
            printf("0\n");
            interpret[0].interpret(factors[0] ^ factors[2]);
            interpret[1].interpret(factors[1] ^ factors[2]);
          }
          else if (fabs(lcont(factors[1], factors[2]).scalar()) >= epsilon) {
            printf("1\n");
            interpret[0].interpret(factors[1] ^ factors[0]);
            interpret[1].interpret(factors[2] ^ factors[0]);
          }
          else if (fabs(lcont(factors[2], factors[0]).scalar()) >= epsilon) {
            printf("2\n");
            interpret[0].interpret(factors[2] ^ factors[1]);
            interpret[1].interpret(factors[0] ^ factors[1]);
          }

          m_point[0][0] = interpret[0].m_point[0][0];
          m_point[0][1] = interpret[0].m_point[0][1];
          m_point[0][2] = interpret[0].m_point[0][2];

          m_point[1][0] = interpret[1].m_point[0][0];
          m_point[1][1] = interpret[1].m_point[0][1];
          m_point[1][2] = interpret[1].m_point[0][2];

          m_vector[0][0] = interpret[0].m_vector[0][0];
          m_vector[0][1] = interpret[0].m_vector[0][1];
          m_vector[0][2] = interpret[0].m_vector[0][2];

          m_vector[1][0] = interpret[1].m_vector[0][0];
          m_vector[1][1] = interpret[1].m_vector[0][1];
          m_vector[1][2] = interpret[1].m_vector[0][2];

          m_vector[2][0] = m_point[1][0] - m_point[0][0];
          m_vector[2][1] = m_point[1][1] - m_point[0][1];
          m_vector[2][2] = m_point[1][2] - m_point[0][2];
        }
        else if ( (X * X).scalar() != 0) {
          m_type |= MVI_REGULUS;
          printf("regulus\n");

          /**
           * scalar0: slant
           * point0 : axis location
           * vector3: main axis direction
           * vector1: axis1 direction
           * vector2: axis2 direction
           */
          VectorXd axis(6), axis1(6), axis2(6);
          int slope = regulusParameters(X, &axis, &axis1, &axis2);

          l3ga mainAxisGA = vectorToNullGA(axis),
            axis1GA = vectorToNullGA(axis1),
            axis2GA = vectorToNullGA(axis2);

          double weight1 = sqrt(axis1GA[GRADE1][L3GA_E01] * axis1GA[GRADE1][L3GA_E01] + axis1GA[GRADE1][L3GA_E02] * axis1GA[GRADE1][L3GA_E02] + axis1GA[GRADE1][L3GA_E03] * axis1GA[GRADE1][L3GA_E03]);
          double weight2 = sqrt(axis2GA[GRADE1][L3GA_E01] * axis2GA[GRADE1][L3GA_E01] + axis2GA[GRADE1][L3GA_E02] * axis2GA[GRADE1][L3GA_E02] + axis2GA[GRADE1][L3GA_E03] * axis2GA[GRADE1][L3GA_E03]);
          std::cout << weight1 << ", " << weight2 << std::endl;

          l3gaLineDirectionOffset(m_point[0], m_vector[3], mainAxisGA);
          l3gaLineDirectionOffset(m_point[1], m_vector[1], axis1GA);
          l3gaLineDirectionOffset(m_point[2], m_vector[2], axis1GA);

          m_scalar[0] = (M_PI / 4.0) * slope;
          m_scalar[1] = weight1;
          m_scalar[2] = weight2;
          m_valid = 1;
        }
        else {
          m_type |= MVI_UNKNOWN;
          printf("grade 3 unknown;\n\tX²:\t%2.2f\n\tf0:\t%s\n\tf1:\t%s\n\tf2:\t%s\n\t#0:\t%d\n", X2, factors[0].string(), factors[1].string(), factors[2].string(), null_vectors);
          m_valid = 0;
        }
        break;
        /* Temporary solution. Interpretation through dualization. */
      case 4: // ******************** regulus pencil, dual line pair, parabolic linear congruence, "[hyperbolic linear congruence], bundle + field"
      case 5: // ******************** regulus bundle, rotation invariants, translation invariants
        printf("dual ");
        tmpInt.interpret(X.dual());
        m_type = tmpInt.m_type;
        m_type |= MVI_DUAL;
        m_valid = tmpInt.m_valid;
        //m_scalar[0] = tmpInt.m_scalar[0];
        m_scalar[1] = tmpInt.m_scalar[1];
        //m_scalar[2] = tmpInt.m_scalar[2];
        if (tmpInt.m_vector[0] != NULL) {
          m_vector[0][0] = tmpInt.m_vector[0][0];
          m_vector[0][1] = tmpInt.m_vector[0][1];
          m_vector[0][2] = tmpInt.m_vector[0][2];
        }
        if (tmpInt.m_vector[1] != NULL) {
          m_vector[1][0] = tmpInt.m_vector[1][0];
          m_vector[1][1] = tmpInt.m_vector[1][1];
          m_vector[1][2] = tmpInt.m_vector[1][2];
        }
        if (tmpInt.m_vector[2] != NULL) {
          m_vector[2][0] = tmpInt.m_vector[2][0];
          m_vector[2][1] = tmpInt.m_vector[2][1];
          m_vector[2][2] = tmpInt.m_vector[2][2];
        }
        if (tmpInt.m_point[0] != NULL) {
          m_point[0][0] = tmpInt.m_point[0][0];
          m_point[0][1] = tmpInt.m_point[0][1];
          m_point[0][2] = tmpInt.m_point[0][2];
        }
        if (tmpInt.m_point[1] != NULL) {
          m_point[1][0] = tmpInt.m_point[1][0];
          m_point[1][1] = tmpInt.m_point[1][1];
          m_point[1][2] = tmpInt.m_point[1][2];
        }
        if (tmpInt.m_point[2] != NULL) {
          m_point[2][0] = tmpInt.m_point[2][0];
          m_point[2][1] = tmpInt.m_point[2][1];
          m_point[2][2] = tmpInt.m_point[2][2];
        }
        break;
      case 6: // ******************** pseudoscalar
        /*********
        Some strange bug occurs when accessing the grade of 6- and 7-blades. Gaigen or GAViewer returns -2147483648.
        I haven't tracked the source of the problem, and I don't know if this is true for >7-blades as well.  Until this is fixed, checks for pseudoscalar is done in the default case.
             - Patrick.
        */
        m_type |= MVI_SPACE;
        //printf("pseudoscalar\n");
        /*
        scalar 0: weight
        scalar 1: orientation
        */
        m_scalar[0] = X[GRADE6][L3GA_I];
        m_valid = 1;
        break;
      default:
        if (X.dual().grade() == 0) {
          m_type |= MVI_SPACE;
          //printf("pseudoscalar\n");
          /*
          scalar 0: weight
          scalar 1: orientation
          */
          m_scalar[0] = X[GRADE6][L3GA_I];
          m_valid;
        }
        else {
          m_type |= MVI_UNKNOWN;
          printf("grade %d object unknown: %s\n",grade, X.string());
          m_valid = 0;
        }
        break;
    }
  }
  // TODO: Add interpretation for versors
  else if (type == GA_VERSOR) {
    l3ga g1, g2, g3, g4;
    g1.takeGrade(X, GRADE1);
    g2.takeGrade(X, GRADE2);
    g3.takeGrade(X, GRADE3);
    g4.takeGrade(X, GRADE4);
    if (X.maxGrade() == 2 && 
        fabs(X.scalar()) >= epsilon &&
        fabs(lcem(g1, g1).scalar()) < epsilon && 
        fabs(lcem(g2, g2).scalar()) >= epsilon) {
      /* For translation, only bivectors of 2 ideal lines are set */
      std::vector<int> translation_coordinates, perspective_coordinates;
      translation_coordinates.push_back(L3GA_E23_E31);
      translation_coordinates.push_back(L3GA_E23_E12);
      translation_coordinates.push_back(L3GA_E31_E12);

      perspective_coordinates.push_back(L3GA_E01_E02);
      perspective_coordinates.push_back(L3GA_E01_E03);
      perspective_coordinates.push_back(L3GA_E02_E03);

      if (only_coordinates_set(X, epsilon, GRADE2, translation_coordinates)) {
        printf("versor: translation\n");
        m_type |= MVI_VERSOR_TRANSLATION;
        m_valid = 1;
      } else if (only_coordinates_set(X, epsilon, GRADE2, perspective_coordinates)) {
        printf("versor: perspective transformation\n");
        m_type |= MVI_VERSOR_PERSPECTIVE_TRANSFORMATION;
        m_valid = 1;
      }
      else {
        printf("versor: directional scaling\n");
        m_type |= MVI_VERSOR_DIRECTIONAL_SCALING;
        m_valid = 1;
      }
    }
    else if (X.maxGrade() == 4 && 
        fabs(X.scalar()) >= epsilon && 
        fabs(lcem(g1, g1).scalar()) < epsilon &&
        fabs(lcem(g2, g2).scalar()) >= epsilon &&
        fabs(lcem(g3, g3).scalar()) < epsilon &&
        fabs(lcem(g4, g4).scalar()) >= epsilon) {
      std::vector<int> rotation_squeeze_coordinates;
      rotation_squeeze_coordinates.push_back(L3GA_E01_E23_E02_E31);
      rotation_squeeze_coordinates.push_back(L3GA_E01_E23_E03_E12);
      rotation_squeeze_coordinates.push_back(L3GA_E02_E31_E03_E12);
      if (only_coordinates_set(X, epsilon, GRADE4, rotation_squeeze_coordinates)) {
        if (X[GRADE4][L3GA_E01_E23_E02_E31] + X[GRADE4][L3GA_E01_E23_E03_E12] + X[GRADE4][L3GA_E02_E31_E03_E12] > 0) {
          printf("versor: rotation\n");
          m_type |= MVI_VERSOR_ROTATION;
          m_valid = 1;
        }
        else {
          printf("versor: lorentz transformation (squeeze)\n");
          m_type |= MVI_VERSOR_SQUEEZE;
          m_valid = 1;
        }
      }
    }
  }
  else {
    m_type |= MVI_UNKNOWN;
    printf("unknown versor\n");
    m_valid = 0;
  }

  return 0;
}


// According to factorization algorithm in Dorst et al., p.535.
double factorize_blade(const l3ga &B, int grade, l3ga (&factors)[7]) {
  if (B.grade() == 6) { 
    return 0; 
  }
  static const int num_basis_blades[] = {1, 6, 15, 20, 15, 6, 1};
  static const l3ga b[6] = {l3ga::e01, l3ga::e23, l3ga::e02, l3ga::e31, l3ga::e03, l3ga::e12};
  l3ga terminal = l3ga(1, 0.0);
  int k = grade, K = 1 << k;
  int basis_blades = num_basis_blades[k];
  double s = lcem(B, B).scalar();
  l3ga Bc = B / s;
  l3ga fi;
  l3ga E;
  l3ga e[6];
  l3ga tmp1, tmp2;
  int i = 0, max_i = 0, j = 0;
  double coordinates[20] = {0};

  if (k < 0) {
    // B is not an n-vector.
    return 0;
  }

  factors[k] = terminal;

  // If B is scalar, then we're quickly finished.
  if (k == 0) {
    factors[0] = B;
  }
  else {
    // Find the basis blade E of B with the largest coordinate.
    for (i = 0; i < basis_blades; ++i) {
      if (fabs(B[K][max_i]) < fabs(B[K][i])) {
        max_i = i;
      }
    }
    coordinates[max_i] = B[K][max_i];
    E.set(K, coordinates);

    // Determine the k basis vectors e[i] that span E.
    for (i = 0; (i < 6) && (j < k); ++i) {
      tmp1.lcem(b[i], E);
      // if b[i] is in E, add b[i] to e.
      if (tmp1.grade() > 0 || fabs(tmp1.scalar()) >= 1e-6) {
        e[j++] = b[i];
      }
    }

    // For all but one of the basis vectors e[i] of E:
    for (i = 0; i < (k - 1); ++i) {
      // Project e[i] onto Bc
      tmp1.scpem(B/s, (B/s).reverse());
      if (fabs(tmp1.scalar()) > 0.0) {
        tmp1.set(GRADE0, (1.0/tmp1.scalar()));
      }
      tmp1.gpem(Bc, tmp1);
      tmp2.lcem(e[i], Bc);
      fi.gpem(tmp2, tmp1);
      // Normalize fi. Add it to the list of factors.
      tmp1.lcem(fi, fi);
      fi.gpem(fi, l3ga(1, 1.0 / sqrt(tmp1.scalar())));
      factors[i] = fi;
      // Update Bc
      tmp1.scpem(fi, fi.reverse());
      tmp1.gpem(fi, l3ga(1.0 / tmp1.scalar()));
      Bc.lcem(fi, Bc);
    }
    factors[k - 1] = Bc;
  }
  return s;
}


int is_parallel(l3ga &a, l3ga &b, double epsilon) {
  // cross_product(a_direction, b_direction) ** 2 < epsilon?
  e3ga x = ((a[GRADE1][L3GA_E02] * b[GRADE1][L3GA_E03]) - (a[GRADE1][L3GA_E03] * b[GRADE1][L3GA_E02])) * e3ga::e1 +
    ((a[GRADE1][L3GA_E03] * b[GRADE1][L3GA_E01]) - (a[GRADE1][L3GA_E01] * b[GRADE1][L3GA_E03])) * e3ga::e2 +
    ((a[GRADE1][L3GA_E01] * b[GRADE1][L3GA_E02]) - (a[GRADE1][L3GA_E02] * b[GRADE1][L3GA_E01])) * e3ga::e3;
  return fabs((x * x).scalar()) < epsilon;
}


bool only_coordinates_set(const l3ga &a, GAIM_FLOAT epsilon, int grade, int coordinate) {
  std::vector<int> coords;
  coords.push_back(coordinate);
  return only_coordinates_set(a, epsilon, grade, coords);
}


bool only_coordinates_set(const l3ga &a, GAIM_FLOAT epsilon, int grade, std::vector<int> coordinates) {
  int i = 0, i_max;
  switch (grade) {
    case GRADE0:
    case GRADE6:
      i_max = 1;
      break;
    case GRADE1:
    case GRADE5:
      i_max = 6;
      break;
    case GRADE2:
    case GRADE4:
      i_max = 15;
      break;
    case GRADE3:
      i_max = 20;
      break;
  }
  for (; i < i_max; ++i) {
    if ((std::find(coordinates.begin(), coordinates.end(), i) == coordinates.end()) &&
        (fabs(a[grade][i]) >= epsilon)) {
      return false;
    }
  }
  return true;
}


p3ga direction(p3ga x) {
  return x[GRADE2][P3GA_E1_E0] * p3ga::e1 + x[GRADE2][P3GA_E2_E0] * p3ga::e2 + x[GRADE2][P3GA_E3_E0] * p3ga::e3;
}


p3ga moment(p3ga x) {
  return (x - (direction(x) * p3ga::e0)) << (p3ga::e1 * p3ga::e2 * p3ga::e3);
}


p3ga cross_product(p3ga x, p3ga y) { 
  return lcem((x * y), -(p3ga::e1 * p3ga::e2 * p3ga::e3));
}


p3ga point_on_line(p3ga x, GAIM_FLOAT t) {
  return (cross_product(moment(x), direction(x)) + (t * direction(x)) + (direction(x) * direction(x) * p3ga::e0)) / (direction(x) * direction(x));
}

// turns an l3ga object to a 6D vector on the {e01, e02, e03, e23, e31, e12}
// basis
VectorXd nullGAToVector(const l3ga &ga)
{
  return (VectorXd(6) <<
          ga[GRADE1][L3GA_E01],
          ga[GRADE1][L3GA_E02],
          ga[GRADE1][L3GA_E03],
          ga[GRADE1][L3GA_E23],
          ga[GRADE1][L3GA_E31],
          ga[GRADE1][L3GA_E12]).finished();
}

// turns a 6D vector on the {e01, e02, e03, e23, e31, e12} basis to an l3ga
// object
l3ga vectorToNullGA(const VectorXd &vec)
{
  return (vec[0] * l3ga::e01) +
         (vec[1] * l3ga::e02) +
         (vec[2] * l3ga::e03) +
         (vec[3] * l3ga::e23) +
         (vec[4] * l3ga::e31) +
         (vec[5] * l3ga::e12);
}

// Returns a copy of A where values very close to 0 are replaced with an actual
// zero
MatrixXd sanitize(const MatrixXd &A)
{
  MatrixXd As(A);
  for (int i = 0; i < As.cols(); ++i) {
    for (int j = 0; j < As.rows(); ++j) {
      if ( (-0.000001 < As(j, i)) && (As(j, i) < 0.000001) ) {
        As(j, i) = 0;
      }
    }
  }

  return As;
}

/* Transforms a vector on on the null basis by a given regulus according to the
 * regulus operator: R[x] = inv(R) x R
 */
VectorXd transformVersor(const l3ga &R, VectorXd vec)
{
  l3ga ga = vectorToNullGA(vec);
  return nullGAToVector(R.inverse() * ga * R);
}

void testTransformation(const MatrixXd &basis, const l3ga &R, const MatrixXd &A)
{
  
  // check if the basis vectors transform properly
  for (int i = 0; i < 6; ++i) {
    VectorXd b = basis.col(i);
    VectorXd bt = A * b;

    l3ga bga = vectorToNullGA(b);
    l3ga bgat = R.inverse() * bga * R;
    VectorXd bgatv = nullGAToVector(bgat);

    // filter out "almost 0" values
    for (int i = 0; i < 6; ++i)
      if ((bgatv[i] <  0.00001) &&
          (bgatv[i] > -0.00001))
        bgatv[i] = 0;
        
    if (bgatv != bt) {
      std::cout << "WARNING: vector transformation incorrect" << std::endl;
      std:: cout << bt.transpose() << std::endl;
      std:: cout << bgatv.transpose() << std::endl << std::endl;
    }
  }

  // construct the metric matrix for the null basis
  MatrixXd M(6, 6);
  M <<
    MatrixXd::Zero(3, 3),     MatrixXd::Identity(3, 3),
    MatrixXd::Identity(3, 3), MatrixXd::Zero(3, 3);

  // symmetry
  // condition: A = M transp(A) M
  if( sanitize(M * A.transpose() * M) != A ) {
    std::cout << "WARNING: transformation not symmetric" << std::endl << std::endl;
    std::cout << "Transformation: " << std::endl;
    std::cout << A << std::endl << std::endl;
  };

  // orthogonality check (Dorst unreleased paper, section 3.7)
  if( sanitize(M * A.transpose() * M * A) != MatrixXd::Identity(6, 6) ) {
    std::cout << "WARNING: transformation not orthogonal" << std::endl;
    std::cout << "M * A.T * M * A = " << std::endl
              << sanitize(M * A.transpose() * M * A) << std::endl << std::endl;
  }
}

/* Geometric Algebra For Computer Science, page 114 */
MatrixXd versorToMatrix(const l3ga &R)
{
  MatrixXd basis(6, 6), transform(6, 6);
  basis << MatrixXd::Identity(6, 6);

  // construct the metric matrix for the null basis
  MatrixXd M(6, 6);
  M <<
    MatrixXd::Zero(3, 3),     MatrixXd::Identity(3, 3),
    MatrixXd::Identity(3, 3), MatrixXd::Zero(3, 3);

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      // take the 0th element because this returns a vector
      transform(j, i) = (basis.row(j) * transformVersor(R, basis.col(i)) )[0];
    }
  }

  transform = sanitize(transform);
  //testTransformation(basis, R, transform);
  
  return transform;
}

int findOddOneOut(int *slope, int *index,
                  MatrixXd &eigenvectors, VectorXd &eigenvalues)
{
  // construct the metric matrix for the null basis
  MatrixXd M(6, 6);
  M <<
    MatrixXd::Zero(3, 3),     MatrixXd::Identity(3, 3),
    MatrixXd::Identity(3, 3), MatrixXd::Zero(3, 3);

  // there are three eigenvectors with corresponding eigenvalue 1
  // either two of them square to 1, or two of them square to -1
  // the odd one out corresponds to the main axis
  std::vector<int> posSquared1, negSquared1,
    posSquaredM1, negSquaredM1;

  // find the squared values for each eigenvector
  for (int i = 0; i < eigenvectors.cols(); ++i) {
    if ( (0.999999 < eigenvalues[i] &&
          eigenvalues[i] < 1.000001) ) {
      if ((eigenvectors.col(i).transpose() * M * eigenvectors.col(i))[0] > 0)
        posSquared1.push_back(i);
      else
        negSquared1.push_back(i);
    }

    else {
      if ((eigenvectors.col(i).transpose() * M * eigenvectors.col(i))[0] > 0)
        posSquaredM1.push_back(i);
      else
        negSquaredM1.push_back(i);
    }
  }

  if (posSquared1.size() == 1) {
    *index = posSquared1[0];
    *slope = 1;
  }

  else if (negSquared1.size() == 1) {
    *index =  negSquared1[0];
    *slope = -1;
  }

  else if (posSquaredM1.size() == 1) {
    *index =  posSquaredM1[0];
    *slope = -1;
  }
  else if (negSquaredM1.size() == 1) {
    *index =  negSquaredM1[0];
    *slope = 1;
  }

  else {
    return 1;
  }

  return 0;
}

bool differentSign(double x, double y)
{
  if ( (x < 0 && y >= 0) || (y < 0 && x >= 0) ) {
    return true;
  } else {
    return false;
  }
}

int findAssociate(int index, MatrixXd &eigenvectors, VectorXd &eigenvalues)
{
  // construct the metric matrix for the null basis
  MatrixXd M(6, 6);
  M <<
    MatrixXd::Zero(3, 3),     MatrixXd::Identity(3, 3),
    MatrixXd::Identity(3, 3), MatrixXd::Zero(3, 3);
  
  std::map<int, int> componentCounts;

  VectorXd axis = eigenvectors.col(index);
  VectorXd squared, added;
 
  VectorXd candidate;

  for (int i = 0; i < eigenvectors.cols(); ++i) {
    // only compare to the eigenvectors with flipped sign
    if (differentSign(eigenvalues[i], eigenvalues[index])) {
      candidate = eigenvectors.col(i);

      added = candidate + axis;
      squared = added.transpose() * M * added;

      // if the added components form a line, check how many grade 1 components it has
      if ((squared[0] < 0.0000001) && (squared[0] > -0.0000001)) {
        componentCounts[i] = 0;
        if (added[0] != 0) componentCounts[i] += 1;
        if (added[1] != 0) componentCounts[i] += 1;
        if (added[2] != 0) componentCounts[i] += 1;
        if (added[3] != 0) componentCounts[i] += 1;
        if (added[4] != 0) componentCounts[i] += 1;
        if (added[5] != 0) componentCounts[i] += 1;
      }
    }
  }
  
  if (componentCounts.empty())
    return -1;
  else {
    int lowestInd, lowestComp = 999;
    for (std::map<int, int>::iterator it = componentCounts.begin();
         it != componentCounts.end(); ++it) {
      if (it->second < lowestComp) {
        lowestInd = it->first;
        lowestComp = it->second;
      }
    }

    return lowestInd;
  }
}

bool closeEnough(double x, double y)
{
  return abs(x - y) < 0.0001;
}

/**
 * Find the secondary axes.
 * index1 and index2: pointers to which the axis indices will be written
 * mainIndex: index of the main axis
 */
void secondaryAxes(int *index1, int *index2, int mainIndex,
                  MatrixXd &eigenvectors, VectorXd &eigenvalues)
{
  int done = 0;

  for (int i = 0; i < eigenvectors.cols(); ++i) {
    if (closeEnough(eigenvalues[i], eigenvalues[mainIndex]) && (i != mainIndex)) {
      switch (done) {
      case 0:
        *index1 = i;
        done += 1;
        break;
      case 1:
        *index2 = i;
        done += 1;
        break;
      default:
        return;
      }
    }
  }
}

void printvector(const l3ga &X)
{
  if (X.grade() != 1) {
    return;
  }

  int components[6] = {L3GA_E01, L3GA_E02, L3GA_E03, L3GA_E23, L3GA_E31, L3GA_E12};
  std::vector<std::string> names;
  names.push_back("e01");
  names.push_back("e02");
  names.push_back("e03");
  names.push_back("e23");
  names.push_back("e31");
  names.push_back("e12");
  bool first = true;

  GAIM_FLOAT c;
  for (int i = 0; i < 6; ++i) {
    c = X[GRADE1][ components[i] ];

    if (c) {
      printf("%.3f*%s", c, names[i].c_str());
      if (first) printf(" + ");
      first = false;
    }
  }
}

void debugprint(const MatrixXd &vectors, const VectorXd &values)
{
  MatrixXd M(6, 6);
  M <<
    MatrixXd::Zero(3, 3),     MatrixXd::Identity(3, 3),
    MatrixXd::Identity(3, 3), MatrixXd::Zero(3, 3);

  std::cout << "Eigenvectors ([index], eigenvalue, square, vector): " << std::endl;
  for (int i = 0; i < vectors.cols(); ++i) {
    std::cout << "[" << i << "]\t" << values[i] << ",\t" << vectors.col(i).transpose() * M * vectors.col(i) << ",\t";
    printvector( vectorToNullGA(vectors.col(i)) );
    printf(",\n");
  }
}

int regulusParameters(const l3ga &X, VectorXd *mainAxis, VectorXd *axis1,
                      VectorXd *axis2)
{
  MatrixXd vectors;
  VectorXd values;

  printf("===== Regulus debug output =====\n");

  MatrixXd transform = versorToMatrix(X);
  std::cout << "Transformation matrix: " << std::endl << transform << std::endl << std::endl;

  // obtain eigenvalues and vectors
  Eigen::EigenSolver<MatrixXd> eigen(transform);
  values  = eigen.eigenvalues().real();
  vectors = eigen.eigenvectors().real();

  debugprint(vectors, values);

  // normalize eigenvectors, unless they square to 0
  for (int i = 0; i < vectors.cols(); ++i) {
    l3ga temp = vectorToNullGA(vectors.col(i));
    if (*(temp * temp)[GRADE0] != 0) {
      vectors.col(i) = nullGAToVector(temp / sqrt( fabs(*(temp * temp)[GRADE0]) ) );
    }
  }

  std::cout << std::endl << "After normalization: " << std::endl;
  debugprint(vectors, values);

  // multiply each vector by -1 if its real component is negative
  // appears to work well for most cases
  for (int i = 0; i < vectors.cols(); ++i) {
    GAIM_FLOAT e01 = vectorToNullGA(vectors.col(i))[GRADE1][L3GA_E01];
    GAIM_FLOAT e02 = vectorToNullGA(vectors.col(i))[GRADE1][L3GA_E02];
    GAIM_FLOAT e03 = vectorToNullGA(vectors.col(i))[GRADE1][L3GA_E03];

    if (e01 < -0.000001) {
      vectors.col(i) *= -1;
    } else if (e01 < 0.000001 && e02 < -0.000001) {
      vectors.col(i) *= -1;
    } else if (e01 < 0.000001 && e02 < 0.000001 && e03 < -0.000001) {
      vectors.col(i) *= -1;
    }
  }

  std::cout << std::endl << "After multiplication by -1: " << std::endl;
  debugprint(vectors, values);

  // find the column indices of the 3 axes
  int slope, mainIndex, index1, index2;
  int status = findOddOneOut(&slope, &mainIndex,
                             vectors, values);

  if (status) {
    // TODO: add logging
    std::cout << "Error: Could not find mainAxis." << std::endl;
  }

  secondaryAxes(&index1, &index2, mainIndex, vectors, values);

  int assocIndexMain = findAssociate(mainIndex, vectors, values);
  int assocIndex1    = findAssociate(index1, vectors, values);
  int assocIndex2    = findAssociate(index2, vectors, values);
  if ((assocIndexMain == -1) || (assocIndex1 == -1) || (assocIndex2 == -1)) {
    std::cout << "Error: Could not find associate." << std::endl;
    assert(false);
  }

  printf("\nMain axis + associate:\t\t[%d] + [%d]\n", mainIndex, assocIndexMain);
  printf("Secondary axis 1 + associate:\t[%d] + [%d]\n", index1, assocIndex1);
  printf("Secondary axis 2 + associate:\t[%d] + [%d]\n", index2, assocIndex2);

  // the axes are proportional to a factor of sqrt(2)
  *mainAxis = (vectors.col(mainIndex) + vectors.col(assocIndexMain)) / sqrt(2.0);
  *axis1 = (vectors.col(index1) + vectors.col(assocIndex1)) / sqrt(2.0);
  *axis2 = (vectors.col(index2) + vectors.col(assocIndex2)) / sqrt(2.0);

  // debug output
  std::cout << "Main axis: ";
  vectorToNullGA(*mainAxis).print();
  std::cout << "Axis 1: ";
  vectorToNullGA(*axis1).print();
  std::cout << "Axis 2: ";
  vectorToNullGA(*axis2).print();

  printf("\n===== End of regulus debug output =====\n\n");
  return slope;
}

// converts an l3ga line to it's offset and direction as arrays of 3 elements
void l3gaLineDirectionOffset(GAIM_FLOAT *offset, GAIM_FLOAT *direction, l3ga line)
{
  // weight
  double weight = sqrt(line[GRADE1][L3GA_E01] * line[GRADE1][L3GA_E01] + line[GRADE1][L3GA_E02] * line[GRADE1][L3GA_E02] + line[GRADE1][L3GA_E03] * line[GRADE1][L3GA_E03]);

  // direction
  direction[0] = line[GRADE1][L3GA_E01] / weight;
  direction[1] = line[GRADE1][L3GA_E02] / weight;
  direction[2] = line[GRADE1][L3GA_E03] / weight;

  // offset
  offset[0] = ((line[GRADE1][L3GA_E12] * direction[1]) - (line[GRADE1][L3GA_E31] * direction[2])) / weight;
  offset[1] = ((line[GRADE1][L3GA_E23] * direction[2]) - (line[GRADE1][L3GA_E12] * direction[0])) / weight;
  offset[2] = ((line[GRADE1][L3GA_E31] * direction[0]) - (line[GRADE1][L3GA_E23] * direction[1])) / weight;
}
