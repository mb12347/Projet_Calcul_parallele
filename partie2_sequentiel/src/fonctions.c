/* fonctions.c
   Données du problème : f, g, h, solution exacte, erreur.

   Validation (Lx = Ly = 1, D = 1) :
     cas 1 : stationnaire,   f = 2*(x - x^2 + y - y^2),       g = 0,                h = 0
     cas 2 : stationnaire,   f = sin(x) + cos(y),             g = sin(x) + cos(y),  h = sin(x) + cos(y)
     cas 3 : instationnaire, f = exp(-(x-Lx/2)^2) * exp(-(y-Ly/2)^2) * cos(pi/2 * t),
                                                               g = 0,                h = 1
*/

#include <math.h>
#include "fonctions.h"

double f(double x, double y, double t, double Lx, double Ly, int cas)
{
   if (cas == 1) {
      return 2.0 * (x - x*x + y - y*y);
   }
   else if (cas == 2) {
      return sin(x) + cos(y);
   }
   else if (cas == 3) {
      return exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos(M_PI/2*t);
   }
   return 0.0;
}

/* Condition de Dirichlet sur Gamma0 (bords bas y = 0 et haut y = Ly).*/
double g(double x, double y, double t, int cas)
{
   if (cas == 1) {
      return 0.0;
   }
   else if (cas == 2) {
      return sin(x) + cos(y);
   }
   else if (cas == 3) {
      return 0.0;
   }
   return 0.0;
}


/* Condition de Dirichlet sur Gamma1 (bords gauche x = 0 et droit x = Lx). */
double h(double x, double y, double t, int cas)
{
   if (cas == 1) {
      return 0.0;
   }
   else if (cas == 2) {
      return sin(x) + cos(y);
   }
   else if (cas == 3) {
      return 1.0;
   }
   return 0.0;
}


/* Solution exacte des cas stationnaires (cas 1 et 2).
   Pas de solution exacte pour le cas 3. */
double solution_exacte(double x, double y, int cas)
{
   if (cas == 1) {
      return x * (1.0 - x) * y * (1.0 - y);
   }
   else if (cas == 2) {
      return sin(x) + cos(y);
   }
   return 0.0;
}


/* Erreur entre la solution calculée et la solution exacte (cas 1 et 2) :
      err = max sur k de | u_k - solution_exacte(x_i, y_j) |
   avec x_i = (i+1)*dx et y_j = (j+1)*dy. */
double erreur_max(double *u, int Nx, int Ny, double dx, double dy, int cas)
{
   double err = 0.0;
   for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
         int k = j*Nx + i;          /* indice vecteur 2d -> 1d */
         double x = (i+1)*dx;       /* coordonnées du point (i,j) */
         double y = (j+1)*dy;
         double e = fabs(u[k] - solution_exacte(x, y, cas));
         if (e > err) err = e;
      }
   }
   return err;
}
