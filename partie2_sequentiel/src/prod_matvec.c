/* prod_matvec.c
   Produit matrice-vecteur y = A x (sans stocker la matrice)
   et second membre b.
   Propriétés (à justifier dans le rapport) :
      - creuse, pentadiagonale (tridiagonale par blocs de taille Nx)
      - symétrique
      - à diagonale strictement dominante (grâce au 1 sur la diagonale)
      - définie positive
   On ne stocke PAS la matrice : les coefficients sont recalculés dans produit_matrice_vecteur.
*/

#include "prod_matvec.h"
#include "fonctions.h"   /* pour f, g, h (second membre) */


/* y = A x, sans stocker A.
(les voisins situés sur le bord du domaine sont déjà dans le second membre) */

void produit_matrice_vecteur(double *x, double *y, int Nx, int Ny,
                             double dt, double dx, double dy, double D)
{
   int n = Nx*Ny;
   double cx = D*dt/(dx*dx);
   double cy = D*dt/(dy*dy);
   double c_diag = 1.0 + 2.0*cx + 2.0*cy;  /* coef_diag */
   for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
         int k = j*Nx + i;  /* indice vecteur 2d -> 1d */
         y[k] = c_diag * x[k];
         if (i > 0)      y[k] += -cx * x[k-1]   ;      /* voisin gauche */
         if (i < Nx-1)   y[k] += -cx * x[k+1];      /* voisin droit */
         if (j > 0)      y[k] += -cy * x[k-Nx];     /* voisin bas */
         if (j < Ny-1)   y[k] += -cy * x[k+Nx];     /* voisin haut */
      }
   }
}

/* Second membre du système linéaire à l'instant t = t^{n+1}.
   Entrée : u = u^n (taille Nx*Ny)
   Sortie : b       (taille Nx*Ny) */
void second_membre(double *b, double *u, double t, int Nx, int Ny, double dx, double dy,
                   double Lx, double Ly, double D, double dt, int cas)
{
   double cx = D*dt/(dx*dx);
   double cy = D*dt/(dy*dy);
   for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
         int k = j*Nx + i;          /* indice vecteur 2d -> 1d */
         double x = (i+1)*dx;       /* coordonnées du point (i,j) */
         double y = (j+1)*dy;
         b[k] = u[k] + dt * f(x, y, t, Lx, Ly, cas);
         if (i == 0)      b[k] += cx * h(0.0, y, t, cas);   /* bord gauche */
         if (i == Nx-1)   b[k] += cx * h(Lx,  y, t, cas);   /* bord droit */
         if (j == 0)      b[k] += cy * g(x, 0.0, t, cas);   /* bord bas */
         if (j == Ny-1)   b[k] += cy * g(x, Ly,  t, cas);   /* bord haut */
      }
   }
}
