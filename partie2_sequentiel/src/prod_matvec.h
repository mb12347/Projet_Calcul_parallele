/* prod_matvec.h
   Le schéma : produit matrice-vecteur y = A x (sans stocker la matrice)
   et second membre b. */

#ifndef PROD_MATVEC_H
#define PROD_MATVEC_H

/* y = A x */
void produit_matrice_vecteur(double *x, double *y, int Nx, int Ny,
                             double dt, double dx, double dy, double D);

/* Second membre du système linéaire à l'instant t = t^{n+1} */
void second_membre(double *b, double *u, double t, int Nx, int Ny, double dx, double dy,
                   double Lx, double Ly, double D, double dt, int cas);

#endif
