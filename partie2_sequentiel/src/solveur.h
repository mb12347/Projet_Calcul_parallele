/* solveur.h
   Produit scalaire, norme et BiCGStab. */

#ifndef SOLVEUR_H
#define SOLVEUR_H

/* (a, b) = somme des a[k]*b[k] */
double produit_scalaire(double *a, double *b, int n);

/* ||v|| = sqrt((v, v)) */
double norme(double *v, int n);

/* Résout A x = b, retourne le nombre d'itérations */
int bicgstab(double *x, double *b, int Nx, int Ny,
             double dt, double dx, double dy, double D,
             double eps, int kmax);

#endif
