/* solveur.c
   Résolution du système linéaire A u^{n+1} = b par BiCGStab.
   Le produit matrice-vecteur est dans prod_matvec.c. */

#include <stdlib.h>
#include <math.h>
#include "solveur.h"
#include "prod_matvec.h"


/* Produit scalaire : (a, b) = somme pour k = 0..n-1 de a[k]*b[k] */
double produit_scalaire(double *a, double *b, int n)
{
   double somme = 0.0;
   for (int k = 0; k < n; k++)
      somme += a[k] * b[k];
   return somme;
}


/* Norme euclidienne : ||v|| = sqrt((v, v)) */
double norme(double *v, int n)
{
   return sqrt(produit_scalaire(v, v, n));
}


/* BiCGStab pour A x = b.
   Contrairement au gradient conjugué, marche aussi si A n'est pas symétrique
   (utile pour les conditions de Robin entre sous-domaines, question 5).
   Entrée : x = point de départ, b = second membre
   Sortie : x = solution
   Retour : nombre d'itérations effectuées
   Critère d'arrêt : ||s|| < eps ou ||r|| < eps (absolu) */
int bicgstab(double *x, double *b, int Nx, int Ny,
             double dt, double dx, double dy, double D,
             double eps, int kmax)
{
   int n = Nx*Ny;
   int nb_iter = 0;

   double *r         = (double*) calloc(n, sizeof(double));
   double *r_chapeau = (double*) calloc(n, sizeof(double));
   double *p         = (double*) calloc(n, sizeof(double));
   double *v         = (double*) calloc(n, sizeof(double));
   double *s         = (double*) calloc(n, sizeof(double));
   double *l         = (double*) calloc(n, sizeof(double));

   double alpha = 1.0, omega = 1.0;
   double rho   = 1.0, rho0  = 1.0;

   /* r = b - A*x */
   produit_matrice_vecteur(x, r, Nx, Ny, dt, dx, dy, D);
   for (int i = 0; i < n; i++) {
      r[i]         = b[i] - r[i];
      r_chapeau[i] = r[i];
   }

   for (int k = 0; k < kmax; k++) {
      nb_iter++;

      rho = produit_scalaire(r_chapeau, r, n);
      if (fabs(rho) < 1e-20) break;

      double beta = (rho / rho0) * (alpha / omega);

      for (int i = 0; i < n; i++)
         p[i] = r[i] + beta * (p[i] - omega * v[i]);

      /* v = A*p */
      produit_matrice_vecteur(p, v, Nx, Ny, dt, dx, dy, D);

      double denom = produit_scalaire(r_chapeau, v, n);
      if (fabs(denom) < 1e-20) break;

      alpha = rho / denom;

      for (int i = 0; i < n; i++)
         s[i] = r[i] - alpha * v[i];

      if (norme(s, n) < eps) {
         for (int i = 0; i < n; i++)
            x[i] += alpha * p[i];
         break;
      }

      /* l = A*s */
      produit_matrice_vecteur(s, l, Nx, Ny, dt, dx, dy, D);

      double l_l = produit_scalaire(l, l, n);
      if (fabs(l_l) < 1e-20) break;

      omega = produit_scalaire(l, s, n) / l_l;

      for (int i = 0; i < n; i++) {
         x[i] += alpha * p[i] + omega * s[i];
         r[i]  = s[i] - omega * l[i];
      }

      if (norme(r, n) < eps) break;

      rho0 = rho;
   }

   free(r); free(r_chapeau); free(p); free(v); free(s); free(l);
   return nb_iter;
}
