/* fonctions.h
   Données du problème : f, g, h, solution exacte, erreur. */

#ifndef FONCTIONS_H
#define FONCTIONS_H

/* Terme source et conditions de bord */
double f(double x, double y, double t, double Lx, double Ly, int cas);
double g(double x, double y, double t, int cas);   /* bords Gamma0 : y = 0 et y = Ly */
double h(double x, double y, double t, int cas);   /* bords Gamma1 : x = 0 et x = Lx */

/* Solution exacte (cas 1 et 2 seulement) */
double solution_exacte(double x, double y, int cas);

/* Erreur max entre u et la solution exacte */
double erreur_max(double *u, int Nx, int Ny, double dx, double dy, int cas);

#endif
