/* main.c
   Partie II - version séquentielle

   Résolution de   du/dt - D Lap(u) = f   sur [0,Lx] x [0,Ly]
                   u = g sur Gamma0,  u = h sur Gamma1
   Différences finies centrées + Euler implicite.
   Système linéaire résolu par BiCGStab, sans stocker la matrice.

   Conventions :
   - Nx x Ny points INTÉRIEURS (les points du bord ne sont pas des inconnues)
   - dx = Lx/(Nx+1),  dy = Ly/(Ny+1)
   - point (i,j), i = 0..Nx-1, j = 0..Ny-1 :  x = (i+1)*dx,  y = (j+1)*dy
   - numéro global k = i + j*Nx, vecteurs de taille N = Nx*Ny
   - Gamma0 = bords y = 0 et y = Ly (côtés "a" de la figure) -> u = g
     Gamma1 = bords x = 0 et x = Lx (côtés "b" de la figure) -> u = h
     (à vérifier sur l'énoncé)

   Fichier de paramètres (valeurs dans cet ordre, sans commentaire dans le fichier) :
      Nx Ny
      Lx Ly
      D
      dt tfinal
      cas          (1, 2 : stationnaires ; 3 : instationnaire)
      eps kmax     (tolérance et nombre max d'itérations de BiCGStab)

   Compilation (depuis partie2_sequentiel) : make
   Exécution   (depuis partie2_sequentiel) : ./sequentiel parametres.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "fonctions.h"
#include "prod_matvec.h"
#include "solveur.h"


/* Lit les paramètres dans le fichier nom_fichier (ordre : voir en tête du fichier).
   Les arguments sont des pointeurs : on écrit dans *Nx, *Ny, ... */
void lecture_parametres(const char *nom_fichier, int *Nx, int *Ny, double *Lx, double *Ly,
                        double *D, double *dt, double *tfinal, int *cas, double *eps, int *kmax)
{
   /* TODO : - fopen en lecture (si fopen renvoie NULL : message d'erreur et exit(1))
             - fscanf de chaque valeur : %d pour les int, %lf pour les double !
             - fclose */
}


/* Écrit la solution dans le fichier nom_fichier : une ligne "x y u" par point intérieur,
   avec une ligne vide après chaque ligne j (format attendu par gnuplot) :
      gnuplot> splot 'solution.dat' w l */
void ecriture_solution(const char *nom_fichier, double *u, int Nx, int Ny, double dx, double dy)
{
   /* TODO : fopen en écriture, double boucle j / i avec fprintf, fclose */
}


int main(int argc, char **argv)
{
   /* Paramètres lus dans le fichier
      (initialisés à 0 : tant que lecture_parametres n'est pas codée, la boucle ne tourne pas) */
   int Nx = 0, Ny = 0, cas = 0, kmax = 0;
   double Lx = 0.0, Ly = 0.0, D = 0.0, dt = 0.0, tfinal = 0.0, eps = 0.0;

   /* Grandeurs calculées */
   int N = 0, nb_pas = 0, n, nb_iter;
   double dx = 0.0, dy = 0.0, t = 0.0;
   double *u, *b;
   clock_t debut, fin;


   /* ---------- 1. Lecture des paramètres ---------- */
   if (argc > 1)
      lecture_parametres(argv[1], &Nx, &Ny, &Lx, &Ly, &D, &dt, &tfinal, &cas, &eps, &kmax);
   else
      lecture_parametres("parametres.txt", &Nx, &Ny, &Lx, &Ly, &D, &dt, &tfinal, &cas, &eps, &kmax);


   /* ---------- 2. Discrétisation ---------- */
   /* TODO : dx = Lx/(Nx+1)
             dy = Ly/(Ny+1)
             N  = Nx*Ny
             nb_pas = nombre de pas de temps pour aller jusqu'à tfinal
                      (arrondir à l'entier le plus proche) */

   printf("Nx = %d, Ny = %d, dt = %g, tfinal = %g, cas = %d\n", Nx, Ny, dt, tfinal, cas);


   /* ---------- 3. Allocations ---------- */
   /* calloc met tout à 0 : condition initiale u^0 = 0 */
   u = (double*) calloc(N, sizeof(double));
   b = (double*) calloc(N, sizeof(double));


   /* ---------- 4. Boucle en temps ---------- */
   /* À chaque pas : on construit b à partir de u^n, puis BiCGStab remplace u^n par u^{n+1}
      (u^n sert de point de départ à BiCGStab).

      Cas stationnaires (1 et 2) : la solution stationnaire est atteinte quand u ne bouge plus.
      Soit on prend tfinal assez grand, soit (mieux) on arrête la boucle quand
      ||u^{n+1} - u^n|| < tolérance (il faut alors garder une copie de u^n). */
   debut = clock();
   for (n = 1; n <= nb_pas; n++) {
      t = n*dt;
      second_membre(b, u, t, Nx, Ny, dx, dy, Lx, Ly, D, dt, cas);
      nb_iter = bicgstab(u, b, Nx, Ny, dt, dx, dy, D, eps, kmax);
      printf("pas %d, t = %f, iterations BiCGStab = %d\n", n, t, nb_iter);
   }
   fin = clock();


   /* ---------- 5. Résultats ---------- */
   printf("Temps de calcul : %f s\n", (double)(fin - debut) / CLOCKS_PER_SEC);

   ecriture_solution("solution.dat", u, Nx, Ny, dx, dy);

   if (cas == 1 || cas == 2)
      printf("Erreur max = %e\n", erreur_max(u, Nx, Ny, dx, dy, cas));


   /* ---------- 6. Libérations ---------- */
   free(u);
   free(b);

   return 0;
}
