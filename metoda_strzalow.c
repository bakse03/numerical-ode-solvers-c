#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dislin.h>
#include "disnrutil.h"
#include "gs_lin.h"
#include "ode.h"

// Rownanie u''+p(x)u'+q(x)u=f(x) x \in (a,b)
//          u(a)=ua,  u(b)=ub



// METODA ROZNICOWA
void rozwiaz_rownanie_ZB_metoda_roznicowa(ZB *q, RYS_KRZYWA_Z_WARTOSCI *rys)
 {
    zainicjuj_ZB_metoda_roznicowa(q);
    rozwiaz_trojdiag_Thomas(q);
    wypisz_rozwiazanie(q);
    zainicjuj_RYS_KRZYWA_Z_WARTOSCI(rys,q);
    rysuj_krzywa_dislin_z_wartosci(rys);
    zwolnij_pamiec_RYS_KRZYWA_Z_WARTOSCI(rys);
    odpusc_pamiec_ZB(q);
 }

void zainicjuj_ZB_metoda_roznicowa(ZB *q)
 {
  int i,j;
  double krokx;   
  q->metoda=1;
  q->a=-1.0; q->b=2.0;
  q->ua=2.0; q->ub=5.0;
  q->p=gee; q->q=ha; q->f=ef;
  q->N=50;
  q->x=dvector(0,q->N+1);
  q->u=dvector(0,q->N+1);
  q->P=dvector(0,q->N+1);
  q->C0=dvector(0,q->N+1);
  q->C1=dvector(0,q->N+1);
  q->C2=dvector(0,q->N+1);
  krokx=(q->b-q->a)/(q->N+1.0);
  for(i=0;i<=q->N+1;i++)
   {
    q->x[i]=q->a+i*krokx;
   }
  q->P[0]=q->ua; q->P[q->N+1]=q->ub;
  for(i=1;i<=q->N;i++)
   {
    q->P[i]=(*q->f)(q->x[i]);
   }   
  q->C0[0]=0.0; q->C1[0]=1.0; q->C2[0]=0.0;
  for(i=1;i<=q->N;i++)
   {
       q->C0[i]=(1.0/krokx-0.5*(*q->p)(q->x[i]))/krokx;       
    q->C1[i]=(*q->q)(q->x[i])-2.0/(krokx*krokx);
       q->C2[i]=(1.0/krokx+0.5*(*q->p)(q->x[i]))/krokx;
   }
  q->C0[q->N+1]=0.0; q->C1[q->N+1]=1.0; q->C2[q->N+1]=0.0;   
 }

void odpusc_pamiec_ZB(ZB *q)
 {
  free_dvector(q->x,0,q->N+1);
  free_dvector(q->u,0,q->N+1);
  free_dvector(q->P,0,q->N+1);
  free_dvector(q->C0,0,q->N+1);
  free_dvector(q->C1,0,q->N+1);
  free_dvector(q->C2,0,q->N+1);
 }
 
void rozwiaz_trojdiag_Thomas(ZB *q)
 {
  int i,j,k;
  double *beta, *gamma;
  beta=dvector(0,q->N+1);
  gamma=dvector(0,q->N+1);
  beta[0]=q->C1[0]; gamma[0]=q->P[0]/q->C1[0];
  for(i=1;i<=q->N+1;i++)
   {
    beta[i]=q->C1[i]-q->C0[i]*q->C2[i-1]/beta[i-1];
    gamma[i]=(q->P[i]-q->C0[i]*gamma[i-1])/beta[i];
   }
  q->u[q->N+1]=gamma[q->N+1];
  for(i=q->N;i>=0;q->u[i]=gamma[i]-q->C2[i]*q->u[i+1]/beta[i],i--);
  free_dvector(beta,0,q->N+1);
  free_dvector(gamma,0,q->N+1);
 }



// METODA STRZALOW

void rozwiaz_rownanie_ZB_metoda_strzalow(ZB *q, RYS_KRZYWA_Z_WARTOSCI *rys)
{
    double sgwiazdka;
    zainicjuj_ZB_metoda_strzalow(q);
    sgwiazdka = bisekcja(q);
    q->s_gwiazdka = sgwiazdka;
    Euler(q);
    wypisz_rozwiazanie(q);
    zainicjuj_RYS_KRZYWA_Z_WARTOSCI(rys,q);
    rysuj_krzywa_dislin_z_wartosci(rys);
    zwolnij_pamiec_RYS_KRZYWA_Z_WARTOSCI(rys);
    odpusc_pamiec_ZB_strzaly(q);
}

void zainicjuj_ZB_metoda_strzalow(ZB *q)
{
    int i;
    double krokx;
    q->metoda = 4;
    q->a = -1.0; q->b = 2.0;
    q->ua = 2.0; q->ub = 5.0;
    q->p = gee; q->q = ha; q->f = ef;
    q->N = 50;
    q->x = dvector(0, q->N+1);
    q->u = dvector(0, q->N+1);
    krokx = (q->b - q->a) / (q->N+1.0);
    for (i = 0; i <= q->N+1; i++) q->x[i] = q->a + i * krokx;
    q->h = krokx;
    q->u[0] = q->ua;
}

double cel(double s, ZB *q)
{
    int i;
    double x = q->a;
    double h = q->h;

    double y1 = q->ua;  // u(a)
    double y2 = s;      // u'(a)

    int nsteps = (int)((q->b - q->a) / h + 0.5);

    for (i = 0; i < nsteps; i++) {
        double dy2 = (*q->f)(x) - (*q->p)(x) * y2 - (*q->q)(x) * y1;

        y1 += h * y2;
        y2 += h * dy2;
        x  += h;
    }

    return y1 - q->ub;
}

void uchwyc_konce(ZB *q)
{
    double s0 = 0.0, ds = 1.0;
    double f1, f2;
    int iter = 0;

    do {
        q->xl = s0 - ds;
        q->xp = s0 + ds;
        f1 = cel(q->xl, q);
        f2 = cel(q->xp, q);
        ds *= 2.0;
        iter++;
    } while (f1 * f2 > 0.0 && iter < 50);
}

double bisekcja(ZB *q)
{
    double a, b, c, fa, fc;
    int it;
    double tol = 1e-8;

    uchwyc_konce(q);
    a = q->xl;
    b = q->xp;
    fa = cel(a, q);

    for (it = 0; it < 60; it++) {
        c = 0.5 * (a + b);
        fc = cel(c, q);

        if (fabs(fc) < tol || 0.5 * (b - a) < tol)
            return c;

        if (fa * fc <= 0.0)
            b = c;
        else {
            a = c;
            fa = fc;
        }
    }
    return 0.5 * (a + b);
}

void Euler(ZB *q)
{
    int i;
    double x = q->a;
    double h = q->h;

    double y1 = q->ua;          // u
    double y2 = q->s_gwiazdka;  // u'

    q->u[0] = y1;
    
    int nsteps = (int)((q->b - q->a) / h + 0.5);

    for (i = 1; i <= nsteps; i++) {
        double dy2 = (*q->f)(x) - (*q->p)(x)*y2 - (*q->q)(x)*y1;

        y1 += h * y2;
        y2 += h * dy2;
        x  += h;

        q->u[i] = y1;
    }
}

void odpusc_pamiec_ZB_strzaly(ZB *q)
{
    if (q->x) free_dvector(q->x,0,q->N+1);
    if (q->u) free_dvector(q->u,0,q->N+1);
}
