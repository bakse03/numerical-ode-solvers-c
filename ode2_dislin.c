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


int main(void) 
{
    ZB q;
    RYS_KRZYWA_Z_WARTOSCI rys;
    rozwiaz_rownanie_ZB_metoda_strzalow(&q,&rys);
    //rozwiaz_rownanie_ZB_metoda_roznicowa(&q,&rys);
    //rozwiaz_rownanie_ZB_metoda_kolokacji(&q,&rys);
    //rozwiaz_rownanie_ZB_metoda_kwadratur_rozniczkowych(&q,&rys);
    //rozwiaz_rownanie_ZB_metoda_Galerkina(&q,&rys);
    return 0;
}

double ge(double x){return sin(x)+2;}
double ha(double x){return 0.1*x;}
double ef(double x){return 2.0+(sin(x)+2.0)*2.0*x+exp(x)*(x*x+1.0);}
double gee(double x){return cos(x);}
double haa(double x){return 0.0;}
double eff(double x){return 2.0+2.0*x*cos(x);}

void zainicjuj_RYS_KRZYWA_Z_WARTOSCI(RYS_KRZYWA_Z_WARTOSCI *q, ZB *s)
 {
 int i;
 float margines;
 double y_min,y_max;
 q->n=s->N+1;
 q->a=(float)s->a; q->b=(float)s->b;
 q->xl=q->a;  q->xp=q->b;
 y_min=-0.1; y_max=0.1;
 for(i=0;i<=q->n;i++)
  {
    if(s->u[i]<y_min)y_min=s->u[i];
    if(s->u[i]>y_max)y_max=s->u[i];
  }
  margines=(y_max-y_min)/20.0;
 q->yd=(float)((int)(y_min-margines))-1.0; q->yg=(float)((int)(y_max+margines))+1.0;
 q->d_x=1.0; q->d_y=1.0;
 switch(s->metoda)
 {
  case 1: strcpy(q->tytul,"Rozwiazanie: metoda roznicowa");               break;
  case 2: strcpy(q->tytul,"Rozwiazanie: metoda kolokacji");               break;
  case 3: strcpy(q->tytul,"Rozwiazanie: metoda kwadratur rozniczkowych"); break;
  case 4: strcpy(q->tytul,"Rozwiazanie: metoda strzalow");                break;
  case 5: strcpy(q->tytul,"Rozwiazanie: metoda Galerkina");               break;
  default:strcpy(q->tytul,"Rozwiazanie ODE");                             break;
 }
 
 q->ticks_x=1; q->ticks_y=1;
 q->grid_x=1; q->grid_y=1;
 q->osie=1; 
 q->xray=fvector(0,q->n);
 q->yray=fvector(0,q->n);
 
 for(i=0;i<=q->n;i++)
  {
    q->xray[i]=s->x[i];
    q->yray[i]=s->u[i];
  }
 }

void zwolnij_pamiec_RYS_KRZYWA_Z_WARTOSCI(RYS_KRZYWA_Z_WARTOSCI *q)
 {
 free_fvector(q->xray,0,q->n);
 free_fvector(q->yray,0,q->n);
 }

void rysuj_krzywa_dislin_z_wartosci(RYS_KRZYWA_Z_WARTOSCI *q)
 {
  int ic;  
  metafl("CONS");
  setpag ("da4l");
  disini();
  pagera();
  complx();
    
  axspos(450,1800);
  axslen(2200,1200);
 
  name("os x","x");
  name("os y","y");

  labdig(1,"xy");
  ticks(q->ticks_x,"x");
  ticks(q->ticks_y,"y");

  titlin(q->tytul,1);
  ic=intrgb(0.95,0.95,0.95);
  axsbgd(ic);
   
   graf(q->xl,q->xp,q->xl,q->d_x,q->yd,q->yg,q->yd,q->d_y);
  
  setrgb(0.7,0.7,0.7);
  if(q->osie)
   {
    xaxgit();
    yaxgit();   
   }
  else grid(q->grid_x,q->grid_y);

  color("fore");
  height(50);
  title();
  color("red");
  curve(q->xray,q->yray,q->n+1); 
  disfin();
}

void wypisz_rozwiazanie(ZB *q)
 {
 FILE *fp;
 int i;
 switch(q->metoda)
  {
   case 1: fp=fopen("wynik_m_roznicowa.dat","wt"); break;
   case 2: fp=fopen("wynik_m_kolokacji.dat","wt"); break;
   case 3: fp=fopen("wynik_m_kwadratur.dat","wt"); break;
   case 4: fp=fopen("wynik_m_strzalow.dat","wt"); break;
   case 5: fp=fopen("wynik_m_Galerkina.dat","wt"); break;
   default:fp=fopen("wynik_ode.dat","wt"); break;
  } 
 fprintf(fp,"%d\n",q->N+2);
 for(i=0;i<=q->N+1;i++)
  {
   fprintf(fp,"%le %le\n",q->x[i],q->u[i]);      
  }
  fclose(fp); 
 }
