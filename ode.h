typedef struct zb
 {
  double a,b,ua,ub;                               // dane zadania
  double (*p)(double),(*q)(double),(*f)(double);  // dane zadania
  int N,K;                            // liczba wezlow - 2, liczba wielomianow Peano
  double *x,*u;                       // wektor wezlow i wartosci
  double *P;
  double *C0,*C1,*C2;	                // dla metody roznicowej
  double **dq_1,**dq_2,**M;           // dla metody kwadratur
  double h,xl,xp,s_gwiazdka;          // dla metody strzalow i Galerkina
  int metoda;
 }ZB;

typedef struct rys_krzywa_z_wartosci
  {
  int n;
  float a,b;
  float xl,yd,xp,yg,d_x,d_y;
  char tytul[80];
  int ticks_x,ticks_y;
  int grid_x,grid_y;
  int osie;
  float *xray, *yray;
  }RYS_KRZYWA_Z_WARTOSCI;  

int main(void);
// DANE DO ZADANIA - WSPOLCZYNNIKI ROWNANIA I PRAWA STRONA
double ge(double x);
double ha(double x);
double ef(double x);
double gee(double x);
double haa(double x);
double eff(double x);
// METODA ROZNICOWA
void rozwiaz_rownanie_ZB_metoda_roznicowa(ZB *q, RYS_KRZYWA_Z_WARTOSCI *rys);
void zainicjuj_ZB_metoda_roznicowa(ZB *q);
void rozwiaz_trojdiag_Thomas(ZB *q);
void odpusc_pamiec_ZB(ZB *q);
// METODA KOLOKACJI
double chi(int j, double x, ZB *q);
double chi_prim(int j, double x, ZB *q);
double chi_bis(int j, double x, ZB *q);
void rozwiaz_rownanie_ZB_metoda_kolokacji(ZB *q, RYS_KRZYWA_Z_WARTOSCI *rys);
void zainicjuj_ZB_metoda_kolokacji(ZB *q);
void odpusc_pamiec_ZB_kolokacja(ZB *q);
// METODA KWADRATUR ROZNICZKOWYCH
double M_prim(int k, ZB *s);
void rozwiaz_rownanie_ZB_metoda_kwadratur_rozniczkowych(ZB *q, RYS_KRZYWA_Z_WARTOSCI *rys);
void zainicjuj_ZB_metoda_kwadratur_rozniczkowych(ZB *q);
void wyznacz_wspolczynniki_kwadratury(ZB *q);
void odpusc_pamiec_ZB_kwadratury(ZB *q);
// METODA STRZALOW
void rozwiaz_rownanie_ZB_metoda_strzalow(ZB *q, RYS_KRZYWA_Z_WARTOSCI *rys);
void zainicjuj_ZB_metoda_strzalow(ZB *q);
double cel(double s, ZB *q);
void uchwyc_konce(ZB *q);
double bisekcja(ZB *q);
void Euler(ZB *q);
void odpusc_pamiec_ZB_strzaly(ZB *q);
//METODA GALERKINA
int iabs(int k);
double C_2(int i, int j, ZB *q);
double Cp_2(int i, int j, ZB *q);
double C_3(int m, int i, int j, ZB *q);
double Chp_3(int m, int i, int j, ZB *q);
void rozwiaz_rownanie_ZB_metoda_Galerkina(ZB *q, RYS_KRZYWA_Z_WARTOSCI *rys);
void zainicjuj_ZB_metoda_Galerkina(ZB *q);
// GRAFIKA
void wypisz_rozwiazanie_kolokacja(ZB *q);
void zainicjuj_RYS_KRZYWA_Z_WARTOSCI(RYS_KRZYWA_Z_WARTOSCI *q, ZB *s);
void zwolnij_pamiec_RYS_KRZYWA_Z_WARTOSCI(RYS_KRZYWA_Z_WARTOSCI *q);
void rysuj_krzywa_dislin_z_wartosci(RYS_KRZYWA_Z_WARTOSCI *q);
// REZULTAT OBLICZEN
void wypisz_rozwiazanie(ZB *q);

