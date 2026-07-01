#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "disnrutil.h"
#include "gs_lin.h"

static double *sa;
static long *ija;

/* NUMERICAL RECIPES */
void linbcg(double *ref_sa, long *ref_ija,long nu, double b[], double x[], int itol, double tol,
	int itmax, int *iter, double *err)
{
	long j,ka;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;
	int wew_iter;
	double wew_err;
	int nuu;
	char bufik[70];
   FILE *fn;

   /* */
    sa=ref_sa;
    ija=ref_ija;
   /* */

	nuu=(int)nu;

	p=dvector(1,nuu);
	pp=dvector(1,nuu);
	r=dvector(1,nuu);
	rr=dvector(1,nuu);
	z=dvector(1,nuu);
	zz=dvector(1,nuu);
	wew_iter=0;
	atimes(nu,x,r,0);
	for (j=1;j<=nu;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	znrm=1.0;
	if (itol == 1) bnrm=snrm(nu,b,itol);
	else if (itol == 2) {
		asolve(nu,b,z,0);
		bnrm=snrm(nu,z,itol);
	}
	else if (itol == 3 || itol == 4) {
		asolve(nu,b,z,0);
		bnrm=snrm(nu,z,itol);
		asolve(nu,r,z,0);
		znrm=snrm(nu,z,itol);
	} else nrerror("illegal itol in linbcg");
	asolve(nu,r,z,0);
	while (wew_iter <= itmax) {
		++(wew_iter);
		zm1nrm=znrm;
		asolve(nu,rr,zz,1);
		for (bknum=0.0,j=1;j<=nu;j++) bknum += z[j]*rr[j];
		if (wew_iter == 1) {
			for (j=1;j<=nu;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=nu;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(nu,p,z,0);
		for (akden=0.0,j=1;j<=nu;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(nu,pp,zz,1);
		for (j=1;j<=nu;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(nu,r,z,0);
		if (itol == 1 || itol == 2) {
			znrm=1.0;
			wew_err=(bnrm!=0.0)?snrm(nu,r,itol)/bnrm:snrm(nu,r,itol)/MIN_VAL;
		} else if (itol == 3 || itol == 4) {
			znrm=snrm(nu,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(nu,p,itol);
				wew_err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				wew_err=(bnrm!=0.0)?znrm/bnrm:znrm/MIN_VAL;
				continue;
			}
			xnrm=snrm(nu,x,itol);
			if (wew_err <= 0.5*xnrm) wew_err /= xnrm;
			else {
				wew_err=(bnrm!=0.0)?znrm/bnrm:znrm/MIN_VAL;
				continue;
			}
		}
     #if !defined(GRAFIKA)
      /*
      sprintf(bufik,"Iteracja: %5d     bˆ¥d = %12.6le\n",wew_iter,wew_err);
		t_put(35,24,bufik,14);
      */
      printf("Iteracja: %5d     bˆ¥d = %12.6le\n",wew_iter,wew_err);
     #endif
     //printf("Iteracja: %5d     bˆ¥d = %12.6le\n",wew_iter,wew_err); getchar();
	if (wew_err <= tol) break;
	}
	*iter=wew_iter;
	*err=wew_err;

	free_dvector(zz,1,nuu);
 free_dvector(z,1,nuu);
	free_dvector(rr,1,nuu);
	free_dvector(r,1,nuu);
	free_dvector(pp,1,nuu);
	free_dvector(p,1,nuu);
}

static void asolve(long nu, double b[], double x[], int itrnsp)
{
	long i;
	for(i=1;i<=nu;i++) x[i]=(sa[i] != 0.0 ? b[i]/sa[i] : b[i]);
}

static void atimes(long nu, double x[], double r[], int itrnsp)
{
	if (itrnsp) dsprstx(sa,ija,x,r,nu);
	else dsprsax(sa,ija,x,r,nu);
}

static void dsprsax(double sa[], long ija[], double x[], double b[], long nu)
{
	long i,k;
	if (ija[1] != nu+2) nrerror("dsprsax: mismatched vector and matrix");
	for (i=1;i<=nu;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++) b[i] += sa[k]*x[ija[k]];
	}
}

static void dsprstx(double sa[], long ija[], double x[], double b[], long nu)
{
	long i,j,k;
	if (ija[1] != nu+2) nrerror("mismatched vector and matrix in dsprstx");
	for (i=1;i<=nu;i++) b[i]=sa[i]*x[i];
	for (i=1;i<=nu;i++) {
		for (k=ija[i];k<=ija[i+1]-1;k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}

double snrm(long nu, double sx[], int itol)
{
	long i,isamax;
	double ans;

	if (itol <= 3) {
		ans = 0.0;
		for (i=1;i<=nu;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=1;
		for (i=1;i<=nu;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}
/* Koniec NUMERICAL RECIPES */
