#include <math.h>
#include "integral.h"

double integral(int m, double aa[], double bb[], double Aizero[], double y, double a, 
				double b, double (*fx)(int,double *,double *,double *,double,double),
				double e[], int ua, int ub)
{
	double x0,x1,x2,f0,f1,f2,re,ae,b1,x;

	b1=0;
	re=e[1];
	if (ub)
		ae=e[2]*180.0/fabs(b-a);
	else
		ae=e[2]*90.0/fabs(b-a);
	if (ua) {
		e[3]=e[4]=0.0;
		x=x0=a;
		f0=(*fx)(m,aa,bb,Aizero,x,y);
	} else {
		x=x0=a=e[5];
		f0=e[6];
	}
	e[5]=x=x2=b;
	e[6]=f2=(*fx)(m,aa,bb,Aizero,x,y);
	e[4] += integralqad(m,aa,bb,Aizero,y,0,fx,e,&x0,&x1,&x2,&f0,&f1,&f2,re,ae,b1);
	if (!ub) {
		if (a < b) {
			b1=b-1.0;
			x0=1.0;
		} else {
			b1=b+1.0;
			x0 = -1.0;
		}
		f0=e[6];
		e[5]=x2=0.0;
		e[6]=f2=0.0;
		ae=e[2]*90.0;
		e[4] -= integralqad(m,aa,bb,Aizero,y,1,fx,e,&x0,&x1,&x2,&f0,&f1,&f2,re,ae,b1);
	}
	return e[4];
}

double integralqad(int m, double aa[], double bb[], double Aizero[], double y, int transf,
 			double (*fx)(int,double *,double *,double *,double,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double re, double ae, double b1)
{
	/* this function is internally used by INTEGRAL */

	double sum,hmin,x,z;

	hmin=fabs((*x0)-(*x2))*re;
	x=(*x1)=((*x0)+(*x2))*0.5;
	if (transf) {
		z=1.0/x;
		x=z+b1;
		(*f1)=(*fx)(m,aa,bb,Aizero,x,y)*z*z;
	} else
		(*f1)=(*fx)(m,aa,bb,Aizero,x,y);
	sum=0.0;
	integralint(m,aa,bb,Aizero,y,transf,fx,e,x0,x1,x2,f0,f1,f2,&sum,re,ae,b1,hmin);
	return sum/180.0;
}

void integralint(int m, double aa[], double bb[], double Aizero[], double y, int transf,
			double (*fx)(int,double *,double *,double *,double,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double *sum, double re, double ae, double b1,
			double hmin)
{
	/* this function is internally used by INTEGRALQAD of INTEGRAL */

	int anew;
	double x3,x4,f3,f4,h,x,z,v,t;

	x4=(*x2);
	(*x2)=(*x1);
	f4=(*f2);
	(*f2)=(*f1);
	anew=1;
	while (anew) {
		anew=0;
		x=(*x1)=((*x0)+(*x2))*0.5;
		if (transf) {
			z=1.0/x;
			x=z+b1;
			(*f1)=(*fx)(m,aa,bb,Aizero,x,y)*z*z;
		} else
			(*f1)=(*fx)(m,aa,bb,Aizero,x,y);
		x=x3=((*x2)+x4)*0.5;
		if (transf) {
			z=1.0/x;
			x=z+b1;
			f3=(*fx)(m,aa,bb,Aizero,x,y)*z*z;
		} else
			f3=(*fx)(m,aa,bb,Aizero,x,y);
		h=x4-(*x0);
		v=(4.0*((*f1)+f3)+2.0*(*f2)+(*f0)+f4)*15.0;
		t=6.0*(*f2)-4.0*((*f1)+f3)+(*f0)+f4;
		if (fabs(t) < fabs(v)*re+ae)
			(*sum) += (v-t)*h;
		else if (fabs(h) < hmin)
			e[3] += 1.0;
		else {
			integralint(m,aa,bb,Aizero,y,transf,fx,e,x0,x1,x2,f0,f1,f2,sum,
							re,ae,b1,hmin);
			*x2=x3;
			*f2=f3;
			anew=1;
		}
		if (!anew) {
			*x0=x4;
			*f0=f4;
		}
	}
}
