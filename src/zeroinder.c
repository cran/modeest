#include <math.h>
#include "main.h"

#define FLT_MIN 1.0e-20

int zeroinder(int n, double aa[], double bb[], double zai[], double vai[], double u, double *x, double *y,
				 double (*fx)(int,double *,double *,double *,double *,double,double),
					double (*dfx)(int,double *,double *,double *,double *,double),
					 double (*tolx)(double))
{
	int ext,extrapolate;
	double b,fb,dfb,a,fa,dfa,c,fc,dfc,d,w,mb,tol,m,p,q;

	b = *x;
	fb=(*fx)(n,aa,bb,zai,vai,*x,u);
	dfb=(*dfx)(n,aa,bb,zai,vai,*x);
	a = *x = *y;
	fa=(*fx)(n,aa,bb,zai,vai,*x,u);
	dfa=(*dfx)(n,aa,bb,zai,vai,*x);
	c=a;
	fc=fa;
	dfc=dfa;
	ext=0;
	extrapolate=1;
	while (extrapolate) {
		if (fabs(fc) < fabs(fb)) {
			a=b;
			fa=fb;
			dfa=dfb;
			b = *x =c;
			fb=fc;
			dfb=dfc;
			c=a;
			fc=fa;
			dfc=dfa;
		}
		tol=(*tolx)(*x);
		m=(c+b)*0.5;
		mb=m-b;
		if (fabs(mb) > tol) {
			if (ext > 2)
				w=mb;
			else {
				if (mb == 0.0)
					tol=0.0;
				else
					if (mb < 0.0) tol = -tol;
				d = (ext == 2) ? dfa : (fb-fa)/(b-a);
				p=fb*d*(b-a);
				q=fa*dfb-fb*d;
				if (p < 0.0) {
					p = -p;
					q = -q;
				}
				w=(p<FLT_MIN || p<=q*tol) ? tol : ((p<mb*q) ? p/q : mb);
			}
			a=b;
			fa=fb;
			dfa=dfb;
			*x = b += w;
			fb=(*fx)(n,aa,bb,zai,vai,*x,u);
			dfb=(*dfx)(n,aa,bb,zai,vai,*x);
			if ((fc >= 0.0) ? (fb >= 0.0) : (fb <= 0.0)) {
				c=a;
				fc=fa;
				dfc=dfa;
				ext=0;
			} else
				ext = (w == mb) ? 0 : ext+1;
		} else
			break;
	}
	*y = c;
	return ((fc >= 0.0) ? (fb <= 0.0) : (fb >= 0.0));
}
