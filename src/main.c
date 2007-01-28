#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gamma.h"
#include "integral.h"
#include "airy.h"
#include "main.h"

#define pi 3.141592653589793
#define c1 1.259921049894873   /* equal to exp(log(2)/3) */
#define c2 1.5874010519681996  /* equal to c1*c1 */
#define sqrt_pi1 1.2533141373155 /* equal to sqrt(pi/2) */
#define sqrt_pi2 5.013256549262  /* equal to 2*sqrt(2*pi)*/

int n = 10;

double sqr(double x)
{
	return x*x;
}

void system_error(char error_message[])
{
	void exit(int);

	printf("%s",error_message);
	exit(1);
}


double *allocate_double_vector(int l, int u)
{
	double *p;

	p = (double *)malloc((unsigned) (u-l+1)*sizeof(double));
	if (!p) system_error("Failure in allocate_double_vector().");
	return p-l;
}


/*void ComputeCoefficients(int n, double a[], double b[])
{
	int i,k,j;
	double sum1,sum2,*c;
	double k1=0.5*sqrt(pi);

	c = allocate_double_vector(0,n);

	a[0]=1.0;
	b[0]=0.0;
	c[0]=1;

	a[1]=7.0/48;
	b[1]=2.0/3;
	c[1]= 3.0/16;

	for (i=2;i<=n;i++) {
		c[i]=-(2*i-3)*(2*i+1)*c[i-1]/(16*i*i*(2*i-1));

		sum1=0;
		sum2=0;

		for (j=0;2*j-1<i;j++) {
			k=2*j-1;
			if (k>0)
				sum2 += a[i-k-1]*beta(3*i-2*k-2,k+1.5)/(gamma(k+1)*pow(2,k+1));
			k=2*j;
			if (k<i)
				sum2 -= a[i-k-1]*beta(3*i-2*k-2,k+1.5)/(gamma(k+1)*pow(2,k+1));
		}
		b[i] =  sum2;

		for (j=0;2*j-1<i;j++) {
			k=2*j-1;
			if (k>0)
				sum1 -= b[i-k]*beta(3*i-2*k-0.5,k+1.5)/(pi*gamma(k+1)*pow(2,k));
			k=2*j;
			if (k<i)
				sum1 += b[i-k]*beta(3*i-2*k-0.5,k+1.5)/(pi*gamma(k+1)*pow(2,k));
		}

		a[i] = c[i] - sum1;
	}

	free_double_vector(c,0,n);
}*/


/*double *a,*b,*zai,*vai;
a = allocate_double_vector(0,n);
b = allocate_double_vector(0,n);
zai = allocate_double_vector(0,n);
vai = allocate_double_vector(0,n);
ComputeCoefficients(n,a,b);
airyzeros(n,0,zai,vai);*/


/*void free_double_vector(double *v, int l, int u)
{
	free((char*) (v+l));
} */


/*double moment_function(int i, double a[], double b[], double zai[], double vai[], double x)
{
	double s;
	s = pow(x,i)*f_Z(20,a,b,zai,vai,x);
	return 2*s;
} */

double tolx(double x)  /* INUTILE ? */
{
	double s=1.0e-10;

	return s;
}


double airyzeros(int n, int d, double zai[], double vai[])
{
	int a,found,i;
	double c,e,r,zaj,zak,vaj,daj,kaj,zz;

	a=((d == 0) || (d == 2));
	r = (d == 0 || d == 3) ? -1.17809724509617 : -3.53429173528852;
	airy(0.0,&zaj,&vaj,&daj,&kaj,&zz,1);
	for (i=1; i<=n; i++) {
		r += 4.71238898038469;
		zz=r*r;
		zaj = (i == 1 && d == 1) ? -1.01879297 :
				((i == 1 && d == 2) ? -1.17371322 :
				pow(r,0.666666666666667)*
				(a ? -(1.0+(5.0/48.0-(5.0/36.0-(77125.0/82944.0-
				(108056875.0/6967296.0-(162375596875.0/334430208.0)/
				zz)/zz)/zz)/zz)/zz) : -(1.0-(7.0/48.0-(35.0/288.0-
				(181223.0/207360.0-(18683371.0/1244160.0-
				(91145884361.0/191102976.0)/zz)/zz)/zz)/zz)/zz)));
		if (d <= 1.0)
			airy(zaj,&vaj,&daj,&c,&e,&zz,0);
		else
			airy(zaj,&c,&e,&vaj,&daj,&zz,0);
		found=(fabs(a ? vaj : daj) < 1.0e-12);
		while (!found) {
			if (a) {
				kaj=vaj/daj;
				zak=zaj-kaj*(1.0+zaj*kaj*kaj);
			} else {
				kaj=daj/(zaj*vaj);
				zak=zaj-kaj*(1.0+kaj*(kaj*zaj+1.0/zaj));
			}
			if (d <= 1)
				airy(zak,&vaj,&daj,&c,&e,&zz,0);
			else
				airy(zak,&c,&e,&vaj,&daj,&zz,0);
			found=(fabs(zak-zaj) < 1.0e-14*fabs(zak) ||
					fabs(a ? vaj : daj) < 1.0e-12);
			zaj=zak;
		}
		vai[i]=(a ? daj : vaj);
		zai[i]=zaj;
	}
	return zai[n];
}


void airyzeros_R(double zai[], double vai[])
{
	int d = 0;
  int a,found,i;
	double c,e,r,zaj,zak,vaj,daj,kaj,zz;

	a=((d == 0) || (d == 2));
	r = (d == 0 || d == 3) ? -1.17809724509617 : -3.53429173528852;
	airy(0.0,&zaj,&vaj,&daj,&kaj,&zz,1);
	for (i=1; i<=n; i++) {
		r += 4.71238898038469;
		zz=r*r;
		zaj = (i == 1 && d == 1) ? -1.01879297 :
				((i == 1 && d == 2) ? -1.17371322 :
				pow(r,0.666666666666667)*
				(a ? -(1.0+(5.0/48.0-(5.0/36.0-(77125.0/82944.0-
				(108056875.0/6967296.0-(162375596875.0/334430208.0)/
				zz)/zz)/zz)/zz)/zz) : -(1.0-(7.0/48.0-(35.0/288.0-
				(181223.0/207360.0-(18683371.0/1244160.0-
				(91145884361.0/191102976.0)/zz)/zz)/zz)/zz)/zz)));
		if (d <= 1.0)
			airy(zaj,&vaj,&daj,&c,&e,&zz,0);
		else
			airy(zaj,&c,&e,&vaj,&daj,&zz,0);
		found=(fabs(a ? vaj : daj) < 1.0e-12);
		while (!found) {
			if (a) {
				kaj=vaj/daj;
				zak=zaj-kaj*(1.0+zaj*kaj*kaj);
			} else {
				kaj=daj/(zaj*vaj);
				zak=zaj-kaj*(1.0+kaj*(kaj*zaj+1.0/zaj));
			}
			if (d <= 1)
				airy(zak,&vaj,&daj,&c,&e,&zz,0);
			else
				airy(zak,&c,&e,&vaj,&daj,&zz,0);
			found=(fabs(zak-zaj) < 1.0e-14*fabs(zak) ||
					fabs(a ? vaj : daj) < 1.0e-12);
			zaj=zak;
		}
		vai[i]=(a ? daj : vaj);
		zai[i]=zaj;
	}
}


double qfunction(int n, double a[], double b[], double zai[], double vai[], double x, double u)
{
	double s;

	s = F_Z(n,a,b,zai,vai,x)-u;

	return s;
}


/*double moment(int i, double a[], double b[], double zai[], double vai[])
{
	double s;
	double e[7];
	
	e[1]=e[2]=1.0e-10;
	e[5]=e[6]=0;

	s=integral1(i,a,b,zai,vai,0,1,moment_function,e,1,0);

	return s;
}*/


double F_Z(int m, double a[], double b[], double zai[], double vai[], double x)
{
	double s;
	double e[7];
	
	e[1]=e[2]=1.0e-10;
	e[5]=e[6]=0;
	
	if (x>=0) {
		if (x<=1)
			s = 0.5+integral1(m,a,b,zai,vai,0,x,f_Z,e,1,1);
		else
			s = 1-integral1(m,a,b,zai,vai,x,x+1,f_Z,e,1,0);
	} else {
		x=-x;
		if (x<=1)
			s = 0.5+integral1(m,a,b,zai,vai,0,x,f_Z,e,1,1);
		else
			s = 1-integral1(m,a,b,zai,vai,x,x+1,f_Z,e,1,0);
	}
	
	return s;
}


/*double F_Z_R(double *x, double *s)*/
void F_Z_R(double *x, double *s, double a[], double b[], double zai[], double vai[])
{
	int nn = 10;
  
  double e[7];
	
	e[1]=e[2]=1.0e-10;
	e[5]=e[6]=0;

	if (*x >= 0) {
		if (*x <= 1)
			*s = 0.5 + integral1(nn,a,b,zai,vai,0,*x,f_Z,e,1,1);
		else
			*s = 1 - integral1(nn,a,b,zai,vai,*x,*x+1,f_Z,e,1,0);
	} else {
		*x = - *x;
		if (*x <= 1)
			*s = 0.5 + integral1(nn,a,b,zai,vai,0,*x,f_Z,e,1,1);
		else
			*s = 1 - integral1(nn,a,b,zai,vai,*x,*x+1,f_Z,e,1,0);
	}
}



double f_Z(int m, double a[], double b[], double zai[], double vai[], double x)
{
	double s,x1,x2,x3,y,y1,y2,y3;
	double e[7];
	
	e[1]=e[2]=1.0e-10;
	e[5]=e[6]=0;

	if (x<0) x =-x;

	x1= integral(m,a,b,zai,x,0,1,p1,e,1,0)/sqrt(2*pi);
	x2= integral2(x,0,1,p2,e,1,0)*2*sqrt(2/pi);
	x3= 2*x-x1+x2;

	y = -x;
	if (y>=-1) {
		y1= integral(m,a,b,zai,y,0,1,p1,e,1,0)/sqrt(2*pi);
		y2= integral2(y,0,1,p2,e,1,0)*2*sqrt(2/pi);
		y3= 2*y-y1+y2;
	}
	else y3= p3(m,zai,vai,y);

	s = x3*y3/2;

	return s;
}


/*void f_Z_R(double *x, double *s, double a[], double b[])*/
void f_Z_R(double *x, double *s, double a[], double b[], double zai[], double vai[])
/*void f_Z_R(double *x, double *s, double *a, double *b)*/
{

	int nn = 10;
  
  double x1,x2,x3,y,y1,y2,y3;
	double e[7];
	
	e[1]=e[2]=1.0e-10;
	e[5]=e[6]=0;

  /*double *a, *b, *zai,*vai;
  
  a = allocate_double_vector(0,nn);
  b = allocate_double_vector(0,nn);
  zai = allocate_double_vector(0,nn);
  vai = allocate_double_vector(0,nn);
  ComputeCoefficients(nn,a,b);
  airyzeros(nn,0,zai,vai); */

	if (*x<0) *x = -*x;

	x1 = integral(nn,a,b,zai,*x,0,1,p1,e,1,0)/sqrt(2*pi);
	x2 = integral2(*x,0,1,p2,e,1,0)*2*sqrt(2/pi);
	x3 = 2*(*x)-x1+x2;

	y = -*x;
	if (y >= -1) {
		y1 = integral(nn,a,b,zai,y,0,1,p1,e,1,0)/sqrt(2*pi);
		y2 = integral2(y,0,1,p2,e,1,0)*2*sqrt(2/pi);
		y3 = 2*y-y1+y2;
	}
	else y3 = p3(nn,zai,vai,y);

	*s = x3*y3/2;
}


/*double quantile_function_R(double *u, double *s)*/
void quantile_function_R(double *u, double *s, double a[], double b[], double zai[], double vai[])
{

	int nn = 10;
  int result;
	double v = 0.5;
	double below,above;
	below=0.0;
	above=5.0;

	if ((*u <= 0) || (*u >= 1))
		system_error("Choose a value strictly between 0 and 1");

	if (*u < 0.5) v = 1-0.5;
	if (*u >= 0.5)
		result = zeroinder(nn, a, b, zai, vai, *u, &below, &above, qfunction, f_Z, tolx);
	else
		result = zeroinder(nn, a, b, zai, vai, v, &below, &above, qfunction, f_Z, tolx);

	*s = (below+above)/2;

	if (*u < 0.5)
		*s = -*s;
}


/*double px(int m, double aa[], double bb[], double Aizero[], double x)
{
	int i;
	double sum,yy,z;

	sum=0;

	if ((0<x) && (x<=1)) {
		sum=-sqrt_pi1;
		for (i=1;i<=m;i++)
			sum+=-sqrt_pi1*aa[i]*exp(3*i*log(x))+bb[i]*exp(3*(i-0.5)*log(x));
	} else {
		if (x>1) {
			sum=0;
			yy=-exp(-1.5*log(x));
			z = sqrt_pi2*exp(-x*sqr(x)/6);
			for (i=1;i<=m;i++)
				sum+= exp(c1*Aizero[i]*x);
			sum = yy+z*sum;
		}
	}
	return sum;
}*/


double p1(int m, double aa[], double bb[], double Aizero[], double x, double y)
{
	int i;
	double sum,yy,z;

	sum=0;

	if ((0<x) && (x<=1)) {
		sum = -sqrt_pi1;
		for (i=1;i<=m;i++)
			sum += -sqrt_pi1*aa[i]*exp(3*i*log(x)) + bb[i]*exp(3*(i-0.5)*log(x));
	} else {
		if (x>1) {
			sum=0;
			yy = -exp(-1.5*log(x));
			z = sqrt_pi2*exp(-x*sqr(x)/6);
			for (i=1;i<=m;i++)
				sum += exp(c1*Aizero[i]*x);
			sum = yy+z*sum;
		}
	}
	
	/*sum = px(m,aa,bb,Aizero,x);*/
  sum = sum*exp(-0.5*x*sqr(2*y+x));

	return sum;
}

double p2(double x, double y)
{
	double s,z0,z1,z2;

	z2 = sqr(x);
	z1 = 2*y+z2; /* modif */
	z0 = 0.5*sqr(z1); /* modif */

	s = (z1*z2 + z0)*exp(-z0*z2);

	return s;
}

double p3(int m, double zai[], double vai[], double x)
{
	int i;
	double s;

	s=0;

	for (i=1;i<=m;i++)
		s += exp(-c1*zai[i]*x)/vai[i];
	s = c2*exp(2*sqr(x)*x/3)*s;

	return s;
}
