
double integral(int m, double aa[], double bb[], double Aizero[], double y, double a, 
				double b, double (*fx)(int,double *,double *,double *,double,double),
				double e[], int ua, int ub);

double integralqad(int m, double aa[], double bb[], double Aizero[], double y, int transf,
 			double (*fx)(int,double *,double *,double *,double,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double re, double ae, double b1);

void integralint(int m, double aa[], double bb[], double Aizero[], double y, int transf,
			double (*fx)(int,double *,double *,double *,double,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double *sum, double re, double ae, double b1,
			double hmin);
			
double integral1(int m, double aa[], double bb[], double Aizero[], double vai[], double a, 
				double b, double (*fx)(int,double *,double *,double *,double *,double),
				double e[], int ua, int ub);

double integralqad1(int m, double aa[], double bb[], double Aizero[], double vai[], int transf,
 			double (*fx)(int,double *,double *,double *,double *,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double re, double ae, double b1);

void integralint1(int m, double aa[], double bb[], double Aizero[], double vai[], int transf,
			double (*fx)(int,double *,double *,double *,double *,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double *sum, double re, double ae, double b1,
			double hmin);
			
double integral2(double y, double a, double b, double (*fx)(double,double), double e[],
					int ua, int ub);

double integralqad2(double y, int transf, double (*fx)(double,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double re, double ae, double b1);

void integralint2(double y,int transf, double (*fx)(double,double), double e[],
			double *x0, double *x1, double *x2, double *f0, double *f1,
			double *f2, double *sum, double re, double ae, double b1,
			double hmin);


