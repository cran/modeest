int n;
void system_error(char error_message[]);
double *allocate_double_vector(int l, int u);
/*void free_double_vector(double *v, int l, int u);*/
double sqr(double x);
/*double px(int m, double aa[], double bb[], double Aizero[], double x);*/
double p1(int m, double aa[], double bb[], double Aizero[], double x, double y);
double p2(double x, double y);
double p3(int m, double zai[], double vai[], double x);
double f_Z(int m, double a[], double b[], double zai[], double vai[], double x);
void f_Z_R(double *x, double *s, double a[], double b[], double zai[], double vai[]);
/*void f_Z_R(double *x, double *s, double a[], double b[]);*/
/*void f_Z_R(double *x, double *s, double *a, double *b);*/
double F_Z(int m, double a[], double b[], double zai[], double vai[], double x);
void F_Z_R(double *x, double *s, double a[], double b[], double zai[], double vai[]);
/*void ComputeCoefficients(double a[], double b[]);*/
/*double moment_function(int i, double a[], double b[], double zai[], double vai[], double x);
double moment(int i, double a[], double b[], double zai[], double vai[]);*/
int zeroinder(int n, double aa[], double bb[], double zai[], double vai[], double u, double *x, double *y,
				 double (*fx)(int,double *,double *,double *,double *,double,double),
					double (*dfx)(int,double *,double *,double *,double *,double),
					 double (*tolx)(double));
double 	qfunction(int n, double a[], double b[], double zai[], double vai[], double x, double u);
void quantile_function_R(double *u, double *s, double a[], double b[], double zai[], double vai[]);
double tolx(double x);
double airyzeros(int n, int d, double zai[], double vai[]);
void airyzeros_R(double zai[], double vai[]);

/*void ComputeCoefficients(int n, double a[], double b[]);*/
/*double quantile_function(int n, double a[], double b[], double zai[], double vai[], double u);*/
