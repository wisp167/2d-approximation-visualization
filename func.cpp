#include "window.h"


void err_output(char* str){
	printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
	str, 6, (double) -1, (double) -1, (double) -1, (double) -1, (double) -1, (double) -1, -1, (double) -1, -1, -1, -1, -1);
}


inline void thread_rows(int n, int p, int kk, int &i1, int &i2){i1 = n*kk; i1/=p; i2=n*(kk+1); i2/=p;}


double get_cpu_time(){
	struct rusage buf;
	getrusage(RUSAGE_THREAD, &buf);
	return buf.ru_utime.tv_sec+buf.ru_utime.tv_usec/1e6;
}

double get_full_time(){
	struct timeval t;
	gettimeofday(&t, 0);
	return t.tv_sec+t.tv_usec/1e6;
}

int read_input(int* temp, double* ttemp, char** argv){
	string ss;
	for(int i = 0; i < 4; ++i){
		if(sscanf(argv[i+1], "%lf", &ttemp[i]) != 1){
			return -1;
		}
	}
	for(int i = 4; i < 6; ++i){
		ss = argv[i+1];
		if((sscanf(argv[i+1], "%d", &temp[i-4]) != 1) || (to_string(temp[i-4]) != ss)){
			return -1;
		}
	}
	ss = argv[9];
	if((sscanf(argv[9], "%d", &temp[2]) != 1) || (to_string(temp[2]) != ss)){
		return -1;
	}
	if(sscanf(argv[10], "%lf", &ttemp[4]) != 1){
		return -1;
	}
	for(int i = 10; i < 12; ++i){
		ss = argv[i+1];
		if((sscanf(argv[i+1], "%d", &temp[i-7]) != 1) || (to_string(temp[i-7]) != ss)){
			return -1;
		}
	}
	return 0;
}


double f(int k, double x, double y){
	switch(k){
		case 0: return 1;
		case 1: return x;
		case 2: return y;
		case 3: return x+y;
		case 4: return sqrt(x*x+y*y);
		case 5: return x*x+y*y;
		case 6: return exp(x*x-y*y);
	}
	return (double)1/(25*(x*x+y*y)+1);
}
/*
inline void init_(int nx, int ny, double* A, size_t* I, double *B, double* x, double *r, double* u, double* v, double *temp1){
	size_t N = (size_t) (nx+1)*(ny+1);
	A = (double*)malloc (sizeof (double) * (get_len_msr(nx,ny)+1));
	I = (size_t*)malloc(sizeof(size_t)*(get_len_msr(nx,ny)+1));
	B = (double*)malloc(sizeof(double) * N);
	x = (double*)malloc(sizeof(double) * N);
	r = (double*)malloc(sizeof(double) * N);
	u = (double*)malloc(sizeof(double) * N);
	v = (double*)malloc(sizeof(double) * N);
	temp1 = (double*) malloc(sizeof(double)*N);
}
*/
void free_vars(double* A, size_t* I, double *B, double* x, double *r, double* u, double* v, double *temp1){
	free(I); free(A); free(B); free(x); free(r); free(u); free(v); free(temp1); //free(arg->x_copy);
	I = nullptr; A = nullptr; B = nullptr; x = nullptr; r = nullptr; u = nullptr; v = nullptr; temp1 = nullptr; 
}
/*
bool init_vars(int nx, int ny, double* A, size_t* I, double *B, double* x, double *r, double* u, double* v, double *temp1){
	size_t N = (size_t) (nx+1)*(ny+1);
	A = (double*)malloc (sizeof (double) * (get_len_msr(nx,ny)+1));
	I = (size_t*)malloc(sizeof(size_t)*(get_len_msr(nx,ny)+1));
	B = (double*)malloc(sizeof(double) * N);
	x = (double*)malloc(sizeof(double) * N);
	//x_copy = (double*) malloc(sizeof(double) * N);
	r = (double*)malloc(sizeof(double) * N);
	u = (double*)malloc(sizeof(double) * N);
	v = (double*)malloc(sizeof(double) * N);
	temp1 = (double*) malloc(sizeof(double)*N);
	return true;
}
*/
void bind_vars(Args* args, double* A, size_t* I, double *B, double* x, double *r, double* u, double* v, double *temp1){
	args->A = A;
	args->I = I;
	args->B = B;
	args->x = x;
	args->u = u;
	args->v = v;
	args->r = r;
	args->temp = temp1;
}
/*
void* thread_f1(void* arg){
	Args *thr = (Args*) arg;
	
	cpu_set_t cpu;
	CPU_ZERO(&cpu);
	int n_cpus = get_nprocs(), p = thr->p;
	int cpu_id = n_cpus-1-(thr->kk%n_cpus);
	CPU_SET(cpu_id, &cpu);
	pthread_t tid = pthread_self();
	pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
	
	int nx = thr->nx, ny = thr->ny, kk = thr->kk, k = thr->k;
	double a = thr->a, b = thr->b, c = thr->c, d = thr->d, eps = thr->eps;
	double hx = (b-a)/nx, hy = (d-c)/ny;
	size_t N = (size_t) (nx+1)*(ny+1);
	size_t l1 = N*kk; l1/=p;
	size_t l2 = N*(kk+1); l2/=p;
	
	while(!thr->finish){
		if(thr->recalc){
			nx = thr->nx, ny = thr->ny, k = thr->k;
			N = (size_t) (nx+1)*(ny+1);
			l1 = N*kk; l1/=p;
			l2 = N*(kk+1); l2/=p;
			hx = (b-a)/nx, hy = (d-c)/ny;
			
			if(!thr->kk){
				if(thr->x != nullptr){
					free(thr->I); free(thr->A); free(thr->B); free(thr->x); free(thr->r); free(thr->u); free(thr->v); free(thr->temp);
				}
				init_(nx, ny, thr);
			}
			reduce_sum_init(p, kk, thr);
			
			
			memset(thr->B+l1, 0, (l2-l1)*sizeof(double));
			if(init_B(thr->B, nx, ny, hx, hy, a, c, k, kk, p)){
				thr->stat = false;
				cout << "FFFFFFFFFFFFF" << endl;
				return nullptr;
			}
			memset(thr->x+l1, 0, (l2-l1)*sizeof(double));
			memset(thr->r+l1, 0, (l2-l1)*sizeof(double));
			memset(thr->v+l1, 0, (l2-l1)*sizeof(double));
			memset(thr->u+l1, 0, (l2-l1)*sizeof(double));
			memset(thr->temp+l1, 0, (l2-l1)*sizeof(double));
			for(int i = 0; i < p; ++i){
				if(i == kk){
					fill_I(nx, ny, thr->I, kk, p);
				}
				reduce_sum(p);
			}
			fill_a(nx, ny, hx, hy, thr->I, thr->A, p, kk);
			reduce_sum(p);
			
			
			if(check_sym(nx, ny, thr->I, thr->A, eps, p, kk)){
				cout << "NOT SYM" << endl;
				thr->stat = false;
				return nullptr;
			}
			thr->t1 = get_full_time();
			thr->it = min_msr_solve(N, thr->A, thr->I, thr->B, thr->x, thr->r, thr->u, thr->v, thr->temp, eps, 50, thr->m, p, kk);
			thr->t1 = get_full_time()-thr->t1;
			reduce_sum(p);
			
			if(thr->it< 0 ){
				cout << "FFFFFFFFFFFFF" << endl;
				return nullptr;
			}
			thr->t2 = get_full_time();
			thr->r1 = r1(N, nx, ny, a, c, hx, hy, thr->x, k, kk, p);
			thr->r2 = r2(N, nx, ny, a, c, hx, hy, thr->x, k, kk, p);
			thr->r3 = r3(N, nx, a, c, hx, hy, thr->x, k, kk, p);
			thr->r4 = r4(N, nx, a, c, hx, hy, thr->x, k, kk, p);
			thr->t2 = get_full_time()-thr->t2;
			
			thr->recalc = false;
			reduce_sum(p);
		}
		reduce_sum_exchange(p, kk, thr);
	}
	cout << "Finished" << endl;
	if(!kk && thr->x != nullptr){
		free(thr->I); free(thr->A); free(thr->B); free(thr->x); free(thr->r); free(thr->u); free(thr->v); free(thr->temp);
	}
	reduce_sum_finishing(p, kk);
	return nullptr;
}
*/

void* thread_f(void* arg){
	Args *thr = (Args*) arg;
	
	cpu_set_t cpu;
	CPU_ZERO(&cpu);
	int n_cpus = get_nprocs(), p = thr->p;
	int cpu_id = n_cpus-1-(thr->kk%n_cpus);
	CPU_SET(cpu_id, &cpu);
	pthread_t tid = pthread_self();
	pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
	
	double* A = thr->A, *B = thr->B;
	size_t* I = thr->I;
	double a = thr->a, b = thr->b, c = thr->c, d = thr->d, eps = thr->eps;
	int nx = thr->nx, ny = thr->ny, kk = thr->kk, k = thr->k;
	double hx = (b-a)/nx, hy = (d-c)/ny;
	size_t N = (size_t) (nx+1)*(ny+1);
	size_t l1 = N*kk; l1/=p;
	size_t l2 = N*(kk+1); l2/=p;
	memset(B+l1, 0, (l2-l1)*sizeof(double));
	double func_mx = 0;
	if(init_B(B, nx, ny, hx, hy, a, c, k, kk, p,thr->p_func, func_mx)){
		thr->stat = false;
		return nullptr;
	}
	thr->func_mx = func_mx;
	memset(thr->x+l1, 0, (l2-l1)*sizeof(double));
	memset(thr->r+l1, 0, (l2-l1)*sizeof(double));
	memset(thr->v+l1, 0, (l2-l1)*sizeof(double));
	memset(thr->u+l1, 0, (l2-l1)*sizeof(double));
	memset(thr->temp+l1, 0, (l2-l1)*sizeof(double));
	for(int i = 0; i < p; ++i){
		if(i == kk){
			fill_I(nx, ny, I, kk, p);
		}
		reduce_sum(p);
	}
	fill_a(nx, ny, hx, hy, I, A, p, kk);
	reduce_sum(p);
	
	if(check_sym(nx, ny, I, A, eps, p, kk)){
		cout << "NOT SYM" << endl;
		thr->stat = false;
		return nullptr;
	}
	thr->t1 = get_full_time();
	thr->it = min_msr_solve(N, A, I, B, thr->x, thr->r, thr->u, thr->v, thr->temp, eps, 50, thr->m, p, kk);
	thr->t1 = get_full_time()-thr->t1;
	reduce_sum(p);
	
	if(thr->it< 0 ){
		return nullptr;
	}
	
	thr->t2 = get_full_time();
	thr->r1 = r1(N, nx, ny, a, c, hx, hy, thr->x, k, kk, p);
	thr->r2 = r2(N, nx, ny, a, c, hx, hy, thr->x, k, kk, p);
	thr->r3 = r3(N, nx, a, c, hx, hy, thr->x, k, kk, p);
	thr->r4 = r4(N, nx, a, c, hx, hy, thr->x, k, kk, p);
	thr->t2 = get_full_time()-thr->t2;
	
	reduce_sum(p);
	return nullptr;
}


int init_B(double* B, int nx, int ny, double hx, double hy, double x0, double y0, int k, int kk, int p, int p_func, double &func_mx){
	size_t N = (size_t) (nx+1)*(ny+1);
	size_t l, l1 = N*kk, l2 = N*(kk+1); l1/=p; l2/=p;
	int i,j;
	double err = 1e307;
	int stat = 0;
	double mx = 0;
	for(l = l1; l < l2; ++l){
		l2ij(nx, i, j, l);
		B[l] = F_ij(nx, ny, hx, hy, x0, y0, i, j, k);
		mx = max(mx, fabs(f(k, x0+i * hx, y0 + j * hy)));
		if(fabs(B[l]) > err){
			stat = 1;
		}
	}
	mx = reduce_sum_mx(p, mx);
	func_mx = mx;
	for(l = l1; l < l2; ++l){
		l2ij(nx, i, j, l);
		B[l] = add_p(nx, ny, hx, hy, x0, y0, i, j, k, mx, p_func);
		if(fabs(B[l]) > err){
			stat = 1;
		}
	}
	return reduce_sum_mx(p, stat);
}


int min_err_matr(size_t n, double *a, size_t *I, double *b, double *x, double *r, double *u, double *v, double *temp, double eps, int max_it, int p, int kk){
	double prec, b_norm, tau, c1, c2;
	int it;
	b_norm = scalar_product(n, b, b, p, kk); //(b,b)
	prec = b_norm*eps*eps;
	matr_mult_vector_msr(n, a, I, x, r, p, kk);
	mult_sub_vector(n, r, b, 1, p, kk); //r_{k+1}
	for(it = 0; it < max_it; ++it){
		apply_precondition(n, a, I, r, v, temp, p, kk);
		matr_mult_vector_msr(n, a, I, v, u, p, kk); // u = Av
		c1 = scalar_product(n, v, r, p, kk);
		c2 = scalar_product(n, v, u, p, kk);
		if(c1 <= prec || c2 <= prec){
			break;
		}
		tau = c1/c2;
		mult_sub_vector(n, x, v, tau, p, kk);//x -= tau*v
		mult_sub_vector(n, r, u, tau, p, kk);
	}
	if(it >= max_it)return -1;
	return it;
}


/*
int min_res_matr(size_t n, double *a, size_t *I, double *b, double *x, double *r, double *u, double *v, double *temp, double eps, int max_it, int p, int kk)
{
	double prec, b_norm, tau, c1, c2;
	int it;
	b_norm = scalar_product(n, b, b, p, kk); //(b,b)
	prec = b_norm*eps*eps;
	matr_mult_vector_msr(n, a, I, x, r, p, kk);
	mult_sub_vector(n, r, b, 1, p, kk); //r_{k+1}
	for(it = 0; it < max_it; ++it){
		apply_precondition(n, a, I, r, v, temp, p, kk);
		matr_mult_vector_msr(n, a, I, v, u, p, kk); // u = AV
		c1 = scalar_product(n, u, r, p, kk);
		c2 = scalar_product(n, u, u, p, kk);
		if(c1 <= prec || c2 <= prec){
			break;
		}
		tau = c1/c2;
		mult_sub_vector(n, x, v, tau, p, kk);//x -= tau*v
		mult_sub_vector(n, r, u, tau, p, kk);
	}
	if(it >= max_it)return -1;
	return it;
}
*/



int min_msr_solve(size_t n, double *a, size_t *I, double *b, double *x, double *r, double *u, double *v, double *temp, double eps, int max_it, int max_steps, int p, int kk){
	int step, ret, its = 0;
	for(step = 0; step < max_steps; ++step){
		ret = min_err_matr(n, a, I, b, x, r, u, v, temp, eps, max_it, p, kk);
		if(ret >= 0){its+=ret; break;}
		its+=max_it;
	}
	if(step >= max_steps){return -1;}
	return its;
}

/*
int min_msr_solve1(size_t n, double *a, size_t *I, double *b, double *x, double *r, double *u, double *v, double *temp, double eps, int max_it, int max_steps, int p, int kk){
	int step, ret, its = 0;
	for(step = 0; step < max_steps; ++step){
		ret = min_res_matr(n, a, I, b, x, r, u, v, temp, eps, max_it, p, kk);
		if(ret >= 0){its+=ret; break;}
		its+=max_it;
	}
	if(step >= max_steps){return -1;}
	return its;
}
*/


double scalar_product(size_t n, double *x, double *y, int p, int kk){
	int i1, i2, i;double s = 0;
	thread_rows(n, p, kk, i1, i2);
	for(i = i1; i < i2; ++i){
		s+=x[i]*y[i];
	}
	s = reduce_sum_det(p, kk, s);
	return s;
}


void mult_sub_vector(size_t n, double *x, double *y, double t, int p, int kk){ //x-=t*y
	int i, i1, i2;
	thread_rows(n, p, kk, i1, i2);
	for(i = i1; i<i2; ++i){
		x[i]-=t*y[i];
	}
	reduce_sum(p);
}


void matr_mult_vector_msr(size_t n, double *a, size_t *I, double *x, double *y, int p, int kk){ //mult matr by vector (par)
	int i, i1, i2, l, J;
	double s;
	thread_rows(n, p, kk, i1, i2);
	for(i = i1; i < i2; ++i){
		s = a[i]*x[i];
		l = I[i+1]-I[i];
		J = I[i];
		for(int j = 0; j < l; ++j){
			s+=a[J+j]*x[I[J+j]];
		}
		y[i] = s;
	}
	reduce_sum(p);
}


void apply_precondition(size_t n, double* a, size_t* I, double *r, double *v, double *temp, int p, int kk){
	int i, j, i1, i2, J, l, ind;
	double s;
	thread_rows(n, p, kk, i1, i2);
	/*
	for(i = i1; i < i2; ++i){
		v[i] = r[i]/a[i];
	}
	reduce_sum(p);
	return;
	*/
	for(i = i1; i < i2; ++i){
		s = r[i];
		J = I[i];
		l = I[i+1]-I[i];
		for(j = 0; j < l; ++j){
			ind = I[J+j];
			if((ind > i1) && (ind < i)){
				s-= a[J+j]*temp[ind];
			}
		}
		temp[i] = s/a[i];
	}
	for(i = i1; i < i2; ++i){
		temp[i]*=a[i];
	}
	for(i = i2-1; i >= i1; --i){
		s = temp[i];
		J = I[i];
		l = I[i+1]-I[i];
		for(j = 0; j < l; ++j){
			ind = I[J+j];
			if((ind > i) && (ind < i2-1)){
				s-= a[J+j]*v[ind];
			}
		}
		v[i] = s/a[i];
	}
	reduce_sum(p);
}


int check_sym(int nx, int ny, size_t *I, double *a, double eps, int p, int kk){
	size_t l1, l2, l;
	int err =0, q2;
	size_t N = (nx+1)*(ny+1);
	l1 = kk*N; l1/=p;
	l2 = N*(kk+1); l2/=p;
	for(l = l1; l < l2; ++l){
		int m = I[l+1]-I[l];
		double *a_offdiag = a+I[l];
		for(int q = 0; q < m; ++q){
			double a_ij = a_offdiag[q];
			int j = I[I[l]+q];
			int m2 = I[j+1]-I[j];
			for(q2 = 0; q2 < m2; q2++){
				if(I[I[j]+q2] == l){
					break;
				}
			}
			if(q2 >= m2){
				cout << I[I[j]+q2-1] << " " << l << " " << nx << " " << ny << endl;
				cout << setprecision(13) << a[I[j]+q2] << " " << a_ij << endl;
				err++;
			}
			else if(fabs(a[I[j]+q2]-a_ij) > eps){
				err++;
			}
		}
	}
	err = reduce_sum(p, err);
	return err;
}

#define F(I, J)(f(k, x0+(I)*hx, y0+(J)*hy)+(int(round(2*I))==nx)*(int(round(2*J))==ny)*0.1*mx*p_func)
double add_p(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, int k, double mx, int p_func){
	double w = hx*hy/192;
	if(i > 0 && i < nx && j > 0 && j < ny){
		return w*(36*F(i, j)+20*(F(i+0.5, j)+F(i, j-0.5)+F(i-0.5, j-0.5) + F(i-0.5, j) + F(i, j+0.5) + F(i+0.5, j+0.5)) +
		4*(F(i+0.5, j-0.5)+F(i-0.5, j-1)+F(i-1, j-0.5)+F(i-0.5, j+0.5)+F(i+0.5,j+1)+F(i+1, j+0.5)) + 
		2*(F(i+1, j)+F(i,j-1)+F(i-1,j-1)+F(i-1, j)+F(i,j+1)+F(i+1,j+1)));
	}
	else if(j == 0 && i > 0 && i < nx){
		return w*(18*F(i, j)+10*(F(i+0.5,j)+F(i-0.5, j))+20*(F(i, j+0.5)+F(i+0.5,j+0.5))+
		4*(F(i-0.5,j+0.5)+F(i+0.5,j+1)+F(i+1,j+0.5))+F(i-1,j)+F(i+1,j)+2*(F(i,j+1)+F(i+1,j+1)));
	}
	else if(j == ny && i > 0 && i < nx){
		return w*(18*F(i,j)+10*(F(i-0.5,j)+F(i+0.5,j))+20*(F(i,j-0.5)+F(i-0.5,j-0.5))+4*(F(i+0.5,j-0.5)+F(i-0.5,j-1)+F(i-1,j-0.5))
		+F(i-1,j)+F(i+1,j)+2*(F(i,j-1)+F(i-1,j-1)));
	}
	else if(i == 0 && j > 0 && j < ny){
		return w*(18*F(i,j)+10*(F(i, j-0.5)+F(i, j+0.5))+20*(F(i+0.5,j)+F(i+0.5,j+0.5))+4*(F(i+0.5,j-0.5)+F(i+0.5,j+1)+F(i+1,j+0.5))+F(i,j-1)+F(i,j+1)
		+2*(F(i+1,j)+F(i+1,j+1)));
	}
	else if(i == nx && j > 0 && j < ny){
		return w*(18*F(i, j)+10*(F(i,j-0.5)+F(i,j+0.5))+20*(F(i-0.5,j)+F(i-0.5,j-0.5))+4*(F(i-0.5,j-1)+F(i-1,j-0.5)+F(i-0.5,j+0.5))+F(i,j-1)+
		F(i,j+1)+2*(F(i-1,j)+F(i-1,j-1)));
	}
	else if(i == 0 && j == 0){
		return w*(12*F(i,j)+10*(F(i+0.5,j)+F(i,j+0.5))+20*F(i+0.5,j+0.5)+4*(F(i+1,j+0.5)+F(i+0.5,j+1))+F(i+1,j)+F(i,j+1)+2*F(i+1,j+1));
	}
	else if(i == nx && j == ny){
		return w*(12*F(i,j)+10*(F(i-0.5,j)+F(i,j-0.5))+20*F(i-0.5,j-0.5)+4*(F(i-0.5,j-1)+F(i-1,j-0.5))+F(i-1,j)+F(i,j-1)+2*F(i-1,j-1));
	}
	else if(i == 0 && j == ny){
		return w*(6*F(i,j)+10*(F(i+0.5,j)+F(i,j-0.5))+4*F(i+0.5,j-0.5)+F(i+1,j)+F(i,j-1));
	}
	else if(i == nx && j == 0){
		return w*(6*F(i,j)+10*(F(i-0.5,j)+F(i,j+0.5))+4*F(i-0.5,j+0.5)+F(i-1,j)+F(i,j+1));
	}
	else{
		abort();
	}
}
#undef F

#define F(I, J)(f(k, x0+(I)*hx, y0+(J)*hy))
double F_ij(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, int k){
	double w = hx*hy/192;
	if(i > 0 && i < nx && j > 0 && j < ny){
		return w*(36*F(i, j)+20*(F(i+0.5, j)+F(i, j-0.5)+F(i-0.5, j-0.5) + F(i-0.5, j) + F(i, j+0.5) + F(i+0.5, j+0.5)) +
		4*(F(i+0.5, j-0.5)+F(i-0.5, j-1)+F(i-1, j-0.5)+F(i-0.5, j+0.5)+F(i+0.5,j+1)+F(i+1, j+0.5)) + 
		2*(F(i+1, j)+F(i,j-1)+F(i-1,j-1)+F(i-1, j)+F(i,j+1)+F(i+1,j+1)));
	}
	else if(j == 0 && i > 0 && i < nx){
		return w*(18*F(i, j)+10*(F(i+0.5,j)+F(i-0.5, j))+20*(F(i, j+0.5)+F(i+0.5,j+0.5))+
		4*(F(i-0.5,j+0.5)+F(i+0.5,j+1)+F(i+1,j+0.5))+F(i-1,j)+F(i+1,j)+2*(F(i,j+1)+F(i+1,j+1)));
	}
	else if(j == ny && i > 0 && i < nx){
		return w*(18*F(i,j)+10*(F(i-0.5,j)+F(i+0.5,j))+20*(F(i,j-0.5)+F(i-0.5,j-0.5))+4*(F(i+0.5,j-0.5)+F(i-0.5,j-1)+F(i-1,j-0.5))
		+F(i-1,j)+F(i+1,j)+2*(F(i,j-1)+F(i-1,j-1)));
	}
	else if(i == 0 && j > 0 && j < ny){
		return w*(18*F(i,j)+10*(F(i, j-0.5)+F(i, j+0.5))+20*(F(i+0.5,j)+F(i+0.5,j+0.5))+4*(F(i+0.5,j-0.5)+F(i+0.5,j+1)+F(i+1,j+0.5))+F(i,j-1)+F(i,j+1)
		+2*(F(i+1,j)+F(i+1,j+1)));
	}
	else if(i == nx && j > 0 && j < ny){
		return w*(18*F(i, j)+10*(F(i,j-0.5)+F(i,j+0.5))+20*(F(i-0.5,j)+F(i-0.5,j-0.5))+4*(F(i-0.5,j-1)+F(i-1,j-0.5)+F(i-0.5,j+0.5))+F(i,j-1)+
		F(i,j+1)+2*(F(i-1,j)+F(i-1,j-1)));
	}
	else if(i == 0 && j == 0){
		return w*(12*F(i,j)+10*(F(i+0.5,j)+F(i,j+0.5))+20*F(i+0.5,j+0.5)+4*(F(i+1,j+0.5)+F(i+0.5,j+1))+F(i+1,j)+F(i,j+1)+2*F(i+1,j+1));
	}
	else if(i == nx && j == ny){
		return w*(12*F(i,j)+10*(F(i-0.5,j)+F(i,j-0.5))+20*F(i-0.5,j-0.5)+4*(F(i-0.5,j-1)+F(i-1,j-0.5))+F(i-1,j)+F(i,j-1)+2*F(i-1,j-1));
	}
	else if(i == 0 && j == ny){
		return w*(6*F(i,j)+10*(F(i+0.5,j)+F(i,j-0.5))+4*F(i+0.5,j-0.5)+F(i+1,j)+F(i,j-1));
	}
	else if(i == nx && j == 0){
		return w*(6*F(i,j)+10*(F(i-0.5,j)+F(i,j+0.5))+4*F(i-0.5,j+0.5)+F(i-1,j)+F(i,j+1));
	}
	else{
		abort();
	}
}
#undef F
