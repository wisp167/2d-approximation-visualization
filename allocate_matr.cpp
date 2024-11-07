#include "window.h"


void ij2l(int nx, int i, int j, size_t &l){l = i+j*(size_t)(nx+1);}
void l2ij(int nx, int &i, int &j, size_t l){j = l/(nx+1); i = l-j*(size_t)(nx+1);}

size_t get_len_msr(int nx, int ny){return (nx+1)*(ny+1)+6*(nx-1)*(ny-1)+4*(2*(nx-1)+2*(ny-1))+2*3+2*2;}
	
/*	
void fill_I(int nx, int ny, int i, int j, int *I_ij){
	int m  =0;
	if(i < nx){ij2l(nx, ny, i+1, j, I_ij[m]);m++;}
	if(j > 0){ij2l(nx, ny, i, j-1, I_ij[m]); m++;}
	if(i > 0 && j > 0){ij2l(nx, ny, I-1, j-1, I_ij[m]); ++m;}
	if(i > 0){ij2l(nx, ny, i-1, j, I_ij[m]); ++m;}
	if(j < ny){ij2l(nx, ny, i, j+1, I_ij[m]); ++m;}
	if(i < nx && j < ny){ij2l(nx, ny, i+1, j+1, I_ij[m]); ++m;}
	return m;
}
*/

#define F(I, J) do{ij2l(nx, I, J, k); if(I_ij){I_ij[m] = k;} m++;} while(0)
int get_off_diag(int nx, int ny, int i, int j, size_t *I_ij = nullptr);
int get_off_diag(int nx, int ny, int i, int j, size_t *I_ij){
	size_t k;
	int m = 0;
	if(i < nx){F(i+1, j);}
	if(j > 0){F(i, j-1);}
	if(i > 0 && j > 0){F(i-1, j-1);}
	if(i > 0){F(i-1, j);}
	if(j < ny){F(i, j+1);}
	if(i < nx && j < ny){F(i+1, j+1);}
	return m;
}
#undef F


/*
size_t get_len_msr_offdiag(int nx, int ny){
	size_t m  =0; int i, j;
	for(i = 0; i <= nx; ++i){
		for(j = 0; j <= ny; ++j){
			m+= get_off_diag(nx, ny, i, j);
		}
	}
	return m;
}


int allocate_msr_matr(int nx, int ny, double **p_a, size_t **p_I){
	size_t diag_len = (nx+1)*(ny+1);
	size_t off_diag = get_len_msr(nx, ny);
	size_t len = diag_len+off_diag+1;
	double *a = nullptr; size_t *I = nullptr;
	a = new double[len]; if(a == nullptr) return 1;
	I = new size_t[len]; if(I == nullptr) return 2;
	*p_a = a; *p_I = I;
	return 0;
}
*/


void fill_I(int nx, int ny, size_t *I, int kk, int p){
	size_t N = (size_t)(nx+1)*(ny+1);
	size_t l;
	static size_t r = N+1;
	int i, j, m = 0;
	if(!kk){r = N+1;}
	size_t l1, l2;
	l1 = kk*N; l1/=p;
	l2 = N*(kk+1); l2/=p;
	for(l = l1; l < l2; ++l){
		l2ij(nx, i, j, l); //l->(i, j)
		I[l] = r;
		m = get_off_diag(nx, ny, i, j, I+r);
		r+=m;
	}
	if(l == N){
		I[l] = r;
	}
}


void fill_a_ij(int nx, int ny, double hx, double hy, int i, int j, double *a_diag, double *a_offdiag){
	double s = (double)hx*hy;
	if(i > 0 && i < nx && j > 0 && j < ny){
		*a_diag = 6*s/12;
		for(int l = 0; l < 6; ++l){
			a_offdiag[l] = s/12;
		}
	}
	else if(j == 0 && i > 0 && i < nx){
		*a_diag = 3*s/12;
		a_offdiag[0] = 1*s/24;
		a_offdiag[1] = 1*s/24;
		a_offdiag[2] = 2*s/24;
		a_offdiag[3] = 2*s/24;
	}
	else if(j == ny && i > 0 && i < nx){
		*a_diag = 3*s/12;
		a_offdiag[0] = 1*s/24;
		a_offdiag[1] = 2*s/24;
		a_offdiag[2] = 2*s/24;
		a_offdiag[3] = 1*s/24;
	}
	else if(i == 0 && j > 0 && j < ny){
		*a_diag = 3*s/12;
		a_offdiag[0] = 2*s/24;
		a_offdiag[1] = 1*s/24;
		a_offdiag[2] = 1*s/24;
		a_offdiag[3] = 2*s/24;
	}
	else if(i == nx && j > 0 && j < ny){
		*a_diag = 3*s/12;
		a_offdiag[0] = 1*s/24;
		a_offdiag[1] = 2*s/24;
		a_offdiag[2] = 2*s/24;
		a_offdiag[3] = 1*s/24;
	}
	else if(i == 0 && j == 0){
		*a_diag = 2*s/12;
		a_offdiag[0] = 1*s/24;
		a_offdiag[1] = 1*s/24;
		a_offdiag[2] = 2*s/24;
	}
	else if(i == nx && j == ny){
		*a_diag = 2*s/12;
		a_offdiag[0] = 1*s/24;
		a_offdiag[1] = 2*s/24;
		a_offdiag[2] = 1*s/24;
	}
	else if(i == 0 && j == ny){
		*a_diag = 1*s/12;
		a_offdiag[0] = 1*s/24;
		a_offdiag[1] = 1*s/24;
	}
	else if(i == nx && j == 0){
		*a_diag = 1*s/12;
		a_offdiag[0] = 1*s/24;
		a_offdiag[1] = 1*s/24;
	}
	else{
		abort();
	}
}


void fill_a(int nx, int ny, double hx, double hy, size_t *I, double *a, int p, int kk){
	size_t l1, l2;
	size_t N = (size_t) (nx+1)*(ny+1);
	l1 = kk*N; l1/=p;
	l2 = N*(kk+1); l2/=p;
	int i, j;
	for(size_t l = l1; l < l2; ++l){
		l2ij(nx, i, j, l);
		double *a_diag = a+l;
		double *a_offdiag = a+I[l];
		fill_a_ij(nx, ny, hx, hy, i, j, a_diag, a_offdiag);
	}
	reduce_sum(p);
}
