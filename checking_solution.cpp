#include "window.h"


inline double det2(double x, double y, double x1, double y1, double x2, double y2){
	return x1*y2+x2*y+x*y1-x1*y-x*y2-x2*y1;
}


double calc_z(double x1, double y1, double f1, double x2, double y2, double f2, double x3, double y3, double f3, double x, double y){
	return f1+((f2-f1)*det2(x1, y1, x, y, x3, y3)-(f3-f1)*det2(x1, y1, x, y, x2, y2))/det2(x1, y1, x2, y2, x3, y3);
}


double f_z(double x_, double y_, double a, double c, int nx, double hx, double hy, double *x){
	int i, j;
	double x1, x2, x3, y1, y2, y3;
	i = int((x_-a)/hx);
	j = int((y_-c)/hy);
	
	x1 = a + i*hx;
	y1 = c + j*hy;
	x3 = a + (i+1)*hx;
	y3 = c + (j+1)*hy;
	size_t l1, l2, l3;
	ij2l(nx, i, j, l1);
	ij2l(nx, i+1, j+1, l3);
	if(y_-y1 > x_-x1){
		x2 = x1;
		y2 = y3;
		ij2l(nx, i, j+1, l2);
	}
	else{
		x2 = x3;
		y2 = y1;
		ij2l(nx, i+1, j, l2);
	}
	return calc_z(x1, y1, x[l1], x2, y2, x[l2], x3, y3, x[l3], x_, y_);
}


double r1(size_t N, int nx, int ny, double a, double c, double hx, double hy, double *x, int k, int kk, int p){
	size_t l, l1 = N*kk, l2 = N*(kk+1); l1/=p; l2/=p;
	double mx = 0;
	int i, j;
	double x_, y_;
	for(l = l1; l < l2; ++l){
		l2ij(nx, i, j, l);
		if((i == nx) || (j == ny)){
			continue;
		}
		x_ = a+hx*i; y_ = c+hy*j;
		mx = max(mx, fabs(f(k, x_ + hx/3, y_ + 2*hy/3)-f_z(x_+hx/3, y_+2*hy/3, a, c, nx, hx, hy, x)));
		mx = max(mx, fabs(f(k, x_ + 2*hx/3, y_ + hy/3)-f_z(x_+2*hx/3, y_+hy/3, a, c, nx, hx, hy, x)));
	}
	return reduce_sum_mx(p, mx);
}


double r2(size_t N, int nx, int ny, double a, double c, double hx,double hy, double *x, int k, int kk, int p){
	size_t l, l1 = N*kk, l2 = N*(kk+1); l1/=p; l2/=p;
	double sum = 0;
	int i, j;
	double x_, y_;
	for(l = l1; l < l2; ++l){
		l2ij(nx, i, j, l);
		if((i == nx) || (j == ny)){
			continue;
		}
		x_ = a+hx*i; y_ = c+hy*j;
		sum += fabs(f(k, x_ + hx/3, y_ + 2*hy/3)-f_z(x_+hx/3, y_+2*hy/3, a, c, nx, hx, hy,x));
		sum += fabs(f(k, x_ + 2*hx/3, y_ + hy/3)-f_z(x_+2*hx/3, y_+hy/3, a, c, nx, hx, hy,x));
	}
	return reduce_sum_det(p, kk, sum)*hx*hy/2;
}


double r3(size_t N, int nx, double a, double c, double hx, double hy, double *x, int k, int kk, int p){
	size_t l, l1 = N*kk, l2 = N*(kk+1); l1/=p; l2/=p;
	double mx = 0;
	int i, j;
	double x_, y_;
	for(l = l1; l < l2; ++l){
		l2ij(nx, i, j, l);
		x_ = a+hx*i; y_ = c+hy*j;
		mx = max(mx, fabs(f(k, x_, y_)-x[l]));
	}
	return reduce_sum_mx(p, mx);
}

double r4(size_t N, int nx, double a, double c, double hx, double hy, double *x, int k, int kk, int p){
	size_t l, l1 = N*kk, l2 = N*(kk+1); l1/=p; l2/=p;
	double sum = 0;
	int i, j;
	double x_, y_;
	for(l = l1; l < l2; ++l){
		l2ij(nx, i, j, l);
		x_ = a+hx*i; y_ = c+hy*j;
		sum += fabs(f(k, x_, y_)-x[l]);
	}
	return reduce_sum_det(p, kk, sum)*hx*hy;
}
