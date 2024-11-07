#include "window.h"

static double *results = nullptr;

bool reduce_sum_finishing(int p, int kk){
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0, t_out = 0;
	pthread_mutex_lock(&m);
	if(kk == p){
		if(t_in != p){
			pthread_mutex_unlock(&m);
			return false;
		}
	}
	t_in++;
	if(t_in >= p+1){
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else{
		while(t_in < p+1){
			pthread_cond_wait(&c_in, &m);
		}
	}
	t_out++;
	if(t_out >= p+1){
		t_in = 0;
		pthread_cond_broadcast(&c_out);
	}
	else{
		while(t_out < p+1){
			pthread_cond_wait(&c_out, &m);
		}
	}
	pthread_mutex_unlock(&m);
	return true;
}
/*
bool reduce_sum_exchange(int p, int kk, Args *a){
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0, t_out = 0;
	static Args* args = new Args[p+1];
	int nx = 0, ny = 0;
	pthread_mutex_lock(&m);
	if(kk == p){
		if(t_in != p){
			pthread_mutex_unlock(&m);
			return false;
		}
		args[kk] = *a;
		if(args[kk] != args[0]){
			a->recalc = true;
		}
		nx = args[0].nx;
		ny = args[0].ny;
		t_in++;
	}
	else{
		args[kk] = *a;
		t_in++;
	}
	if(t_in >= p+1){
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else{
		while(t_in < p+1){
			pthread_cond_wait(&c_in, &m);
		}
	}
	if(kk != p){
		if(*a != args[p]){
			a->recalc = true;
			a->nx = args[p].nx;
			a->ny = args[p].ny;
			a->k = args[p].k;
			a->finish = args[p].finish;
		}
	}
	else if(a->recalc){
		if(args[0].x != nullptr){
			a->f_min = a->f_max = args[0].x[0];
			size_t N = (size_t)(nx+1)*(ny+1);
			for(size_t i = 0; i < N; ++i){
				a->x[i] = args[0].x[i];
				a->f_min = min(a->f_min, a->x[i]);
				a->f_max = max(a->f_max, a->x[i]);
			}
		}
		else{
			a->x = nullptr;
		}
	}
	t_out++;
	if(t_out >= p+1){
		t_in = 0;
		pthread_cond_broadcast(&c_out);
	}
	else{
		while(t_out < p+1){
			pthread_cond_wait(&c_out, &m);
		}
	}
	pthread_mutex_unlock(&m);
	return a->finish;
}
*/


void reduce_sum_init(int p, int kk, Args* thr){
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0, t_out = 0;
	static double* A = nullptr, *B = nullptr, *x = nullptr, *u = nullptr, *v = nullptr, *r = nullptr, *temp = nullptr;
	static size_t* I = nullptr;
	if(p <= 1){
		return;
	}
	pthread_mutex_lock(&m);
	if(!kk){
		A = thr->A;
		B = thr->B;
		x = thr->x;
		u = thr->u;
		v = thr->v;
		r = thr->r;
		temp = thr->temp;
		I = thr->I;
	}
	t_in++;
	if(t_in >= p){
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else{
		while(t_in < p){
			pthread_cond_wait(&c_in, &m);
		}
	}
	if(kk){
		thr->A = A;
		thr->B = B;
		thr->x = x;
		thr->u = u;
		thr->v = v;
		thr->r = r;
		thr->temp = temp;
		thr->I = I;
	}
	t_out++;
	if(t_out >= p){
		t_in = 0;
		pthread_cond_broadcast(&c_out);
	}
	else{
		while(t_out < p){
			pthread_cond_wait(&c_out, &m);
		}
	}
	pthread_mutex_unlock(&m);
	return;
}


void reduce_sum(int p){
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0, t_out = 0;
	if(p <= 1){
		return;
	}
	pthread_mutex_lock(&m);
	t_in++;
	if(t_in >= p){
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else{
		while(t_in < p){
			pthread_cond_wait(&c_in, &m);
		}
	}
	t_out++;
	if(t_out >= p){
		t_in = 0;
		pthread_cond_broadcast(&c_out);
	}
	else{
		while(t_out < p){
			pthread_cond_wait(&c_out, &m);
		}
	}
	pthread_mutex_unlock(&m);
	return;
}


double reduce_sum(int p, double val){
	static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0, t_out = 0;
	static double ans = 0;
	if(p <= 1){
		return val;
	}
	pthread_mutex_lock(&mu);
	t_in++;
	if(val > ans){
		ans += val;
	}
	if(t_in >= p){
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else{
		while(t_in < p){
			pthread_cond_wait(&c_in, &mu);
		}
	}
	val = ans;
	t_out++;
	if(t_out >= p){
		t_in = 0;
		ans = 0;
		pthread_cond_broadcast(&c_out);
	}
	else{
		while(t_out < p){
			pthread_cond_wait(&c_out, &mu);
		}
	}
	pthread_mutex_unlock(&mu);
	return val;
}


double reduce_sum_mx(int p, double val){
	static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0, t_out = 0;
	static double ans = 0;
	if(p <= 1){
		return val;
	}
	pthread_mutex_lock(&mu);
	t_in++;
	if(val > ans){
		ans = val;
	}
	if(t_in >= p){
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else{
		while(t_in < p){
			pthread_cond_wait(&c_in, &mu);
		}
	}
	val = ans;
	t_out++;
	if(t_out >= p){
		t_in = 0;
		ans = 0;
		pthread_cond_broadcast(&c_out);
	}
	else{
		while(t_out < p){
			pthread_cond_wait(&c_out, &mu);
		}
	}
	pthread_mutex_unlock(&mu);
	return val;
}


int init_reduce_sum(int p){
	results = new double[p];
	if(results == nullptr) return -1;
	return 0;
}


int delete_reduce_sum(){
	if(results == nullptr){return -1;}
	delete[] results;
	return 0;
}


double reduce_sum_det(int p, int k, double s){
	double sum = 0; int l;
	results[k] = s;
	reduce_sum(p);
	for(l = 0; l < p; ++l){
		sum+=results[l];
	}
	reduce_sum(p);
	return sum;
}
