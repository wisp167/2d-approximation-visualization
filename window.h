#ifndef WINDOW_H
#define WINDOW_H
#pragma once

#include <QtWidgets/QtWidgets>
#include <fenv.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <pthread.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/sysinfo.h>
using namespace std;

struct Args{
	double a, b, c, d, eps, t1 = 0, t2 = 0, f_min, f_max;
	int nx, ny, kk, k, m, p, it = 0, p_func = 0, func_mx;
	char* name = nullptr;
	double* A = nullptr, *B = nullptr, *x = nullptr, *u = nullptr, *v = nullptr, *r = nullptr, *temp = nullptr, *x_copy = nullptr;
	size_t* I = nullptr;
	bool stat = true;
	bool calculating = false;
	bool update_req = false;
	//bool recalc = false;
	//bool calculations_finished = false;
	bool finish = false;
	double r1 = -1, r2 = -1, r3 = -1, r4 = -1;
	/*void operator=(const Args &arg){
		nx = arg.nx;
		ny = arg.ny;
		k = arg.k;
		p_func = arg.p_func;
	}
	*/
	bool operator!=(const Args &arg){
		if(nx != arg.nx
		|| ny != arg.ny
		|| k != arg.k
		|| p_func != arg.p_func){
			return true;
		}
		return false;
	}
	Args() = default;
	~Args() = default;
};

class foo : public QMainWindow
{
    Q_OBJECT
signals:
    void signalName();
private:
    void closeEvent(QCloseEvent *bar);
    bool stat = false;
public slots:
	void change_stat();
};

class Window : public QWidget
{
  Q_OBJECT

private:
  Args data;
  int func_id, method;
  const char *f_name;
  const char *min_text = "min =";
  const char *max_text = "max =";
  const char *scale_text = "scale =";
  const char *nx_text = "nx =";
  const char *ny_text = "ny =";
  const char *my_text = "my =";
  const char *mx_text = "mx =";
  const char *method_text = "mode =";
  const char *p_text = "p =";
  const char *calc_text = "Waiting for the end of calculations";
  const char *ending_message = "PREPARING FOR CLOSING, double tap please";
  char *f_min;
  char *f_max;
  char *sscale;
  char *nnx;
  char *nny;
  char *nmx;
  char *nmy;
  char *mmethod;
  char *pp_;
  double *x;
  double scale, output_scale;
  int s;
  double a_, b_, c_, d_, eps;
  double ag, bg, cg, dg;
  double delta_y;
  double f_mx, f_mn, true_mn, true_mx;
  bool ending = false, not_init = true;
  size_t curr_N;
  int nx, ny, mx, my, m, p, pp;
  int nx_upd, ny_upd, pp_upd, k_upd;
  double (*f) (double,double);
  double (*df) (double, double);
  ~Window() = default;
signals:
   void signalName1();
public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;
  int parse_command_line (int argc, char *argv[]);
  void init(int nx, int ny);
  QPointF l2g (double x_loc, double y_loc, double scale);
  void paint_triangle(QPainter *painter, const QPointF &p1, const QPointF &p2, const QPointF &p3, double val, double y_min, double y_max);
  void update_func();
  void print_text();
  void update_x();
  double g2lx(double x, double scale);
  double g2ly(double y, double scale);
public slots:
  void preparation_for_finishing();
  void change_func ();
  void increase_scale();
  void decrease_scale();
  void change_mode();
  void increase_n();
  void decrease_n();
  void increase_p();
  void decrease_p();
  void increase_m();
  void decrease_m();

protected:
  void simple();
  void swap();
  void paintEvent (QPaintEvent *event);
};

void* thread_window_f(void* arg);
double add_p(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, int k, double mx, int p_func);
void free_vars(double* A, size_t* I, double *B, double* x, double *r, double* u, double* v, double *temp1);
bool init_vars(int nx, int ny, double* A, size_t* I, double *B, double* x, double *r, double* u, double* v, double *temp1);
void bind_vars(Args* arg, double* A, size_t* I, double *B, double* x, double *r, double* u, double* v, double *temp1);
void err_output(char* str);
double get_cpu_time();
double get_full_time();
int read_input(int* temp, double* ttemp, char** argv);
double f(int k, double x, double y);
void* thread_f(void* arg);
void* thread_f1(void* arg);
double scalar_product(size_t n, double *x, double *y, int p, int kk);
int min_err_matr(size_t n, double *a, size_t *I, double *b, double *x, double *r, double *u, double *v, double *temp, double eps, int max_it, int p, int kk);
int min_msr_solve(size_t n, double *a, size_t *I, double *b, double *x, double *r, double *u, double *v, double *temp, double eps, int max_it, int max_steps, int p, int kk);
void mult_sub_vector(size_t n, double *x, double *y, double t, int p, int kk);
void matr_mult_vector_msr(size_t n, double *a, size_t *I, double *x, double *y, int p, int kk);
void apply_precondition(size_t n, double* a, size_t* I, double *v, double *r, double *temp, int p, int kk);
int check_sym(int nx, int ny, size_t *I, double *a, double eps, int p, int kk);
double F_ij(int nx, int ny, double hx, double hy, double x0, double y0, int i, int j, int k);
int init_B(double* B, int nx, int ny, double hx, double hy, double x0, double y0, int k, int kk, int p, int p_func, double &func_mx);
int min_msr_solve1(size_t n, double *a, size_t *I, double *b, double *x, double *r, double *u, double *v, double* temp, double eps, int max_it, int max_steps, int p, int kk);


void ij2l(int nx, int i, int j, size_t &l);
void l2ij(int nx, int &i, int &j, size_t l);
void fill_a(int nx, int ny, double hx, double hy, size_t *I, double *a, int p, int kk);
size_t get_len_msr(int nx, int ny);
void fill_I(int nx, int ny, size_t *I, int kk, int p);
void reduce_sum(int p);
double reduce_sum_mx(int p, double val);
double reduce_sum(int p, double val);
int init_reduce_sum(int p);
int delete_reduce_sum();
double reduce_sum_det(int p, int k, double s);
bool reduce_sum_finishing(int p, int kk);
void reduce_sum_init(int p, int kk, Args* thr);
bool reduce_sum_exchange(int p, int kk, Args *a);

double f_z(double x_, double y_, double a, double c, int nx, double hx, double hy, double *x);
double r1(size_t N, int nx, int ny, double a, double c, double hx, double hy, double *x, int k, int kk, int p);
double r2(size_t N, int nx, int ny, double a, double c, double hx, double hy, double *x, int k, int kk, int p);
double r3(size_t N, int nx, double a, double c, double hx, double hy, double *x, int k, int kk, int p);
double r4(size_t N, int nx, double a, double c, double hx, double hy, double *x, int k, int kk, int p);
using namespace std;

double f_0 (double x, double y);
double f_1 (double x, double y);
double f_2 (double x, double y);
double f_3 (double x, double y);
double f_4 (double x, double y);
double f_5 (double x, double y);
double f_6 (double x, double y);
double f_7 (double x, double y);

/*
double df_0 (double x, double y);
double df_1 (double x, double y);
double df_2 (double x, double y);
double df_3 (double x, double y);
double df_4 (double x, double y);
double df_5 (double x, double y);
double df_6 (double x, double y);
*/

#endif
