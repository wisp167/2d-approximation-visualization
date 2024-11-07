
#include <QPainter>
#include <stdio.h>

#include "window.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_C -10
#define DEFAULT_D 10
#define DEFAULT_N 10
#define L2G(X,Y) (l2g ((X), (Y), scale))

Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a_ = DEFAULT_A;
  b_ = DEFAULT_B;
  c_ = DEFAULT_C;
  d_ = DEFAULT_D;
  nx = DEFAULT_N;
  ny = DEFAULT_N;
  scale = output_scale = 1;
  s = 0;
  method = 0;
  func_id = 0;
  curr_N = 0;
  delta_y = 0;
  p = 0;
  pp = 0;
  x = nullptr;
  change_func ();
}
QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc < 13)
    return -1;

  if (   sscanf (argv[1], "%lf", &a_) != 1
      || sscanf (argv[2], "%lf", &b_) != 1
      || sscanf (argv[3], "%lf", &c_) != 1
      || sscanf (argv[4], "%lf", &d_) != 1
      || sscanf (argv[5], "%d", &nx) != 1
      || sscanf (argv[6], "%d", &ny) != 1
      || sscanf (argv[7], "%d", &mx) != 1
      || sscanf (argv[8], "%d", &my) != 1
      || sscanf (argv[9], "%d", &func_id) != 1
      || sscanf (argv[10], "%lf", &eps) != 1
      || sscanf (argv[11], "%d", &m) != 1
      || sscanf (argv[12], "%d", &p) != 1
      || b_ - a_ < 1.e-6
      || d_ - c_ < 1.e-6
      || mx <= 1
      || my <= 1
      || func_id < 0
      || func_id > 8
      || nx <= 1
      || ny <= 1)
    return -2;
  bg = b_;
  cg = c_;
  ag = a_;
  dg = d_;
  data.nx = nx;
  data.ny = ny;
  data.p = p;
  data.k = func_id;
  data.kk = p;
  nx_upd = nx;
  ny_upd = ny;
  k_upd = func_id;
  pp_upd = 0;
  //update_x();
  update_func();
  return 0;
}

void Window::preparation_for_finishing(){
	data.finish = true;
	update();
}
inline int ipow(int base, int exp)
{
    int result = 1;
    while(true){
        if(exp & 1){
            result *= base;
		}
        exp >>= 1;
        if (!exp){
            break;
		}
        base *= base;
	}
    return result;
}
void Window::increase_scale(){
	if(s < 16){
		s++;
		double temp = (bg-ag)/ipow(2, s), temp1 = (bg+ag)/2;
		b_ = temp1+temp;
		a_ = temp1-temp;
		temp = (bg-ag)/ipow(2, s), temp1 = (bg+ag)/2;;
		d_ = temp1+temp;
		c_ = temp1-temp;
		output_scale*=2;
	}
	update();
}
void Window::decrease_scale(){
	if(s){
		s--;
		double temp = (bg-ag)/ipow(2, s), temp1 = (bg+ag)/2;
		b_ = temp1+temp;
		a_ = temp1-temp;
		temp = (bg-ag)/ipow(2, s), temp1 = (bg+ag)/2;;
		d_ = temp1+temp;
		c_ = temp1-temp;
		output_scale/=2;
	}
	if(!s){
		a_ = ag;
		b_ = bg;
		c_ = cg;
		d_ = dg;
	}
	update();
}
void Window::increase_n(){
	nx_upd*=2;
	ny_upd*=2;
	update();
}
void Window::decrease_n(){
	if(nx > 3){
		nx_upd/=2;
	}
	if(ny > 3){
		ny_upd/=2;
	}
	update();
}
void Window::increase_m(){
	mx*=2;
	my*=2;
	update();
}
void Window::decrease_m(){
	if(mx > 3){
		mx/=2;
	}
	if(my > 3){
		my/=2;
	}
	update();
}
void Window::increase_p(){
	pp_upd++;
	update();
}
void Window::decrease_p(){
	pp_upd--;
	update();
}
void Window::change_mode(){
	method = (method+1)%3;
	update();
}

void* thread_window_f(void* arg){
	Args *thr = (Args*) arg;
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_mutex_t m1 = PTHREAD_MUTEX_INITIALIZER;
	static Args vars_to_copy;
	static double *A = nullptr, *B = nullptr, *x = nullptr, *x_copy = nullptr, *u = nullptr, *v =nullptr, *r = nullptr, *temp1 = nullptr;
	static size_t *I = nullptr;
	static bool recalc = false, calculations_finished = false, finish = false;
	int p = thr->p;
	if(thr->kk == p){ // graphic
		if(A == nullptr){
			vars_to_copy.nx = thr->nx;
			vars_to_copy.ny = thr->ny;
			vars_to_copy.k = thr->k-1;
			vars_to_copy.p_func = thr->p_func-1;
		}
		if(thr->finish){
			finish = true;
			if(x_copy != nullptr){
				free(x_copy);
				x_copy = nullptr;
			}
			return nullptr;
		}
		
		
		if(calculations_finished){
			if(x_copy != nullptr){
				free(x_copy);
				x_copy = nullptr;
			}
			size_t N = (size_t) (vars_to_copy.nx+1)*(vars_to_copy.ny+1);
			x_copy = (double*)malloc(sizeof(double) * N);
			for(size_t i = 0; i < N; ++i){
				x_copy[i] = x[i];
			}
			thr->func_mx = vars_to_copy.func_mx;
			thr->x_copy = x_copy;
			thr->calculating = false;
			thr->update_req=true;
			calculations_finished = false;
		}
		if((*thr != vars_to_copy) && (!thr->calculating)){
			vars_to_copy.nx = thr->nx;
			vars_to_copy.ny = thr->ny;
			vars_to_copy.k = thr->k;
			vars_to_copy.p_func = thr->p_func;
			thr->calculating = true;
			recalc = true;
		}
	}
	else{ // calc
		bool enter = false;
		int counter = 0;
		while(!finish){
			//cout << recalc << endl;
			pthread_mutex_lock(&m1);
			if(recalc){
				enter = true;
			}
			pthread_mutex_unlock(&m1);
			if(enter){
				enter = false;
				reduce_sum(p);
				thr->nx = vars_to_copy.nx;
				thr->ny = vars_to_copy.ny;
				thr->k = vars_to_copy.k;
				thr->p_func = vars_to_copy.p_func;
				size_t N = (size_t) (thr->nx+1)*(thr->ny+1);
				if(!thr->kk){
					if(A != nullptr){
						free_vars(A, I, B, x,r, u, v, temp1);
					}
					//init_vars(thr->nx, thr->ny, A, I, B, x,r, u, v, temp1);
					A = (double*)malloc (sizeof (double) * (get_len_msr(thr->nx,thr->ny)+1));
					I = (size_t*)malloc(sizeof(size_t)*(get_len_msr(thr->nx,thr->ny)+1));
					B = (double*)malloc(sizeof(double) * N);
					x = (double*)malloc(sizeof(double) * N);
					r = (double*)malloc(sizeof(double) * N);
					u = (double*)malloc(sizeof(double) * N);
					v = (double*)malloc(sizeof(double) * N);
					temp1 = (double*) malloc(sizeof(double)*N);
				}
				//bind_vars(thr, A, I, B, x,r, u, v, temp1);
				reduce_sum(p);
				
				bind_vars(thr, A, I, B, x,r, u, v, temp1);
				reduce_sum(p);
				
				thread_f(arg);
				reduce_sum(p);
				if(!thr->kk){
					pthread_mutex_lock(&m);
					recalc = false;
					vars_to_copy.func_mx = thr->func_mx;
					calculations_finished = true;
					pthread_mutex_unlock(&m);
					if(counter)
						printf ("%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
	thr->name, 6, thr->r1, thr->r2, thr->r3, thr->r4, thr->t1, thr->t2, thr->it, thr->eps, thr->k, thr->nx, thr->ny, thr->p);
					++counter;
				}
				reduce_sum(p);
			}
		}
		thr->finish = true;
		reduce_sum(p);
		if(!thr->kk){
			free(I); free(A); free(B); free(r); free(x); free(u); free(v); free(temp1);
		}
		reduce_sum_finishing(p, thr->kk);
		cout << "PROC ENDING" << endl;
	}
	return nullptr;
}

/// change current function for drawing
void Window::change_func ()
{
  k_upd = (k_upd+1)%8;
  update();
}
void Window::update_func()
{
  func_id = data.k;
  switch (func_id)
    {
      case 0:
        f_name = "k = 0, f (x) = 1";
        f = f_0;
        //df = df_0;
        break;
      case 1:
        f_name = "k = 1, f (x) = x";
        f = f_1;
        //df = df_1;
        break;
      case 2:
        f_name = "k = 2, f (x) = y";
        f = f_2;
        //df = df_2;
        break;
      case 3:
        f_name = "k = 3, f (x) = x+y";
        f = f_3;
        //df = df_3;
        break;
      case 4:
        f_name = "k = 4, f (x) = sqrt(x^2+y^2)";
        f = f_4;
        //df = df_4;
        break;
      case 5:
        f_name = "k = 5, f (x) = x^2+y^2";
        f = f_5;
        //df = df_5;
        break;
      case 6:
        f_name = "k = 6, f (x) = exp(x^2-y^2)";
        f = f_6;
        //df = df_6;
        break;
      case 7:
		f_name = "k = 7, f (x) = 1/(25*(x^2+y^2)+1)";
		f = f_7;
    }
    /*
    double temp;
    double hx = (b_-a_)/nx, hy = (d_-c_)/ny;
    true_mn = true_mx = f(a_, c_);
    for(int i = 0; i < nx; ++i){
		for(int j  =0; j < ny; ++j){
			temp = f(a_+i*hx, c_ + j*hy);
			true_mn = min(temp, true_mn);
			true_mx = max(temp, true_mx);
		}
	}
	*/
}
/*
void Window::update_func(){
  double x1; int m;
  for (x1 = a_, m = 0; m < nx; ++m, x1 = a_+(m*(b_-a_))/(nx-1))
  {
      d[m] = df(x1, x1);
      f_ar[m] = f(x1, x1);
  }
}
void Window::update_x(){
  for (int m = 0; m < nx; ++m)
  {
      x[m] = a_+((m*(b_-a_))/(nx-1));
  }
}
*/
/*
void Window::init(int nx, int ny){
	size_t N = (size_t)(nx+1)*(ny+1);
	double *x1 = (double*) malloc(sizeof(double)*N);
	if(x != nullptr){for(size_t i = 0; i < curr_N; ++i){x1[i] = x[i];} free(x);}
	data.x = x1;
	x = x1;
}
*/
QPointF Window::l2g (double x_loc, double y_loc, double scale)
{
  double l = (a_*(scale+1)+b_*(scale-1))/(2*scale);
  double l1 = (c_*(scale+1)+d_*(scale-1))/(2*scale);
  double x_gl = (x_loc - l) / (b_ - a_) * width () * scale;
  double y_gl = (y_loc - l1) / (d_ - c_) * height () * scale;
  return QPointF (x_gl, height() - y_gl);
}
double Window::g2lx(double x, double scale){
	double c2 = (a_*(scale+1)+b_*(scale-1))/(2*scale);
	double c1 = ((a_+b_)/2+(b_-a_)/(2*scale)-c2)/width();
	return c1*x+c2;
}
double Window::g2ly(double y, double scale){
	double c2 = (c_*(scale+1)+d_*(scale-1))/(2*scale);
	double c1 = ((d_+c_)/2+(d_-c_)/(2*scale)-c2)/height();
	return c1*y+c2;
}
void Window::print_text(){
  QPainter painter (this);
  painter.setPen ("yellow");
  f_min = (char*) malloc(sizeof(char)*100);
  f_max = (char*) malloc(sizeof(char)*100);
  sscale = (char*) malloc(sizeof(char)*100);
  nnx = (char*) malloc(sizeof(char)*100);
  nny = (char*) malloc(sizeof(char)*100);
  nmx = (char*) malloc(sizeof(char)*100);
  nmy = (char*) malloc(sizeof(char)*100);
  mmethod = (char*) malloc(sizeof(char)*100);
  pp_ = (char*) malloc(sizeof(char)*100);
  snprintf(sscale, 100, "%lf", output_scale);
  snprintf(f_min, 100, "%e", f_mn);
  snprintf(f_max, 100, "%e", f_mx);
  snprintf(nnx, 100, "%d", nx);
  snprintf(nny, 100, "%d", ny);
  snprintf(nmx, 100, "%d", mx);
  snprintf(nmy, 100, "%d", my);
  snprintf(pp_, 100, "%d", pp);
  snprintf(mmethod, 100, "%d", method);
  painter.drawText (20, 20, f_name);
  painter.drawText (20, 40, min_text);
  painter.drawText(60, 40, f_min);
  painter.drawText(20, 60, max_text);
  painter.drawText(60, 60, f_max);
  painter.drawText(20, 80, scale_text);
  painter.drawText(70, 80, sscale);
  painter.drawText(20, 100, nx_text);
  painter.drawText(50, 100, nnx);
  painter.drawText(20, 120, ny_text);
  painter.drawText(50, 120, nny);
  painter.drawText(20, 140, mx_text);
  painter.drawText(50, 140, nmx);
  painter.drawText(20, 160, my_text);
  painter.drawText(50, 160, nmy);
  painter.drawText(20, 180, method_text);
  painter.drawText(70, 180, mmethod);
  painter.drawText(20, 200, p_text);
  painter.drawText(40, 200, pp_);
  if(data.calculating){
	  painter.drawText(250, 220, calc_text);
  }
  if(data.finish){
	  painter.drawText(300, 220, ending_message);
  }
  free(f_min); free(f_max); free(sscale); free(nnx); free(nny); free(nmx); free(nmy); free(mmethod); free(pp_);
}

void Window::paint_triangle(QPainter *painter, const QPointF &p1, const QPointF &p2, const QPointF &p3, double val, double y_min, double y_max)
{
        QPolygonF polygon;
        polygon << p1 << p2 << p3;
        double k1, k2;
        if(fabs(y_min-y_max) < eps){
			k1 = 0;
			k2 = 1;
		}
		else{
			k1 = (val-y_min)/(y_max-y_min);
			k2 = (y_max - val)/(y_max-y_min);
		}
        painter->setBrush(QBrush(QColor(min(max(k1, (double)0), (double)1)*255, 0, min(max(k2, (double)0), (double)1)*255)));
        painter->drawPolygon(polygon);
        //Q_UNUSED(option);
        //Q_UNUSED(widget);
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{
  QPainter painter (this);
  double x1, y1, x2, y2, tem;
  //QPen pen_black(Qt::black, 0, Qt::SolidLine); 
  //QPen pen_red(Qt::red, 0, Qt::SolidLine); 
  //painter.setPen (pen_black);
  /*if(!ending){
	ending = reduce_sum_exchange(p, p, &data);
	if(data.recalc){
		data.recalc = false;
		not_init = false;
		nx = data.nx;
		ny = data.ny;
		x = data.x;
		f_mn = data.f_min;
		f_mx = data.f_max;
		update_func();
		curr_N = (size_t)(nx+1)*(ny+1);
	}
	if(not_init){
		return;
	}
  }
  else if(reduce_sum_finishing(p, p)){
	emit signalName1();close();
	return;//close();
  }
  */
  //bool redraw = false;
  if(data.finish){
	  thread_window_f(&data);
	  while(!reduce_sum_finishing(p, p)){
		  print_text();
	  }
	  //free(data.x); free(data.x_copy);
	  QCoreApplication::quit();
	  return;
	  emit signalName1();
	  close();
	  return;
  }
  if(!not_init){
	  thread_window_f(&data);
  }
  else{
	  while(!data.update_req){
		  thread_window_f(&data);
	  }
	  not_init =false;
  }
  if(data.update_req){ 
	  nx = data.nx;
	  ny = data.ny;
	  if(func_id != data.k){
		 update_func(); 
	  }
	  pp = data.p_func;
	  x = data.x_copy;
	  data.update_req = false;
	  //redraw = true;
  }
  if(!data.calculating){
	  data.nx = nx_upd;
	  data.ny = ny_upd;
	  data.k = k_upd;
	  data.p_func = pp_upd;
  }
  /*
  if(!redraw){
	  cout << "redra" << endl;
	  if(data.calculating){
		  cout << "redraw" << endl;
		  painter.setPen("black");
		  painter.drawText(250, 220, calc_text);
	  }
	  update();
  }
  */
  
  double y_max = fabs(f(a_, c_)), y_min = fabs(f(a_, c_));
  //double ym_max = f(a_, c_), ym_min = f(a_, c_);
  //size_t N = (size_t) (nx+1)*(ny+1);
  //size_t l;
  double hx = (bg-ag)/nx, hy = (dg-cg)/ny;
  
  
  QPointF p1, p2, p3;
  if(!method){
	  double hmx = (b_-a_)/mx, hmy = (d_-c_)/my;
	  /*
	  for(int i = 0; i <= nx; ++i){
		  for(int j = 0; j <= ny; ++j){
			  x1 = a_ + i*hx; y1 = c_ + j*hy;
			  tem = f(x1, y1);
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  //tem = f(x2, y2);
			  //y_min = min(y_min, tem); y_max = max(y_max, tem);
		  }
	  }
	  */
	  for(int i = 0; i < mx; ++i){
		  for(int j = 0; j < my; ++j){
			  x1 = a_ + i*hmx; y1 = c_ + j*hmy;
			  tem = fabs(f(x1 + 2*hmx/3, y1 + hmy/3));
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  tem = fabs(f(x1 + hmx/3, y1 + 2*hmy/3));
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
		  }
	  }
	  for(int i = 0; i < mx; ++i){
		for(int j = 0; j < my; ++j){
			  x1 = a_ + i*hmx; y1 = c_ + j*hmy;
			  p1 = L2G(x1, y1);
			  p2 = L2G(x1 + hmx, y1 + hmy);
			  p3 = L2G(x1 + hmx, y1);
			  //cout << i << " " << j << "  " << x1 << " " << y1 << "    " << p1.x() << " " << p1.y() << "  " << p2.x() << " " << p2.y() << "  " << p3.x() << " " << p3.y() << endl;
			  tem = fabs(f(x1 + 2*hmx/3, y1 + hmy/3));
			  paint_triangle(&painter, p1, p2, p3, tem, y_min, y_max);
			  p3 = L2G(x1, y1 + hmy);
			  tem = fabs(f(x1 + hmx/3, y1 + 2*hmy/3));
			  paint_triangle(&painter, p1, p2, p3, tem, y_min, y_max);
		}
	  }
  }
  else if(method == 2){
	  double hmx = (b_-a_)/mx, hmy = (d_-c_)/my;
	  y_min = fabs(f(a_ + 2*hmx/3, c_ + hmy/3) - f_z(a_ + 2*hmx/3, c_ + hmy/3, ag, cg, nx, hx, hy, x));
	  y_max = fabs(f(a_ + 2*hmx/3, c_ + hmy/3) - f_z(a_ + 2*hmx/3, c_ + hmy/3, ag, cg, nx, hx, hy, x));
	  /*
	  for(int i = 0; i <= nx; ++i){
		  for(int j = 0; j <= ny; ++j){
			  ij2l(nx, i, j, l);
			  x1 = a_ + i*hx; y1 = c_ + j*hy;
			  //x1 = a_ + i*hx + hx/3; y1 = c_ + j*hy + 2*hy/3;
			  //x2 = a_ + i*hx + 2*hx/3; y2 = c_ + j*hmy + hmy/3;
			  tem = f(x1, y1) - x[l];
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  //tem = f(x2, y2) - f_z(x2, y2, a_, c_, nx, hx, hy, x);
			  //y_min = min(y_min, tem); y_max = max(y_max, tem);
		  }
	  }
	  */
	  for(int i = 0; i < mx; ++i){
		  for(int j = 0; j < my; ++j){
			  x1 = a_ + i*hmx; y1 = c_ + j*hmy;
			  //x1 = a_ + i*hx + hx/3; y1 = c_ + j*hy + 2*hy/3;
			  //x2 = a_ + i*hx + 2*hx/3; y2 = c_ + j*hmy + hmy/3;
			  tem = fabs(f(x1 + 2*hmx/3, y1 + hmy/3) - f_z(x1 + 2*hmx/3, y1 + hmy/3, ag, cg, nx, hx, hy, x));
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  tem = fabs(f(x1 + hmx/3, y1 + 2*hmy/3) - f_z(x1 + hmx/3, y1 + 2*hmy/3, ag, cg, nx, hx, hy, x));
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  
		  }
	  }
	  for(int i = 0; i < mx; ++i){
		for(int j = 0; j < my; ++j){
			  x1 = a_ + i*hmx; y1 = c_ + j*hmy;
			  p1 = L2G(x1, y1);
			  p2 = L2G(x1 + hmx, y1 + hmy);
			  p3 = L2G(x1 + hmx, y1);
			  //cout << i << " " << j << "  " << x1 << " " << y1 << "    " << p1.x() << " " << p1.y() << "  " << p2.x() << " " << p2.y() << "  " << p3.x() << " " << p3.y() << endl;
			  tem = fabs(f(x1 + 2*hmx/3, y1 + hmy/3) - f_z(x1 + 2*hmx/3, y1 + hmy/3, ag, cg, nx, hx, hy, x));
			  paint_triangle(&painter, p1, p2, p3, tem, y_min, y_max);
			  p3 = L2G(x1, y1 + hmy);
			  tem = fabs(f(x1 + hmx/3, y1 + 2*hmy/3) - f_z(x1 + hmx/3, y1 + 2*hmy/3, ag, cg, nx, hx, hy, x));
			  paint_triangle(&painter, p1, p2, p3, tem, y_min, y_max);
		}
	  }
  }
  else{
	  double hmx = (b_-a_)/mx, hmy = (d_-c_)/my;
	  /*
	  for(size_t i = 0; i < (size_t)(nx+1)*(ny+1); ++i){
		  tem =  x[i];
		  y_min = min(y_min, tem); y_max = max(y_max, tem);
	  }
	  */
	  /*
	  for(int i = 0; i <= nx; ++i){
		  for(int j = 0; j <= ny; ++j){
			  ij2l(nx, i, j, l);
			  x1 = a_ + i*hx; y1 = c_ + j*hy;
			  //x1 = a_ + i*hx + hx/3; y1 = c_ + j*hy + 2*hy/3;
			  //x2 = a_ + i*hx + 2*hx/3; y2 = c_ + j*hmy + hmy/3;
			  tem = x[l];
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  
			  
			  x1 = a_ + i*hmx + hmx/3; y1 = c_ + j*hmy + 2*hmy/3;
			  x2 = a_ + i*hmx + 2*hmx/3; y2 = c_ + j*hmy + hmy/3;
			  tem = f_z(x1, y1, a_, c_, nx, hx, hy, x);
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  tem = f_z(x2, y2, a_, c_, nx, hx, hy, x);
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  
		  }
	  }
	  */
	  for(int i = 0; i < mx; ++i){
		  for(int j = 0; j < my; ++j){
			  x1 = a_ + i*hmx + hmx/3; y1 = c_ + j*hmy + 2*hmy/3;
			  x2 = a_ + i*hmx + 2*hmx/3; y2 = c_ + j*hmy + hmy/3;
			  tem = fabs(f_z(x1, y1, ag, cg, nx, hx, hy, x));
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
			  tem = fabs(f_z(x2, y2, ag, cg, nx, hx, hy, x));
			  y_min = min(y_min, tem); y_max = max(y_max, tem);
		  }
	  }
	  for(int i = 0; i < mx; ++i){
		for(int j = 0; j < my; ++j){
			  x1 = a_ + i*hmx; y1 = c_ + j*hmy;
			  p1 = L2G(x1, y1);
			  p2 = L2G(x1 + hmx, y1 + hmy);
			  p3 = L2G(x1 + hmx, y1);
			  //cout << i << " " << j << "  " << x1 << " " << y1 << "    " << p1.x() << " " << p1.y() << "  " << p2.x() << " " << p2.y() << "  " << p3.x() << " " << p3.y() << endl;
			  tem = fabs(f_z(x1 + 2*hmx/3, y1 + hmy/3, ag, cg, nx, hx, hy, x));
			  paint_triangle(&painter, p1, p2, p3, tem, y_min, y_max);
			  p3 = L2G(x1, y1 + hmy);
			  tem = fabs(f_z(x1 + hmx/3, y1 + 2*hmy/3, ag, cg, nx, hx, hy, x));
			  paint_triangle(&painter, p1, p2, p3, tem, y_min, y_max);
		}
	  }
  }
  
  
  
  
  
  //cout << y_min << " " << y_max << " " << hx << " " << hy << " " << a_ << " " << b_ << " " << width() << endl;
  f_mn = y_min, f_mx = y_max;
  print_text();
  
  update();
  
}
