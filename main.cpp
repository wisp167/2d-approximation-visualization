
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QCloseEvent>
#include <QObject>

#include "window.h"

void foo::change_stat(){
	stat = true;
	close();
}
void foo::closeEvent(QCloseEvent *bar)
{
    // Do something
    if(stat){
		cout << "CLOSING APP" << endl;
		bar->accept();
		return;
	}
    cout << "PREPARING FOR CLOSING APP" << endl;
    emit signalName();
    bar->ignore();
}

int main (int argc, char *argv[])
{
  feenableexcept(FE_ALL_EXCEPT ^ FE_INEXACT);
  QApplication app (argc, argv);
  foo *window = new foo;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window);
  QAction *action;

  if (graph_area->parse_command_line (argc, argv))
    {
      QMessageBox::warning (0, "Wrong input arguments!", 
                            "Wrong input arguments!");
      return -1;
    }
	
	
	
  int* temp = (int*)malloc (sizeof (int) * 5);
  double* ttemp = (double*) malloc(sizeof(double)*5);
  if(read_input(temp, ttemp, argv)){
  	err_output(argv[0]);
  	cout << read_input(temp, ttemp, argv) << endl;
  	free(temp); free(ttemp);
  	return -1;
  }
  double a, b, c, d, eps;
  int nx, ny, k, m, p;
  a = ttemp[0]; b = ttemp[1]; c = ttemp[2]; d = ttemp[3]; eps = ttemp[4];
  nx = temp[0]; ny = temp[1]; k = temp[2]; m = temp[3]; p = temp[4];
  free(temp); free(ttemp);
  
  
  
  init_reduce_sum(p);
  //size_t N = (size_t)(nx+1)*(ny+1);
  //static double *A = nullptr, *B = nullptr, *x = nullptr, *x_copy = nullptr, *u = nullptr, *v =nullptr, *r = nullptr, *temp1 = nullptr;
  //size_t *I = nullptr;
  Args *args = nullptr;
  pthread_t *tid = nullptr;
  try{
		args = new Args[p];
		tid = new pthread_t[p];
		/*A = (double*)malloc (sizeof (double) * (get_len_msr(nx,ny)+1));
		I = (size_t*)malloc(sizeof(size_t)*(get_len_msr(nx,ny)+1));
		B = (double*)malloc(sizeof(double) * N);
		x = (double*)malloc(sizeof(double) * N);
		x_copy = (double*) malloc(sizeof(double) * N);
		r = (double*)malloc(sizeof(double) * N);
		u = (double*)malloc(sizeof(double) * N);
		v = (double*)malloc(sizeof(double) * N);
		temp1 = (double*) malloc(sizeof(double)*N);
		*/
  }
  catch(std::bad_alloc&){
	err_output(argv[0]);
	delete_reduce_sum();
 	return -3;
  }
  
  for(int i = 0; i < p; ++i){
	    args[i].name = argv[0];
		args[i].p = p;
		args[i].k = k;
		args[i].kk = i;
		args[i].a = a;
		args[i].b = b;
		args[i].c = c;
		args[i].d = d;
		args[i].nx = nx;
		args[i].ny = ny;
		args[i].m = m;
		args[i].eps = eps;
		args[i].stat = true;
		/*
		args[i].A = A;
		args[i].I = I;
		args[i].B = B;
		args[i].x = x;
		args[i].u = u;
		args[i].v = v;
		args[i].r = r;
		args[i].x_copy = x_copy;
		args[i].temp = temp1;
		*/
  }
  for(int i = 0; i < p; ++i){
	if(pthread_create(tid+i, 0, thread_window_f, args+i)){
			abort();
	}
  }
  if(!args[0].finish){
	  action = tool_bar->addAction ("&Change function", graph_area, SLOT (change_func ()));
	  action->setShortcut (QString ("0"));

	  action = tool_bar->addAction ("&Change mode", graph_area, SLOT (change_mode ()));
	  action->setShortcut (QString ("1"));
	  
	  action = tool_bar->addAction ("&Increase scale", graph_area, SLOT (increase_scale ()));
	  action->setShortcut (QString ("2"));
	 
	  action = tool_bar->addAction ("&Decrease scale", graph_area, SLOT (decrease_scale ()));
	  action->setShortcut (QString ("3"));

	  action = tool_bar->addAction ("&Increase number of points", graph_area, SLOT (increase_n ()));
	  action->setShortcut (QString ("4"));
	  
	  action = tool_bar->addAction ("&Decrease number of points", graph_area, SLOT (decrease_n ()));
	  action->setShortcut (QString ("5"));
	  
	  action = tool_bar->addAction ("&+delta", graph_area, SLOT (increase_p ()));
	  action->setShortcut (QString ("6"));
	  
	  action = tool_bar->addAction ("&-delta", graph_area, SLOT (decrease_p ()));
	  action->setShortcut (QString ("7"));
	  
	  action = tool_bar->addAction ("&Increase number of points m", graph_area, SLOT (increase_m ()));
	  action->setShortcut (QString ("8"));
	  
	  action = tool_bar->addAction ("&Decrease number of points m", graph_area, SLOT (decrease_m ()));
	  action->setShortcut (QString ("9"));
	  
	  QObject::connect(window, SIGNAL(signalName()), graph_area, SLOT(preparation_for_finishing()));
	  QObject::connect(graph_area, SIGNAL(signalName1()), window, SLOT(change_stat()));

	  tool_bar->setMaximumHeight (30);

	  window->setMenuBar (tool_bar);
	  window->setCentralWidget (graph_area);
	  window->setWindowTitle ("Graph");

	  window->show ();
	  app.exec ();
	  //delete window;
  }
  for(int i  = 0; i <p; ++i){
  	pthread_join(tid[i], 0);
  }
  delete_reduce_sum();
  delete[] args; delete[] tid;
  delete window;
  //free(I); free(A); free(B); free(x); free(r); free(u); free(v); free(temp1); free(x_copy);
  return 0;
}
