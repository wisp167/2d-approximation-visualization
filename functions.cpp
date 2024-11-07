#include <stdio.h>
#include <window.h>

double f_0 (double x, double y)
{
  x*=10;
  y*=10;
  return 1;
}
double f_1 (double x, double y)
{
  y*=10;
  return x;
}
double f_2 (double x, double y)
{
  x*=10;
  return y;
}
double f_3 (double x, double y)
{
  return x+y;
}
double f_4 (double x, double y)
{
  return sqrt(x*x+y*y);
}
double f_5 (double x, double y)
{
  return x*x+y*y;
}
double f_6 (double x, double y)
{
  return exp(x*x-y*y);
}
double f_7(double x, double y){
  return 1/(25*(x*x+y*y)+1);
}
/*
double df_0 (double x)
{
  x*=10;
  return 0;
}
double df_1 (double x)
{
  x*=10;
  return 1;
}
double df_2 (double x)
{
  return 2*x;
}
double df_3 (double x)
{
  return 3*x*x;
}
double df_4 (double x)
{
  return 4*x*x*x;
}
double df_5 (double x)
{
  return exp(x);
}
double df_6 (double x)
{
  return (50*x)/((25*x*x+1)*(25*x*x+1));
}
*/
