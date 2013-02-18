/**************************************************************************
 This program calculates the ground state of the GPE in the random
 potential created by laser speckles with square aperturte in 2D.
 Numerical calculations are based on the imaginary-time propagation technique.
 Calculated wavefunction is real.
 The wavefunction is normalized to unity.
 **************************************************************************/

# include <omp.h>

# include <cstdlib>
# include <cstdio>
# include <cstring>

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../NR_MAT.h"
using namespace std;

/////////////////////////////////

//include "nrutil_nr.h"
//redo it !!!!
//try to introduce "FF" numbers


typedef const NRMat<DP> Mat_I_DP;
typedef NRMat<DP> Mat_DP, Mat_O_DP, Mat_IO_DP;
////////////////////////////////






// all parameters are in the file "input2d.h"
#include "input2d.h"
typedef NRMat3d<double> Mat3D_DP;
typedef double DP;

//this file
void speckle_pattern_2d_square(Mat3D_DP &potential, int rstart);
DP k0_fraction(Mat3D_DP &wavefunction);

void normalize_wavefunction_2d(Mat3D_DP &wavefunction);
void initialize_wavefunction_2d(Mat3D_DP &potential, Mat3D_DP &wavefunction);
void Energy(Mat3D_DP &potential, Mat3D_DP &wavefunction, DP *energy, DP *mu,
    DP *ekin);
void calculate_time_step_2d(Mat3D_DP &potential, Mat3D_DP &wavefunction, DP dt);
DP Discrepancy(Mat3D_DP &potential, Mat3D_DP &wavefunction, DP mu);
DP condensate_Bogoliubov_thermodynamic(void);

void initialization_of_wavefunction(int r, Mat3D_DP &potential,
    Mat3D_DP &wavefunction);

DP denmat_2d_last(Mat3D_DP &wavefunction);
void denmat_2d(Mat3D_DP &wavefunction, float dm[]);
DP denmat_2d(Mat3D_DP &wavefunction, int R);
void momentum_distribution_2d(char filename[], Mat3D_DP &func);

void write_function_to_file(char filename[], Mat3D_DP &function);
void write_density_to_file(char filename[], Mat3D_DP &wavefunction);

void read_wavefunction_from_file(char filename[], Mat3D_DP &wavefunction);
void read_wavefunction_from_matlab_output(char filename[],
    Mat3D_DP &wavefunction);

void output_filename(int r, char fname[]);
void matlab_density_filename(int r, char fname[]);
void density_filename(int r, char fname[]);
void potential_filename(int r, char fname[]);

void etot_filename(int r, char fname[]);
void mu_filename(int r, char fname[]);
void ekin_filename(int r, char fname[]);
void k0_filename(int r, char fname[]);

//*********************************************************
// TODO: begin of main

int main() {
  Mat3D_DP wavefunction(1, N_POINTS, N_POINTS);
  Mat3D_DP potential(1, N_POINTS, N_POINTS);

  int realizations[] = RR;
  int i, nr = sizeof(realizations) / sizeof(realizations[0]);

  DP h, cprev, condensate;
  DP dt, eprev, en, mu, ekin, discr;
  DP dprev;

  unsigned long it;
  int rstart, r;
  int counter, count_step;

  FILE *f_out;
  char fname[100], fnameout[100];

  char buf[10];

  printf("sp2dsq_par.c\n");
  printf("\n");
  //////////
  /// TODO:
  /// Undefined "VS"
  ///
#define VS ATOM_NUMBER
  ///
  gcvt(VS, 10, buf);
  printf(" VS = %s", buf);

  gcvt(ATOM_NUMBER, 10, buf);
  printf(" N = %s", buf);

  gcvt(SCATTERING_LENGTH, 10, buf);
  printf(" A = %s", buf);

  printf("\n");

  omp_set_num_threads(nr);
  //omp_set_num_threads(2);

#pragma omp parallel for private(wavefunction,potential,h,cprev,condensate,dt,eprev,en,mu,ekin,dprev,discr,it,rstart,r,counter,count_step,f_out,fname,fnameout)

  for(i = 0; i < nr; i++) {
    r = realizations[i];
    rstart = -r;

    Mat3D_DP wavefunction(1, N_POINTS, N_POINTS);
    Mat3D_DP potential(1, N_POINTS, N_POINTS);

    // create random potential

#pragma omp critical
    {
      speckle_pattern_2d_square(potential, rstart);
    }
    printf("  r = %d\n", r);
    potential_filename(r, fname);
    //write_function_to_file( fname, potential );
    initialization_of_wavefunction(r, potential, wavefunction);

    //h = LENGTH/N_POINTS;
    h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
    dt = h * h / 2;

    output_filename(r, fnameout);

    if (INITIAL_STATE == PREVIOUS) {
      f_out = fopen(fnameout, "a");
      fprintf(f_out, "\n    time step = %e\n", dt);
      fclose(f_out);
    } else {
      f_out = fopen(fnameout, "w");

      if (SCATTERING_LENGTH > 0.) {
        condensate = condensate_Bogoliubov_thermodynamic();
        fprintf(f_out, "    condensate = %f\n\n", condensate);
      }
      fprintf(f_out, " it    E    mu   Ekin   condensate   Discrepancy\n");
      fprintf(f_out, "\n    time step = %e\n\n", dt);
      fclose(f_out);
    }

    Energy(potential, wavefunction, &en, &mu, &ekin);
    //condensate = denmat_2d_last( wavefunction );
    condensate = k0_fraction(wavefunction);
    discr = Discrepancy(potential, wavefunction, mu);

    eprev = en;
    cprev = condensate;
    dprev = discr;

    //imaginary time propagation
    it = 0;
    f_out = fopen(fnameout, "a");
    fprintf(f_out, "%i  %e  %e  %e  %e  %e\n", it, en, mu, ekin, condensate,
        discr);
    fclose(f_out);

    counter = 0;
    count_step = 0;

    for(it = 1; counter < 10 && count_step < 50; it++) {
      calculate_time_step_2d(potential, wavefunction, dt);

      if (!(it % 1000)) {
        Energy(potential, wavefunction, &en, &mu, &ekin);

        if ((en - eprev) > 1.e-6) {
          density_filename(r, fname);
          read_wavefunction_from_file(fname, wavefunction);

          dt /= 2;

          count_step++;

          f_out = fopen(fnameout, "a");
          fprintf(f_out, "\n    time step = %e", dt);
          fprintf(f_out, " ( de = %e )\n\n", eprev - en);
          fclose(f_out);
        } else {
          eprev = en;

          density_filename(r, fname);
          write_density_to_file(fname, wavefunction);

          //condensate = denmat_2d_last( wavefunction );
          condensate = k0_fraction(wavefunction);
          discr = Discrepancy(potential, wavefunction, mu);

          f_out = fopen(fnameout, "a");
          fprintf(f_out, "%i  %e  %e  %e  %e  %e\n", it, en, mu, ekin,
              condensate, discr);
          fclose(f_out);

          //if( (condensate-cprev) > TOLERANCE )
          if ((discr - dprev) > TOLERANCE) {
            dt /= 2;

            count_step++;

            f_out = fopen(fnameout, "a");
            fprintf(f_out, "\n    time step = %e", dt);
            //fprintf(f_out," ( dc = %e )\n\n",cprev-condensate);
            fprintf(f_out, " ( dd = %e )\n\n", dprev - discr);
            fclose(f_out);
          }

          //if( fabs(cprev-condensate) < TOLERANCE )  counter++;
          if (fabs(dprev - discr) < TOLERANCE)
            counter++;
          else {
            counter = 0;
            //cprev = condensate;
            dprev = discr;
          }
        }
      }
    }
    etot_filename(r, fname);
    f_out = fopen(fname, "a");
    fprintf(f_out, "%f %e\n", VS, en);
    fclose(f_out);

    mu_filename(r, fname);
    f_out = fopen(fname, "a");
    fprintf(f_out, "%f %e\n", VS, mu);
    fclose(f_out);

    ekin_filename(r, fname);
    f_out = fopen(fname, "a");
    fprintf(f_out, "%f %e\n", VS, ekin);
    fclose(f_out);

    k0_filename(r, fname);
    f_out = fopen(fname, "a");
    fprintf(f_out, "%f %e\n", VS, condensate);
    fclose(f_out);
  }

  return 0;
}

//*************************************************************************
void initialization_of_wavefunction(int r, Mat3D_DP &potential,
    Mat3D_DP &wavefunction)
//*************************************************************************
{
  char fname[100];
  if (INITIAL_STATE == PREVIOUS) {
    density_filename(r, fname);
    read_wavefunction_from_file(fname, wavefunction);
  } else if (INITIAL_STATE == MATLAB) {
    matlab_density_filename(r, fname);
    read_wavefunction_from_matlab_output(fname, wavefunction);

    density_filename(r, fname);
    write_density_to_file(fname, wavefunction);
  } else {
    initialize_wavefunction_2d(potential, wavefunction);

    density_filename(r, fname);
    write_density_to_file(fname, wavefunction);
  }

}

//*************************************************************************
void etot_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "etot");
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-N_");
  gcvt(ATOM_NUMBER, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-A_");
  gcvt(SCATTERING_LENGTH, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".dat");
}

//*************************************************************************
void mu_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "mu");
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-N_");
  gcvt(ATOM_NUMBER, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-A_");
  gcvt(SCATTERING_LENGTH, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".dat");
}

//*************************************************************************
void ekin_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "ekin");
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-N_");
  gcvt(ATOM_NUMBER, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-A_");
  gcvt(SCATTERING_LENGTH, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".dat");
}

//*************************************************************************
void k0_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "k0");
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-N_");
  gcvt(ATOM_NUMBER, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-A_");
  gcvt(SCATTERING_LENGTH, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".dat");
}

//*************************************************************************
void output_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "sp2dsq_exp");
  strcat(fname, "-VS_");
  gcvt(VS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-N_");
  gcvt(ATOM_NUMBER, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-A_");
  gcvt(SCATTERING_LENGTH, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".out");

}

//*************************************************************************
void matlab_density_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "density");
  strcat(fname, "-VS_");
  gcvt(VS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".dat");

}

//*************************************************************************
void density_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "sp2dsq_den");
  strcat(fname, "-VS_");
  gcvt(VS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-N_");
  gcvt(ATOM_NUMBER, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-A_");
  gcvt(SCATTERING_LENGTH, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".dat");

}

//*************************************************************************
void potential_filename(int r, char fname[])
//*************************************************************************
    {
  char buf[10];

  strcpy(fname, "sp2dsq");
  strcat(fname, "-VS_");
  gcvt(VS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-LC_");
  gcvt(CORRELATION_LENGTH, 10, buf);
  strcat(fname, buf);
  /*strcat(fname,"-L_");
   gcvt(LENGTH,10,buf);
   strcat(fname,buf);*/
  strcat(fname, "-RHO_");
  gcvt(DENSITY, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-NP_");
  gcvt(N_POINTS, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-N_");
  gcvt(ATOM_NUMBER, 10, buf);
  strcat(fname, buf);
  strcat(fname, "-r_");
  gcvt(r, 10, buf);
  strcat(fname, buf);
  strcat(fname, ".dat");

}

//*************************************************************************
void write_wavefunction_to_file(char filename[], Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 filename[] - the name of the file where to write
 wavefunction[1][n][n] - wavefunction

 Output:

 file in the format:
 x y wavefunction

 **************************************************************************/
{
  float h;
  int i, j;
  FILE *f_out;
  f_out = fopen(filename, "w");
  //h = LENGTH/N_POINTS;
  h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  for(i = 0; i < N_POINTS; i++) {
    for(j = 0; j < N_POINTS; j++)
      fprintf(f_out, "%e %e %e\n", h * i, h * j, wavefunction[0][i][j]);
    fprintf(f_out, "%e %e %e\n", h * i, h * N_POINTS, wavefunction[0][i][0]);
  }
  for(j = 0; j < N_POINTS; j++)
    fprintf(f_out, "%e %e %e\n", h * N_POINTS, h * j, wavefunction[0][0][j]);
  fprintf(f_out, "%e %e %e\n", h * N_POINTS, h * N_POINTS,
      wavefunction[0][0][0]);
  fclose(f_out);
}

//*************************************************************************
void write_density_to_file(char filename[], Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 filename[] - the name of the file where to write
 wavefunction[1][n][n] - wavefunction

 Output:

 file in the format:
 x y wavefunction^2

 **************************************************************************/
{
  float h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  int i, j;
  FILE *f_out;
  f_out = fopen(filename, "w");
  for(i = 0; i < N_POINTS; i++) {
    for(j = 0; j < N_POINTS; j++)
      fprintf(f_out, "%e %e %e\n", h * i, h * j,
          wavefunction[0][i][j] * wavefunction[0][i][j]);

    fprintf(f_out, "%e %e %e\n", h * i, h * N_POINTS,
        wavefunction[0][i][0] * wavefunction[0][i][0]);
  }
  for(j = 0; j < N_POINTS; j++)
    fprintf(f_out, "%e %e %e\n", h * N_POINTS, h * j,
        wavefunction[0][0][j] * wavefunction[0][0][j]);
  fprintf(f_out, "%e %e %e\n", h * N_POINTS, h * N_POINTS,
      wavefunction[0][0][0] * wavefunction[0][0][0]);
  fclose(f_out);
}

//*************************************************************************
void read_wavefunction_from_file_c(char filename[], Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 filename[] - the name of the file where to write
 wavefunction[1][n][n] - wavefunction

 Output:

 file in the format:
 x y wavefunction^2

 **************************************************************************/
{
  float den;
  int i, j;
  FILE *f_out;
  f_out = fopen(filename, "r");
  for(i = 0; i < N_POINTS; i++) {
    for(j = 0; j < N_POINTS; j++) {
      fscanf(f_out, "%e", &den);
      fscanf(f_out, "%e", &den);
      fscanf(f_out, "%e", &den);
      wavefunction[0][i][j] = sqrt(den);
    }
    fscanf(f_out, "%e", &den);
    fscanf(f_out, "%e", &den);
    fscanf(f_out, "%e", &den);
  }
  fclose(f_out);
}

//*************************************************************************
void read_wavefunction_from_file(char filename[], Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 filename[] - the name of the file where to write
 wavefunction[1][n][n] - wavefunction

 Output:

 file in the format:
 x y wavefunction^2

 **************************************************************************/
{
  DP den;
  int i, j;
  ifstream f_out;
  f_out.open(filename, ios::in);
  for(i = 0; i < N_POINTS; i++) {
    for(j = 0; j < N_POINTS; j++) {
      f_out >> den;
      f_out >> den;
      f_out >> den;
      wavefunction[0][i][j] = sqrt(den);
    }
    f_out >> den;
    f_out >> den;
    f_out >> den;
  }
  f_out.close();
}

//*************************************************************************
void read_wavefunction_from_matlab_output(char filename[],
    Mat3D_DP &wavefunction)
    /**************************************************************************
     Input:

     filename[] - the name of the file where to write
     wavefunction[1][n][n] - wavefunction

     Output:

     file in the format:
     x y wavefunction^2

     **************************************************************************/
    {
  float den;
  int i, j;
  FILE *f_out;
  f_out = fopen(filename, "r");
  for(i = 0; i < N_POINTS; i++) {
    for(j = 0; j < N_POINTS; j++) {
      fscanf(f_out, "%e", &den);
      fscanf(f_out, "%e", &den);
      fscanf(f_out, "%e", &den);
      wavefunction[0][i][j] = sqrt(den);
    }
  }
  fclose(f_out);
}

//*************************************************************************
void write_function_to_file(char filename[], Mat3D_DP &function)
/**************************************************************************
 Input:

 filename[] - the name of the file where to write
 wavefunction[1][n][n] - wavefunction

 Output:

 file in the format:
 x y wavefunction

 **************************************************************************/
{
  float h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  int i, j;
  FILE *f_out;
  f_out = fopen(filename, "w");

  for(i = 0; i < N_POINTS; i++) {
    for(j = 0; j < N_POINTS; j++)
      fprintf(f_out, "%e %e %e\n", h * i, h * j, function[0][i][j]);
    fprintf(f_out, "%e %e %e\n", h * i, h * N_POINTS, function[0][i][0]);
  }

  for(j = 0; j < N_POINTS; j++)
    fprintf(f_out, "%e %e %e\n", h * N_POINTS, h * j, function[0][0][j]);
  fprintf(f_out, "%e %e %e\n", h * N_POINTS, h * N_POINTS, function[0][0][0]);
  fclose(f_out);
}

//*************************************************************************
DP condensate_Bogoliubov_thermodynamic(void)
//*************************************************************************
    {
  DP en, condensate;

  en = 2. / M_PI * sqrt(2. / M_PI) * SCATTERING_LENGTH * CORRELATION_LENGTH
      * CORRELATION_LENGTH * ATOM_NUMBER / (ATOM_NUMBER / DENSITY);
  condensate = 2 * sqrt(1. + en) * atanf(1. / sqrt(1 + en));
  en = sqrt(en);
  condensate -= 2 * en * atanf(1. / en);
  condensate -= en * en
      * (logf(en) - logf(1. + en * en) + 0.5 * logf(2. + en * en));
  en *= VS * (ATOM_NUMBER / DENSITY)
      / (2. * sqrt(2. * M_PI) * SCATTERING_LENGTH * ATOM_NUMBER);
  condensate *= en * en;
  condensate *= 0.5;
  condensate = 1. - condensate;

  return condensate;
}



//TODO: Laplace is applied here.
//*************************************************************************
void Energy(Mat3D_DP &potential, Mat3D_DP &wavefunction, DP *energy, DP *mu,
    DP *ekin)
    /**************************************************************************
     Input:

     wavefunction[n][n][n] - wavefunction

     Output:

     ground-state energy

     **************************************************************************/
    {
  DP h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  DP nonlin = sqrt(2 * M_PI) * SCATTERING_LENGTH * (ATOM_NUMBER - 1);
  DP Ekin, Epot, Eint, laplace, wf2;

  int i, j, k, i1;

  //potential energy

  Epot = 0.;
  Eint = 0.;

  for(i = 0; i <= 0; i++)
    for(j = 0; j < N_POINTS; j++)
      for(k = 0; k < N_POINTS; k++) {
        wf2 = wavefunction[i][j][k] * wavefunction[i][j][k];
        Epot += potential[i][j][k] * wf2;
        Eint += nonlin * wf2 * wf2;
      }

  Epot *= h * h;
  Eint *= h * h;

  //kinetic energy

  Ekin = 0.;

  for(i = 0; i <= 0; i++)
    for(j = 0; j < N_POINTS; j++)
      for(k = 0; k < N_POINTS; k++) {
        i1 = j - 1;
        if (i1 == -1)
          i1 = N_POINTS - 1;
        laplace = wavefunction[i][i1][k];

        i1 = j + 1;
        if (i1 == N_POINTS)
          i1 = 0;
        laplace += wavefunction[i][i1][k];

        i1 = k - 1;
        if (i1 == -1)
          i1 = N_POINTS - 1;
        laplace += wavefunction[i][j][i1];

        i1 = k + 1;
        if (i1 == N_POINTS)
          i1 = 0;
        laplace += wavefunction[i][j][i1];

        laplace -= 4 * wavefunction[i][j][k];

        Ekin += wavefunction[i][j][k] * laplace;
      }

  Ekin *= (-0.5);

  *energy = Ekin + Epot + Eint;
  *mu = Ekin + Epot + 2 * Eint;
  *ekin = Ekin;
}



// TODO: Laplace is applied here
//*************************************************************************
void calculate_time_step_2d(Mat3D_DP &potential, Mat3D_DP &wavefunction, DP dt)
/**************************************************************************
 Input:

 potential[n][n] - external potential
 wavefunction[1][n][n] - wavefunction(t)

 Output:

 wavefunction[1][n][n] - wavefunction(t+dt)

 **************************************************************************/
{
  Mat3D_DP wavefunction_prev(1, N_POINTS, N_POINTS);
  int i, j, k, i1;
  DP h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  DP nonlin = 2 * sqrt(2 * M_PI) * SCATTERING_LENGTH * (ATOM_NUMBER - 1);
  DP laplace;

  for(i = 0; i <= 0; i++)
    for(j = 0; j < N_POINTS; j++)
      for(k = 0; k < N_POINTS; k++)
        wavefunction_prev[i][j][k] = wavefunction[i][j][k];

  for(i = 0; i <= 0; i++)
    for(j = 0; j < N_POINTS; j++)
      for(k = 0; k < N_POINTS; k++) {
        i1 = j - 1;
        if (i1 == -1)
          i1 = N_POINTS - 1;
        laplace = wavefunction_prev[i][i1][k];

        i1 = j + 1;
        if (i1 == N_POINTS)
          i1 = 0;
        laplace += wavefunction_prev[i][i1][k];

        i1 = k - 1;
        if (i1 == -1)
          i1 = N_POINTS - 1;
        laplace += wavefunction_prev[i][j][i1];

        i1 = k + 1;
        if (i1 == N_POINTS)
          i1 = 0;
        laplace += wavefunction_prev[i][j][i1];

        laplace -= 4 * wavefunction[i][j][k];

        wavefunction[i][j][k] -= dt
            * (-laplace / (2 * h * h)
                + nonlin * wavefunction[i][j][k] * wavefunction[i][j][k]
                    * wavefunction[i][j][k]
                + potential[i][j][k] * wavefunction[i][j][k]);
      }

  normalize_wavefunction_2d(wavefunction);
}



//TODO: Laplace is applied here; matrix size = N_POINTS
//*************************************************************************
DP Discrepancy(Mat3D_DP &potential, Mat3D_DP &wavefunction, DP mu)
/**************************************************************************
 Input:

 potential[n][n] - external potential
 wavefunction[1][n][n] - wavefunction(t)

 Output:

 wavefunction[1][n][n] - wavefunction(t+dt)

 **************************************************************************/
{
  int i, j, k, i1;
  DP a, norm;
  DP h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  DP nonlin = 2 * sqrt(2 * M_PI) * SCATTERING_LENGTH * (ATOM_NUMBER - 1);
  DP laplace;

  norm = 0.;

  for(i = 0; i <= 0; i++)
    for(j = 0; j < N_POINTS; j++)
      for(k = 0; k < N_POINTS; k++) {
        i1 = j - 1;
        if (i1 == -1)
          i1 = N_POINTS - 1;
        laplace = wavefunction[i][i1][k];

        i1 = j + 1;
        if (i1 == N_POINTS)
          i1 = 0;
        laplace += wavefunction[i][i1][k];

        i1 = k - 1;
        if (i1 == -1)
          i1 = N_POINTS - 1;
        laplace += wavefunction[i][j][i1];

        i1 = k + 1;
        if (i1 == N_POINTS)  i1 = 0;
        laplace += wavefunction[i][j][i1];
        laplace -= 4 * wavefunction[i][j][k];
        laplace /= (h * h);

        a = -laplace / 2
            + nonlin * wavefunction[i][j][k] * wavefunction[i][j][k]
                * wavefunction[i][j][k]
            + (potential[i][j][k] - mu) * wavefunction[i][j][k];

        norm += a * a;
      }
  return sqrt(norm);
}

//*************************************************************************
void initialize_wavefunction_2d(Mat3D_DP &potential, Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 potential[1][n][n] - external potential

 Output:

 wavefunction[1][n][n]
 initial guess for the wavefunction taken as a square root of the inverted potential

 **************************************************************************/
{
  DP max;
  int i, j;
  max = 0.;
  for(i = 0; i < N_POINTS; i++)
    for(j = 0; j < N_POINTS; j++)
      if (potential[0][i][j] > max)
        max = potential[0][i][j];
  max = 0.;
  if (max == 0.)
    for(i = 0; i < N_POINTS; i++)
      for(j = 0; j < N_POINTS; j++)
        wavefunction[0][i][j] = 1.;
  else
    for(i = 0; i < N_POINTS; i++)
      for(j = 0; j < N_POINTS; j++)
        wavefunction[0][i][j] = sqrt(max - potential[0][i][j]);
  normalize_wavefunction_2d(wavefunction);
}



// TODO: This is used for GP groundstate.
//*************************************************************************
void normalize_wavefunction_2d(Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 Output:

 normalized wavefunction
 **************************************************************************/
{
  float h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  DP norm;
  int i, j, k;

  norm = 0.;

  for(i = 0; i <= 0; i++)
    for(j = 0; j < N_POINTS; j++)
      for(k = 0; k < N_POINTS; k++)
        norm += wavefunction[i][j][k] * wavefunction[i][j][k];

  norm *= h * h;

  norm = sqrt(norm);

  for(i = 0; i <= 0; i++)
    for(j = 0; j < N_POINTS; j++)
      for(k = 0; k < N_POINTS; k++)
        wavefunction[i][j][k] /= norm;

}

//*************************************************************************
DP denmat_2d(Mat3D_DP &wavefunction, int R)
/**************************************************************************
 Input:

 Output:

 normalized wavefunction
 **************************************************************************/
{
  DP denmat;
  int i, j, k;
  int jpR, jmR, kpR, kmR;

  denmat = 0.;

  for(i = 1; i <= 1; i++) {
    for(j = 1; j <= N_POINTS; j++) {
      jpR = j + R;
      if (jpR > N_POINTS)
        jpR -= N_POINTS;

      jmR = j - R;
      if (jmR < 1)
        jmR += N_POINTS;

      for(k = 1; k <= N_POINTS; k++) {
        kpR = k + R;
        if (kpR > N_POINTS)
          kpR -= N_POINTS;

        kmR = k - R;
        if (kmR < 1)
          kmR += N_POINTS;

        denmat += wavefunction[i][j][k]
            * (wavefunction[i][jpR][k] + wavefunction[i][jmR][k]
                + wavefunction[i][j][kpR] + wavefunction[i][j][kmR]);
      }
    }
  }

  denmat /= 4.;
  denmat /= (N_POINTS * N_POINTS);
  denmat *= (ATOM_NUMBER / DENSITY);

  return (denmat);
}

//*************************************************************************
void denmat_2d(Mat3D_DP &wavefunction, float dm[])
/**************************************************************************
 Input:

 Output:
 **************************************************************************/
{
  Vec_INT counts(N_POINTS / 2);
  int i, j;
  int di, dj, ipdi, jpdj, R;

  for(R = 0; R < N_POINTS / 2; R++)
    dm[R + 1] = 0.;

  for(R = 0; R < N_POINTS / 2; R++)
    counts[R + 1] = 0;

  for(i = 1; i <= N_POINTS; i++)
    for(j = 1; j <= N_POINTS; j++) {
      for(di = -N_POINTS / 2 + 1; di < N_POINTS / 2; di++) {
        ipdi = i + di;
        if (ipdi > N_POINTS)
          ipdi -= N_POINTS;
        if (ipdi < 1)
          ipdi += N_POINTS;

        for(dj = -N_POINTS / 2 + 1; dj <= N_POINTS; dj++) {
          jpdj = j + dj;
          if (jpdj > N_POINTS)
            jpdj -= N_POINTS;
          if (jpdj < 1)
            jpdj += N_POINTS;

          R = (int) round(sqrt((DP) di * di + dj * dj));

          if (R < N_POINTS / 2) {
            counts[R + 1]++;
            dm[R + 1] += wavefunction[1][i][j] * wavefunction[1][ipdi][jpdj];
          }
        }
      }
    }

  for(R = 0; R < N_POINTS / 2; R++)
    dm[R + 1] /= counts[R + 1];

  for(R = 0; R < N_POINTS / 2; R++)
    dm[R + 1] *= (ATOM_NUMBER / DENSITY) / (N_POINTS * N_POINTS);
}

//*************************************************************************
void momentum_distribution_2d(char filename[], Mat3D_DP &func)
/**************************************************************************
 Input:

 filename[] - the name of the file where to write
 wavefunction[1][n][n] - wavefunction

 Output:

 file in the format:
 x y wavefunction

 **************************************************************************/
{
  Mat_DP speq(1, 2 * N_POINTS);
  float h, f2;
  int i, j;

  FILE *fout;

  h = 2 * M_PI / sqrt(ATOM_NUMBER / DENSITY);

  NR::rlft3(func, speq, 1);

  //printf(" %e\n",func[1][1][1]/N_POINTS/N_POINTS);
  fout = fopen(filename, "w");

  for(i = 0; i <= (N_POINTS / 2); i++) {
    for(j = 0; j < (N_POINTS / 2); j++) {
      f2 = func[1][i + 1][2 * j + 1] * func[1][i + 1][2 * j + 1]
          + func[1][i + 1][2 * j + 2] * func[1][i + 1][2 * j + 2];

      f2 *= (ATOM_NUMBER / DENSITY) / N_POINTS / N_POINTS / N_POINTS / N_POINTS;
      fprintf(fout, "%e %e %e\n", h * i, h * j, f2);
      if (i || j) fprintf(fout, "%e %e %e\n", -h * i, -h * j, f2);
    }
  }

  for(i = -N_POINTS / 2 + 1; i <= -1; i++) {
    for(j = 0; j < (N_POINTS / 2); j++) {
      f2 = func[1][N_POINTS + i + 1][2 * j + 1]
          * func[1][N_POINTS + i + 1][2 * j + 1]
          + func[1][N_POINTS + i + 1][2 * j + 2]
              * func[1][N_POINTS + i + 1][2 * j + 2];
      f2 *= (ATOM_NUMBER / DENSITY) / N_POINTS / N_POINTS / N_POINTS / N_POINTS;
      fprintf(fout, "%e %e %e\n", h * i, h * j, f2);
      if (j) fprintf(fout, "%e %e %e\n", -h * i, -h * j, f2);
    }
  }

  // Nyquist critical frequency
  j = N_POINTS / 2;

  for(i = 0; i <= (N_POINTS / 2); i++) {
    f2 = speq[1][2 * i + 1] * speq[1][2 * i + 1]
        + speq[1][2 * i + 2] * speq[1][2 * i + 2];

    f2 *= (ATOM_NUMBER / DENSITY) / N_POINTS / N_POINTS / N_POINTS / N_POINTS;
    fprintf(fout, "%e %e %e\n", h * i, h * j, f2);
    fprintf(fout, "%e %e %e\n", -h * i, -h * j, f2);
  }

  for(i = -N_POINTS / 2 + 1; i <= -1; i++) {
    f2 = speq[1][2 * (N_POINTS + i) + 1] * speq[1][2 * (N_POINTS + i) + 1]
        + speq[1][2 * (N_POINTS + i) + 2] * speq[1][2 * (N_POINTS + i) + 2];

    f2 *= (ATOM_NUMBER / DENSITY) / N_POINTS / N_POINTS / N_POINTS / N_POINTS;
    fprintf(fout, "%e %e %e\n", h * i, h * j, f2);
    fprintf(fout, "%e %e %e\n", -h * i, -h * j, f2);
  }
  fclose(fout);
  NR::rlft3(func, speq, -1);
}

//*************************************************************************
void speckle_pattern_2d_square(Mat3D_DP &potential, int rstart)
/**************************************************************************
 Input:

 n - number of spatial points (must be power of 2)
 L - size of the system
 Lc - correlation length
 Vs - mean value of the speckle potential
 rstart - starting value for the random number generator (must be negative)

 Output:

 The values of the speckle potential will be in the two-dimensional array
 potential[1][n][n].
 **************************************************************************/
{
  const int dim = 2;
  Vec_INT nn(dim);
  Vec_DP field(2 * N_POINTS * N_POINTS);

  unsigned long i, j, imin, imax;
  DP norm;
  int idum = rstart;

  nn[0] = nn[1] = N_POINTS;

  for(i = 0; i < (2 * N_POINTS * N_POINTS); i++)
    field[i] = NR::gasdev2(idum);

  NR::fourn(field, nn, 1);

  /*imin = (unsigned long)ceil(L/(2*Lc));
   imax = (unsigned long)floor(n-L/(2*Lc));*/

  imin = (int) ceil(sqrt(ATOM_NUMBER / DENSITY) / (2 * CORRELATION_LENGTH));
  imax = (int) floor(
      N_POINTS - sqrt(ATOM_NUMBER / DENSITY) / (2 * CORRELATION_LENGTH));

  //for( i=imin ; i <= imax ; i++ )  field[2*i]=field[2*i+1]=0.;

  for(i = 0; i < imin; i++)
    for(j = imin; j <= imax; j++)
      field[2 * (i * N_POINTS + j)] = field[2 * (i * N_POINTS + j) + 1] = 0.;

  for(i = imin; i <= imax; i++)
    for(j = 0; j < N_POINTS; j++)
      field[2 * (i * N_POINTS + j)] = field[2 * (i * N_POINTS + j) + 1] = 0.;

  for(i = imax + 1; i < N_POINTS; i++)
    for(j = imin; j <= imax; j++)
      field[2 * (i * N_POINTS + j)] = field[2 * (i * N_POINTS + j) + 1] = 0.;

  NR::fourn(field, nn, -1);

  norm = 0.;

  for(i = 0; i < (N_POINTS * N_POINTS); i++)
    norm += (field[2 * i] * field[2 * i] + field[2 * i + 1] * field[2 * i + 1]);

  norm /= (N_POINTS * N_POINTS);

  for(i = 0; i < N_POINTS; i++)
    for(j = 0; j < N_POINTS; j++)
      potential[0][i][j] = VS
          * (field[2 * (i * N_POINTS + j)] * field[2 * (i * N_POINTS + j)]
              + field[2 * (i * N_POINTS + j) + 1]
                  * field[2 * (i * N_POINTS + j) + 1]) / norm;
}

//*************************************************************************
DP denmat_2d_last(Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 Output:
 **************************************************************************/
{
  DP dm;
  unsigned long counts;
  int i, j;
  int di, dj, ipdi, jpdj, R2;

  R2 = N_POINTS / 2 - 1;
  R2 *= R2;

  dm = 0.;

  counts = 0;

  for(di = -N_POINTS / 2 + 1; di < N_POINTS / 2; di++) {
    dj = (int) floor(sqrt(R2 - di * di));

    for(i = 0; i < N_POINTS; i++) {
      ipdi = i + di;
      if (ipdi >= N_POINTS)
        ipdi -= N_POINTS;
      else if (ipdi < 0)
        ipdi += N_POINTS;

      for(j = 0; j < N_POINTS; j++) {
        jpdj = j + dj;
        if (jpdj >= N_POINTS)
          jpdj -= N_POINTS;

        dm += wavefunction[0][i][j] * wavefunction[0][ipdi][jpdj];

        if (dj) {
          jpdj = j - dj;
          if (jpdj < 0)
            jpdj += N_POINTS;

          dm += wavefunction[0][i][j] * wavefunction[0][ipdi][jpdj];
        }
      }
    }

    counts += N_POINTS * N_POINTS;

    if (dj)
      counts += N_POINTS * N_POINTS;
  }

  dm /= counts;

  dm *= (ATOM_NUMBER / DENSITY);

  return dm;
}

//*************************************************************************
DP k0_fraction(Mat3D_DP &wavefunction)
/**************************************************************************
 Input:

 Output:
 **************************************************************************/
{
  DP sum;
  DP h = sqrt(ATOM_NUMBER / DENSITY) / N_POINTS;
  int i, j;

  sum = 0.;

  for(i = 0; i < N_POINTS; i++)
    for(j = 0; j < N_POINTS; j++)
      sum += wavefunction[0][i][j];

  return (sum * sum * h * h / (N_POINTS * N_POINTS));
}
