#include "Laplacian2D.h"
#include <iostream>
#include <fstream>
#include <chrono>


using namespace Eigen;
using namespace std;

int main()
{
  int Nx = 200;
  int Ny = 2;
  double xmin = 0.;
  double xmax = 0.01;
  double ymin = 0.;
  double ymax = 0.005;
  double a = 1./(1500*1000);
  double deltaT = 0.01;
  double tfinal = 10;
  Eigen::VectorXd CI;
  string CL_bas = "Neumann"; // "Neumann" , "Dirichlet"
  string CL_haut = "Neumann";
  string CL_gauche = "Dirichlet";
  string CL_droite = "Neumann";
  double Val_CL_bas = 0; //Flux si CL_bas == "Neumann", Température si CL_bas == "Dirichlet"
  double Val_CL_haut = 0;
  double Val_CL_gauche = 2900;
  double Val_CL_droite = 0;
  int nb_iterations = int(ceil(tfinal/deltaT));
  string Equation = "EC_ClassiqueP";

  CI.resize(Nx*Ny);
  for(int j=0; j < Ny; j++)
    {
      for(int i=0; i < Nx; i++)
  	    {
  	      CI(i + j*Nx) = 293.;
  	    }
    }

  Laplacian2D *Lap;

  int cas = 0;
  if (Equation == "EC_ClassiqueM")
  {
    cas = 1;
  }
  if (Equation == "EC_ClassiqueP")
  {
    cas = 2;
  }

  switch(cas)
  {
    case 1:
      Lap = new EC_ClassiqueM();
      Lap->Initialize(xmin,xmax,ymin,ymax,Nx,Ny,a,deltaT);
      Lap->InitializeCI(CI);
      Lap->InitializeCL(CL_bas, CL_haut, CL_gauche, CL_droite, Val_CL_bas, Val_CL_haut, Val_CL_gauche, Val_CL_droite);
      Lap->InitializeMatrix();
      Lap->DirectSolver(nb_iterations);
      break;

    case 2:
      Lap = new EC_ClassiqueP();
      Lap->Initialize(xmin,xmax,ymin,ymax,Nx,Ny,a,deltaT);
      Lap->InitializeCI(CI);
      Lap->InitializeCL(CL_bas, CL_haut, CL_gauche, CL_droite, Val_CL_bas, Val_CL_haut, Val_CL_gauche, Val_CL_droite);
      Lap->InitializeMatrix();
      Lap->DirectSolver(nb_iterations);
      break;

    default:
      std::cout << "Ce choix n'est pas disponible" << std::endl;
      exit(0);
  }
  return 0;
  }

// Pas trop pratique -> Faire vite def un prod mat/Vect pour Eul_exp et utiliser le solveur pour Eul_imp -> voir comment on va faire en propre aussi




  /*
  auto start = chrono::high_resolution_clock::now();
  Lap.DirectSolver();
  auto finish = chrono::high_resolution_clock::now();
  Lap.SaveSol("sol_dir.dat");
  err = Lap.GetError();
  double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
  cout << "l'erreur pour le solveur direct est : " << err << " et il a mis "<< t << " microsecondes a s'effectuer" << endl;


  start = chrono::high_resolution_clock::now();
  Lap.IterativeSolver();
  finish = chrono::high_resolution_clock::now();
  Lap.SaveSol("sol_ite.dat");
  err = Lap.GetError();
  t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
  cout << "l'erreur pour le solveur itératif est : " << err << " et il a mis "<< t << " microsecondes a s'effectuer"<< endl;



  */
