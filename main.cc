#include "Datafile.h"
#include "Laplacian2D.h"
#include <iostream>
#include <fstream>
#include <chrono>



using namespace Eigen;
using namespace std;

int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Please, enter the name of your data file." << endl;
      abort();
    }

  // ----------------------- Fichier de données --------------------------------
  DataFile data_file(argv[1]);

  // Lecture du fichier de données
  data_file.ReadDataFile();

  // int Nx = 200;
  // int Ny = 2;
  // double xmin = 0.;
  // double xmax = 0.01;
  // double ymin = 0.;
  // double ymax = 0.005;
  // double a = 1./(1500*1000);
  // double deltaT = 0.01;
  // double tfinal = 10;
  // Eigen::VectorXd CI;
  // string CL_bas = "Neumann"; // "Neumann" , "Dirichlet"
  // string CL_haut = "Neumann";
  // string CL_gauche = "Dirichlet";
  // string CL_droite = "Neumann";
  // double Val_CL_bas = 0; //Flux si CL_bas == "Neumann", Température si CL_bas == "Dirichlet"
  // double Val_CL_haut = 0;
  // double Val_CL_gauche = 2900;
  // double Val_CL_droite = 0;
  // string Equation = "EC_ClassiqueP";

  double Tfinal, deltaT;
  Tfinal = data_file.Get_T_final();
  deltaT = data_file.Get_deltaT();
  int nb_iterations = int(ceil(Tfinal/deltaT));

  Laplacian2D *Lap;


  // Lap = new EC_PyrolyseMC();
  // Lap->Initialize(data_file);
  // Lap->Advance(nb_iterations);


  if (data_file.Get_eq() == "EC_ClassiqueM")
    {
      Lap = new EC_ClassiqueM();
      Lap->Initialize(data_file);
      Lap->InitializeMatrix();
      Lap->IterativeSolver(nb_iterations);
    }
  if (data_file.Get_eq() == "EC_ClassiqueP")
     {
       Lap = new EC_ClassiqueP();
       Lap->Initialize(data_file);
       Lap->InitializeMatrix();
       cout<<"bonjour"<<endl;
       auto start = chrono::high_resolution_clock::now();
       Lap->IterativeSolver(nb_iterations);
       auto finish = chrono::high_resolution_clock::now();
       double t = chrono::duration_cast<chrono::microseconds>(finish-start).count();
       cout << "Le prog a mis " << t*0.000001 << " secondes a s'effectuer" << endl;
     }
  if (data_file.Get_eq() == "EC_PyrolyseMC")
     {
       Lap = new EC_PyrolyseMC();
       Lap->Initialize(data_file);
       Lap->IterativeSolver(nb_iterations);
     }



  //   // On aurait pu faire juste avec les ifs mais bon:
  //   switch(cas)
  //   {
  //     case 1:
  //       Lap = new EC_ClassiqueM();
  //       Lap->Initialize(data_file);
  //       Lap->InitializeMatrix();
  //       Lap->DirectSolver(nb_iterations);
  //       break;

  //     case 2:
  //       Lap = new EC_ClassiqueP();
  //       Lap->Initialize(data_file);
  //       Lap->InitializeMatrix();
  //       Lap->DirectSolver(nb_iterations);
  //       break;

  //   case 3:
  //     Lap = new EC_Pyrolyse();
  //     Lap->Initialize(data_file);
  //     Lap->Advance(nb_iterations);
  //     break;

  //     default:
  //       std::cout << "Ce choix n'est pas disponible" << std::endl;
  //       exit(0);
  //   }
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
