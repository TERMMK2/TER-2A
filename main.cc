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

  double Tfinal, deltaT;
  Tfinal = data_file.Get_T_final();
  deltaT = data_file.Get_deltaT();
  int nb_iterations = int(ceil(Tfinal/deltaT));

  Laplacian2D *Lap;

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
    if (data_file.Get_Schema() == "Implicite")
    {
      Lap->IterativeSolver(nb_iterations);
    }
    else if (data_file.Get_Schema() == "Explicite")
    {
      Lap->Advance(nb_iterations);
    }
  }
  if (data_file.Get_eq() == "EC_PyrolyseMV")
  {
    Lap = new EC_PyrolyseMV();
    Lap->Initialize(data_file);
    if (data_file.Get_Schema() == "Implicite")
    {
      Lap->IterativeSolver(nb_iterations);
    }
    else if (data_file.Get_Schema() == "Explicite")
    {
      // cout << "Explicite n'est pas encore codé tête de gland" << endl;
      // abort();
      Lap->Advance(nb_iterations);
    }
  }
  delete Lap;
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
