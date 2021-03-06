#include "Sparse"
#include "Dense"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <memory>
#include "Datafile.h"

class Laplacian2D // pas fini de modifier
{
 protected: // Les attributs de la classe

  double _x_min, _x_max, _y_min, _y_max, _h_x, _h_y, _a, _deltaT;
  int _Nx, _Ny;
  Eigen::SparseMatrix<double> _LapMat; // matrice creuse du laplacien
  Eigen::VectorXd _f; // vecteur source _f qui prend les données de _sol(i) pour calculer _sol(i+1)
  Eigen::VectorXd _sol; // vecteur solution U
  std::string _CL_bas, _CL_haut, _CL_gauche, _CL_droite;
  double _Val_CL_bas, _Val_CL_haut, _Val_CL_gauche, _Val_CL_droite;

  std::string _Solveur;

  std::string _save_all_file, _save_points_file, _restart_file;
  int _number_saved_points;
  std::vector<std::vector<double>> _saved_points;



 public: // Méthodes et opérateurs de la classe
  Laplacian2D();
  // Constructeur : Initialiser _x_min, _x_max, _y_min; _y_max; _N; _h; _LapMat; _x; _y et _sol.
  virtual ~Laplacian2D();

  virtual void Initialize(DataFile datafile);
  virtual void InitializeMatrix();
  void UpdateCL (int num_it);
  virtual void DirectSolver(int nb_iterations) = 0;   // Résout le système _LapMat * _sol = _f avec un solveur direct.
  virtual void IterativeSolver(int nb_iterations) = 0;   // Résout le système _LapMat * _sol = _f avec un solveur itératif.
  virtual void SaveSol(std::string name_file); // Écrit les solutions dans le fichier "name_file".
  void ConditionsLimites(int num_it);
  virtual void Advance(int nb_iterations)=0;
};

class EC_ClassiqueM : public Laplacian2D //Première version avec un matrice identique quelles que soient les conditions aux bords
{
 public:
  void DirectSolver(int nb_iterations);
  void IterativeSolver(int nb_iterations);
  void ConditionsLimites(int num_it);
  inline void Advance(int nb_iteration) {};
};

class EC_ClassiqueP : public Laplacian2D //Seconde version avec une matrice qui dépend des conditions aux bords
{
 public:
  void InitializeMatrix();
  void DirectSolver(int nb_iterations);
  void IterativeSolver(int nb_iterations);
  void ConditionsLimites(int num_it);
  inline void Advance(int nb_iteration) {};
};

class EC_PyrolyseMC : public Laplacian2D //Schéma équation correction à matériau constant
{
 protected:
  double _FS, _FN, _FE, _FO;
  double _Cp, _Lambda, _A, _Ta;
  double _rho_v, _rho_p;
  Eigen::VectorXd _RhoTilde; //valeur provisoire de rho
  Eigen::VectorXd _sol_T; //solution en température
  Eigen::VectorXd _sol_R; //solution en masse volumique

 public:
  void Initialize(DataFile datafile);
  void InitializeMatrix();
  void SaveSol(int iteration);
  void Flux_Cal(int i, int j); //Calcul des flux Nord/Sud/Est/Ouest à lambda constant
  void Rho_Cal_P(); //Calcul de _RhoTilde (prédiction)
  void Rho_Cal_C(); //Calcul de _sol_R (correction)
  void T_Cal(); //Calcul de T si Cp est CONSTANT
  void Advance(int nb_iterations);
  void IterativeSolver(int nb_iterations);
  void ConditionsLimites(int num_it);

  inline void DirectSolver(int nb_iterations){};

};

class EC_PyrolyseMV : public EC_PyrolyseMC
{
 private:
  Eigen::VectorXd _lambdaMV,_CpMV;
  double _Cpp,_Cpv;

 public:
  void Initialize (DataFile datafile);
  void lambda_Cal();
  void Cp_Cal();
  void Advance(int nb_iterations); //version explicite
  void T_Cal(); //Calcul de T si Cp est variable
  void Flux_Cal(int i, int j); // Lambda variable
  void InitializeMatrix();
  void IterativeSolver(int nb_iterations);
  void ConditionsLimites(int num_it);
  void Newton(double epsilon);
};
