#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "Laplacian2D.h"

using namespace Eigen;
using namespace std ;


//Constructeur :
Laplacian2D::Laplacian2D()
{}
//Destructeur :
Laplacian2D::~Laplacian2D()
{}

void Laplacian2D::Initialize(double x_min, double x_max, double y_min, double y_max, int Nx, int Ny, double a, double deltaT)
{
  _x_min = x_min;
  _y_min = y_min;
  _x_max = x_max;
  _y_max = y_max;
  _Nx = Nx;
  _Ny = Ny;
  _a = a;
  _deltaT = deltaT;
  _h_y = (y_min-y_max)/(Ny+1.);
  _h_x = (x_max-x_min)/(Nx+1.);
}

void Laplacian2D::InitializeCI(Eigen::VectorXd CI)
{
  _sol = CI;
}

void Laplacian2D::InitializeCL(std::string CL_bas, std::string CL_haut, std::string CL_gauche, std::string CL_droite, double Val_CL_bas, double Val_CL_haut, double Val_CL_gauche, double Val_CL_droite)
{
  _CL_bas = CL_bas;
  _CL_haut = CL_haut;
  _CL_gauche = CL_gauche;
  _CL_droite = CL_droite;
  _Val_CL_bas = Val_CL_bas;
  _Val_CL_haut = Val_CL_haut;
  _Val_CL_gauche = Val_CL_gauche;
  _Val_CL_droite = Val_CL_droite;
}

void Laplacian2D::UpdateCL(int num_it)
{
  double t = num_it*_deltaT;
  if (t <= 50.)
  {
    _Val_CL_gauche = 10000.*t;
  }
  else
  {
    _Val_CL_gauche = -9000.*(t-50.) + 500000.;
  }
}

void EC_ClassiqueM::InitializeMatrix()
{
  system("rm -Rf EC_ClassiqueM");
  system("mkdir -p ./EC_ClassiqueM");

  _x.resize(_Nx);
  _y.resize(_Ny);
  for (int j =0; j < _Nx ; j++)
  {
    _x(j) = _x_min + j*_h_x;
  }

  for (int i = 0; i < _Ny ; i++)
  {
    _y(i) = _y_min + i*_h_y;
  }

  _LapMat.resize(_Nx*_Ny,_Nx*_Ny);

  double alpha = 1 + 2*_a*_deltaT/(_h_x*_h_x) + 2*_a*_deltaT/(_h_y*_h_y);
  double beta = -_a*_deltaT/(_h_x*_h_x);
  double gamma = -_a*_deltaT/(_h_y*_h_y);

  vector<Triplet<double>> liste_elem;

  for (int i = 0 ; i<_Nx*_Ny ; i++)
  {
    liste_elem.push_back({i,i,alpha});
  }

  for (int i = 0 ; i<_Ny*_Nx-1; i++)
  {
    if ((i+1)%_Nx!=0)
    {
      liste_elem.push_back({i,i+1,beta});
      liste_elem.push_back({i+1,i,beta});
    }
  }
  for (int i = 0 ; i<_Nx*(_Ny-1) ; i++)
  {
    liste_elem.push_back({i,_Nx+i,gamma});
    liste_elem.push_back({_Nx+i,i,gamma});
  }

  _LapMat.setFromTriplets(liste_elem.begin(), liste_elem.end());
}

void EC_ClassiqueP::InitializeMatrix()
{
  system("rm -Rf EC_ClassiqueP");
  system("mkdir -p ./EC_ClassiqueP");

  _x.resize(_Nx);
  _y.resize(_Ny);
  for (int j =0; j < _Nx ; j++)
  {
    _x(j) = _x_min + j*_h_x;
  }

  for (int i = 0; i < _Ny ; i++)
  {
    _y(i) = _y_min + i*_h_y;
  }

  _LapMat.resize(_Nx*_Ny,_Nx*_Ny);

  double alpha = 1 + 2*_a*_deltaT/(_h_x*_h_x) + 2*_a*_deltaT/(_h_y*_h_y);
  double beta = -_a*_deltaT/(_h_x*_h_x);
  double gamma = -_a*_deltaT/(_h_y*_h_y);

  vector<Triplet<double>> liste_elem;

  for (int i = 0 ; i<_Nx*_Ny ; i++)
  {
    liste_elem.push_back({i,i,alpha});
  }

  for (int i = 0 ; i<_Ny*_Nx-1; i++)
  {
    if ((i+1)%_Nx!=0)
    {
      liste_elem.push_back({i,i+1,beta});
      liste_elem.push_back({i+1,i,beta});
    }
  }
  for (int i = 0 ; i<_Nx*(_Ny-1) ; i++)
  {
    liste_elem.push_back({i,_Nx+i,gamma});
    liste_elem.push_back({_Nx+i,i,gamma});
  }

  _LapMat.setFromTriplets(liste_elem.begin(), liste_elem.end());


  if (_CL_gauche == "Neumann" or _CL_gauche == "Neumann_non_constant")
  {
    for (int i = 0 ; i < _Ny; i++)
    {
      _LapMat.coeffRef(_Nx*i,_Nx*i) += beta;  //Bord gauche
    }
  }

  if (_CL_droite == "Neumann")
  {
    for (int i = 0 ; i < _Ny; i++)
    {
      _LapMat.coeffRef(_Nx*(i+1) - 1, _Nx*(i+1) - 1) += beta; //Bord droite
    }
  }

  if (_CL_haut == "Neumann")
  {
    for (int i = 0; i < _Nx ; i++)
    {
      _LapMat.coeffRef(i,i) += gamma; //Bord haut
    }
  }

  if (_CL_bas == "Neumann")
  {
    for (int i = 0; i < _Nx ; i++)
    {
      _LapMat.coeffRef((_Ny - 1)* _Nx + i , (_Ny - 1)* _Nx + i) += gamma; //Bord bas
    }
  }
}

void EC_ClassiqueM::DirectSolver (int nb_iterations)
{
  SimplicialLLT <SparseMatrix<double> > solver;
  solver.compute(_LapMat);

  for( int i=0 ; i<=nb_iterations ; i++)
    {
      EC_ClassiqueM::SaveSol("EC_ClassiqueM/sol_it_"+to_string(i)+".vtk");
      EC_ClassiqueM::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
        {
          _f(j) = _sol(j);
        }
      _sol = solver.solve(_f);
    }
}

void EC_ClassiqueP::DirectSolver (int nb_iterations)
{
  SimplicialLLT <SparseMatrix<double> > solver;
  solver.compute(_LapMat);

  for( int i=0 ; i<=nb_iterations ; i++)
    {
      EC_ClassiqueP::SaveSol("EC_ClassiqueP/sol_it_"+to_string(i)+".vtk");
      EC_ClassiqueP::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
        {
          _f(j) = _sol(j);
        }
      _sol = solver.solve(_f);
    }
}

void EC_ClassiqueM::IterativeSolver (int nb_iterations)
{
  ConjugateGradient <SparseMatrix<double> > solver;
  solver.compute(_LapMat);

  for( int i=0 ; i<=nb_iterations ; i++)
    {
      EC_ClassiqueM::SaveSol("EC_ClassiqueM/sol_it_"+to_string(i)+".vtk");
      EC_ClassiqueM::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
        {
          _f(j) = _sol(j);
        }
      _sol = solver.solve(_f);
    }
}

void EC_ClassiqueP::IterativeSolver (int nb_iterations)
{
  ConjugateGradient <SparseMatrix<double> > solver;
  solver.compute(_LapMat);

  for( int i=0 ; i<=nb_iterations ; i++)
    {
      EC_ClassiqueP::SaveSol("EC_ClassiqueP/sol_it_"+to_string(i)+".vtk");
      EC_ClassiqueP::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
        {
          _f(j) = _sol(j);
        }
      _sol = solver.solve(_f);
    }
}

void Laplacian2D::SaveSol(string name_file)
{
  using namespace std;
  ofstream mon_flux;
  mon_flux.open(name_file, ios::out);
  mon_flux << "# vtk DataFile Version 3.0" << endl
	   << "cell" << endl
	   << "ASCII" << endl
	   << "DATASET STRUCTURED_POINTS" << endl
	   << "DIMENSIONS " << _Nx << " " << _Ny << " 1" << endl
	   << "ORIGIN 0 0 0" << endl
	   << "SPACING " + to_string((_x_max-_x_min)/_Nx)+ " " + to_string((_y_max-_y_min)/_Ny) +" 1" << endl
	   << "POINT_DATA " << _Nx*_Ny << endl
	   << "SCALARS sample_scalars double" << endl
	   << "LOOKUP_TABLE default" << endl;

  for(int i=_Ny-1; i>=0; i--)
    {
      for(int j=0; j<_Nx; j++)
	{
	  mon_flux << _sol(j + i*_Nx) << " ";
	}
      mon_flux << endl;
    }

  mon_flux.close();
}

void EC_ClassiqueM::ConditionsLimites(int num_it)
{
  Eigen :: VectorXd temp(_Nx*_Ny);
  double gamma = -_a*_deltaT/(_h_y*_h_y);
  double beta = -_a*_deltaT/(_h_x*_h_x);

  for (int i =0; i<_Nx*_Ny ; i++)
  {
    temp(i) = _sol(i);
  }

  if (_CL_haut == "Neumann")
  {
    for (int j = 0; j < _Nx ; j++) //Condition de flux en haut
    {
      _sol(j) = _sol(j)-gamma*temp(j) - gamma*_Val_CL_haut*_h_y;
    }
  }

  if (_CL_bas == "Neumann") //Condition de flux en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(_Nx*(_Ny -1)+ j) = _sol(_Nx*(_Ny -1)+ j)-gamma*temp(_Nx*(_Ny -1)+ j) - gamma*_Val_CL_bas*_h_y;
    }
  }

  if (_CL_gauche == "Neumann") //Condition de flux à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)-beta*temp(i*_Nx) - beta*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_gauche == "Neumann_non_constant") //Condition de flux à gauche
  {
    Laplacian2D::UpdateCL(num_it);
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)-beta*temp(i*_Nx) - beta*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_droite == "Neumann") //Condition de flux à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol((i+1)*_Nx - 1) = _sol((i+1)*_Nx - 1)-beta*temp((i+1)*_Nx - 1) - beta*_Val_CL_droite*_h_x;
    }
  }

  if (_CL_haut == "Dirichlet") //Condition de température en haut
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(j) = _sol(j)-gamma*_Val_CL_haut;
    }
  }
  if (_CL_bas == "Dirichlet") //Condition de température en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(_Nx*(_Ny -1)+ j) = _sol(_Nx*(_Ny -1)+ j)-gamma*_Val_CL_bas;
    }
  }

  if (_CL_gauche == "Dirichlet")  //Condition de température à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)-beta*_Val_CL_gauche;
    }
  }
  if (_CL_droite == "Dirichlet") //Condition de température à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol((i+1)*_Nx - 1) = _sol((i+1)*_Nx - 1)-beta*_Val_CL_droite;
    }
  }
}

void EC_ClassiqueP::ConditionsLimites(int num_it)
{
  double gamma = -_a*_deltaT/(_h_y*_h_y);
  double beta = -_a*_deltaT/(_h_x*_h_x);

  if (_CL_haut == "Dirichlet") //Condition de température en haut
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(j) = _sol(j)-gamma*_Val_CL_haut;
    }
  }
  if (_CL_bas == "Dirichlet") //Condition de température en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(_Nx*(_Ny -1)+ j) = _sol(_Nx*(_Ny -1)+ j)-gamma*_Val_CL_bas;
    }
  }

  if (_CL_gauche == "Dirichlet")  //Condition de température à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)-beta*_Val_CL_gauche;
    }
  }

  if (_CL_droite == "Dirichlet") //Condition de température à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol((i+1)*_Nx - 1) = _sol((i+1)*_Nx - 1)-beta*_Val_CL_droite;
    }
  }

  if (_CL_haut == "Neumann") //Condition de flux en haut
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(j) = _sol(j)-gamma*_Val_CL_haut*_h_y;
    }
  }
  if (_CL_bas == "Neumann") //Condition de flux en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(_Nx*(_Ny -1)+ j) = _sol(_Nx*(_Ny -1)+ j)-gamma*_Val_CL_bas*_h_y;
    }
  }

  if (_CL_gauche == "Neumann")  //Condition de flux à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)-beta*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_gauche == "Neumann_non_constant")  //Condition de flux à gauche
  {
    Laplacian2D::UpdateCL(num_it);
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)-beta*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_droite == "Neumann") //Condition de flux à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol((i+1)*_Nx - 1) = _sol((i+1)*_Nx - 1)-beta*_Val_CL_droite*_h_x;
    }
  }
}
