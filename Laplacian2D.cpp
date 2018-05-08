#include "Laplacian2D.h"

using namespace Eigen;
using namespace std ;


//Constructeur :
Laplacian2D::Laplacian2D()
{}
//Destructeur :
Laplacian2D::~Laplacian2D()
{}

void Laplacian2D::Initialize(DataFile data_file)
{
  _CL_droite = data_file.Get_CL_droite();
  _CL_gauche = data_file.Get_CL_gauche();
  _CL_haut = data_file.Get_CL_haut();
  _CL_bas = data_file.Get_CL_bas();
  _Val_CL_droite = data_file.Get_Val_CL_droite();
  _Val_CL_gauche = data_file.Get_Val_CL_gauche();
  _Val_CL_haut = data_file.Get_Val_CL_haut();
  _Val_CL_bas = data_file.Get_Val_CL_bas();

  _Nx = data_file.Get_N_x();
  _Ny = data_file.Get_N_y();
  _x_min = data_file.Get_x_min();
  _x_max = data_file.Get_x_max();
  _y_min = data_file.Get_y_min();
  _y_max = data_file.Get_y_max();
  _deltaT = data_file.Get_deltaT();
  //_T_final = data_file.Get_T_final(); -> inutile ici
  _Solveur = data_file.Get_Solveur();
  //Schema = data_file.Get_Schema(); -> pour l'instant ne sert à rien


  _save_all_file = data_file.Get_save_all_file();

  _save_points_file = data_file.Get_save_points_file();
  _number_saved_points = data_file.Get_number_saved_points();
  _saved_points = data_file.Get_saved_points();

  _restart_file = data_file.Get_restart_file();


  _h_y = (_y_max-_y_min)/(_Ny+1.);
  _h_x = (_x_max-_x_min)/(_Nx+1.);


  //-----------Def de la CI--------------
  if (_restart_file == "non")
    {
      _sol.resize(_Nx*_Ny);
      for(int j=0; j < _Ny; j++)
	{
	  for(int i=0; i < _Nx; i++)
  	    {
  	      _sol(i + j*_Nx) = data_file.Get_CI();
  	    }
	}

      if ((data_file.Get_eq() == "EC_ClassiqueM") or (data_file.Get_eq() == "EC_ClassiqueP")) //Virer ça et mettre _a et ce truc dans un initialize pour ces 2 classes filles.
	{
	  _a = data_file.Get_lambda()/(data_file.Get_rho()*data_file.Get_Cp());
	}


    }
  else
    {//Mettre ici le truc pour le restart file
    }
  //---------------------------------------

  //-----Creation des fichiers de sauvegarde------
  //Je pense que ça pourrait être pas mal de faire un truc pour vérifier si on a déja enregistrer quelque chose sous le même nom et demander à l'utilisateur si il est sur de vouloir suprimer ce truc. (à voir)


  if (_save_all_file !="non")
    {
      system(("rm -Rf "+_save_all_file).c_str());
      system(("mkdir -p ./"+_save_all_file).c_str());
    }


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

void Laplacian2D::InitializeMatrix()
{
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

  Laplacian2D::InitializeMatrix();

  double beta = -_a*_deltaT/(_h_x*_h_x);
  double gamma = -_a*_deltaT/(_h_y*_h_y);

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

  ofstream* flux_pts(new ofstream);

  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      flux_pts->open(_save_points_file+".txt", ios::out);
    }

  for( int i=0 ; i<=nb_iterations ; i++)
    {

      // Systeme de sauvegarde de points :---------------------------------
      if (_save_all_file != "non")
	{
	  EC_ClassiqueM::SaveSol(_save_all_file+"/sol_it_"+to_string(i)+".vtk");
	}


      if (_save_points_file != "non")
	{
	  *flux_pts<<i*_deltaT<<" ";
	  //char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      int pos = floor(_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y);
	      *flux_pts<<_sol(pos)<<" ";
	    }
	  *flux_pts<<endl;
	}
      //-------------------------------------------------------------------

      EC_ClassiqueM::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
        {
          _f(j) = _sol(j);
        }
      _sol = solver.solve(_f);
    }


  if(_save_points_file != "non")
    {
      //On referme les flux qu'on a ouvert
      flux_pts->close();
    }

  delete flux_pts;
}

void EC_ClassiqueP::DirectSolver (int nb_iterations)
{
  SimplicialLLT <SparseMatrix<double> > solver;
  solver.compute(_LapMat);

  ofstream* flux_pts(new ofstream);

  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      flux_pts->open(_save_points_file+".txt", ios::out);
    }

  for( int i=0 ; i<=nb_iterations ; i++)
    {

      // Systeme de sauvegarde de points :---------------------------------
      if (_save_all_file != "non")
	{
	  EC_ClassiqueP::SaveSol(_save_all_file+"/sol_it_"+to_string(i)+".vtk");
	}


       if (_save_points_file != "non")
	 {
	   *flux_pts<<i*_deltaT<<" ";
	  //char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      int pos = floor(_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y);
	      *flux_pts<<_sol(pos)<<" ";
	    }
	  *flux_pts<<endl;
	 }
      //-------------------------------------------------------------------

      EC_ClassiqueP::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
        {
          _f(j) = _sol(j);
        }
      _sol = solver.solve(_f);
    }


  if(_save_points_file != "non")
    {
      //On referme les flux qu'on a ouvert
      flux_pts->close();
    }

  delete flux_pts;

}

void EC_ClassiqueM::IterativeSolver (int nb_iterations)
{
  ConjugateGradient <SparseMatrix<double> > solver;
  solver.compute(_LapMat);

  ofstream* flux_pts(new ofstream);

  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      flux_pts->open(_save_points_file+".txt", ios::out);
    }

  for( int i=0 ; i<=nb_iterations ; i++)
    {

      // Systeme de sauvegarde de points :---------------------------------
      if (_save_all_file != "non")
	{
	  EC_ClassiqueM::SaveSol(_save_all_file+"/sol_it_"+to_string(i)+".vtk");
	}


      if (_save_points_file != "non")
	{
	  *flux_pts<<i*_deltaT<<" ";
	  //char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      int pos = floor(_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y);
	      *flux_pts<<_sol(pos)<<" ";
	    }
	  *flux_pts<<endl;
	}
      //-------------------------------------------------------------------

      EC_ClassiqueM::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
        {
          _f(j) = _sol(j);
        }
      _sol = solver.solve(_f);
    }


  if(_save_points_file != "non")
    {
      flux_pts->close();
    }

  delete flux_pts;

}

void EC_ClassiqueP::IterativeSolver (int nb_iterations)
{
  ConjugateGradient <SparseMatrix<double> > solver;
  solver.compute(_LapMat);

  ofstream* flux_pts(new ofstream);

  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      flux_pts->open(_save_points_file+".txt", ios::out);
    }

  for( int i=0 ; i<=nb_iterations ; i++)
    {

      // Systeme de sauvegarde de points :---------------------------------
      if (_save_all_file != "non")
	{
	  EC_ClassiqueP::SaveSol(_save_all_file+"/sol_it_"+to_string(i)+".vtk");
	}


      if (_save_points_file != "non")
	{
	  *flux_pts<<i*_deltaT<<" ";
	  //char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      int pos = floor(_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y);
	      *flux_pts<<_sol(pos)<<" ";
	    }
	  *flux_pts<<endl;
	}

      //-------------------------------------------------------------------

      EC_ClassiqueP::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
	{
	  _f(j) = _sol(j);
	}
      _sol = solver.solve(_f);
    }


  if(_save_points_file != "non")
    {
      flux_pts->close();
    }
  delete flux_pts;

}

void Laplacian2D::SaveSol(string name_file)
{


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
	  //mon_flux.write("bonjour",7);

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

//-------------------------------------------------------------------------------------------------------------------------
//Début de EC_Pyrolyse:
//-------------------------------------------------------------------------------------------------------------------------

void EC_PyrolyseMC::Initialize(DataFile data_file)
{
  Laplacian2D::Initialize(data_file);

  _rho_v = data_file.Get_rhov(); //1500.
  _rho_p = data_file.Get_rhop();//1000.
  _Cp = data_file.Get_Cp();//1500
  _Lambda = data_file.Get_lambda();
  _A = data_file.Get_Aexp();//1000.
  _Ta = data_file.Get_Ta();//6000.



  _sol_T.resize(_Nx*_Ny);
  _sol_R.resize(_Nx*_Ny);
  _RhoTilde.resize(_Nx*_Ny);
  for(int j=0; j < _Ny; j++)
    {
      for(int i=0; i < _Nx; i++)
	{
	  _sol_T(i + j*_Nx) = data_file.Get_CI();
	  _sol_R(i+j*_Nx) = _rho_v;
	}
    }

  cout << "valeur de _h_y " << _h_y << endl;
  cout << "valeur de _h_x " << _h_x << endl;
  cout << "cfl :" << _Lambda*_deltaT/(_rho_v*_Cp)*(1/(_h_x*_h_x) + 1/(_h_y*_h_y)) << endl;


}

void EC_PyrolyseMC::InitializeMatrix()
{
  _LapMat.resize(_Nx*_Ny,_Nx*_Ny);

  Eigen::VectorXd alpha;
  Eigen::VectorXd beta;
  Eigen::VectorXd gamma;

  alpha.resize(_Nx*_Ny);
  beta.resize(_Nx*_Ny);
  gamma.resize(_Nx*_Ny);



  for (int i=0;i<_Nx*_Ny;i++)
    {
      alpha[i] = 1 + 2*(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_x*_h_x) + 2*(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_y*_h_y);
      beta[i] = -(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_x*_h_x);
      gamma[i] = -(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_y*_h_y);
    }


  vector<Triplet<double>> liste_elem;

  for (int i = 0 ; i<_Nx*_Ny ; i++)
  {
    liste_elem.push_back({i,i,alpha[i]});
  }

  for (int i = 0 ; i<_Ny*_Nx-1; i++)
  {
    if ((i+1)%_Nx!=0)
    {
      liste_elem.push_back({i,i+1,beta[i]});
      liste_elem.push_back({i+1,i,beta[i+1]});
    }
  }
  for (int i = 0 ; i<_Nx*(_Ny-1) ; i++)
  {
    liste_elem.push_back({i,_Nx+i,gamma[i]});
    liste_elem.push_back({_Nx+i,i,gamma[_Nx+i]});
  }

  _LapMat.setFromTriplets(liste_elem.begin(), liste_elem.end());


  if (_CL_gauche == "Neumann" or _CL_gauche == "Neumann_non_constant")
   {
     for (int i = 0 ; i < _Ny; i++)
     {
       _LapMat.coeffRef(_Nx*i,_Nx*i) += beta[_Nx*i];  //Bord gauche
     }
   }

   if (_CL_droite == "Neumann")
   {
     for (int i = 0 ; i < _Ny; i++)
     {
       _LapMat.coeffRef(_Nx*(i+1) - 1, _Nx*(i+1) - 1) += beta[_Nx*(i+1) - 1]; //Bord droite
     }
   }


   if (_CL_haut == "Neumann")
   {
     for (int i = 0; i < _Nx ; i++)
     {
       _LapMat.coeffRef(i,i) += gamma[i]; //Bord haut
     }
   }

   if (_CL_bas == "Neumann")
   {
     for (int i = 0; i < _Nx ; i++)
     {
       _LapMat.coeffRef((_Ny - 1)* _Nx + i , (_Ny - 1)* _Nx + i) += gamma[(_Ny - 1)* _Nx + i]; //Bord bas
     }
   }
}

void EC_PyrolyseMC::Flux_Cal(int i, int j)
{
  _FN = 0; _FS = 0; _FE = 0; _FO = 0;

  double Tg,Td,Tb,Th;

  //Modifier les valeurs limites en fonction des conditions
  if (_CL_gauche == "Neumann" or _CL_gauche == "Neumann_non_constant")
    Tg = _sol_T[_Nx*i+0] + _h_x*_Val_CL_gauche;
  else
    Tg = _Val_CL_gauche;

  if (_CL_droite == "Neumann")
    Td = _sol_T[_Nx*i+(_Nx-1)] + _h_x*_Val_CL_droite;
  else
    Td = _Val_CL_droite;

  if (_CL_haut == "Neumann")
    Th = _sol_T[_Nx*(_Ny-1)+j] + _h_y*_Val_CL_haut;
  else
    Th = _Val_CL_haut;;

  if (_CL_bas == "Neumann")
    Tb = _sol_T[_Nx*0+j] + _h_y*_Val_CL_bas;
  else
    Tb = _Val_CL_bas;


  if((i>0) and (i<_Ny-1))
    {
      if((j>0) and (j<_Nx-1))
	{
	  _FN = _Lambda*(_sol_T[_Nx*(i+1)+j] - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*(i-1)+j])/_h_y;
	  _FE = _Lambda*(_sol_T[_Nx*i+(j+1)] - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*i+(j-1)])/_h_x;
	}
      else if (j==0)
	{
	  _FN = _Lambda*(_sol_T[_Nx*(i+1)+j] - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*(i-1)+j])/_h_y;
	  _FE = _Lambda*(_sol_T[_Nx*i+(j+1)] - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - Tg)/_h_x;
	}
      else if (j==_Nx-1)
	{
	  _FN = _Lambda*(_sol_T[_Nx*(i+1)+j] - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*(i-1)+j])/_h_y;
	  _FE = _Lambda*(Td - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*i+(j-1)])/_h_x;
	}
    }
  else if (i==0)
    {
      if((j>0) and (j<_Nx-1))
	{
	  _FN = _Lambda*(_sol_T[_Nx*(i+1)+j] - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - Tb)/_h_y;
	  _FE = _Lambda*(_sol_T[_Nx*i+(j+1)] - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*i+(j-1)])/_h_x;
	}
      else if (j==0)
	{
	  _FN = _Lambda*(_sol_T[_Nx*(i+1)+j] - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - Tb)/_h_y;
	  _FE = _Lambda*(_sol_T[_Nx*i+(j+1)] - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - Tg)/_h_x;
	}
      else if (j==_Nx-1)
	{
	  _FN = _Lambda*(_sol_T[_Nx*(i+1)+j] - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - Tb)/_h_y;
	  _FE = _Lambda*(Td - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*i+(j-1)])/_h_x;
	}
    }
  else if (i==_Ny-1)
    {
      if((j>0) and (j<_Nx-1))
	{
	  _FN = _Lambda*(Th - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*(i-1)+j])/_h_y;
	  _FE = _Lambda*(_sol_T[_Nx*i+(j+1)] - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*i+(j-1)])/_h_x;
	}
      else if (j==0)
	{
	  _FN = _Lambda*(Th - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*(i-1)+j])/_h_y;
	  _FE = _Lambda*(_sol_T[_Nx*i+(j+1)] - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] -  Tg)/_h_x;
	}
      else if (j==_Nx-1)
	{
	  _FN = _Lambda*(Th - _sol_T[_Nx*i+j])/_h_y;
	  _FS = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*(i-1)+j])/_h_y;
	  _FE = _Lambda*(Td - _sol_T[_Nx*i+j])/_h_x;
	  _FO = _Lambda*(_sol_T[_Nx*i+j] - _sol_T[_Nx*i+(j-1)])/_h_x;
	}
    }
}

void EC_PyrolyseMC::Rho_Cal_P()
{
  double k;
  // k = _A*_rho_p/(_rho_v - _rho_p);
  k = -_A*_rho_v/(_rho_v - _rho_p);

  Eigen::VectorXd C;
  C.resize(_Nx*_Ny);
  for (int i=0; i<_Nx*_Ny; i++)
    {
      //Euler Exp
      // C[i] = k*exp(-_Ta/_sol_T[i]);
      // _RhoTilde[i] = (_deltaT*C[i]+1)*_sol_R[i] - _deltaT*C[i]*_rho_v;

      //sol exacte
      C[i] = exp(-_Ta/_sol_T[i]);
      // _RhoTilde[i] = _sol_R[i]*((exp(k*_deltaT)-1)*C[i]+1) - k*_rho_v*C[i]*_deltaT;
      _RhoTilde[i] = (_sol_R[i]-_rho_p)*exp(k*C[i]*_deltaT)  + _rho_p;
    }
}

void EC_PyrolyseMC::Rho_Cal_C()
{
  double k;
  // k = _A*_rho_p/(_rho_v - _rho_p);
  k = -_A*_rho_v/(_rho_v - _rho_p);

  Eigen::VectorXd C;
  C.resize(_Nx*_Ny);
  for (int i=0; i<_Nx*_Ny; i++)
    {
      //Euler Exp
      // C[i] = k*exp(-_Ta/_sol_T[i]);
      // _sol_R[i] = (_deltaT*C[i]+1)*_sol_R[i] - _deltaT*C[i]*_rho_v;

      //sol exacte
      C[i] = exp(-_Ta/_sol_T[i]);
      // _sol_R[i] = _sol_R[i]*((exp(k*_deltaT)-1)*C[i]+1) - k*_rho_v*C[i]*_deltaT;
      _sol_R[i] = (_sol_R[i]-_rho_p)*exp(k*C[i]*_deltaT) + _rho_p;
    }
}

void EC_PyrolyseMC::T_Cal()
{
  Eigen::VectorXd Tbis;
  Tbis.resize(_Nx*_Ny);
  for (int i =0; i<_Nx*_Ny ; i++)
  {
    Tbis(i) = _sol_T(i);
  }

  for (int i=0; i<_Ny; i++)
    {
      for (int j=0; j<_Nx; j++)
	{
	  EC_PyrolyseMC::Flux_Cal(i,j);
	  Tbis[i*_Nx+j] = (_deltaT*((_FE-_FO)/_h_x + (_FN-_FS)/_h_y))/(_RhoTilde[i*_Nx+j]*_Cp) + Tbis[i*_Nx+j];
    //Tbis[i*_Nx+j] = (_deltaT*((_FE-_FO)/_h_x + (_FN-_FS)/_h_y))/(1000*_Cp) + Tbis[i*_Nx+j];
	}
    }

  for (int i =0; i<_Nx*_Ny ; i++)
  {
    _sol_T(i) = Tbis(i);
  }



}

void EC_PyrolyseMC::SaveSol(int iteration)
{
  ofstream mon_flux;
  string name_file1 = _save_all_file+"/solT_it_"+to_string(iteration)+".vtk";
  mon_flux.open(name_file1, ios::out);
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
	  mon_flux << _sol_T[j + i*_Nx] << " ";
	}
      mon_flux << endl;
    }

  mon_flux.close();

  string name_file2 = _save_all_file+"/solR_it_"+to_string(iteration)+".vtk";
  mon_flux.open(name_file2, ios::out);
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
  	  mon_flux << _sol_R(j + i*_Nx) << " ";
  	}
      mon_flux << endl;
    }

  mon_flux.close();
}

void EC_PyrolyseMC::Advance(int nb_iterations)
{

  ofstream* flux_pts(new ofstream);


  if(_save_points_file != "non")
  {
    //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
    flux_pts->open(_save_points_file+".txt", ios::out);

  }

  for( int i=0 ; i<=nb_iterations ; i++)
  {
    // Systeme de sauvegarde de points :---------------------------------
    if (_save_all_file != "non")
    {
      EC_PyrolyseMC::SaveSol(i);
    }

    if (_save_points_file != "non")
    {
      *flux_pts<<i*_deltaT<<" ";
      //char* truc = new char;
      for (int j=0; j<_number_saved_points; j++)
      {
        int pos = floor(_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y);
        *flux_pts<<_sol_T(pos)<<" "<<_sol_R(pos)<<" ";
      }
      *flux_pts<<endl;
    }

    Rho_Cal_P();
    if (_CL_gauche == "Neumann_non_constant")
    {
      Laplacian2D::UpdateCL(i+1);
    }
    T_Cal();
    Rho_Cal_C();


    //Barre de chargement
    int i_barre;
    int p = floor((((double)i)/((double)nb_iterations))*100);
    printf( "[" );
    for(i_barre=0;i_barre<=p;i_barre+=2) printf( "*" );
    for (;i_barre<100; i_barre+=2 ) printf( "-" );
    printf( "] %3d %%", p );

    for(i_barre=0;i_barre<59;++i_barre) printf( "%c", 8 );

    fflush(stdout );

  }

  if(_save_points_file != "non")
  {
    flux_pts->close();
  }

  printf( "\n" );
  delete flux_pts;

}

void EC_PyrolyseMC::IterativeSolver (int nb_iterations)
{
  BiCGSTAB <SparseMatrix<double> > solver;

  ofstream* flux_pts(new ofstream);

  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      flux_pts->open(_save_points_file+".txt", ios::out);
    }

  for( int i=0 ; i<=nb_iterations ; i++)
    {

      // Systeme de sauvegarde de points :---------------------------------
      if (_save_all_file != "non")
	{
	  EC_PyrolyseMC::SaveSol(i);
	}

      if (_save_points_file != "non")
	{
	  *flux_pts<<i*_deltaT<<" ";
	  //char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      int pos = floor(_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y);
	      *flux_pts<<_sol_T(pos)<<" "<<_sol_R(pos)<<" ";
	    }
	  *flux_pts<<endl;
	}
      //-------------------------------------------------------------------

      Rho_Cal_P();
      EC_PyrolyseMC::InitializeMatrix();
      solver.compute(_LapMat);
      EC_PyrolyseMC::ConditionsLimites(i);
      _f.resize(_Nx*_Ny);
      for (int j =0; j<_Nx*_Ny ; j++)
	{
	  _f(j) = _sol_T(j);
	}
      _sol_T = solver.solve(_f);
      VectorXd r = _LapMat*_sol_T - _f;
      Rho_Cal_C();
      // cout << "Rhotilde[14] = " << _RhoTilde[14] << " , T[14] = " << _sol_T[14] << " , Rho[14] = " << _sol_R[14] <<endl;

      //barre_de_chargement

      int i_barre;
      int p = floor((((double)i)/((double)nb_iterations))*100);
      printf( "[" );
      for(i_barre=0;i_barre<=p;i_barre+=2) printf( "*" );
      for (;i_barre<100; i_barre+=2 ) printf( "-" );
      printf( "] %3d %%", p );

      for(i_barre=0;i_barre<59;++i_barre) printf( "%c", 8 );

      fflush(stdout );
    }

  if(_save_points_file != "non")
    {
      //On referme les flux qu'on a ouvert
      flux_pts->close();
    }


  delete flux_pts;
  printf("\n");

}

void EC_PyrolyseMC::ConditionsLimites(int num_it)
{
  Eigen::VectorXd alpha;
  Eigen::VectorXd beta;
  Eigen::VectorXd gamma;

  alpha.resize(_Nx*_Ny);
  beta.resize(_Nx*_Ny);
  gamma.resize(_Nx*_Ny);

  for (int i=0;i<_Nx*_Ny;i++)
  {
    alpha[i] = 1 + 2*(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_x*_h_x) + 2*(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_y*_h_y);
    beta[i] = -(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_x*_h_x);
    gamma[i] = -(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_y*_h_y);
  }

  if (_CL_haut == "Dirichlet") //Condition de température en haut
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(j) = _sol_T(j)-gamma(j)*_Val_CL_haut;
    }
  }
  if (_CL_bas == "Dirichlet") //Condition de température en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(_Nx*(_Ny -1)+ j) = _sol_T(_Nx*(_Ny -1)+ j)-gamma(_Nx*(_Ny -1)+ j)*_Val_CL_bas;
    }
  }

  if (_CL_gauche == "Dirichlet")  //Condition de température à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T(i*_Nx) = _sol_T(i*_Nx)-beta(i*_Nx)*_Val_CL_gauche;
    }
  }
  if (_CL_droite == "Dirichlet") //Condition de température à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T((i+1)*_Nx - 1) = _sol_T((i+1)*_Nx - 1)-beta((i+1)*_Nx - 1)*_Val_CL_droite;
    }
  }

  if (_CL_haut == "Neumann") //Condition de flux en haut
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(j) = _sol_T(j)-gamma(j)*_Val_CL_haut*_h_y;
    }
  }
  if (_CL_bas == "Neumann") //Condition de flux en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(_Nx*(_Ny -1)+ j) = _sol_T(_Nx*(_Ny -1)+ j)-gamma(_Nx*(_Ny -1)+ j)*_Val_CL_bas*_h_y;
    }
  }

  if (_CL_gauche == "Neumann")  //Condition de flux à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T(i*_Nx) = _sol_T(i*_Nx)-beta(i*_Nx)*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_gauche == "Neumann_non_constant")  //Condition de flux à gauche
  {
    Laplacian2D::UpdateCL(num_it);
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T(i*_Nx) = _sol_T(i*_Nx)-beta(i*_Nx)*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_droite == "Neumann") //Condition de flux à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T((i+1)*_Nx - 1) = _sol_T((i+1)*_Nx - 1)-beta((i+1)*_Nx - 1)*_Val_CL_droite*_h_x;
    }
  }
}

//-----------------------------------------------------------------------------------
//Pyrolyse Materiau variable
//-----------------------------------------------------------------------------------

void EC_PyrolyseMV::Initialize(DataFile data_file)
{
  EC_PyrolyseMC::Initialize(data_file);

  _lambdaMV.resize(_Nx*_Ny);
  _Cpp = data_file.Get_Cpp();
  _Cpv = data_file.Get_Cpv();

  // Si on fait un truc en explicite mettre la cfl ici.
}

void EC_PyrolyseMV::lambda_Cal() // On calcule Lambda à partir de _RhoTilde et de _sol_T
{
  double lv; //Lambda vierge
  double lp; //Lambda pyrolysé
  double ksi;
  for (int i=0; i<_Nx*_Ny; i++)
  {
    if (_sol_T(i) < 1000)
    {
      lv = -0.5/728.*_sol_T(i) + 1;
      lp = 1./728.*_sol_T(i) + 1;
    }
    else
    {
      lv = 0.5;
      lp = 2.;
    }
    ksi = (_rho_v - _RhoTilde(i))/(_rho_v - _rho_p);
    _lambdaMV(i) = (1-ksi)*lv + ksi*lp;
  }
}

void EC_PyrolyseMV::Cp_Cal() // On calcule Lambda à partir de _RhoTilde et de _sol_T
{
  for (int i = 0; i < _Nx*_Ny; i++)
  {
    double ksi = (_rho_v - _RhoTilde(i))/(_rho_v - _rho_p);
    _CpMV(i) = (1-ksi)*_Cpv + ksi*_Cpp;
  }
}

void EC_PyrolyseMV::IterativeSolver(int nb_iterations)
{
  BiCGSTAB <SparseMatrix<double> > solver;

  ofstream* flux_pts(new ofstream);

  if(_save_points_file != "non")
  {
    //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
    flux_pts->open(_save_points_file+".txt", ios::out);
  }

  for( int i=0 ; i<=nb_iterations ; i++)
  {
    // Systeme de sauvegarde de points :---------------------------------
    if (_save_all_file != "non")
    {
      EC_PyrolyseMC::SaveSol(i);
    }

    if (_save_points_file != "non")
    {
      *flux_pts<<i*_deltaT<<" ";
      //char* truc = new char;
      for (int j=0; j<_number_saved_points; j++)
      {
        int pos = floor(_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y);
        *flux_pts<<_sol_T(pos)<<" "<<_sol_R(pos)<<" ";
      }
      *flux_pts<<endl;
    }
    //-------------------------------------------------------------------

    Rho_Cal_P();
    Cp_Cal();
    EC_PyrolyseMV::InitializeMatrix();
    solver.compute(_LapMat);
    EC_PyrolyseMV::ConditionsLimites(i);
    _f.resize(_Nx*_Ny);
    for (int j =0; j<_Nx*_Ny ; j++)
    {
      _f(j) = _sol_T(j);
    }


    //---------------------------Newton-----------------------------------



  }
}

void EC_PyrolyseMV::ConditionsLimites(int num_it)
{
  Eigen::VectorXd alpha;
  Eigen::VectorXd beta;
  Eigen::VectorXd gamma;

  alpha.resize(_Nx*_Ny);
  beta.resize(_Nx*_Ny);
  gamma.resize(_Nx*_Ny);

  for (int i=0;i<_Nx*_Ny;i++)
  {
    _Lambda = _lambdaMV(i);
    alpha[i] = 1 + 2*(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_x*_h_x) + 2*(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_y*_h_y);
    beta[i] = -(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_x*_h_x);
    gamma[i] = -(_Lambda/(_RhoTilde[i]*_Cp))*_deltaT/(_h_y*_h_y);
  }

  if (_CL_haut == "Dirichlet") //Condition de température en haut
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(j) = _sol_T(j)-gamma(j)*_Val_CL_haut;
    }
  }
  if (_CL_bas == "Dirichlet") //Condition de température en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(_Nx*(_Ny -1)+ j) = _sol_T(_Nx*(_Ny -1)+ j)-gamma(_Nx*(_Ny -1)+ j)*_Val_CL_bas;
    }
  }

  if (_CL_gauche == "Dirichlet")  //Condition de température à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T(i*_Nx) = _sol_T(i*_Nx)-beta(i*_Nx)*_Val_CL_gauche;
    }
  }
  if (_CL_droite == "Dirichlet") //Condition de température à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T((i+1)*_Nx - 1) = _sol_T((i+1)*_Nx - 1)-beta((i+1)*_Nx - 1)*_Val_CL_droite;
    }
  }

  if (_CL_haut == "Neumann") //Condition de flux en haut
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(j) = _sol_T(j)-gamma(j)*_Val_CL_haut*_h_y;
    }
  }
  if (_CL_bas == "Neumann") //Condition de flux en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol_T(_Nx*(_Ny -1)+ j) = _sol_T(_Nx*(_Ny -1)+ j)-gamma(_Nx*(_Ny -1)+ j)*_Val_CL_bas*_h_y;
    }
  }

  if (_CL_gauche == "Neumann")  //Condition de flux à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T(i*_Nx) = _sol_T(i*_Nx)-beta(i*_Nx)*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_gauche == "Neumann_non_constant")  //Condition de flux à gauche
  {
    Laplacian2D::UpdateCL(num_it);
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T(i*_Nx) = _sol_T(i*_Nx)-beta(i*_Nx)*_Val_CL_gauche*_h_x;
    }
  }

  if (_CL_droite == "Neumann") //Condition de flux à droite
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol_T((i+1)*_Nx - 1) = _sol_T((i+1)*_Nx - 1)-beta((i+1)*_Nx - 1)*_Val_CL_droite*_h_x;
    }
  }
}

void EC_PyrolyseMV::InitializeMatrix()
{
  _LapMat.resize(_Nx*_Ny,_Nx*_Ny);

  Eigen::VectorXd alpha;
  Eigen::VectorXd beta_g;
  Eigen::VectorXd beta_d;
  Eigen::VectorXd gamma_h;
  Eigen::VectorXd gamma_b;

  alpha.resize(_Nx*_Ny);
  beta_g.resize(_Nx*_Ny);
  beta_d.resize(_Nx*_Ny);
  gamma_h.resize(_Nx*_Ny);
  gamma_b.resize(_Nx*_Ny);



  for (int i=0;i<_Nx*_Ny;i++)
    {
      _Lambda = _lambdaMV(i);
      double mh = _Lambda; //Moyenne des lambda à l'interface Nord
      double mb = _Lambda; //Interface Sud
      double mg = _Lambda; //Interface Ouest
      double md = _Lambda; //Interface Est

      if (i >= _Nx) //On n'est pas sur le bord haut
      {
        mh = 0.5*(_Lambda+_lambdaMV(i-_Nx));
      }
      if (i <= _Nx*_Ny-1 - _Nx) //On n'est pas sur le bord bas
      {
        mb = 0.5*(_Lambda+_lambdaMV(i+_Nx));
      }
      if (i%_Nx != 0)
      {
        mg = 0.5*(_Lambda+_lambdaMV(i-1));
      }
      if (i%_Nx != _Nx-1)
      {
        md = 0.5*(_Lambda+_lambdaMV(i+1));
      }
      alpha[i] = _RhoTilde(i)*_CpMV(i)/_deltaT -(mh+mb)/(_h_y*_h_y) - (mg+md)/(_h_x*_h_x);
      beta_g[i] = mg/(_h_x*_h_x);
      beta_d[i] = md/(_h_x*_h_x);
      gamma_h[i] = mh/(_h_y*_h_y);
      gamma_b[i] = mb/(_h_y*_h_y);
    }


  vector<Triplet<double>> liste_elem;

  for (int i = 0 ; i<_Nx*_Ny ; i++)
  {
    liste_elem.push_back({i,i,alpha[i]});
  }

  for (int i = 0 ; i<_Ny*_Nx-1; i++)
  {
    if ((i+1)%_Nx!=0)
    {
      liste_elem.push_back({i,i+1,beta_d[i]});
      liste_elem.push_back({i+1,i,beta_g[i+1]});
    }
  }
  for (int i = 0 ; i<_Nx*(_Ny-1) ; i++)
  {
    liste_elem.push_back({i,_Nx+i,gamma_b[i]});
    liste_elem.push_back({_Nx+i,i,gamma_h[_Nx+i]});
  }

  _LapMat.setFromTriplets(liste_elem.begin(), liste_elem.end());


  if (_CL_gauche == "Neumann" or _CL_gauche == "Neumann_non_constant")
   {
     for (int i = 0 ; i < _Ny; i++)
     {
       _LapMat.coeffRef(_Nx*i,_Nx*i) += beta_g[_Nx*i];  //Bord gauche
     }
   }

   if (_CL_droite == "Neumann")
   {
     for (int i = 0 ; i < _Ny; i++)
     {
       _LapMat.coeffRef(_Nx*(i+1) - 1, _Nx*(i+1) - 1) += beta_d[_Nx*(i+1) - 1]; //Bord droite
     }
   }


   if (_CL_haut == "Neumann")
   {
     for (int i = 0; i < _Nx ; i++)
     {
       _LapMat.coeffRef(i,i) += gamma_h[i]; //Bord haut
     }
   }

   if (_CL_bas == "Neumann")
   {
     for (int i = 0; i < _Nx ; i++)
     {
       _LapMat.coeffRef((_Ny - 1)* _Nx + i , (_Ny - 1)* _Nx + i) += gamma_b[(_Ny - 1)* _Nx + i]; //Bord bas
     }
   }
}
