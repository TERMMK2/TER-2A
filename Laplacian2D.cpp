#include "Laplacian2D.h"

// Demander à Motte ou il a foutu la CL chelou de la prof


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

      if ((data_file.Get_eq() == "EC_ClassiqueM") or (data_file.Get_eq() == "EC_ClassiqueP"))
	{
	  _a = data_file.Get_lambda()/(data_file.Get_rho()*data_file.Get_Cp());
	}
      //else if : mettre ici les CI pour lambda et rho et def Cp.

    }
  else
    {//Mettre ici le truc pour le restart file
    }
  //---------------------------------------

  //-----Creation des fichiers de sauvegarde------
  

  if (_save_all_file !="non")
    {
      system(("rm -Rf "+_save_all_file).c_str());
      system(("mkdir -p ./"+_save_all_file).c_str());
    }

  if (_save_points_file !="non")
    {
      system(("rm -Rf "+_save_points_file).c_str());
      system(("mkdir -p ./"+_save_points_file).c_str());
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



void EC_ClassiqueM::InitializeMatrix()
{
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

// Faire quelque chose pour fusionner les solveur c'est moche comme ça


void EC_ClassiqueM::DirectSolver (int nb_iterations) //il reste un petit truc à régler la dedans mais (normalement) ça devrait être rapide
{
  SimplicialLLT <SparseMatrix<double> > solver;
  solver.compute(_LapMat);
  
  vector< shared_ptr<ofstream> > mes_flux;
  
  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  shared_ptr<ofstream> flux(new ofstream);
	  flux->open(_save_points_file+"/point_"+to_string(i), ios::out);
	  mes_flux.push_back(flux);
	}
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
	  char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      double truc_b1, truc_b2;
	      truc_b1 = i*_deltaT;
	      int pos = floor((_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y));
	      truc_b2 = _sol(pos) ;
	      sprintf(truc, "%f  %f", truc_b1, truc_b2);
	      mes_flux[j]->write(truc,16);
	      mes_flux[j]->write("\n",1);
	      //Bon alors je vais expliquer un peu le bordel que c'est :
	      //En gros je suis obligé de faire des pointeurs pour ouvrir un certain nombre de fichier dynamiquement 
	      //Et à cause de ça je peux pas faire << comme d'hab
	      //Donc j'utilise write (qui marche avec les pointeurs mais qui est plus chiant à utiliser) et c'est pour ça que c'est dégueu
	      //Mais en gros ça marche comme ça (d'après ce que j'ai pigé) :
	      //Tu utilise sprintf pour transformer tes nombres en des chaine de charactères puis write pour l'écrire dans ton fichier
	      //Voilà. (Si vous avez besoin d'y retoucher demander moi avant svp)
	    }
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
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  mes_flux[i]->close();
	}
    }
  

}

void EC_ClassiqueP::DirectSolver (int nb_iterations)
{
  SimplicialLLT <SparseMatrix<double> > solver;
  solver.compute(_LapMat);
  
  vector< shared_ptr<ofstream> > mes_flux;
  
  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  shared_ptr<ofstream> flux(new ofstream);
	  flux->open(_save_points_file+"/point_"+to_string(i), ios::out);
	  mes_flux.push_back(flux);
	}
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
	  char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      double truc_b1, truc_b2;
	      truc_b1 = i*_deltaT;
	      int pos = floor((_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y));
	      truc_b2 = _sol(pos) ;
	      sprintf(truc, "%f  %f", truc_b1, truc_b2);
	      mes_flux[j]->write(truc,16);
	      mes_flux[j]->write("\n",1);
	      //Bon alors je vais expliquer un peu le bordel que c'est :
	      //En gros je suis obligé de faire des pointeurs pour ouvrir un certain nombre de fichier dynamiquement 
	      //Et à cause de ça je peux pas faire << comme d'hab
	      //Donc j'utilise write (qui marche avec les pointeurs mais qui est plus chiant à utiliser) et c'est pour ça que c'est dégueu
	      //Mais en gros ça marche comme ça (d'après ce que j'ai pigé) :
	      //Tu utilise sprintf pour transformer tes nombres en des chaine de charactères puis write pour l'écrire dans ton fichier
	      //Voilà. (Si vous avez besoin d'y retoucher demander moi avant svp)
	    }
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
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  mes_flux[i]->close();
	}
    }

}


void EC_ClassiqueM::IterativeSolver (int nb_iterations)
{
  ConjugateGradient <SparseMatrix<double> > solver;
  solver.compute(_LapMat);
  
  vector< shared_ptr<ofstream> > mes_flux;
  
  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  shared_ptr<ofstream> flux(new ofstream);
	  flux->open(_save_points_file+"/point_"+to_string(i), ios::out);
	  mes_flux.push_back(flux);
	}
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
	  char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      double truc_b1, truc_b2;
	      truc_b1 = i*_deltaT;
	      int pos = floor((_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y));
	      truc_b2 = _sol(pos) ;
	      sprintf(truc, "%f  %f", truc_b1, truc_b2);
	      mes_flux[j]->write(truc,16);
	      mes_flux[j]->write("\n",1);
	      //Bon alors je vais expliquer un peu le bordel que c'est :
	      //En gros je suis obligé de faire des pointeurs pour ouvrir un certain nombre de fichier dynamiquement 
	      //Et à cause de ça je peux pas faire << comme d'hab
	      //Donc j'utilise write (qui marche avec les pointeurs mais qui est plus chiant à utiliser) et c'est pour ça que c'est dégueu
	      //Mais en gros ça marche comme ça (d'après ce que j'ai pigé) :
	      //Tu utilise sprintf pour transformer tes nombres en des chaine de charactères puis write pour l'écrire dans ton fichier
	      //Voilà. (Si vous avez besoin d'y retoucher demander moi avant svp)
	    }
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
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  mes_flux[i]->close();
	}
    }

}

void EC_ClassiqueP::IterativeSolver (int nb_iterations)
{
  ConjugateGradient <SparseMatrix<double> > solver;
  solver.compute(_LapMat);
  
  vector< shared_ptr<ofstream> > mes_flux;
  
  if(_save_points_file != "non")
    {
      //Si on sauvegarde des points en particulier, on initialise l'ouverture des fichiers ici.
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  shared_ptr<ofstream> flux(new ofstream);
	  flux->open(_save_points_file+"/point_"+to_string(i), ios::out);
	  mes_flux.push_back(flux);
	}
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
	  char* truc = new char;
	  for (int j=0; j<_number_saved_points; j++)
	    {
	      double truc_b1, truc_b2;
	      truc_b1 = i*_deltaT;
	      int pos = floor((_saved_points[j][0]/_h_x) + _Nx*floor(_saved_points[j][1]/_h_y));
	      truc_b2 = _sol(pos) ;
	      sprintf(truc, "%f  %f", truc_b1, truc_b2);
	      mes_flux[j]->write(truc,16);
	      mes_flux[j]->write("\n",1);
	      //Bon alors je vais expliquer un peu le bordel que c'est :
	      //En gros je suis obligé de faire des pointeurs pour ouvrir un certain nombre de fichier dynamiquement 
	      //Et à cause de ça je peux pas faire << comme d'hab
	      //Donc j'utilise write (qui marche avec les pointeurs mais qui est plus chiant à utiliser) et c'est pour ça que c'est dégueu
	      //Mais en gros ça marche comme ça (d'après ce que j'ai pigé) :
	      //Tu utilise sprintf pour transformer tes nombres en des chaine de charactères puis write pour l'écrire dans ton fichier
	      //Voilà. (Si vous avez besoin d'y retoucher demander moi avant svp)
	    }
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
      for (int i=0; i<_number_saved_points; i++)
	{ 
	  mes_flux[i]->close();
	}
    }

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
      _sol(j) = _sol(j)-gamma*temp(j) + gamma*_Val_CL_haut*_h_y;
    }
  }

  if (_CL_bas == "Neumann") //Condition de flux en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(_Nx*(_Ny -1)+ j) = _sol(_Nx*(_Ny -1)+ j)-gamma*temp(_Nx*(_Ny -1)+ j) + gamma*_Val_CL_bas*_h_y;
    }
  }

  if (_CL_gauche == "Neumann") //Condition de flux à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)-beta*temp(i*_Nx) + beta*_Val_CL_gauche*_h_x;
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
      _sol((i+1)*_Nx - 1) = _sol((i+1)*_Nx - 1)-beta*temp((i+1)*_Nx - 1) + beta*_Val_CL_droite*_h_x;
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
      _sol(j) = _sol(j)+gamma*_Val_CL_haut*_h_y;
    }
  }
  if (_CL_bas == "Neumann") //Condition de flux en bas
  {
    for (int j = 0; j < _Nx ; j++)
    {
      _sol(_Nx*(_Ny -1)+ j) = _sol(_Nx*(_Ny -1)+ j)+gamma*_Val_CL_bas*_h_y;
    }
  }

  if (_CL_gauche == "Neumann")  //Condition de flux à gauche
  {
    for (int i = 0; i < _Ny; i++)
    {
      _sol(i*_Nx) = _sol(i*_Nx)+beta*_Val_CL_gauche*_h_x;
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
      _sol((i+1)*_Nx - 1) = _sol((i+1)*_Nx - 1)+beta*_Val_CL_droite*_h_x;
    }
  }
}

