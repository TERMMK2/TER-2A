#ifndef FILE_DATA_FILE_CPP //ça compile et tout mais j'ai pas l'impression que tout marche parfaitement il faudrait reprendre les trucs 

#include "Datafile.h"
#include <fstream>
#include <iostream>


using namespace std;
using namespace Eigen;

// Constructeur

DataFile::DataFile(string file_name)
  : _file_name(file_name), _if_CL_droite(false), _if_CL_gauche(false), _if_CL_haut(false), _if_CL_bas(false),_if_Val_CL_droite(false), _if_Val_CL_gauche(false), _if_Val_CL_haut(false), _if_Val_CL_bas(false), _if_CI(false), _if_eq(false), _if_N_x(false), _if_N_y(false), _if_x_min(false), _if_x_max(false), _if_y_min(false), _if_y_max(false), _if_deltaT(false), _if_T_final(false), _if_Solveur(false), _if_Schema(false), _if_save_all_file(false), _if_save_points_file(false), _if_number_saved_points(false), _if_saved_points(false), _if_restart_file(false)
{}


//Lecture du fichier de donnée

void DataFile::ReadDataFile()
{
  ifstream data_file(_file_name.data());
  if (!data_file.is_open())
    {
      cout << "Unable to open file " << _file_name << endl;
      abort();
    }
  else
    {
      cout << "--------------------------------------------------" << endl;
      cout << "Reading data file " << _file_name << endl;
    }
  
  string file_line; 
  
  string truc;
  size_t taille;
  char carac;

  // Parcourir le fichier pour rechercher les paramètres
  while (!data_file.eof())
    {
      getline(data_file, file_line);
      cout << file_line <<endl;
      if (file_line.find("CL droite") != std::string::npos)
	{
	  data_file >> _CL_droite; _if_CL_droite = true;
	  data_file >> truc;
	  taille = truc.size();
	  data_file.seekg(-(taille),ios::cur);
	  data_file.get(carac);
	  cout <<"carac :" <<carac<<"|  bool :"<<((carac == '1') or (carac == '2') or (carac == '3') or (carac == '4') or (carac == '5') or (carac == '6') or  (carac == '7') or (carac == '8') or (carac == '9') or (carac == '0'))<< endl;
	  data_file.seekg(-1,ios::cur);
	  if ((carac == '1') or (carac == '2') or (carac == '3') or (carac == '4') or (carac == '5') or (carac == '6') or  (carac == '7') or (carac == '8') or (carac == '9') or (carac == '0'))
	    {
	      _Val_CL_droite = atoi(truc.c_str());
	      _if_Val_CL_droite = true;
	      data_file.seekg(taille,ios::cur);
	    }
	  else
	    {
	      data_file.seekg(-taille,ios::cur);
	    }
	  
	}


      if (file_line.find("CL gauche") != std::string::npos)
	{ 
	  data_file >> _CL_gauche; _if_CL_gauche = true;
	  data_file >> truc;
	  taille = truc.size();
	  data_file.seekg(-(taille),ios::cur);
	  data_file.get(carac);
	  data_file.seekg(-1,ios::cur);
	  if ((carac == '1') or (carac == '2') or (carac == '3') or (carac == '4') or (carac == '5') or (carac == '6') or  (carac == '7') or (carac == '8') or (carac = '9') or (carac == '0'))
	    {
	      _Val_CL_gauche = atoi(truc.c_str());
	      _if_Val_CL_gauche = true;
	      data_file.seekg(taille,ios::cur);
	    }
	  else
	    {data_file.seekg(-taille,ios::cur);}
	}


      if (file_line.find("CL haut") != std::string::npos)
	{
	  data_file >> _CL_haut; _if_CL_haut = true;
	  data_file >> truc;
	  taille = truc.size();
	  data_file.seekg(-(taille),ios::cur);
	  data_file.get(carac);
	  data_file.seekg(-1,ios::cur);
	  if ((carac == '1') or (carac == '2') or (carac == '3') or (carac == '4') or (carac == '5') or (carac == '6') or  (carac == '7') or (carac == '8') or (carac = '9') or (carac == '0'))
	    {
	      _Val_CL_haut = atoi(truc.c_str());
	      _if_Val_CL_haut = true;
	      data_file.seekg(taille,ios::cur);
	    }
	  else
	    {data_file.seekg(-taille,ios::cur);}
	}


      if (file_line.find("CL bas") != std::string::npos)
	{
	  data_file >> _CL_bas; _if_CL_bas = true;
	  data_file >> truc;
	  taille = truc.size();
	  data_file.seekg(-(taille),ios::cur);
	  data_file.get(carac);
	  data_file.seekg(-1,ios::cur);
	  if ((carac == '1') or (carac == '2') or (carac == '3') or (carac == '4') or (carac == '5') or (carac == '6') or  (carac == '7') or (carac == '8') or (carac = '9') or (carac == '0'))
	    {
	      _Val_CL_bas = atoi(truc.c_str());
	      _if_Val_CL_bas = true;
	      data_file.seekg(taille,ios::cur);
	    }
	  else
	    {data_file.seekg(-taille,ios::cur);}
	}

      
      if (file_line.find("CI") != std::string::npos)
	{
	  data_file >> _CI; _if_CI = true;
	}
      
      if (file_line.find("Equation") != std::string::npos)
      	{
      	  data_file >> _eq; _if_eq = true;
      	}
      

      if (file_line.find("N_x") != std::string::npos)
	{
	  data_file >> _N_x; _if_N_x = true;
	}
      if (file_line.find("N_y") != std::string::npos)
	{
	  data_file >> _N_y; _if_N_y = true;
	}
      if (file_line.find("x_min") != std::string::npos)
	{
	  data_file >> _x_min; _if_x_min = true;
	}
      if (file_line.find("x_max") != std::string::npos)
	{
	  data_file >> _x_max; _if_x_max = true;
	}
      if (file_line.find("y_min") != std::string::npos)
	{
	  data_file >> _y_min; _if_y_min = true;
	}
      if (file_line.find("y_max") != std::string::npos)
	{
	  data_file >> _y_max; _if_y_max = true;
	}
      if (file_line.find("delta t") != std::string::npos)
	{
	  data_file >> _deltaT; _if_deltaT = true;
	}
      if (file_line.find("T final") != std::string::npos)
	{
	  data_file >> _T_final; _if_T_final = true;
	}

      
      if (file_line.find("Solveur") != std::string::npos)
	{
	  data_file >> _Solveur; _if_Solveur = true;
	}
      if (file_line.find("Schema") != std::string::npos)
	{
	  data_file >> _Schema; _if_Schema = true;
	}
      

      if (file_line.find("Fichier Paraview") != std::string::npos)
	{
	  data_file >> _save_all_file; _if_save_all_file = true;
	}
      
      if (file_line.find("Fichier sauvegarde points") != std::string::npos)
	{
	  data_file >> _save_points_file; _if_save_points_file = true;
	}
      if (file_line.find("Nombre de points à sauver") != std::string::npos)
	{
	  data_file >> _number_saved_points; _if_number_saved_points = true;
	  _saved_points.resize(_number_saved_points);
	  for (int i = 0 ; i<_number_saved_points; i++) 
	    { _saved_points[i].resize(2);  }
	}
      if (file_line.find("Points à sauver") != std::string::npos)
	{
	  int i =0;
	  bool test = true;
	  
	  while ((i<2*_number_saved_points) and (test)) 
	    {
	      data_file >> truc;
	      taille = truc.size();
	      data_file.seekg(-(taille),ios::cur);
	      data_file.get(carac);
	      data_file.seekg(-1,ios::cur);
	      if ((carac == '1') or (carac == '2') or (carac == '3') or (carac == '4') or (carac == '5') or (carac == '6') or  (carac == '7') or (carac == '8') or (carac == '9') or (carac == '0'))
		{  
		  _saved_points[i/2][i%2] = atof(truc.c_str());
		  data_file.seekg(taille,ios::cur);
		}
	      else 
		{ 
		  test = false;
		  data_file.seekg(-(taille),ios::cur);
		}
	      i++;
	    }
	  if (test)
	      _if_saved_points = true;
	}  

      if(file_line.find("Fichier de redemarrage") != std::string::npos)
	{
	  data_file >> _restart_file; _if_restart_file = true;
	}


    }
  



  // Initialisation par défaut des paramètres non fixés dans le fichier
  // Un message prévient l'utilisateur

  if ((!_if_CL_droite) or ((_CL_droite != "Dirichlet") and (_CL_droite != "Neumann")))
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention la valeur par défaut pour la CL droite (Neumann) est utilisé." << endl;
      _CL_droite = "Neumann";
    }
  if ((!_if_CL_gauche) or ((_CL_gauche != "Dirichlet") and (_CL_gauche != "Neumann")))
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention la valeur par défaut pour la CL gauche (Dirichlet) est utilisé." << endl;
      _CL_gauche = "Dirichlet";
    }
  if ((!_if_CL_haut) or ((_CL_haut != "Dirichlet") and (_CL_haut != "Neumann")))
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention la valeur par défaut pour la CL haut (Neumann) est utilisé." << endl;
      _CL_haut = "Neumann";
    }
  if ((!_if_CL_bas) or ((_CL_bas != "Dirichlet") and (_CL_bas != "Neumann")))
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention la valeur par défaut pour la CL bas (Neumann) est utilisé." << endl;
      _CL_bas = "Neumann";
    }
  

  
  if (!_if_Val_CL_droite)
    {
      if (_CL_droite =="Neumann")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL droite  pour "<<_CL_droite<<" (0) est utilisée."  << endl;
	  _Val_CL_droite = 0;
	}
      if (_CL_droite == "Dirichlet")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL droite pour "<<_CL_droite<<" (2900) est utilisée."  << endl;
	  _Val_CL_droite = 2900;
	}
    }
  if (!_if_Val_CL_gauche)
    {
      if (_CL_gauche =="Neumann")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL gauche pour "<<_CL_gauche<<" (0) est utilisée."  << endl;
	  _Val_CL_gauche = 0;
	}
      if (_CL_gauche == "Dirichlet")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL gauche pour "<<_CL_gauche<<" (2900) est utilisée."  << endl;
	  _Val_CL_gauche = 2900;
	}
    }
  if (!_if_Val_CL_haut)
    {
      if (_CL_haut =="Neumann")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL haut pour  "<<_CL_haut<<" (0) est utilisée."  << endl;
	  _Val_CL_haut = 0;
	}
      if (_CL_haut == "Dirichlet")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL haut pour "<<_CL_haut<<" (2900) est utilisée."  << endl;
	  _Val_CL_haut = 2900;
	}
    }
  if (!_if_Val_CL_bas)
    {
      if (_CL_bas =="Neumann")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL bas pour "<<_CL_bas<<" (0) est utilisée."  << endl;
	  _Val_CL_bas = 0;
	}
      if (_CL_bas == "Dirichlet")
	{
	  cout << "---------------------------------------------------" << endl;
	  cout << "Attention, la valeur par défaut de la CL bas pour "<<_CL_bas<<" (2900) est utilisée."  << endl;
	  _Val_CL_bas = 2900;
	}
    }


  
  if (!_if_CI)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, la valeur par défaut pour la CI (293) est utilisé. " << endl;
      _CI = 293;      
    }

  if ((!_if_eq) or ((_eq != "EC_ClassiqueM") and (_eq != "EC_ClassiqueP")))
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, vous n'avez pas entré d'équation ou celle que vous avez choisie n'existe pas." <<endl;
      cout <<"L'equation par défaut (EC_ClassiqueP) va donc être utilisée." << endl;
      _eq = "EC_ClassiqueP";
    }


  
  if (!_if_N_x)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le N_x par défaut est utilisé (100).----" << endl;
      _N_x = 100;
    }
  if (!_if_N_y)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le N_y par défaut est utilisé (100).----" << endl;
      _N_y = 100;
    }
  if (!_if_x_min)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le x_min par défaut est utilisé (0).----" << endl;
      _x_min = 0;
    }
  if (!_if_x_max)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le x_max par défaut est utilisé (1).----" << endl;
      _x_max = 1;
    }
  if (!_if_y_min)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le y_min par défaut est utilisé (0).----" << endl;
      _y_min = 0;
    }
  if (!_if_y_max)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le y_max par défaut est utilisé (1).----" << endl;
      _y_max = 100;
    }
  if (!_if_deltaT)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le deltaT par défaut est utilisé (0.01)." << endl;
      _deltaT = 0.01;
    }
  if (!_if_T_final)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le T_final par défaut est utilisé (10).-" << endl;
      _T_final = 10;
    }

  
  if (!_if_Solveur)
    {
      cout << "--------------------------------------------------" << endl;
      cout << "Attention, le solveur par défaut (BiCGStab) est utilisé." << endl;
      _Solveur = "BiCGStab";
    }

  
  if (!_if_Schema)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, le schéma n'est pas défini (pas grave y'en a pas pour l'instant ;)" << endl;
      _Schema = "tkt_bro";
    }
  
  if (!_if_save_all_file)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, la solution ne sera pas sauvegardé dans sa totalité" << endl;
      _save_all_file = "non";
    }

  if (!_if_save_points_file)
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, il manque un ou plusieurs points à sauvegarder" << endl;
      cout << "Nous n'allons donc pas sauvegarder de points de points" <<endl;
      _save_points_file = "non";
    }

  if ((!_if_number_saved_points) and (_if_save_points_file))
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, vous n'avez pas entré de nombre de points à sauvegarder " << endl;
      cout << "Nous n'allons donc pas sauvegarder de points-------" <<endl;
      _save_points_file = "non";
    }
  
  if ((_save_points_file != "non") and (!_if_saved_points))
    {
      cout << "---------------------------------------------------" << endl;
      cout << "Attention, vous n'avez pas entré de points à sauvegarder (ou vous en avez oublier un)" << endl;
      cout << "Nous n'allons donc pas sauvegarder de points-------" << endl;
      _save_points_file = "non";
    }


  if ((_save_points_file == "non") and (_save_all_file == "non"))
    {
      // Potentiellement ajouter des trucs qui t'affiche ton erreur et tout avec exit(1) -> a voir !
      cout << "---------------------------------------------------" <<endl;
      cout << "--------------------Attention----------------------" <<endl;
      cout << "---------------------------------------------------" <<endl;
      cout << "---------------------------------------------------" <<endl<<endl;
      cout << "----Erreur critique, auto-destruction en cours-----" <<endl<<endl;
      cout << "---------------------------------------------------" <<endl;
      cout << "---------------------------------------------------" <<endl<<endl<<endl;
      
      cout << "Non je déconne, mais en vrai t'as pas mis de fichier de sauvegarde donctu risque pas de voir si ça a marché ou non (bon sauf si tu regarde juste si ta partie du code marche sans erreur, mais tu peux être sûr qu'il n'y a pas d'erreur si tu regarde pas la solution non ?)" <<endl<<endl;
      cout << "Bref, je vais pas faire un stop au programme parce que je sais pas le faire mais sache que si ton PC fait des calculs en se moment, ils servent à rien."<<endl;
      cout << "Sur ce mon message est fini, bonne journée à toi, jeune utilisateur quelque peu étourdi"<<endl;
      cout << "P.S. : saucissons" <<endl;
      cout << "P.P.S. : 2" <<endl;
      cout << "P.P.P.S. : 50" <<endl;
      cout << "P.P.P.P.S. : 1729"<<endl;
      cout << "P.P.P.P.P.S. : 635318657" <<endl;
      cout << "P.P.P.P.P.P.S. : Saura tu trouver de quelle suite il s'agit ?" << endl;
      cout << "P.P.P.P.P.P.P.S. : Ainsi que le terme suivant ?" <<endl;
      cout << "P.P.P.P.P.P.P.P.S. : En vrai te fatique pas c'est un problème ouvert de math de calculer le suivant, alors à moins que tu as 2 ans devant toi et un super calculateur tu trouvera surrement pas sans une idée de génie."<<endl;
      cout << "P.P.P.P.P.P.P.P.P.S. : J'ai perdu ma soirée à essayer de le trouver avant de m'en rendre compte. :'(" <<endl;
      cout << "P.P.P.P.P.P.P.P.P.P.S. : Mais j'ai appris plein de truc, ce qui est cool." <<endl;
      cout << "P.P.P.P.P.P.P.P.P.P.P.S : Sur ce je vais me coucher en essayant de faire marcher git" <<endl;
    }


  if( !_if_restart_file )
    {
      cout << "Attention, pas de fichier de redémarage--------------"<<endl;
      _restart_file = "non";
    }
}


	
#define FILE_DATA_FILE_CPP
#endif
