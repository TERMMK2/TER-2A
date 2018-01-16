#ifndef FILE_DATA_FILE_H  // Pas a jour du tout

#include <string>
#include <iostream>
#include <Eigen>

class DataFile
{
 private:
  //Nom du fichier de donnée d'entrée
  const std::string _file_name;

  // Variables à utiliser // Ajouter un truc pour le systeme d'adaptation du maillage
  std::string _CL_droite;
  std::string _CL_haut;
  std::string _CL_bas;
  std::string _CL_gauche;
  double _Val_CL_droite;
  double _Val_CL_haut;
  double _Val_CL_bas;
  double _Val_CL_gauche;

  double _CI; // Voir si besoin de changer ça ou pas (CI temperature constante)
  
  std::string _eq;
  
  int _N_x;
  int _N_y;
  double _x_min;
  double _x_max;
  double _y_min;
  double _y_max;
  double _deltaT;
  double _T_final;

  std::string _Solveur;
  std::string _Schema; // A voir si vraiment utile !

  std::string _save_all_file;  // Nom du fichier qui contiendra la solution pour tout points
  std::string _save_points_file; // Nom du ou des fichiers qui contiendra l'évolution de T en fonction du temps pour 1 unique point ( si plusieurs points selectionnés les fichiers s'apellerons tous comme _save_points_file + le numéro du point
  int _number_saved_points;
  std::vector<std::vector<double>> _saved_points;
  
  
  


  // Pour savoir si l'utilisateur a donner les paramètres 
  // ou si il faut utiliser les paramètres par défaut.
  bool _if_CL_droite;
  bool _if_Cl_gauche;
  bool _if_Cl_haut;
  bool _if_CL_bas;
  bool _if_Val_CL_droite;
  bool _if_Val_Cl_gauche;
  bool _if_Val_Cl_haut;
  bool _if_Val_CL_bas;
  bool _if_CI;
  bool _if_eq;
  bool _if_N_x;
  bool _if_N_y;
  bool _if_x_min;
  bool _if_x_max;
  bool _if_y_min;
  bool _if_y_max;
  bool _if_deltaT;
  bool _if_T_final;

  bool _if_Solveur;
  bool _if_Schema; // a voir si utile ou pas pas fait ça et la suite

  bool _if_save_all_file;

  bool _if_save_points_file;
  bool _if_number_saved_points;
  bool _if_saved_points;


 public:

  // Constructeur
  DataFile(std::string file_name);

  // Lecture du fichier de donnée // Pas fait encore ça !
  void ReadDataFile();

  inline int Get_n() {return _n; };
  inline Eigen::MatrixXd Get_A() {return _A; };
  inline std::string Get_Solveur() {return _Solveur;};
  inline std::string Get_prec() {return _prec;};
  inline std::string Get_methode_prec() {return _methode_prec;};
  inline double Get_arg_prec() {return _arg_prec;};
  inline int Get_m() {return _m;};
  inline double Get_eps() {return _eps;};
  inline int Get_kmax() {return _kmax;};
  inline bool Get_check() {return _check;}; 
  inline std::string Get_save_r() {return _save_r;};
  inline std::string Get_save_sol() {return _save_sol;};

};

#define FILE_DATA_FILE_H
#endif
