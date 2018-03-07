#include "Sparse"
#include "Dense"

class Laplacian2D // pas fini de modifier
{
  protected: // Les attributs de la classe
    double _x_min, _x_max, _y_min, _y_max, _h_x, _h_y, _a, _deltaT;
    int _Nx, _Ny;
    Eigen::SparseMatrix<double> _LapMat; // matrice creuse du laplacien
    Eigen::VectorXd _x, _y; // points de discretisation x et y
    Eigen::VectorXd _f; // vecteur source _f qui prend les données de _sol(i) pour calculer _sol(i+1)
    Eigen::VectorXd _sol; // vecteur solution U
    std::string _CL_bas, _CL_haut, _CL_gauche, _CL_droite;
    double _Val_CL_bas, _Val_CL_haut, _Val_CL_gauche, _Val_CL_droite;

  public: // Méthodes et opérateurs de la classe
    Laplacian2D();
    // Constructeur : Initialiser _x_min, _x_max, _y_min; _y_max; _N; _h; _LapMat; _x; _y et _sol.
    virtual ~Laplacian2D();

    void Initialize(double x_min, double x_max, double y_min, double y_max, int Nx, int Ny, double a, double deltaT);
    void InitializeCL(std::string CL_bas, std::string CL_haut, std::string CL_gauche, std::string CL_droite, double Val_CL_bas, double Val_CL_haut, double Val_CL_gauche, double Val_CL_droite);
    void UpdateCL(int num_it);
    virtual void InitializeMatrix() = 0;
    void InitializeCI(Eigen::VectorXd CI);
    virtual void DirectSolver(int nb_iterations) = 0;   // Résout le système _LapMat * _sol = _f avec un solveur direct.
    virtual void IterativeSolver(int nb_iterations) = 0;   // Résout le système _LapMat * _sol = _f avec un solveur itératif.
    void SaveSol(std::string name_file); // Écrit les solutions dans le fichier "name_file".
    virtual void ConditionsLimites(int num_it) = 0;
  };

class EC_ClassiqueM : public Laplacian2D //Première version avec un matrice identique quelles que soient les conditions aux bords
{
  public:
    void InitializeMatrix();
    void DirectSolver(int nb_iterations);
    void IterativeSolver(int nb_iterations);
    void ConditionsLimites(int num_it);
};

class EC_ClassiqueP : public Laplacian2D //Seconde version avec une matrice qui dépend des conditions aux bords
{
  public:
    void InitializeMatrix();
    void DirectSolver(int nb_iterations);
    void IterativeSolver(int nb_iterations);
    void ConditionsLimites(int num_it);
};
