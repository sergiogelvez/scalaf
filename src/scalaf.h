// estas son las constantes físicas necesarias para los cálculos.
// algunas se redefinen por parámetros de entrada del programa.
#define density 2500.0
#define gravity 9.8
#define heatCapacity 840.0
#define emisivity 0.9
#define SBConst 0.0000000568
#define time_delta 1

// Esta es la estructura de cada mapCell
typedef struct {
  double thickness;
  double temperature;
  double altitude;
  int isVent;
  double yield;
  double viscosity;
  short exits;
  double inboundV;
  double outboundV;
  double inboundQ;
} mapCell;

// estructura de los puntos de los cráteres
typedef struct {
  int x;
  int y;
} point2D;

// estructura de las condiciones iniciales
typedef struct {
  int maxRows;
  int maxColumns;
  double cellWidth;
  double eruptionRate;
  double eruptionTemperature;
  double deltat;
  int timeSteps;
} initialConditions;

// prototipos de las funciones principales
void FuncionPrincipal(int filas, int columnas, mapCell *A, mapCell *C);
int leerArchivoTexto_Matriz(char *path, int filas, int columnas,
                            mapCell *matriz);
void preFuncion(int, int, const mapCell *, mapCell *);
void postFuncion(int, int, const mapCell *A, mapCell *C);
int leerArchivoPuntos(char *, int, point2D *);
int colocarCrateres(mapCell *, const point2D *, int, int, int);

// funciones de utilidades, prototipos
int limpiarPath(char[], char[]);
int obtenerPath(char[]);
void imprimirMatrizPantalla(int, int, const mapCell *, int);
void imprimirMatrizPantalla_2(int, int, const mapCell *, int);
int prepararVisualizacionGNUPlot(int, char *, int, int, mapCell *, int, double,
                                 double, double);
int prepararVisualizacionGNUPlot_2(int, char *, int, int, mapCell *, int,
                                   double, double, double);

// funciones de cálculo de valores físicos
double visc(double);
double yield(double);

// funciones de prueba o incompletas, prototipos
int generarAnimacionGNUPlot(char[], int);
void testAnimacion(void);