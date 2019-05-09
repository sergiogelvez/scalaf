// estas son las constantes físicas necesarias para los cálculos.
// algunas se redefinen por parámetros de entrada del programa.
#define density 2500.0
#define gravity 9.8
#define cellWidth 10.0
#define heatCapacity 840.0
#define emisivity 0.9
#define SBConst 0.0000000568
#define extrudeRate 130 // 0.2
#define extrudeTemp 1300.0
#define dt 1

// Esta es la estructura de cada celda
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
} celda;

// estructura de los puntos de los cráteres
typedef struct {
  int x;
  int y;
} point2D;

// estructura de las condiciones iniciales
typedef struct {
  int tFilas;
  int tColumnas;
  double anchoCelda;
  double eRate;
  double eTemp;
  double deltat;
  int nPasos;
} cIni;

// prototipos de las funciones principales
void FuncionPrincipal(int filas, int columnas, celda *A, celda *C);
int leerArchivoTexto_Matriz(char *path, int filas, int columnas, celda *matriz);
void preFuncion(int, int, const celda *, celda *);
void postFuncion(int, int, const celda *A, celda *C);
int leerArchivoPuntos(char *, int, point2D *);
int colocarCrateres(celda *, const point2D *, int, int, int);

// funciones de utilidades, prototipos
int limpiarPath(char[], char[]);
int obtenerPath(char[]);
void imprimirMatrizPantalla(int, int, const celda *, int);
void imprimirMatrizPantalla_2(int, int, const celda *, int);
int prepararVisualizacionGNUPlot(int, char *, int, int, celda *, int, double,
                                 double, double);
int prepararVisualizacionGNUPlot_2(int, char *, int, int, celda *, int, double,
                                   double, double);

// funciones de cálculo de valores físicos
double visc(double);
double yield(double);

// funciones de prueba o incompletas, prototipos
int generarAnimacionGNUPlot(char[], int);
void testAnimacion(void);