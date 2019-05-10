/*
Código secuencial para la simulación de flujos de lava volcánica
                        -Sergio Augusto Gélvez Cortés
*/

#include "scalaf.h"
#include "math.h"
#include "string.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Función que calcula la viscosidad a partir de la temperatura,
// tomada del artículo de Miyamoto y Sasaki
double visc(double temperature) {
  return pow(10, (20.0 * (exp(-0.001835 * (temperature - 273.0)))));
}

// Función que calcula la tensión cortante a partir de la temperatura,
// tomada del artículo de Miyamoto y Sasaki
double yield(double temperature) {
  return pow(10, (11.67 - 0.0089 * (temperature - 273.0)));
}

// Función para colocar los cráteres en la matriz,
// toma los datos de un punto en 2D y revisa que estén en el rango correcto.
int placeCraters(celda *A, const point2D *P, int totalRows, int totalColumns,
                 int totalCraters) {

  int craterColumn, craterRow;
  for (int i = 0; i < totalCraters; ++i) {
    craterRow = P[i].x;
    craterColumn = P[i].y;
    if ((craterRow > -1) && (craterRow < totalRows)) {
      if ((craterColumn > -1) && (craterColumn < totalColumns)) {
        A[craterRow * totalColumns + craterColumn].isVent = 1;
      }
    }
  }
}

// Función para leer puntos en 2D desde un archivo de texto plano.
// X para Filas, Y para Columnas.
int readCratersPositionFile(char *path, int numberOfCraters,
                            point2D *craterPositions) {
  FILE *cratersPositionFile;
  char lineBuffer[1000];
  char *token;
  int error;

  if (cratersPositionFile = fopen(path, "r")) {
    printf("\nLeyendo posición de los crátetes...");
    int i = 0;
    while ((fgets(lineBuffer, 1000, cratersPositionFile) != NULL) &&
           (i < numberOfCraters)) {
      token = strtok(lineBuffer, ",");
      craterPositions[i].x = atol(token);
      token = strtok(NULL, ",");
      craterPositions[i].y = atol(token);
      i += 1;
    }
    printf("\nArchivo leido correctamente.");
    error = 0;

  } else {
    printf("\nERROR: Archivo con posición de los cráteres no encontrado!");
    error = 1;
  }
  fclose(cratersPositionFile);
  return error;
}

// función para leer las altitudes de los puntos de la matriz de un archivo de
// texto plano (las altitudes vienen directamente en un arreglo de dos
// dimensiones que representa la posición de dichas altitudes en la matriz)
int leerArchivoTexto_Matriz(char *path, int filas, int columnas,
                            celda *matriz) {
  size_t sizep;
  char buff[8192];
  char *tok;
  long int i, j;
  int error;
  double elem;
  FILE *datosAltitud;
  sizep = strlen(path);
  if (sizep != 0) {
    printf("\nIntentando leer el archivo de altitudes...\n", path);
    printf("Archivo a leer: %s\n", path);
    datosAltitud = fopen(path, "r");
    if (datosAltitud != NULL) {
      printf("El archivo existe! \n");
      i = 0;
      while (fgets(buff, 1000, datosAltitud) != NULL) {
        // Toma una linea completa y convierte en token
        j = 0;
        tok = strtok(buff, ",");
        while (tok) {
          if (j < columnas) {
            elem = atof(tok);
            matriz[i * columnas + j].altitude = elem;
            // acá se aprovecha para llenar todo de 0
            // en las demás casillas.
            matriz[i * columnas + j].thickness = 0;
            matriz[i * columnas + j].isVent = 0;
            matriz[i * columnas + j].temperature = 273.0;
            matriz[i * columnas + j].yield = 0.0;
            matriz[i * columnas + j].viscosity = 0.0;
            matriz[i * columnas + j].exits = 0;
            j += 1;
          }
          tok = strtok(NULL, ",");
        }
        if (j >= columnas) {
          i += 1;
        }
      }
      printf("Archivo leido correctamente.\n");
      error = 0;
    } else {
      printf("\n** ERROR: El archivo mno existe, o hubo problemas con el "
             "archivo **\n");
      error = 1;
    }
    // cerrar el archivo
    fclose(datosAltitud);
  }
  return (error);
}

// función para imprimir fila a fila, columna por columna los valores de
// altitud, temperatura, grosor de la capa de lava, existencia o no de crater
// y otros de cada celda de la matriz. Es esencialmente una función para
// pruebas.
void imprimirMatrizPantalla(int filas, int columnas, const celda *matriz,
                            int cifrasSignif) {
  /*	función para imprimir toda la matriz en pantalla, para propositos de
   *	de depuración. */
  int i, j;
  double elem = 0.0;
  printf("\n");
  for (i = 0; i < filas; ++i) {
    printf("\nFila %d\n", i);
    for (j = 0; j < columnas; j++) {
      elem = matriz[columnas * i + j].altitude;
      printf("COL %d - alt: %4.2lf,", j, elem);
      elem = matriz[columnas * i + j].temperature;
      printf("temp: %4.2lf,", elem);
      elem = matriz[columnas * i + j].thickness;
      printf("tck: %4.2lf,", elem);
      elem = matriz[columnas * i + j].viscosity;
      printf("visc: %4.2lf,", elem);
      elem = matriz[columnas * i + j].yield;
      printf("yield: %4.2lf\t", elem);
      elem = matriz[columnas * i + j].exits;
      printf("exits: %4.2lf\t", elem);
      if (matriz[i * columnas + j].isVent == 1) {
        printf(" cráter!! ");
      }
      // se pueden agregar acá más propiedades si se requiere.
    }
    printf("\n");
  }
}

// Esta es la versión de imprimir pantalla que sirve con la matriz agrandada
// directamente
void imprimirMatrizPantalla_2(int filas, int columnas, const celda *matriz,
                              int cifrasSignif) {
  /*	función para imprimir toda la matriz en pantalla, para propositos de
   *	de depuración. */
  int i, j;
  double elem = 0.0;
  printf("\n");
  for (i = 1; i < (filas - 1); ++i) {
    printf("\nFila %d\n", i - 1);
    for (j = 1; j < (columnas - 1); j++) {
      elem = matriz[columnas * i + j].altitude;
      printf("COL %d - alt: %4.2lf,", j - 1, elem);
      elem = matriz[columnas * i + j].temperature;
      printf("temp: %4.2lf,", elem);
      elem = matriz[columnas * i + j].thickness;
      printf("tck: %4.8lf,", elem);
      elem = matriz[columnas * i + j].viscosity;
      printf("visc: %4.2lf,", elem);
      elem = matriz[columnas * i + j].yield;
      printf("yield: %4.2lf\t", elem);
      elem = matriz[columnas * i + j].exits;
      printf("exits: %4.2lf\t", elem);
      if (matriz[i * columnas + j].isVent == 1) {
        printf(" cráter!! ");
      }
      // se pueden agregar acá más propiedades si se requiere.
    }
    printf("\n");
  }
}

/*  La función prefunción "agranda" la matriz, colocando una fila y una
        columna al principio y al final (la matriz tiene dimensiones (MAX_ROWS
        +2)*(MAX_COLS+2). De antemano me disculpo por la cantidad de veces que
        imprimo la matriz. */
void preFuncion(int filas, int columnas, const celda *A, celda *C) {
  int i, j, c, f;
  printf("\nAgregando filas y columnas extras en los bordes...\n");
  // cargar matriz a memoria
  // y acá debería agrandar la matriz
  // las celdas extras tienen altitud 100000 metros, grosor de capas 0,
  // temperatura 0
  celda *B;
  B = (celda *)malloc((filas + 2) * (columnas + 2) * sizeof(celda));
  // crear elementos para la matriz agrandada. B[filas+2][columnas+2];
  f = 0;
  for (i = 0; i < filas + 2; ++i) {
    if (!(i == 0 || i == (filas + 1))) {
      c = 0;
      for (j = 0; j < columnas + 2; ++j) {
        if (!(j == 0 || j == (columnas + 1))) {
          // acá van los valores no en los bordes
          B[(columnas + 2) * i + j].altitude = A[(columnas)*f + c].altitude;
          B[(columnas + 2) * i + j].thickness = A[(columnas)*f + c].thickness;
          B[(columnas + 2) * i + j].temperature =
              A[(columnas)*f + c].temperature;
          B[(columnas + 2) * i + j].isVent = A[(columnas)*f + c].isVent;
          B[(columnas + 2) * i + j].yield = 0;
          B[(columnas + 2) * i + j].viscosity = 0;
          B[(columnas + 2) * i + j].exits = 0;
          B[(columnas + 2) * i + j].inboundV = 0;
          B[(columnas + 2) * i + j].outboundV = 0;
          B[(columnas + 2) * i + j].inboundQ = 0;
          c += 1;
          // llenar los valores originales
          // ya que acá se trate de los posiciones diferentes a los
          // nuevos bordes.
        } else {
          // crear las filas y columnas extras
          // en la fila y columna 0, y en la fila y columna
          // filas y columnas, o sea los máximos.
          // esos valores son 0 y altitud muy grande, más grande
          // que la altitud del everest.
          // en este caso son solo los bordes de las "columnas"
          B[(columnas + 2) * i + j].altitude = 100000;
          B[(columnas + 2) * i + j].thickness = 0;
          B[(columnas + 2) * i + j].temperature = 0;
          B[(columnas + 2) * i + j].isVent = 0;
          B[(columnas + 2) * i + j].yield = 0;
          B[(columnas + 2) * i + j].viscosity = 0;
          B[(columnas + 2) * i + j].exits = 0;
          B[(columnas + 2) * i + j].inboundV = 0;
          B[(columnas + 2) * i + j].outboundV = 0;
          B[(columnas + 2) * i + j].inboundQ = 0;
        }
      }
      f += 1;
    } else {
      for (j = 0; j < columnas + 2; ++j) {
        // crear las filas y columnas extras
        // en la fila y columna 0, y en la fila y columna
        // filas y columnas, o sea los máximos.
        // esos valores son 0 y altitud muy grande, más grande
        // que la altitud del everest.
        // en este caso son solo los bordes de las "filas"
        B[(columnas + 2) * i + j].altitude = 100000;
        B[(columnas + 2) * i + j].thickness = 0;
        B[(columnas + 2) * i + j].temperature = 0;
        B[(columnas + 2) * i + j].isVent = 0;
        B[(columnas + 2) * i + j].yield = 0;
        B[(columnas + 2) * i + j].viscosity = 0;
        B[(columnas + 2) * i + j].exits = 0;
        B[(columnas + 2) * i + j].inboundV = 0;
        B[(columnas + 2) * i + j].outboundV = 0;
        B[(columnas + 2) * i + j].inboundQ = 0;
      }
    }
  }
  // acá se copia la nueva matriz a la matrix pasada por referencia
  // esa será el resultado.
  memcpy(C, B, (filas + 2) * (columnas + 2) * sizeof(celda));
  printf("Salir de la función agrandar...\n");
  // En C se guardaron los resultados.
  // Llamamos la función principal, la que hace los cálculos de la
  // simulación:
  // FuncionPrincipal(MAX_ROWS+1, MAX_COLS+1, B, C);
}

/* La función postfunción "reduce" la matriz, eliminando una fila y una
columna al principio y al final (la matriz tiene dimensiones (MAX_ROWS
+2)*(MAX_COLS+2) y el resultado MAX_ROWS*MAX_COLS. */
void postFuncion(int filas, int columnas, const celda *A, celda *C) {
  // veamos si entra a la funcion
  int i, j;
  printf("eliminando filas y columnas extras de la matriz...\n");
  celda *B;
  int n_columnas = columnas - 2;
  int n_filas = filas - 2;
  B = (celda *)malloc(n_filas * n_columnas * sizeof(celda));
  int f, c;
  c = 0;
  f = 0;
  // reducir la matriz
  for (i = 0; i < filas; i++) {
    if (!(i == 0 || i >= (filas - 1))) {
      c = 0;
      for (j = 0; j < columnas; j++) {
        if (!(j == 0 || j >= (columnas - 1))) {
          B[f * n_columnas + c].altitude = A[i * columnas + j].altitude;
          B[f * n_columnas + c].thickness = A[i * columnas + j].thickness;
          B[f * n_columnas + c].temperature = A[i * columnas + j].temperature;
          B[f * n_columnas + c].isVent = A[i * columnas + j].isVent;
          B[f * n_columnas + c].yield = A[i * columnas + j].yield;
          B[f * n_columnas + c].viscosity = A[i * columnas + j].viscosity;
          B[f * n_columnas + c].exits = A[i * columnas + j].exits;
          B[f * n_columnas + c].inboundV = A[i * columnas + j].inboundV;
          B[f * n_columnas + c].outboundV = A[i * columnas + j].outboundV;
          c += 1;
        }
      }
      f += 1;
    }
  }
  // copiar la matriz reducida a la pasada por parámetro.
  memcpy(C, B, n_filas * n_columnas * sizeof(celda));
}

// declaración de la variable global que guarda las condiciones iniciales
// no es muy recomendable, pero no hay tanto tiempo.
cIni c0;

// Esta es la función donde se calcula todo.
// en la versión CUDA es sustituida.
void FuncionPrincipal(int filas, int columnas, celda *A, celda *C) {
  // falta hacer una función visualizar que ignore las columnas extras
  // Asumimos al iniciar la función que la matriz ya viene aumentada.
  int i, j, l, m;
  int flujos = 0;
  double deltaV = 0.0, deltaT = 0.0, deltaH = 0.0, maxV = 0.0;
  double delta_t, deltaQ = 0.0, deltaQ_rad = 0.0, deltaQ_flu = 0.0;
  double Q_base = 0.0;
  delta_t = c0.deltat;
  double Aref = 0.0, Href = 0.0;
  double Acomp = 0.0, Hcomp = 0.0, Hcrit = 0.0, Hcrit2 = 0.0;
  double cArea = c0.anchoCelda * c0.anchoCelda;
  double alfa = 0.0;
  // primer ciclo que es solo para inicializar
  for (i = 1; i < filas - 1; i++) {
    for (j = 1; j < columnas - 1; j++) {
      // primero calcular la viscosidad y el yield, e inicializar los flujos
      A[i * columnas + j].viscosity = visc(A[i * columnas + j].temperature);
      A[i * columnas + j].yield = yield(A[i * columnas + j].temperature);
      A[i * columnas + j].inboundV = 0.0;
      A[i * columnas + j].inboundQ = 0.0;
      A[i * columnas + j].outboundV = 0.0;
      A[i * columnas + j].exits = 0;
    }
  }
  // ciclo para evaluar la cantidad de salidas que tiene cada celda
  for (i = 1; i < filas - 1; i++) {
    for (j = 1; j < columnas - 1; j++) {
      Href = A[i * (columnas) + j].thickness;
      Aref = A[i * (columnas) + j].altitude;
      flujos = 0;
      for (l = -1; l < 2; l++) {
        for (m = -1; m < 2; m++) {
          if (!(m == 0 && l == 0)) {
            int ni = i + l;
            int nj = j + m;
            Acomp = A[(ni) * (columnas) + (nj)].altitude;
            Hcomp = A[(ni) * (columnas) + (nj)].thickness;
            // alfa = atan((Acomp-Aref)/c0.anchoCelda);
            Hcrit =
                fabs((A[(ni) * (columnas) + (nj)].yield *
                      sqrt((Acomp - Aref) * (Acomp - Aref) +
                           c0.anchoCelda * c0.anchoCelda)) /
                     (density * gravity * ((Acomp - Aref) - (Hcomp - Href))));
            if ((Hcomp > Hcrit) && (Hcrit > 1e-8)) {
              // mas depuración
              // como esta operacion se repite es mejor hacerla una sola vez
              // calcular el valor del volumen que sale
              if (fabs(Hcomp + Acomp) >
                  fabs(Href + Aref)) { // aca ya asumi que es plano
                A[ni * (columnas) + nj].exits += 1;
              }
            }
          }
        }
      }
    }
  }

  // acá solo se calculan los flujos, en el primer ciclo.  Como uno de los
  // los supuestos de los autómatas celulares es que los estados solo se
  // actualizan al final, se necesita un segundo ciclo.  Esto es crucial para el
  // mapeo.
  for (i = 1; i < filas - 1; i++) {
    for (j = 1; j < columnas - 1; j++) {
      double deltaQ_flu_vent = 0.0;
      if (A[i * columnas + j].isVent == 1) {
        // variables temporales para almacenar el delta de volumen y de
        // temperatura. si hay un crater se aumenta el thickness en un valor
        // igual a la tasa de erupción sobre el area de la celda por el delta de
        // tiempo
        deltaV = (c0.eRate) * delta_t;
        A[i * columnas + j].inboundV += (deltaV);
        deltaQ_flu_vent = (deltaV * c0.eTemp) * heatCapacity * density;
        A[i * columnas + j].inboundQ += deltaQ_flu_vent;
        // printf("\nCrater %d,%d inbound %lf ", i-1, j-1, deltaV);
      }
      deltaV = 0.0;
      deltaH = 0.0;
      // depuracion, recordar que las celda se cuentan descartando las filas y
      // columnas extras printf("\n\n\nDatos de la celda %d,%d : Altitud: %lf -
      // Grosor: %lf - Temperatura: %lf - Viscosidad: %lf - Esfuerzo: %lf\n",
      // i-1, j-1, A[i*columnas+j].altitude, A[i*columnas+j].thickness,
      // A[i*columnas+j].temperature, A[i*columnas+j].viscosity,
      // A[i*columnas+j].yield); acá arranco el proceso de detectar cuantos
      // flujos salen de una celda y de esta manera puedo saber como dividir el
      // flujo cuando necesite sumarlos en el paso siguiente. reviso las celdas
      // adyacentes:
      Href = A[i * (columnas) + j].thickness;
      Aref = A[i * (columnas) + j].altitude;
      // printf("\n\nEl valor de flujos que sale de %d,%d es %d
      // \n",i-1,j-1,A[i*(columnas) + j].exits); crear el ciclo que recorre los
      // adyacentes y calcula el grosor critico en los flujos, agregar el deltaV
      // correcto. y el correspondiente deltaQ comparar a ver cuantas salidas
      // tiene cada celda. aca pongo salidas de depuración printf("Analizando
      // flujos de lava relacionados con la celda %d,%d que tiene valores de
      // altitud %lf y grosor de capa %lf:\n", i-1, j-1, Aref, Href); aca
      // calculo los flujos como tal, vamos a calcular directamente, sin
      // provisión de orden en los flujos, agregar el deltaV correcto. y el
      // correspondiente deltaQ
      for (l = -1; l < 2; l++) {
        for (m = -1; m < 2; m++) {
          if (!(m == 0 && l == 0)) {
            int ni = i + l;
            int nj = j + m;
            Acomp = A[(ni) * (columnas) + (nj)].altitude;
            Hcomp = A[(ni) * (columnas) + (nj)].thickness;
            deltaH = Hcomp - Href;
            // alfa = atan((Acomp-Aref)/c0.anchoCelda);
            alfa = 0;
            // printf("\nComparando celda %d,%d con %d,%d (ref)  Angulo %lf,
            // Acomp: %lf - Aref: %lf - Cos alfa %lf, Sin alfa %lf -- yield %lf
            // - visc %lf - gravedad %lf densidad %lf - ancho %lf", ni-1, nj-1,
            // i-1, j-1, alfa/M_PI, Acomp, Aref, cos(alfa), sin(alfa),
            // A[(ni)*(columnas)+(nj)].yield, A[(ni)*(columnas)+(nj)].viscosity,
            // gravity, density, c0.anchoCelda);
            Hcrit = fabs((A[(ni) * (columnas) + (nj)].yield *
                          sqrt((Acomp - Aref) * (Acomp - Aref) +
                               c0.anchoCelda * c0.anchoCelda)) /
                         (density * gravity * ((Acomp - Aref) - (deltaH))));
            // Hcrit =
            // ((A[(ni)*(columnas)+(nj)].yield)/((density*gravity)*(sin(alfa)-((Hcomp-Href)/c0.anchoCelda)*cos(alfa))));
            // esta version funcionaría perfecto para topografías planas? o por
            // el angulo ya esta incluido? hallar el grosor critico para cada
            // celda adyacente. Hcrit=(Hcomp-Href)/c0.anchoCelda; salidas de
            // depuración que hacen falta
            // printf("Comparando con celda %d,%d que tiene valores de altitud
            // %lf y grosor de capa %lf: -- grosor critico %lf \n", ni, nj,
            // Acomp, Hcomp, Hcrit);
            // printf("\nEl Hcomp es %lf, el Href es %lf, el ancho es %lf, El Hc
            // es %lf, salidas de la celda de comp %d, celda %d,%d", Hcomp,
            // Href, c0.anchoCelda, Hcrit, A[ni*(columnas)+nj].exits, ni-1,
            // nj-1);
            if ((Hcomp > Hcrit) && (Hcrit > 1e-8)) {
              // mas depuración
              // como esta operacion se repite es mejor hacerla una sola vez
              double h_hc = Hcomp / Hcrit;
              // calcular el valor del volumen que sale
              if ((fabs(Hcomp + Acomp) > fabs(Href + Aref)) &&
                  (A[ni * (columnas) + nj].exits > 0)) {

                // aca ya asumi que es plano
                // deltaV =
                // (1.0/A[ni*(columnas)+nj].exits)*((A[ni*columnas+nj].yield*Hcrit*Hcrit*c0.anchoCelda)/(3*A[ni*columnas+nj].viscosity))*(h_hc*h_hc*h_hc-1.5*h_hc*h_hc+0.5)*(delta_t);
                // estoy jugando
                deltaV = (1.0 / A[ni * (columnas) + nj].exits) *
                         ((A[ni * columnas + nj].yield * Hcrit * Hcrit *
                           c0.anchoCelda) /
                          (3 * A[ni * columnas + nj].viscosity)) *
                         (h_hc * h_hc * h_hc - 1.5 * h_hc * h_hc + 0.5) *
                         (delta_t);
                maxV = (deltaH * cArea) / (2 * A[ni * (columnas) + nj].exits);
                if (maxV < deltaV) {
                  // luz, fuego, destrucción
                  deltaV = maxV;
                }
              } else {
                deltaV = 0.0;
              }
              // printf("\nEl valor de dv que sale de %d,%d es
              // %lf",ni-1,nj-1,deltaV); El volumen cedido se suma al volumen
              // total a restar a la celda Y tambien al volumen a sumar.
              // machetazos horribles
              A[(ni) * (columnas) + (nj)].outboundV += (deltaV);
              A[(i) * (columnas) + (j)].inboundV += (deltaV);
              A[(i) * (columnas) + (j)].inboundQ +=
                  deltaV * A[ni * columnas + nj].temperature * density *
                  heatCapacity;
              // se incluye un estimado del calor transferido, temp del origen
              // por delta
              // A[(i)*(columnas)+(j)].inboundQ += (deltaV *
              // A[ni*columnas+nj].temperature) * heatCapacity * density;
              // A[(ni)*(columnas)+(nj)].inboundQ -= (deltaV *
              // A[ni*columnas+nj].temperature) * heatCapacity * density;
              // reviso el número de flujos
              // flujos += 1;
              // calcular el inboundQ por flujos
            }
          }
        }
      }
      // aca se agregan los flujos en caso de ser un crater
      // acá sumamos los volumenes y convertimos en grosor
    }
  }
  // segundo ciclo, consolidamos los flujos, calculamos nuevos grosores
  // y temperaturas.
  // luego agregamos crateres y calculamos la temperatura perdida por radiacion

  for (i = 1; i < filas - 1; i++) {
    for (j = 1; j < columnas - 1; j++) {
      deltaT = 0.0;
      deltaQ = 0.0;
      deltaQ_rad = 0.0;
      deltaQ_flu = 0.0;
      Q_base = 0.0;
      double deltaQ_flu_in = 0.0, deltaQ_flu_out = 0.0;
      double thickness_0 = 0.0, temperature_0 = 0.0;
      // Solo necesito calcular el valor de T teniendo en cuenta el calor
      // y el valor de thickness teniendo en cuenta el volumen
      // balance de volumenes, ojo.
      thickness_0 = A[i * columnas + j].thickness;
      temperature_0 = A[i * columnas + j].temperature;
      Q_base = thickness_0 * temperature_0 * density * heatCapacity;
      // Cuando el grosor el negligible con relación al area, no hay perdida de
      // calor if (A[i*columnas+j].thickness > 1e-8) {
      if (A[i * columnas + j].thickness > 1e-4) {
        // A[i*columnas+j].temperature = (A[i*columnas+j].inboundQ +
        // (thickness_0 * cArea *
        // A[i*columnas+j].temperature))/(A[i*columnas+j].thickness * cArea);
        deltaQ_rad =
            (-1.0) * SBConst * (cArea)*emisivity * delta_t *
            (A[i * columnas + j].temperature * A[i * columnas + j].temperature *
             A[i * columnas + j].temperature * A[i * columnas + j].temperature);
      } else {
        deltaQ_rad = 0;
      }
      A[i * columnas + j].thickness = thickness_0 +
                                      (A[i * columnas + j].inboundV / (cArea)) -
                                      (A[i * columnas + j].outboundV / (cArea));
      deltaQ_flu_in = A[i * columnas + j].inboundQ;
      deltaQ_flu_out = A[i * columnas + j].outboundV * temperature_0 * density *
                       heatCapacity;
      deltaQ_flu = deltaQ_flu_in - deltaQ_flu_out;
      // revisar por 0 en el thickness y no hacer la operacion en ese caso
      // calculo de la perdida, un nuevo deltaQ
      // printf("\nPara Celda %d,%d antes de crateres el valor de");
      // printf("\nPara Celda %d,%d antes de crateres el valor de grosor es %lf,
      // el de calor por flujo %lf, el de calor perdido por radiacion es %lf, el
      // de calor base es %lf, el de temperatura es %lf \n", i-1, j-1,
      // A[i*columnas + j].thickness, deltaQ_flu, deltaQ_rad, Q_base,
      // A[i*columnas + j].temperature); printf("\nPara Celda %d,%d antes de
      // crateres el valor de calor de salida es %lf, de calor de entrada es
      // %lf", i-1, j-1, deltaQ_flu_out, deltaQ_flu_in); mensajes de depuración.
      // Acá se cálcula si es un crater o no, y con eso se cálcula
      // un nuevo grosor.
      // se revisa si en la celda hay un crater
      deltaQ = Q_base + deltaQ_flu + deltaQ_rad;
      if (A[i * columnas + j].thickness > 1e-8) {
        A[i * columnas + j].temperature =
            deltaQ /
            (density * heatCapacity * cArea * A[i * columnas + j].thickness);
      } else {
        A[i * columnas + j].temperature = 273.0;
      }
      // printf("\nPara Celda %d,%d despues de crateres el valor de grosor es
      // %lf, el de calor por flujo %lf, el de calor perdido por radiacion es
      // %lf, el de calor base es %lf, el de temperatura es %lf, el calor
      // agregado por crater es bastante", i-1, j-1, A[i*columnas +
      // j].thickness, deltaQ_flu, deltaQ_rad, Q_base, A[i*columnas +
      // j].temperature); printf("\nEl la celda %d,%d el valor de grosor es %lf,
      // el de inbound es %lf y el de outboud es %lf.  El valor de calor es %lf
      // - la temperatura anterior es %lf, el grosor anterior es %lf", i-1, j-1,
      // A[i*columnas+j].thickness, A[i*columnas+j].inboundV,
      // A[i*columnas+j].outboundV, A[i*columnas+j].inboundQ, temperature_0,
      // thickness_0); printf("\nParametros cArea %lf, heatCapacity %lf, density
      // %lf temp %lf", cArea, heatCapacity, density,
      // A[i*columnas+j].temperature); printf("\nDeltaT %lf - Grosor: %lf -
      // deltaQ_rad %lf - inboundQ %lf\n\n", deltaT,
      // A[(i)*(columnas)+(j)].thickness, deltaQ_rad, A[i*columnas+j].inboundQ);
    }
  }
  memcpy(C, A, filas * columnas * sizeof(celda));
}

// nota, falta implementar las cifras significativas
int prepararVisualizacionGNUPlot(int secuencia, char *path, int filas,
                                 int columnas, celda *matriz, int cifrasSignif,
                                 double w, double x0, double y0) {
  // Esta función genera los dos archivos necesarios para producir una imagen
  // en GNU plot.  Los archivos son:
  // 1. datafile.dat (o variantes) que contiene los datos de altitud,
  // temperatura y grosor de la capa, teniendo en cuenta que altitud ya tiene en
  // cuenta el de la capa, pero este se incluye para dar mas opciones
  // 2. archivo de comandos, que incluye lo necesario para poder configurar
  // el área de dibujo, las escalas, las leyendas, etc.
  // Los parámetros son matriz, que es la matriz completa con los resultados
  // nombreArchivo, que es el nombre base de los archivos de datos a producir
  // y secuencia, en el caso en el que se produzcan varios archivos para
  // generar una animación, nos dará el orden.  El resultado es un plot de gnu
  // plot mas una imagen png exportada a partir de eso, que llevará como nombre
  // nombreArchivo + secuencia.png
  FILE *datosAltitud;
  FILE *encabezado;
  size_t sizep;
  int i, j, flag;
  double xcoord, ycoord, zcoord, temp;
  sizep = strlen(path);
  long int cont;
  char newpath[500], tsec[20], pathenca[550], command[700], f_path[1024];
  // printf("\nNumero de caracteres de la ruta del archivo %d\n",(int)sizep);
  if (sizep != 0) {
    strcpy(newpath, path);
    sprintf(tsec, "%d", secuencia);
    strcat(newpath, tsec);
    printf("Escribiendo archivo: %s\n", newpath);
    datosAltitud = fopen(newpath, "w");
    if (datosAltitud != NULL) {
      printf("Guardando en archivo...\n");
      // Escribiendo los datos
      cont = 0;
      for (j = 0; j < columnas; ++j) { // aca debo cambiar esto, altitude +
                                       // thickness da la altura, osea la z
        for (i = 0; i < filas; ++i) {
          xcoord = x0 + i * w;
          ycoord = y0 + j * w;
          zcoord = matriz[columnas * i + j].thickness +
                   matriz[columnas * i + j].altitude;
          temp = matriz[columnas * i + j].temperature;
          fprintf(datosAltitud, "%6.3lf %6.3lf %6.3lf %6.3lf\n", xcoord, ycoord,
                  zcoord, temp);
          cont += 1;
        }
        // para las isolineas
        fprintf(datosAltitud, "\n");
      }
      printf("Se escribieron %ld elementos.\n\n", cont);
      fclose(datosAltitud);
    } else {
      printf("***\nError al intentar escribir el archivo de datos.\n***\n\n");
    }
    strcpy(pathenca, newpath);
    strcat(pathenca, "_enca");
    printf("Guardando archivo de encabezado %s \n", pathenca);
    encabezado = fopen(pathenca, "w");
    if (encabezado != NULL) {
      // escribiendo el encabezado, pero en GNUPLOT es diferente
      // se muestran opciones para png y eps
      // fprintf(encabezado,"set terminal png size 400,300 enhanced font
      // \"Helvetica,20\"\n\n",newpath);
      fprintf(encabezado, "set terminal png size 1366,720 enhanced\n\n",
              newpath);
      fprintf(encabezado, "set output '%s.png'\n", newpath);
      // dejo estas líneas comentadas por si necesito configurar mejor
      // la visualización.
      // fprintf(encabezado,"set term postscript eps enhanced color\n");
      // fprintf(encabezado,"set output '%s.eps'\n",newpath);
      // fprintf(encabezado,"set autoscale y #set	autoscale	# scale
      // axes automatically\n"); fprintf(encabezado,"unset label	# remove
      // any previous labels\n"); fprintf(encabezado,"set encoding utf8\n\n");
      // fprintf(encabezado,"set isosamples 100\n");
      // fprintf(encabezado,"set samples 100\n\n");
      // fprintf(encabezado,"unset key\n");
      // fprintf(encabezado,"unset grid\n");
      // fprintf(encabezado,"unset border\n");
      // fprintf(encabezado,"unset tics\n");
      fprintf(encabezado, "set xrange [0:20]\n");
      fprintf(encabezado, "set yrange [0:20]\n");
      fprintf(encabezado, "set zrange [0:20]\n\n");
      fprintf(encabezado, "set colorbox\n");
      // fprintf(encabezado,"unset surface\n");
      fprintf(encabezado, "set pm3d\n");
      fprintf(encabezado, "set title \"%d\"\n\n", secuencia);
      fprintf(encabezado, "splot \"%s\" using 1:2:3:4 with pm3d\n", newpath);
      fprintf(encabezado, "reset\n");
      fclose(encabezado);
    } else {
      printf(
          "\n***\nError al intentar escribir el archivo de ecabezado.\n***\n");
    }
    // se ejecuta el gnu plot para generar el archivo de imagen
    limpiarPath(pathenca, f_path);
    strcpy(command, "gnuplot ");
    strcat(command, f_path);
    printf("Ejecutando comando: %s\n", command);
    flag = system(command);
    printf("Mensaje del comando %d\n", flag);
  }
}

// Esta es una copia de la funcion anterior que ignora los espacios extras
// de la matriz aumentada
int prepararVisualizacionGNUPlot_2(int secuencia, char *path, int filas,
                                   int columnas, celda *matriz,
                                   int cifrasSignif, double w, double x0,
                                   double y0) {
  // Esta función genera los dos archivos necesarios para producir una imagen
  // en GNU plot.  Los archivos son:
  // 1. datafile.dat (o variantes) que contiene los datos de altitud,
  // temperatura y grosor de la capa, teniendo en cuenta que altitud ya tiene en
  // cuenta el de la capa, pero este se incluye para dar mas opciones
  // 2. archivo de comandos, que incluye lo necesario para poder configurar
  // el área de dibujo, las escalas, las leyendas, etc.
  // Los parámetros son matriz, que es la matriz completa con los resultados
  // nombreArchivo, que es el nombre base de los archivos de datos a producir
  // y secuencia, en el caso en el que se produzcan varios archivos para
  // generar una animación, nos dará el orden.  El resultado es un plot de gnu
  // plot mas una imagen png exportada a partir de eso, que llevará como nombre
  // nombreArchivo + secuencia.png
  FILE *datosAltitud;
  FILE *encabezado;
  size_t sizep;
  int i, j, flag;
  double xcoord, ycoord, zcoord, temp;
  sizep = strlen(path);
  long int cont;
  char newpath[500], tsec[20], pathenca[550], command[700], f_path[1024];
  // printf("\nNumero de caracteres de la ruta del archivo %d\n",(int)sizep);
  if (sizep != 0) {
    strcpy(newpath, path);
    sprintf(tsec, "%d", secuencia);
    strcat(newpath, tsec);
    printf("Escribiendo archivo: %s\n", newpath);
    datosAltitud = fopen(newpath, "w");
    if (datosAltitud != NULL) {
      printf("Guardando en archivo...\n");
      // Escribiendo los datos
      cont = 0;
      for (j = 1; j < columnas - 1; ++j) { // aca debo cambiar esto, altitude +
                                           // thickness da la altura, osea la z
        for (i = 1; i < filas - 1; ++i) {
          xcoord = x0 + (i - 1) * w;
          ycoord = y0 + (j - 1) * w;
          zcoord = matriz[columnas * i + j].thickness +
                   matriz[columnas * i + j].altitude;
          temp = matriz[columnas * i + j].temperature;
          fprintf(datosAltitud, "%6.3lf %6.3lf %6.8lf %6.8lf\n", xcoord, ycoord,
                  zcoord, temp);
          cont += 1;
        }
        // para las isolineas
        fprintf(datosAltitud, "\n");
      }
      printf("Se escribieron %ld elementos.\n\n", cont);
      fclose(datosAltitud);
    } else {
      printf("***\nError al intentar escribir el archivo de datos.\n***\n\n");
    }
    strcpy(pathenca, newpath);
    strcat(pathenca, "_enca");
    printf("Guardando archivo de encabezado %s \n", pathenca);
    encabezado = fopen(pathenca, "w");
    if (encabezado != NULL) {
      // escribiendo el encabezado, pero en GNUPLOT es diferente
      // se muestran opciones para png y eps
      // fprintf(encabezado,"set terminal png size 400,300 enhanced font
      // \"Helvetica,20\"\n\n",newpath);
      fprintf(encabezado, "set terminal png size 1366,720 enhanced\n\n",
              newpath);
      fprintf(encabezado, "set output '%s.png'\n", newpath);
      // dejo estas líneas comentadas por si necesito configurar mejor
      // la visualización.
      // fprintf(encabezado,"set term postscript eps enhanced color\n");
      // fprintf(encabezado,"set output '%s.eps'\n",newpath);
      // fprintf(encabezado,"set autoscale y #set	autoscale	# scale
      // axes automatically\n"); fprintf(encabezado,"unset label	# remove
      // any previous labels\n"); fprintf(encabezado,"set encoding utf8\n\n");
      // fprintf(encabezado,"set isosamples 100\n");
      // fprintf(encabezado,"set samples 100\n\n");
      // fprintf(encabezado,"unset key\n");
      fprintf(encabezado, "set grid\n");
      // fprintf(encabezado,"unset border\n");
      // fprintf(encabezado,"unset tics\n");
      fprintf(encabezado, "set xrange [0:%d]\n", filas - 2);
      fprintf(encabezado, "set yrange [0:%d]\n", filas - 2);
      fprintf(encabezado, "set zrange [0:3]\n\n");
      fprintf(encabezado, "set colorbox\n");
      // fprintf(encabezado,"unset surface\n");
      fprintf(encabezado, "set pm3d\n");
      fprintf(encabezado, "set title \"%d\"\n\n", secuencia);
      fprintf(encabezado, "splot \"%s\" using 1:2:3:4 with pm3d\n", newpath);
      fprintf(encabezado, "reset\n");
      fclose(encabezado);
    } else {
      printf(
          "\n***\nError al intentar escribir el archivo de ecabezado.\n***\n");
    }
    // se ejecuta el gnu plot para generar el archivo de imagen
    limpiarPath(pathenca, f_path);
    strcpy(command, "gnuplot ");
    strcat(command, f_path);
    printf("Ejecutando comando: %s\n", command);
    flag = system(command);
    printf("Mensaje del comando %d\n", flag);
  }
}

/*
 * Función para obtener el path actual, y de esta manera generar
 * las imagenes en el directorio actual.
 */
int obtenerPath(char path[]) {
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
    fprintf(stdout, "Current working dir: %s\n", cwd);
    strcpy(path, cwd);
    return 0;
  } else {
    perror("getcwd() error");
    return 1;
  }
}

// función para colocar \ donde sean necesarios y de esta manera
// poder usar paths y nombres de archivos con espacios.
int limpiarPath(char path[], char spath[]) {
  char r_path[1024], f_path[1024];
  int i, j;
  int lg = 0;
  strcpy(r_path, path);
  lg = strlen(r_path);
  j = 0;
  for (i = 0; i < lg; ++i) {
    if (r_path[i] == ' ') {
      f_path[j] = '\\';
      f_path[j + 1] = ' ';
      j += 2;
    } else {
      f_path[j] = r_path[i];
      j += 1;
    }
  }
  f_path[j] = '\0';
  strcpy(spath, f_path);
}

int generarAnimacionGNUPlot(char nombreArchivo[], int secuencia) {
  /* ----------------------------------------------------------
   * En construcción.
   * "mencoder mf://*.%s -mf w=1366:h=720:fps=5:type=png -ovc copy -oac copy -o
   * %S", "png", nombre.avi Función para generar un video a partir de imágenes
   * fijas.
   * ---------------------------------------------------------- */
}

// Acá va la función main.
int main(int argc, char *argv[]) {
  int i, flag = 0;
  celda *testPoint, *resultPoint, *resultCalc, *resultPoint2;
  point2D *crateres;
  int fila = 0;
  int columna = 0;
  int puntosCrater = 0;
  char s_path[1024];
  char a_path[1024];
  char etiqueta[1024];
  char path[1024];
  // acá revisamos los parámetros
  /*if(argc<21) {
          printf("Numero de parámetros incorrecto\n");
          printf("./scalaf -t temperatura_erupcion -v velocidad_erupcion -w
  ancho_celda -s archivo_salidas -a archivo_altitudes -r filas -c columnas -p
  puntos_crater -\n\n"); exit(1);
  }*/
  // ojo aca se asume que los parametros se pasan correctamente,
  // o sea se revisa la existencia de estos, no sou congruencia.
  for (i = 0; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 't': {
        // temperatura de extruccion
        c0.eTemp = atof(argv[++i]);
        // c0.eTemp = 273.00+atof(argv[++i]);
        break;
      }
      case 'v': {
        // veloicidad de extruccion
        c0.eRate = atof(argv[++i]);
        break;
      }
      case 'w': {
        // ancho de las celdas, recordar que son cuadradas
        c0.anchoCelda = atof(argv[++i]);
        break;
      }
      case 's': {
        // archivo de los cráteres
        strcpy(s_path, argv[++i]);
        break;
      }
      case 'a': {
        // archivo de las altitudes
        strcpy(a_path, argv[++i]);
        break;
      }
      case 'r': {
        // numero de filas
        c0.tFilas = atol(argv[++i]);
        break;
      }
      case 'c': {
        // numero de columnas
        c0.tColumnas = atol(argv[++i]);
        break;
      }
      case 'p': {
        // numero de crateres
        puntosCrater = atol(argv[++i]);
        break;
      }
      case 'e': {
        // etiqueta de salida
        strcpy(etiqueta, argv[++i]);
        break;
      }
      case 'n': {
        // etiqueta de salida
        c0.nPasos = atol(argv[++i]);
        break;
      }
      }
    }
  }
  // dt esta de momento hardcoded, lo cambiaré más adelante
  c0.deltat = dt;
  // prueba de los parámetros ingresados
  // printf("\n\nparametros leidos filas=%d columnas=%d ancho=%lf velocidad=%lf
  // temperatura=%lf crateres=%d \n", c0.tFilas, c0.tColumnas, c0.anchoCelda,
  // c0.eRate, c0.eTemp, puntosCrater); leer datos de topografía ahora se tienen
  // las filas y columnas, se pueden crear el puntero que representa la matriz
  // tambien se pueden crear los punteros que representan las variantes
  // agrandadas y reducidas de la matriz
  // falta: liberar la memoria en cada paso.
  testPoint = (celda *)malloc(c0.tFilas * c0.tColumnas * sizeof(celda));
  resultPoint =
      (celda *)malloc((c0.tFilas + 2) * (c0.tColumnas + 2) * sizeof(celda));
  resultCalc =
      (celda *)malloc((c0.tFilas + 2) * (c0.tColumnas + 2) * sizeof(celda));
  resultPoint2 = (celda *)malloc(c0.tFilas * c0.tColumnas * sizeof(celda));
  leerArchivoTexto_Matriz(a_path, c0.tFilas, c0.tColumnas, testPoint);
  // al tener el número de crateres se puede crear el puntero que representa el
  // conjunto de estos.
  crateres = (point2D *)malloc(puntosCrater * sizeof(point2D));
  readCratersPositionFile(s_path, puntosCrater, crateres);
  // codigo de prueba:
  /*printf("Test de los puntos de salida de lava\n\n");
  for (i = 0; i < puntosCrater; ++i) {
          printf("Punto %d, x=%d y=%d\n", i, crateres[i].x, crateres[i].y );
  }*/
  placeCraters(testPoint, crateres, c0.tFilas, c0.tColumnas, puntosCrater);
  // printf("\n\n\n Primera matriz \n");
  // imprimirMatrizPantalla(c0.tFilas , c0.tColumnas, testPoint, 3);
  /*printf("\n test de la viscosidad\n\n");
  for (i=0; i< 2000; i+=100) {
          printf("%lf %lf\n",i+500.0, visc(i+500.0));
  }
  printf("\n test del yield \n\n");
  for (i=0; i< 2000; i+=100) {
          printf("%lf %lf\n",i+500.0, yield(i+500.0));
  }*/
  // fin codigo de prueba;
  // place-holders de las funciones del flujo de agrandar reducir
  preFuncion(c0.tFilas, c0.tColumnas, testPoint, resultPoint);
  for (i = 0; i < c0.nPasos; i++) {
    printf("\n\nPaso de Tiempo %d: \n\n", i);
    FuncionPrincipal(c0.tFilas + 2, c0.tColumnas + 2, resultPoint, resultCalc);
    // imprimirMatrizPantalla_2(c0.tFilas+2, c0.tColumnas+2, resultCalc, 3);
    flag = obtenerPath(path);
    strcat(path, "/");
    strcat(path, etiqueta);
    strcat(path, "_");
    // poner el path
    if (!(flag)) {
      if (i % 5 == 0) {
        prepararVisualizacionGNUPlot_2(i, path, c0.tFilas + 2, c0.tColumnas + 2,
                                       resultCalc, 3, c0.anchoCelda, 0, 0);
      }
    } else {
      printf("Problemas con el path\n");
    }
    memcpy(resultPoint, resultCalc,
           (c0.tFilas + 2) * (c0.tColumnas + 2) * sizeof(celda));
  }
  postFuncion(c0.tFilas + 2, c0.tColumnas + 2, resultCalc, resultPoint2);
  // place-holder de la funcion de escribir gnuplot
  // flag = prepararVisualizacionGNUPlot(1, "/home/sergio/salida1", MAX_ROWS,
  // MAX_COLS, test, 3, 1, 0, 0); printf("%d",flag); fin del programa
  return 0;
}
