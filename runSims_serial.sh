#!/bin/bash

# Acá asigno las variables, para que sea más facil editar los valores
NUM_SIMS=10
ARCH_EXEC=scalaf_serial_imp
ARCH_ALTITUD=altitudesMil.csv
ARCH_CRATER=aCraterMil
TEMP_ERUP=1500
TASA_ERUP=0.15
WIDTH=1
FILAS=1024
COLS=1024
RAIZ=serial_imp_3600_fullmesh
NUM_CRATER=1
PASOS_TIEMPO=3600
HILOS_BLOQUE=171
BLOQUES=$(($FILAS/$HILOS_BLOQUE))

clear
echo "numero de hilos $HILOS_BLOQUE, numero de bloques $BLOQUES"
echo "Iniciando iteraciones"
for i in `seq 0 $((${NUM_SIMS}-1))`;
	do
	# Nombre de salida, con el nombre correcto
	NOMBRE=${RAIZ}_EXP_${i}
	ARCH_SALIDA=${NOMBRE}
	echo ""
	echo "Iteración $i"
	echo "$ARCH_SALIDA es el nombre de la salida"
	ARCH_ERROR=${NOMBRE}_ERR
	echo "$ARCH_ERROR es el nombre de la salida"
	echo "time ./$ARCH_EXEC -t $TEMP_ERUP -v $TASA_ERUP -w $WIDTH -s $ARCH_CRATER -a $ARCH_ALTITUD -r $FILAS -c $COLS -p $NUM_CRATER -e $NOMBRE -n $PASOS_TIEMPO > $ARCH_SALIDA 2 > $ARCH_ERROR"
	time ./$ARCH_EXEC -t $TEMP_ERUP -v $TASA_ERUP -w $WIDTH -s $ARCH_CRATER -a $ARCH_ALTITUD -r $FILAS -c $COLS -p $NUM_CRATER -e $NOMBRE -n $PASOS_TIEMPO > $ARCH_SALIDA
	# acá tambien tiene que ir el comando de comprimir el resultado
	tar -zcvpf ${NOMBRE}.tgz ${NOMBRE}_*
	mkdir /sergio_tmp/backups
    mv ${NOMBRE}.tgz /sergio_tmp/backups
	echo "Archivo copiado en backups"
	rm ${NOMBRE}_*
	echo "Archivos limpiados"
	# y luego copiar el archivo en otra carpeta
	#time ./$ARCH_EXEC -t $TEMP_ERUP -v $TASA_ERUP -w $WIDTH -s $ARCH_CRATER -a $ARCH_ALTITUD -r $FILAS -c $COLS -p $NUM_CRATER -e $NOMBRE -n $NUM_CRATER -b $HILOS_BLOQUE > $ARCH_SALIDA 2 > $ARCH_ERROR
	done
echo ""
echo "Fin de las iteraciones"
