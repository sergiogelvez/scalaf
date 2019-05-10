#!/bin/bash

# Simulation Parameters
eruption_temperature=1500
eruption_rate=100
cell_width=1
map_rows=1024
map_columns=1024
number_of_craters=1
time_steps=10

# Config Parameters
ARCH_EXEC=build/scalaf
ARCH_OUTPUT=output/
ARCH_ALTITUD=altitudesMil.csv
ARCH_CRATER=crater_location
BASENAME=${map_rows}x${map_columns}_EXP
number_of_repetitions=1

HILOS_BLOQUE=171
BLOQUES=$(($map_rows/$HILOS_BLOQUE))

clear
printf "Número de hilos $HILOS_BLOQUE, número de bloques $BLOQUES"
printf "\nIniciando iteraciones"
for i in `seq 0 $((${number_of_repetitions} -1))`;
	do
	ITERATION_NAME=${BASENAME}_${i}
	ARCH_SALIDA=${ITERATION_NAME}
	ARCH_ERROR=${ITERATION_NAME}_ERR
	printf "\nIteración $i"
	printf "\nOUTPUT_FILE: $ARCH_SALIDA"
	printf "\nERROR_FILE: $ARCH_ERROR"
	
	printf "\ntime ./$ARCH_EXEC -t $eruption_temperature -v $eruption_rate -w $cell_width -s $ARCH_CRATER -a $ARCH_ALTITUD -r $map_rows -c $map_columns -p $number_of_craters -e $ITERATION_NAME -n $time_steps > $ARCH_SALIDA 2 > $ARCH_ERROR"
	time ./$ARCH_EXEC -t $eruption_temperature -v $eruption_rate -w $cell_width -s $ARCH_CRATER -a $ARCH_ALTITUD -r $map_rows -c $map_columns -p $number_of_craters -e $ITERATION_NAME -n $time_steps > $ARCH_SALIDA
	
	# Comprimir resultados y mover a directorio de salida
	tar -zcvpf ${ITERATION_NAME}.tgz ${ITERATION_NAME}_*
	mkdir -p ${ARCH_OUTPUT}
    mv ${ITERATION_NAME}* ${ARCH_OUTPUT}
	printf "\nArchivo copiado en ${ARCH_OUTPUT}"
	rm ${ITERATION_NAME}_*
	printf "\nArchivos limpiados"
	done

printf "\nFin de las iteraciones"
