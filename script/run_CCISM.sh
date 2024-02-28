#!/bin/bash

mkdir -p ../../zenodo_repository/CCISM
for lib in p0{07,08,13,14,16,20,21}t p0{26,35}t p009t{1,2}; do 
	CCISM -i ../../zenodo_repository/cellSNP/${lib} -o ../../zenodo_repository/CCISM/${lib} -m 0 --estimate_power --thetaN 1.e-4
	for frac in 05 1 2 5; do 
		lib2=${lib}_${frac}
		echo ${lib2}
		CCISM -i ../../zenodo_repository/cellSNP/downsampled/${lib2} -o ../../zenodo_repository/CCISM/${lib2} -m 0 --estimate_power --thetaN 1.e-4
	done
done
