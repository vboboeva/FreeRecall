# Free Recall

Codes with simulations for the data published in

> V. Boboeva, A. Pezzotta, C. Clopath, *Free-recall scaling laws and short-term memory effects in a latching attractor network*

## Codes


### Main simulations

1. `DetailedLR`: Figs. 1,2,3 and 4 of the main text, S4-S6 of SI 

3. `SimplifiedLR`: Figs. 5 and 6A, B, C, F, G of main text

4. `TwoOrthSets`: Fig. 6D, E, F, G of main text, S8 of SI

5. `TwoListsFFR`: Fig. S7 of SI


### Mean field

`MeanField`: folder containing code for producing Figs. S1-S3 of SI


## Compile and run

```bash
cd DetailedLR  # or any of the directories 1-4
cd main
make
./main <trial number>  <pattern set number>
```
The code uses the files in the input folder to produce files in the output folder, that can be analyzed to produce the figures. 

**Note**: working with `gcc-7` compiler. You can test other compilers by changing `CC` in `Makefile`.