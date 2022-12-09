# reads-compression-optimization
Optimizing read ordering for a better compression efficiency.

## Goal
New gen sequencing requires a large disk space in order to store the read output. In this project we try to reorder the reads contained in an output file to have similar reads close by.
This repository contains a few functions that utilize different strategies in order to match reads.
Those functions are based on two strategies : either a k-mer based sorting strategy or a dimension reduction method
Stratégies
- Fonction baseline de sorting (à trouver) pour comparer le reste
- Fonctions à tester
