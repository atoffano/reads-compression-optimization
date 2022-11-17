### Projet ALG 2022 -- Aide au codage

> V1 (26/10/2022)


Vous pouvez utiliser ces suggestions pour faciliter la mesure du temps d'execution d'une partie de votre code et l'affichage d'une barre de progression lors de différents calculs.

#### Mesure du temps d'execution d'une portion de code: 

timer.py :

```python
import time


# This class allows us to monitor the time spent by a function
class Timer():
    # Function called right before the execution of the function
    def __enter__(self):
        self.t1 = time.perf_counter()
        return self

    # Function called right after the function
    def __exit__(self, type, value, traceback):
        self.t2 = time.perf_counter()
        self.t = self.t2 - self.t1

    # Function that prints on the shell the time spent by the instructions
    def print(self, template: str = "{}"):
        print(template.format(round(self.t, 2)))

if __name__ == "__main__":
    # Exemple how to use this class
    with Timer() as total_time:  # time all instructions in the ’with’ statements
        for i in range(5):
            time.sleep(0.471)
    total_time.print("Durée = {} secondes")
```

#### Affichage d'une barre de progression: 

progress_bar.py :

```python
import sys
def update_progress(progress):
    barLength = 50  # Modify this to change the length of the progress bar
    status = ""
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100,2), status)
    sys.stderr.write(text)
    sys.stderr.flush()
    
if __name__ == "__main__":
    n = 10000000
    for i in range(n):
        if i%1000 == 0: # Avoid to update at each step (time consuming)
            update_progress(i/float(n))
            # your real code here
    update_progress(1) # Ends the line & writes "Done"
```

**Exemple (avec importation et utilisation de barre de progression et de calcul de temps)**

```python
import time
from timer import Timer	# should be in the same directory 
from progress_bar import update_progress # should be in the same directory 

if __name__ == "__main__":
    n = 60
    with Timer() as total_time:  # time all instructions in the ’with’ statements
        for i in range(n):
            update_progress(i/float(60))
            time.sleep(0.1) # your code here
        update_progress(1)
    total_time.print("Time spent = {} seconds")
```



