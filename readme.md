Dans ce répertoire, vous trouvez tout le matériel nécessaire pour le TP de l'ES problèmes inverse 2020/2021.

## I - Installation de Julia


Pour gagner du temps, il serait souhaitable que vous l'ayez installé sur votre machine. C'est assez simple et rapide. Vous trouverez toutes les instructions ici :

https://julialang.org/

## II - Jupyter

Ensuite, ce n'est pas obligatoire mais ce sera plus confortable, nous allons travailler avec des notebooks que vous pouvez éditer grâce à Jupyter notebook

que vous connaissez déjà.

Pour que Julia puisse être utilisé dans Jupyter, il y a une petite manip à faire (suivre les étapes  2 et 3 du lien suivant)

https://datatofish.com/add-julia-to-jupyter/

## III - Packages

Dernière étape (qui pourra être faite mardi mais il est préférable de le faire en avance), installer les packages utiles.

Il suffit de copier-coller les lignes suivantes dans une console julia ou bien directement dans votre notebook.

```
using Pkg  #chargement de l'utilitaire de package

Pkg.add("Plots")
Pkg.add("Images")
Pkg.add("TestImages")
Pkg.add("Printf")
Pkg.add("StatsBase")
Pkg.add("LaTeXStrings")
Pkg.add("Optim")
Pkg.add("Dierckx")
Pkg.add("DSP")
Pkg.add("Dierckx")
```

## IV - Télécharger les notebooks et les lancer (à faire en début de TP)
Télécharger le contenu du dépôt et lancer jupyter.
Pour ce dernier point on peut executer `julia` depuis un terminal et taper:
```
using IJulia
notebook()
```


