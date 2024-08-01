# Color_RNA

##Beta Version##

This package uses web scraping to generate custom colored images of RNA secondary structures by using forna[1].
Currently supports 5 colors, Red, Blue, Green, Yellow and Pink. 
Further color support to be added in the future.

How to use:
Suitable chromedriver.exe should be downloaded from https://chromedriver.chromium.org/downloads
We recommend adding chromedriver.exe to the working path but you can place it elsewhere and supply its path (see instructions below).
Internet connection and Chrome browser are required

The function that you are going to use is called "create_image".
It's inputs are:
##
sequence (str): The sequence to be drawn
Accepts characters:
"A", "a", "T", "t", "G", "g", "C", "c", "U", "u"
According to nitrogenous bases.

##
structure (str): The secondary structure to be drawn
Accepts characters:
"." - Unbound base
"(" - Bound base - must be accompanied by ")" later in the sequence
")" - Bound base - first ")" is bound to first "(", etc. 
"[", "]" - Pseudoknots connection, works as "(" and ")"

##
colors_string (str): Position of wanted colors in the image.
For images with no colors, supply an empty string. 
Accepts charcters:
"r" - Red
"g" - Green
"b" - Blue
"p" - Pink
"y" - yellow
" " - White


## 
Optional: image_path (str): The path and image name in which the user wants the image to be saved.
Default: "RNA_image.png"
If not supplied, the image will be saved in the working path as "RNA_image.png".

##
Optional: pseudoknots_strength (str): the strength of the pseudoknots.
If not supplied, the default is weak pseudoknots.
Else, strong pseudoknots will be enabled.

## 
Optional: driver_path (str): The path in which "chromedriver.exe" is present.
If not supplied, the function expects chromedriver.exe to be in working path.
If the user wants, they can save chromedriver.exe elsewhere and supply the path.

Example of use:
After installation and chromedriver download, import package and use function as follows:
from color_rna import Color_RNA

In this example, we wish to draw:
*sequence - "ATGCCGTA"
*structue - "(......)"
*colors - positions 1,2 are red, positions 3-6 are yellow and position 6,7 are green
*image to be saved at working path with "Example.png" as its name.
*pseudoknots bonds are weak (default)
*chromedriver.exe is saved on the desktop and not in working path.
The code should look like:
Color_RNA.create_image("ATGCCGTA", "(......)", "rryyyygg", image_path = "Exapmle.png", driver_path = r"C:\Users\Username\Desktop\chromedriver")

[1] Kerpedjiev P, Hammer S, Hofacker IL (2015). Forna (force-directed RNA): Simple and effective online RNA secondary structure diagrams. Bioinformatics 31(20):3377-9.