# Map-generator-for-matlab
This program can be used to create large terrain maps for matlab, good for fantasy RPGs.


This program creates a continent or island for RPG games using a simple
form of the Plurigaussian simulation. 
The resulting png is composed by pixels, with different colors for the 
various terrain types. Each pixel represents an area that the appropriate
terrain is dominant. This program -doesn't- place cities or rivers.
There is also matrix "VAL" with values representing the terrain type for 
each pixel. The png (or matrix) can be used with different programs to 
trace or convert them to hex maps.
Personally, I use the image with the Hexographer program 
(http://www.hexographer.com) to "convert Underlay" and create hex maps.

The program creates a "base" map first as a trend and uses this to create
a series of more detailed maps that it then joins together using some
basic (and crude) smoothing techniques and corrections.

ATTENTION: This program is crudely made, and it suits me and my needs. It
requires time and RAM to run and it is no way optimized.
For example, if you want to change any of the parameters you need to get
in the code and change them as this program is not a function with
inputs, nargin etc.

DEPENDENCIES: This program requires the statistics toolbox and image 
processing toolbox of Matlab.
OCTAVE: This program can be used with Octave but it takes signficantly longer than the times presented here. 

This program uses the exponential covariance so that every point on the
grid has some correlation (even if very small) to all other points. The
exponential model is the simplest to implement but probably not the best.
Bibliography suggests that 2D random-walk correlation would be the best
fit but this program was not created specifically to create terrain maps. 

If you want to use sparse matrices to save on space in order to create
trully huge maps, then I suggest you use the spherical covariance model
and triple the correlation lenghts. 

EXAMPLE 1: You can run the program as is and get five maps of a big temperate
island, without tundras or deserts. 240 x 240 pixels, or 480 x 480 miles.
It takes my computer 3-4GB of RAM and 14 minutes to make the necessary 
covariance matrix and once the matrix is calculated then it takes 3-4 
mins to make the simulations.
Time varies depending on how many "Attempts" are required to get a good
simulation.

EXAMPLE 2: You can change SyntheM value to SyntheM=10, change North and
South to North=1 and South=1, change corl to corl=2 and have 5 large
continents, comparable to Europe in size with cold climate on top and warm
climate at the bottom. The resulting map will be 800x800 pixels or 2400 x
2400 miles.
It takes my computer 3-4GB of RAM and about 2 hours to finish the 5
continents.
Time varies depending on how many "Attempts" are required to get a good
simulation.

PARAMETERS:
Since the program was made for me, it's not a function with inputs, so to
change anything you have to change it directly from the code.
Explanations of the parameters are in comments.
