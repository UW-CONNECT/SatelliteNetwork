# Satellite_Tracking
 Satellite tracking from TLEs in Python

## track_singleTxRx.py 
txlatlong is the latitude and longitude coordinates of the ground station of interest. 
satnames corresponds to satellites that we are tracking. 
Running this script shows a timeplot of next passes (according to minElevationAngle)

## heatmap.m 
Mostly similar to the implementation of heatmap.m -- plots relative 
connectivity of rx points about a certain centerpoint. 

## track_singleGS.py 
Plots connectivity of a single groundstation for multiple satellites. 

## track_twoGS.py 
Plots simultaneous connectivity of two groundstations for multiple satellites. 
 
## TODO
* Accept multiple Tx/Rx locations in latitude/longitude representations 
* What is the number of Rx/gateways needed to have a certain duration of connectivity for a given Tx? 
	* Minimum number of Rx we need to cover an area? 
* Determine the maximum coverage of a given transmitter and set of Rx locations 

## Example Scripts
Example from the web to base following code off of. 
 
sat_tracker_tle.py

## Other notes 
How MATLAB access between satellites was determined: https://www.mathworks.com/help/aerotbx/ug/matlabshared.satellitescenario.access.html#mw_d2108208-f31e-4cbd-b8b5-63b0ff995dce [algorithms section]

To update TLEs automatically without webscraping: https://stackoverflow.com/questions/34459285/can-scraping-be-applied-to-this-page-which-is-actively-recalculating
(actually just use built-in functions from skyfield)