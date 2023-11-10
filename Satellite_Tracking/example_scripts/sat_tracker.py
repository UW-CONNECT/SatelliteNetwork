'''
This is example code obtained from https://medium.com/@meetgandhi586/python-in-space-building-a-satellite-tracker-with-python-9ff43b24f0bf

TLE from https://celestrak.org/NORAD/elements/stations.txt
'''
import ephem 
import matplotlib.pyplot as plt 
from datetime import datetime 


# Define a function to track a satellite's position
def track_satellite(satellite_name, observer_location, duration_hours):
    observer = ephem.Observer()
    observer.lon, observer.lat = observer_location # Set your observer's longitude and latitude
    t = datetime.utcnow()
    print(t)

    for _ in range(duration_hours):
        satellite = ephem.readtle(
        satellite_name,
        "1 25544U 98067A   23312.57012176  .00025167  00000+0  44286-3 0  9997",
        "2 25544  51.6447 333.3637 0001209  90.7991  16.8031 15.50299604424183",
        )

        satellite.compute(t)
        print(f"Satellite: {satellite_name} - Altitude: {satellite.alt}, Azimuth: {satellite.az}")
        #print(satellite.alt)
        t += ephem.hour

# Example usage
track_satellite("ISS (ZARYA)", ("1.5", "52.3"), 5) # Track the ISS for 5 hour