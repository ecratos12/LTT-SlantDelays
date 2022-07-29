# LTT-SlantDelays
This is a slant delay forward model based on ray tracing for determination of the GNSS signal path through the Numerical Weather Prediction model grid.

## LTT - "Least Travel Time"
Concept of Least-Travel time operator is based on search for the signal path that provides the shortest possible travel time from the transmitter to the receiver. 

The propagation is stricted to two-dimensional plane. According to this assumption, Snell's Law equations were implemented and solved within the operator.

See  Clive D. Rodgers' book "Inverse methods for atmospheric sounding", section 9.4, for verbose description of the approach.

# The input
The input to the LTT algorithm consists of a two-dimensional array of refractivity as a function of radius and angular separation (polar coordinates counterpart), as interpolated to the two-dimensional plane from the three-dimensional weather data.

## Weather data
Refractivity solely depends on the atmospheric state, as a function of temperature, specific humidity, pressure, and non-gaseous constituents. The atmospheric state is provided by Numerical Weather Prediction model solutions.

Our slant delay model uses forecasted states from [OpenIFS model](https://confluence.ecmwf.int/display/OIFS/OpenIFS+Home)
