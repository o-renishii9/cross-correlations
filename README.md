cross-correlations
==================

Time delays at SdS phases using cross correlations
This is a preliminary code, it needs to be optimised and some other functions need to be integrated as well.
It also performs a linear regression between the predicted travel times for each of the waveforms 
(using the epicentral distance and cubic spline interpolation to obtain prem travel time predictions) and 
the measured by cross correlation.
What is important to observe is the general scatter between these two sets of dt.
