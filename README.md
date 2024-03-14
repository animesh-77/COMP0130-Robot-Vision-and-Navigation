Scripts and report for CW1 for COMP0130 Robot-Vision-and-Navigation taught at UCL 2023-24.

Here are some feedback that we recieved regarding the submission.

* Position solution is generally good, but affected by GNSS outliers.
Your north velocity solution is good, but your east velocity is completely wrong; it seems to be the same as your north velocity.
Heading is noisy, but unbiased
* Generally good. Cascading a GNSS KF into the integration KF can cause stability issues if R is not big enough in the integration filter
* No GNSS outlier detection in your code
* In  your heading filter, “predictedRateOfChange” needs to be multiplied by the time interval.