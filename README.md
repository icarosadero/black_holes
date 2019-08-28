# black_holes
Code to solve Kerr and Schwarzschild Black Hole geodesics in the equator of the metrics. This code can be used to find the orbits of a time-like particle near the black holes.
For 'round' solutions in the Kerr Black Hole, use the Python script `main.py` to optimize for the smallest variance in radius.
The metrics ae solved standard ODE numerical methods. For the case of the Schwarzschild Black Hole, the method used is Runge-Kutta and for the Kerr Black Hole it is employed the Dormand Prince Integrator 853.
