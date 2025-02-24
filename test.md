### Methodology: Poisson Spike Train Generation  

The **simple_spike_generator** function simulates neuronal firing using a Poisson process. It generates a binary spike train over a duration \( T \) with time steps of size \( dt \). The firing rate \( r \) can be a constant (homogeneous Poisson process) or a time-varying array (inhomogeneous Poisson process).  

The function initializes a binary array of length \( T/dt \), representing spike occurrences. At each time step, a spike is generated probabilistically based on \( r \cdot dt \). If a spike occurs, the corresponding index is set to one. An optional refractory period prevents immediate subsequent spikes, consisting of an absolute refractory period \( abs\_ref \) and a relative refractory period \( rel\_ref\_mean \), sampled from an exponential distribution.  

The output is a binary numpy array indicating spike occurrences. This method effectively models neuronal firing patterns, accommodating both constant and dynamic firing rates while incorporating biologically realistic refractory effects.