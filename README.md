# Approximate Bayesian computation (ABC) approximation of the potassium channel current parameters of the classical Hodgkin-Huxley membrane potential model.

In 1952, the first biophysical computational model of the membrane potential of a neuron emerged from two Nobel prize winning scientist Hodgkin and Huxley at Cambridge University.
The membrane potential describes the activity of a neuron. Their computational model captured the shape of an action potential in a squid gaint axon. The 
pair was able to accomplish this feat with their invention of the voltage patch clamp.  This allowed them to hold the membrane potential (voltage) of the squid neuron constant 
and measure the change in current through different ion channels, such as the sodium and potassium channels. As their names suggest, the sodium channel allows sodium ions to 
flow into the neuron and the potassium channel allows potassium ions to follow out of the neuron. The flow of ions changes the potential across the membrane of the neuron and 
is the cause of the action potential. Using these measurements, they built their famous model.
	
The development of the Hodgkin-Huxley model is even more impressive considering they fit and developed their model with out modern day computers. In fact, they solved their 
differential equations with a mechanical calculator and developed the differential equations by laying graph paper over the data they collected with the proposed differential 
equations solutions drawn on the paper. Since the parameters of the ion channels were not fit using present day computational approaches, I fit the dynamics of the potassium 
ion channel data collected by Hodgkin and Huxley in their original paper using an approximate Bayesian computation sequential Monte-Carlo (ABC-SMC) approach.
