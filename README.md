# PCM
Pothole Cascade Model - simulating effects of storage in depressions in the Canadian Prairies

The Pothole Cascade Model is a conceptual model of the behaviour of depressions in the Canadian Prairies.
It is very much based on the work of Shaw et al., 2012. The behavior of PCM is described in Shook et al., 2013.
The PCM algorithms have also been incorporated in a CRHM (Cold Regions Hydrological Modelling) platform model of Smith Creek (Pomeroy et al., 2014).

The PCM is a command-line program written in Fortran 95. Using the program requires compilation and construction of input forcing files. There is no form of GUI.  

When PCM is executed, it applies a single flux of water (input or output) to sets of depressions. Simulation of a sequence of operations will require a script (such as bash on Linux) or calling the program from another language.


References

Pomeroy, J.W., Shook, K., Fang, X., Dumanski, S., Westbrook, C., Brown, T., 2014. Improving and Testing the Prairie Hydrological Model at Smith Creek Research Basin (No. 14). Centre for Hydrology, Saskatoon, Saskatchewan.

Shaw, D.A., Pietroniro, A., Martz, L.W., 2012. Topographic analysis for the prairie pothole region of Western Canada. Hydrological Processes 27, 3105–3114. https://doi.org/10.1002/hyp.9409

Shook, K., Pomeroy, J.W., Spence, C., Boychuk, L., 2013. Storage dynamics simulations in prairie wetland hydrology models: evaluation and parameterization. Hydrological Processes 27, 1875–1889. https://doi.org/10.1002/hyp.9867

