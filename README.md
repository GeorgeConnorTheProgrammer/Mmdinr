# MMDINR
GPU-based modeling of material degradation in a nuclear reactor with Cluster Dynamics and Mean Field Rate Theory

## Mean Field Rate Theory (MFRT)
### USAGE: 
  ```
  g++ mfrt.cpp
  ./a.out
  ```

### CONFIGURATION:
  Edit values in the `config.json` file. 
 - `total_time_seconds` == total time you want the simulation to run
 - `dt_seconds` == size of the model's increments
 - `temperature_kelvin` == temperature of the system in Kelvin
 - `K_0_exp` == defines the defect production rate. `K0` is calculated as 10 ^ `K_0_exp`
 - `C_s_exp` == defines the sink strength. `Cs` is calculated as 10 ^ `C_s_exp`
 - `sample_interval` == how often to take data points from the model and output them to the .csv file
### DESCRIPTION:
This program models the rate of change of the concentration of interstitials and vacancies in a material, using the following equations: ![MFRT Equations](https://github.com/GeorgeConnorTheProgrammer/Mmdinr/assets/148592312/2d116231-c031-4122-a44b-e0581b6d63d3)

The first equation refers to the rate of change of the concentration of vacancies and the second refers to the concentration of interstitials. 

The K0 term refers to the defect production rate. This is some predetermined value that comes from the circumstances of the system. 

The second term in both equations refers to the rate of defects being removed from interstitials and vacancies meeting (out of place atoms moving to empty spaces in the structure). This term is determined by the product of some constant predetermined rate (Kiv), the concentration of interstitials (Ci), and the concentration of vacancies (Cv). 

The third term in both equations refers to the rate of defects (vacancies or interstitials) being removed by sinks. This rate is determined as the product of the rate of point defects being absorbed by sinks (Kvs or Kis), the concentration of point defects (Cv or Ci), and the concentration of sinks (Cs).

Notice that both the second and third term contribute to a lower rate of increasing defect concentration.

### DEFINITIONS:
_These definitions are meant to provide a basic understanding of the program, and do not go in depth._
- **interstitial**: Atom out of place
- **vacancy**: Empty spot left by an atom out of place
- **point defect**: Any "disruption" to the crystal structure. In this case, either an interstitial or a vacancy.
- **sink**: Some force/preexisting microstructure that causes the concentration of point defects to go down. For example, the tendency for an "air bubble" to rise to the surface of a metal and disappear could be considered a sink.
=======