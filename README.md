 ## <p align="center">CCWSIM: An Efficient and Fast Wavelet-Based CCSIM for Categorical Characterization in Large-Scale Geological Domains
 </p>


<p align="justify"> We propose a new implementation of Multiple-point Simulation (MPS) that combines CCSIM with the Discrete Wavelet Transform (DWT) to rapidly simulate grids with millions of cells. Beyond its computational efficiency, the proposed methodology prioritizes the maintenance of spatial continuity and adherence to data conditioning principles to a considerable degree. 
<p align="justify">For swift testing and easy reproduction of the results, please clone the repository and execute the "CCWSIM_runme.m" file. Ensure that all required functions and MAT files are in the designated path. For a comprehensive guide, please refer to the Help file.
</p>


#### -Schematic illustration of the identification and reconstruction of best-matched patterns using the CCW score map in CCWSIM algorithm.

<p align="center">
  <img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/CCWSIM_Algorithm.png" width=90% height=90%>
</p>

## Experiments
#### -(a) The channel training image and (b) three unconditional realizations obtained using the proposed method (average CPU time: 6.52 seconds).

<p align="center">
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/UnconditionalSim.png" width=90% height=90%>
</p>

#### -Conditional simulation using the proposed method: (a) three distinct conditional simulations and (b) the ensemble average map obtained from 100 conditional realizations using the proposed method (average CPU time: 2.42 seconds).
<p align="center">
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/ConditionalSim.png" width=90% height=90%>
</p>

#### -(a)The variograms and (b) the probability of connection for the channel training image and 100 realizations obtained from the proposed method.
<p align="center">
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/Validation.png" width=90% height=90%>
</p>

#### -Conditional realizations of the Gange delta TI obtained using the proposed method (average CPU time: 7.7 seconds). The blue and red circles mark the positions of the conditional data points
<p align="center">
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/GangDelatSim.png" width=90% height=90%>
</p>

#### -Analysis of distance (ANODI): Between and within-realization variability ratio for three resolution levels obtained from the proposed method and MS-CCSIM.
<p align="center">
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/ANODIST.png" width=70% height=70%>
</p>

#### -Conditional simulation obtained through the proposed method and MS-CCSIM algorithm using different multi-scale, DWT levels and conditional hard data.
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/Benchmark.png" width=90% height=90%>

#### -Conditional simulation obtained from the proposed method and MS-CCSIM algorithm using stone wall training image.
<p align="center">
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/StoneWallSim.png" width=90% height=90%>
</p>

#### -Average CPU times for five different simulations grid size using TI in Fig. 6a at different DWT levels (1, 2 and 3).
<p align="center">
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/CpuTime.png" width=70% height=70%>
</p>

#### -Reference image and (b) one  simulation obtained using 4000 conditional data through CCWSIM (average CPU time 1.06 sec).
<img src="https://github.com/MBS1984/CCWSIM/blob/main/Figures/DenseSim.png" width=90% height=90%>














