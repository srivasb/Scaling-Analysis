# Scaling-Analysis
Study the performance of a quantum transport parallel code based on NEGF formalism + tight bindings, that explores the electronic properties of Josephson Junctions heterostructures.

## Step 1

When studying the curren-phase-relation (CPR) of a Josephson Junction (JJ), super-noraml-super (SNS), the intensive computations are due to the operations that have to be performed in matrices of dimmension equal the junction number of sites $N$, i.e. $ N \times N$. However, as a first step in deciding which parameters might be improve to increase the performace of the parallel computation, it is explored the relation between the runtime and the number of independend k-points calculations that are needed to compute the CPR. 

> Labels:
  - number of k-points: Chunk (Chunk Size)
  - Runtime: Time of execution of the parallel code + serial code

> Metric Quantities to Study the performance: Fixed the following parameters and only change the Chunk size"
  - Junction Lenght = $127.15\,nm$ (this is the length of the N-regin o f the junction.
  - k-points Infinite Mass boundary conditions (Width >> Length)
  - pair potential: $\Delta_0=0.0012\,eV$
  - vargamma_SC (hopping between N-region and S-leads ): $\vargamm_{SC}=2.97*0.67$ 
  - $EF = 0.350\,eV$
  - $T = 0.04\,k$
  - Valley Zeeman SOC ($\lambda_{VZ}$): $\lambda_{VZ}=0.002\, eV$
  - Inplane magnetic field = 0
  - Rashba SOC = 0

> Results:

  - Runtime (h)
    
  ![Alternative Text](images/step1_chunksize_vs_runtime_1.pdf)
  
  - Throughtput: $T=\frac{\text{Chunk Size}}{h}$
