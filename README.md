# Scaling-Analysis
Study the performance of a quantum transport parallel code based on NEGF formalism + tight bindings, that explores the electronic properties of Josephson Junctions heterostructures.

## Step 1

When studying the curren-phase-relation (CPR) of a Josephson Junction (JJ), the intensive computations are due to the operations that have to be performed in matrices of dimmension equal the junction number of sites $N$, i.e. $ N\times N$. However, as a first step in deciding which parameters might be improve to increase the performace of the parallel computation, it is explored the relation between the runtime and the number of independend k-points calculations that are needed to compute the CPR. 

> Labels:
  - number of k-points: Chunk (Chunk Size)
  - Runtime: Time of execution of the parallel code + serial code

  
