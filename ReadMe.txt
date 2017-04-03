This folder include 6 .m files :

1 main_QST_fixedpoint
  this is the main program.The fuction of this program is using the FP_ADMM algorithm to 	  reconstruct the density matrix by numerical simulaiton in incremental measurement rates . 

2 generate_rho_outlier
  this fuction is used to generate the numerical density matrix of rank R and the outlier.
  when do the simulation without outlier, set kk and tt in main_QST_fixedpoint to 0. then 
  outlier = 0.

3 generate_AA_pauli
  this fuction is used to generate the complete Pauli measuement set. 
  the measurement matrix is a subset of AA.
  When do the simulation by different measurement sets, please change this function to
  7 generate_AA_stokes or 8 generate_AA_platonic

4 generate_A_y
  this function is used to generate A and y from AA according the measurment
  rate eta. It's a simulatoin of the process of measurement.

5 Robust_Quamtum_fixpoint
  this is the funtion of reconstructe the density matrix by FP_ADMM algorithm. X_hat is the 
  reconstructed matrix, S is the reconstructed outlier. If the simulation without outlier, set
  it to zero.
  When do the simulation by different algorithm, please change this function to
  9 Robust_Quantum_ADMM3

6 project2Hermitian_singular_shrink
  this function is used to project the matrix rho1 to a low rank and
  Conjugate symmetric matrix rho2. This function is called by Robust_Quantum_fixpoint.

7 generate_AA_stokes
  this fuction is used to generate the complete Stokes measuement set. 

8 generate_AA_platonic
  this fuction is used to generate the complete platonic solid measuement set. 

9 Robust_Quantum_ADMM3
  this is the funtion of reconstructe the density matrix by ADMM algorithm. 