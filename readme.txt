This repository contains the Matlab implementations of the NEO-GGIW-PMBM tracker
this version divides the target into 2 subobject. Note that the method can be extended to situations where the number of subobject is greater than 2 (an algorithm for determining the number of subobject is needed in this situation). Remember change the second parameter of the function NEOpartition, the third parameter of the function initShape and so on.
More details on GGIW-PMBM can be found in the paper (accepted for publication in IEEE Transactions on Aerospace and Electric Systems)

Granstrom, Karl, Maryam Fatemi, and Lennart Svensson. "Poisson multi-Bernoulli conjugate prior for multiple extended object estimation." arXiv preprint arXiv:1605.06311 (2016).

Full text is available at https://arxiv.org/abs/1605.06311

The filters are evaluated using the generalised optimal subpattern-assignment (GOSPA) integrated with the Gaussian Wasserstein distance and We add the newly proposed extension error into GOSPA to evaluate shape errors better.

GOSPA

[C] A. S. Rahmathullah, A. F. García-Fern‡ndez, and L. Svensson, Generalized optimal sub-pattern assignment metric, in 20th International
Conference on Information Fusion, 2017.

[D] Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM
 
extension error

PMBM Filter for Non-ellipsoidal Extended Object Tracking using Sub-object Model.

=================================================================

- main.m runs the NEO-GGIW-PMBM tracker

- data association algorithm: Stochastic Optimization and NEOpartition. 

Details of Stochastic Optimization can be found in the paper

Granstršm, Karl, Lennart Svensson, Stephan Reuter, Yuxuan Xia, and Maryam Fatemi. "Likelihood-based data association for extended object tracking using sampling methods." IEEE Transactions on intelligent vehicles 3, no. 1 (2017): 30-45.

NEOpartition

PMBM Filter for Non-ellipsoidal Extended Object Tracking using Sub-object Model.
