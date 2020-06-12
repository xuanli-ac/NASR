# Non-negative Adaptive Sparse Representation (NASR)
>Matlab code for calculating non-negative and sparse functional connectivity on fMRI data

The (non-negative) adaptive sparse representation is a technique to calculate the functional connectivity among nodes, thus constructing the brain functional network (represented by an association matrix). For each node, the NASR method computes its associations with all the other nodes simultaneously, rather than pairwise associations. This is achieved by using the time series of a node as a probe sample and the time series of all the other nodes as a dictionary to learn the sparse encoding vector, which represents the associations (i.e. functional connectivity).

As a sparse representation method, NASR adopts a trace norm as its regularizer. It adaptively balances the sparsity of l1-norm and the grouping-effect of l2-norm, suitable for processing highly-correlated data. 

Given a normalized fMRI data matrix (T by N), where T is the number of time points and N is the number of nodes, this code returns an association matrix (N by N).

## Description
* main.m - main function to perform the calculation
* cal_assocmatrix.m - calculate the association matrices given the fMRI data of M subjects
* solve_nntl.m - solve the optimization of the non-negative trace lasso problem of NASR, by using the alternating direction method
* get_init_lambda.m - estimate the upper bound of the regularizer parameter lambda where it achieves the most sparse solution for each sample (node) in the data matrix
* ROISignals - folder consists of a sample normalized fMRI data for test, Dn

## Usage
1. use get_init_lambda.m to determine an appropriate range of values of lambda for each sample

1. run main.m 
```
% Input:
% Dn: m*1 cell, data matrix of normalized fMRI series of m subjects
% lamlist: set the values of lambda to be tested

% Output:  
% M: m*1 cell, association matrices of m subjects, saved in a folder named 'Results_M'
% Para: obj&reconstruction error, for results control
```

## Reference
- Li, X., and Wang, H. (2015). **[Identification of functional networks in resting state fMRI data using adaptivesparse representation and affinity propagation clustering.](https://www.frontiersin.org/articles/10.3389/fnins.2015.00383/full)** *Frontiers in Neuroscience*, vol. 9, pp. 383.1-383.16.

- Li, X., Hu, Z., and Wang, H. (Oct. 2016). **[Overlapping community structure detection of brain functionalnetwork using non-negative matrix factorization.](https://link.springer.com/chapter/10.1007/978-3-319-46675-0_16)** *In: International Conference on Neural Information Processing(ICONIP 2016)*, pp. 140-147.

- Lu, C., Feng, J., Lin, Z., and Yan, S. (2013). **Correlation adaptive subspace segmentation by trace LASSO.** *In: Computer Vision (ICCV), 2013 IEEE International Conference on (Sydney, NSW:IEEE)*, pp. 1345-1352.

- Wang, J., Lu, C., Wang, M., Li, P., Yan, S., and Hu, X. (2014). **Robust face recognition via adaptive sparse representation.** *IEEE Trans. Cybern.*, vol. 44, pp. 2368-2378.

- Grave, E., Obozinski, G.R., and Bach, F.R. (2011). **Trace lasso: a trace norm regularization for correlated designs.** *In: Advances in Neural Information Processing Systems (Granada)*, pp. 2187-2195.

