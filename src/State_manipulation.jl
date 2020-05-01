"""
  Permute_rc(M,u,d)

  Return the matric `M` where row `u` has been exchanged with row `d`,
  and column `u`  has been exchanged with column `d`.
"""
function Permute_rc(M,u,d)
    temp = M[u,:];
    M[u,:] = M[d,:];
    M[d,:] = temp;
    temp = M[:,u];
    M[:,u] = M[:,d];
    M[:,d] = temp;
    return M;
end

"""
  Reduce_gamma(M, N_partition, first_index)

  If `M` is a Dirac correlation matrix, returns the reduced density matrix of the system
  `[first_index, first_index+1,...,first_index+N_partition-1]`.
"""
function Reduce_gamma(M, N_partition, first_index)
   N_f = div(size(M,1),2);
   first_index = first_index-1;
   periodic_dimension = max((N_partition.+first_index-N_f),0)
   dim_UL = N_partition-periodic_dimension;

   redgamma = zeros(Complex{Float64}, N_partition*2, N_partition*2);
   #Copy the upper left left part of the correlation matrix
   redgamma[1:dim_UL,1:dim_UL] = M[(1:dim_UL).+first_index,(1:dim_UL).+first_index];
   redgamma[(1:dim_UL).+N_partition,1:dim_UL] = M[(1:dim_UL).+N_f.+first_index,(1:dim_UL).+first_index];
   redgamma[1:dim_UL,(1:dim_UL).+N_partition] = M[(1:dim_UL).+first_index,(1:dim_UL).+N_f.+first_index];
   redgamma[(1:dim_UL).+N_partition,(1:dim_UL).+N_partition] = M[(1:dim_UL).+N_f.+first_index,(1:dim_UL).+N_f.+first_index];

   if (periodic_dimension>0)
     redgamma[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension))] = M[1:periodic_dimension,1:periodic_dimension];
     redgamma[1:dim_UL,(dim_UL.+(1:periodic_dimension))] = M[(first_index.+(1:dim_UL)),1:periodic_dimension];
     redgamma[(dim_UL.+(1:periodic_dimension)),1:dim_UL] = M[1:periodic_dimension,(first_index.+(1:dim_UL))];

     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,(dim_UL.+(1:periodic_dimension))] = M[(1:periodic_dimension).+N_f,1:periodic_dimension];
     redgamma[(1:dim_UL).+N_partition,(dim_UL.+(1:periodic_dimension))] = M[(first_index.+(1:dim_UL)).+N_f,1:periodic_dimension];
     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,1:dim_UL] = M[(1:periodic_dimension).+N_f,(first_index.+(1:dim_UL))];

     redgamma[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension)).+N_partition] = M[1:periodic_dimension,(1:periodic_dimension).+N_f];
     redgamma[1:dim_UL,(dim_UL.+(1:periodic_dimension)).+N_partition] = M[(first_index.+(1:dim_UL)),(1:periodic_dimension).+N_f];
     redgamma[(dim_UL.+(1:periodic_dimension)),(1:dim_UL).+N_partition] = M[1:periodic_dimension,(first_index.+(1:dim_UL)).+N_f];

     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,(dim_UL.+(1:periodic_dimension)).+N_partition] = M[(1:periodic_dimension).+N_f,(1:periodic_dimension).+N_f];
     redgamma[(1:dim_UL).+N_partition,(dim_UL.+(1:periodic_dimension)).+N_partition] = M[(first_index.+(1:dim_UL)).+N_f,(1:periodic_dimension).+N_f];
     redgamma[(dim_UL.+(1:periodic_dimension)).+N_partition,(1:dim_UL).+N_partition] = M[(1:periodic_dimension).+N_f,(first_index.+(1:dim_UL)).+N_f];
   end

   return redgamma
end

"""
  Inject_gamma(gamma, injection, first_index)

  If `gamma` is a Dirac correlation matrix, returns a new gamma such that
  `Reduce_gamma(gamma, size(injection,1), first_index) = injection`.
"""
function Inject_gamma(gamma, injection, first_index)
 dim_gamma     = div(size(gamma, 1),2);
 dim_injection = div(size(injection, 1), 2);

 first_index = first_index-1;
 periodic_dimension = max((dim_injection+first_index-dim_gamma),0)
 dim_UL = dim_injection-periodic_dimension;

 #Injecto la parte Z nei 4 riquadri
 gamma[(1:dim_UL).+first_index,(1:dim_UL).+first_index]                      = injection[(1:dim_UL),(1:dim_UL)]
 gamma[(1:dim_UL).+(first_index+dim_gamma),(1:dim_UL).+first_index]            = injection[(1:dim_UL).+dim_injection,(1:dim_UL)]
 gamma[(1:dim_UL).+first_index,(1:dim_UL).+first_index.+dim_gamma]            = injection[(1:dim_UL),(1:dim_UL).+dim_injection]
 gamma[(1:dim_UL).+first_index.+dim_gamma,(1:dim_UL).+first_index.+dim_gamma]  = injection[(1:dim_UL).+dim_injection,(1:dim_UL).+dim_injection]


 if (periodic_dimension>0)
   #Injecto A,B,C  per ogni riquadro
   gamma[1:periodic_dimension,1:periodic_dimension] = injection[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension))];
   gamma[(first_index.+(1:dim_UL)),1:periodic_dimension] = injection[1:dim_UL,(dim_UL.+(1:periodic_dimension))];
   gamma[1:periodic_dimension,(first_index.+(1:dim_UL))] = injection[(dim_UL.+(1:periodic_dimension)),1:dim_UL];

   gamma[(1:periodic_dimension).+dim_gamma,1:periodic_dimension] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,(dim_UL.+(1:periodic_dimension))];
   gamma[(first_index.+(1:dim_UL)).+dim_gamma,1:periodic_dimension] = injection[(1:dim_UL).+dim_injection,(dim_UL.+(1:periodic_dimension))];
   gamma[(1:periodic_dimension).+dim_gamma,(first_index.+(1:dim_UL))] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,1:dim_UL];

   gamma[1:periodic_dimension,(1:periodic_dimension).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)),(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[(first_index.+(1:dim_UL)),(1:periodic_dimension).+dim_gamma] = injection[1:dim_UL,(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[1:periodic_dimension,(first_index.+(1:dim_UL)).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)),(1:dim_UL).+dim_injection];

   gamma[(1:periodic_dimension).+dim_gamma,(1:periodic_dimension).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[(first_index.+(1:dim_UL)).+dim_gamma,(1:periodic_dimension).+dim_gamma] = injection[(1:dim_UL).+dim_injection,(dim_UL.+(1:periodic_dimension)).+dim_injection];
   gamma[(1:periodic_dimension).+dim_gamma,(first_index.+(1:dim_UL)).+dim_gamma] = injection[(dim_UL.+(1:periodic_dimension)).+dim_injection,(1:dim_UL).+dim_injection];
 end

 return gamma
end

"""
  Project_diagonals(M4,off_diagonals)

  Returns a 4-blocks matrix, in wich in each block only
  the first `off_diagonals` diagonals off diagonal are mantained
  the rest is set to 0.
  If `off_diagonals=0` then it mantains only the diagonal of each block
"""
function Project_diagonals(M4,off_diagonals)
 N_f = convert(Int64, size(M4,1)/2.)
 M_finale = zeros(Complex{Float64}, 2*N_f, 2*N_f)

 for iiter=1:N_f
   M_finale[iiter, iiter]          = M4[iiter, iiter]
   M_finale[iiter+N_f, iiter]      = M4[iiter+N_f, iiter]
   M_finale[iiter, iiter+N_f]      = M4[iiter, iiter+N_f]
   M_finale[iiter+N_f, iiter+N_f]  = M4[iiter+N_f, iiter+N_f]
   for jiter=1:off_diagonals
     M_finale[iiter, mod(iiter+jiter-1,N_f)+1]          = M4[iiter, mod(iiter+jiter-1,N_f)+1]
     M_finale[iiter+N_f, mod(iiter+jiter-1,N_f)+1]      = M4[iiter+N_f, mod(iiter+jiter-1,N_f)+1]
     M_finale[iiter, mod(iiter+jiter-1,N_f)+1+N_f]      = M4[iiter, mod(iiter+jiter-1,N_f)+1+N_f]
     M_finale[iiter+N_f, mod(iiter+jiter-1,N_f)+1+N_f]  = M4[iiter+N_f, mod(iiter+jiter-1,N_f)+1+N_f]

     M_finale[iiter, mod(iiter-jiter-1,N_f)+1]          = M4[iiter, mod(iiter-jiter-1,N_f)+1]
     M_finale[iiter+N_f, mod(iiter-jiter-1,N_f)+1]      = M4[iiter+N_f, mod(iiter-jiter-1,N_f)+1]
     M_finale[iiter, mod(iiter-jiter-1,N_f)+1+N_f]      = M4[iiter, mod(iiter-jiter-1,N_f)+1+N_f]
     M_finale[iiter+N_f, mod(iiter-jiter-1,N_f)+1+N_f]  = M4[iiter+N_f, mod(iiter-jiter-1,N_f)+1+N_f]
   end
 end
 return M_finale
end


"""
  Build_GDE(g,U)

  If `g` is a Dirac correlation matrix, and `U` is a fermionic transormation that diagonalise
  an f.q.h. `H`, `Build_GDE(g,U)` returns the Gaussian Diagonal Ensemble of G with respect to `H`.
"""
function Build_GDE(g,U)
  return (U*Project_diagonals(U'*g*U,0)*U');;
end
