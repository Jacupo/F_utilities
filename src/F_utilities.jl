module F_utilities;
using PyPlot;
using LinearAlgebra;

function Print_matrix(title, matrix)
 figure(title)
 pcolormesh(matrix)
 colorbar()
 ylim(size(matrix,1),0)
end

function Build_Omega(N)
#Build the matrix omega of dimension 2*N, that is for N sites.
  Ω                     = zeros(Complex{Float64}, 2*N, 2*N);
  ν                     = (1/(sqrt(2)));
  Ω[1:N,1:N]            = Diagonal(ν*ones(N));
  Ω[1:N,(1:N).+N]       = Diagonal(ν*ones(N));
  Ω[(1:N).+N,1:N]       = Diagonal(im*ν*ones(N));
  Ω[(1:N).+N,(1:N).+N]  = Diagonal(-im*ν*ones(N));
  return Ω;
end

function Build_FxxTxp(N)
 FxxTxp = zeros(Int64, 2*N, 2*N);
 for iiter=1:N
   FxxTxp[2*iiter-1,iiter]    = 1;
   FxxTxp[2*iiter, iiter+N]   = 1;
 end
 return FxxTxp;
end

function Build_FxpTxx(N)
 FxpTxx = zeros(Int64, 2*N, 2*N)
 for iiter=1:N
   FxpTxx[iiter,2*iiter-1]    = 1;
   FxpTxx[iiter+N, 2*iiter]   = 1;
 end
 return FxpTxx;
end

function Diag_real_skew(M, rand_perturbation::Int64=0)
 N = div(size(M,1),2);

 #Random perturbation before forcing skew symmetrisation
 if (rand_perturbation != 0)
   if (rand_perturbation == 1)
     random_M = 1*rand(2*N,2*N)*eps();
     random_M = (random_M-random_M')/2.;
     M += random_M;
   end
 end

 M = real((M-M')/2.); #Force skew-symmetry
 #Random pertubation after the skew symmetrization
 if (rand_perturbation != 0)
   if (rand_perturbation == 2)    #Perturb the diagonal elements (loose perfect skew-symmetry)
     M += diagm(rand(2*N)*eps())
   end
   if (rand_perturbation == 3)  #Perturb the whole matrix (loose perfect skew-symmetry)
     random_M = 1*rand(2*N,2*N)*eps();
     random_M = (random_M-random_M')/2.;
     M += random_M;
   end
 end

 Schur_object   = LinearAlgebra.schur(M);

 Schur_ort_i    = Schur_object.vectors;
 Schur_blocks_i = Schur_object.Schur;

 Schur_adjust = zeros(Int64, 2*N, 2*N);
 for iiter=1:N
   if (Schur_blocks_i[2*iiter,2*iiter-1] >= 0.)
     Schur_adjust[2*iiter-1, 2*iiter] = 1;
     Schur_adjust[2*iiter, 2*iiter-1] = 1;
   else
     Schur_adjust[2*iiter-1,2*iiter-1] = 1;
     Schur_adjust[2*iiter, 2*iiter]    = 1;
   end
 end

 M_temp   = Schur_adjust*Schur_blocks_i*Schur_adjust';
 O_temp   = (Schur_ort_i*Schur_adjust);

 #Sort the blocks, λ_1>=λ_2>=...>=λ_N with λ_1 the coefficient in the upper left block
 not_sorted = true;
 while not_sorted
     not_sorted = false;
     for jiter=2:(N)
       if (M_temp[2*(jiter-1)-1, 2*(jiter-1)] <= M_temp[2*jiter-1, 2*jiter])
         not_sorted = true;
         temp = M_temp[2*(jiter-1)-1, 2*(jiter-1)];
         M_temp[2*(jiter-1)-1, 2*(jiter-1)] = M_temp[2*jiter-1, 2*jiter];
         M_temp[2*(jiter-1), 2*(jiter-1)-1] = -M_temp[2*(jiter-1)-1, 2*(jiter-1)];
         M_temp[2*jiter-1, 2*jiter] = temp;
         M_temp[2*jiter, 2*jiter-1] = -temp;
         temp1 = O_temp[:,2*(jiter-1)];
         O_temp[:,2*(jiter-1)] = O_temp[:,2*jiter];
         O_temp[:,2*jiter] = temp1;
         temp1 = O_temp[:,2*(jiter-1)-1];
         O_temp[:,2*(jiter-1)-1] = O_temp[:,2*jiter-1];
         O_temp[:,2*jiter-1] = temp1;
       end
     end
     N = N-1;
 end
   M_f = Schur_blocks_i;#M_temp;
   O_f = Schur_ort_i;#real.(O_temp);
   M_f = M_temp;
   O_f = real.(O_temp);

 return M_f, O_f;
end


function Diag_h(M,rand_perturbation::Int64=0)
    N = div(size(M,1),2);

    F_xptxx = Build_FxpTxx(N);

    Ω  = Build_Omega(N);
    M_temp = real(-im*Ω*M*Ω')
    # M_temp += 0.01*rand(2*N,2*N);
    M_temp = (M_temp-M_temp')/2.;
    #M_temp[1,1] += eps();
    M_temp, O = Diag_real_skew(M_temp,rand_perturbation)
    M_temp = F_xptxx*M_temp*(F_xptxx')
    M_temp = im*Ω'*M_temp*(Ω)

    M_f = M_temp;
    U_f = Ω'*O*(F_xptxx')*Ω;

    return U_f'*M*U_f, U_f;#real(M_f), U_f;
end
export Diag_h;

function Diag_gamma(Γ,rand_perturbation::Int64=0)
 Γ = (Γ+Γ')/2.;
 γ, U = Diag_h(Γ-0.5*I,rand_perturbation);

 return U'*Γ*U,U;#real(γ+0.5*eye(size(Γ,1))),U
end


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

function Eigenvalues_of_rho(M)
   N = convert(Int64, size(M,1)/2);

   evor = ones(Float64, 2^N)

   D,U = Diag_gamma((M+M')/2.)



   for iiter=1:2^N
       index = iiter-1;
       for jiter=1:N
           evor[iiter] = evor[iiter]*round((mod(index-1,2)*D[jiter,jiter]+(1-mod(index-1,2))*(D[jiter+N,jiter+N])),18)
           index -= mod(index,2);
           index = index/2;
       end
   end

   return evor;
end

function VN_entropy(M)
   N = convert(Int64, size(M,1));

   D,U = LinearAlgebra.eigen((M+M')/2.);

   S = 0;
   for iiter=1:N
       nu = abs(round.(D[iiter];digits=14))
       if (nu != 0 && nu != 1)
           S -= log(nu)*nu;
       end
   end

   return S;
end


function Purity(M)
   N_f = convert(Int64, size(M,1)/2.)
   M[1,1] += eps()
   D,U = Diag_gamma(M);

   purity = 1;

   for iiter=1:N_f
       purity = real(purity*(2*(D[iiter,iiter]-1)*D[iiter,iiter]+1));
   end

   return purity
end

function Contour(Γ)
  N = div(size(Γ,1),2);
  D,U = Diag_gamma((Γ+Γ')/2.);

  p = zeros(Float64,N,N);

  for i=1:N
    for k=1:N
      dand = real(U[i,k]*conj(U[i,k]));
      ndda = real(U[i+N,k+N]*conj(U[i+N,k+N]));
      dada = real(U[i,k+N]*conj(U[i,k+N]));
      ndnd = real(U[i+N,k]*conj(U[i+N,k]));
      p[i,k] = 0.5*(dand+ndda+dada+ndnd);
    end
  end

  Ent_Cont = zeros(Float64, N);
  for i=1:N
    for k=1:N
      ν = real(D[k,k]);
      if (ν<0.)
        print()
        ν = 0;
      end
      if (ν>1.)
        ν = 1;
      end
      if (ν!=0.0 && ν!=1.0)
        Ent_Cont[i] -= p[i,k]*(ν*log(ν)+(1-ν)*log(1-ν));
      end
    end
  end

  return Ent_Cont;
end

function GS_gamma(D,U)
   N = div(size(D,1),2);

   Gamma_diag_base = zeros(Complex{Float64}, 2*N, 2*N);
   # for iiter=1:N
   #     Gamma_diag_base[iiter+N, iiter+N] = 1;
   # end
   # Gamma = U*Gamma_diag_base*U';
   #
   #


  for index=1:N
    if real(D[index,index])<=0
      Gamma_diag_base[index+N,index+N] = 1;
    end
    if real(D[index+N,index+N])<=0
      Gamma_diag_base[index,index] = 1;
    end
  end
  Gamma = U*Gamma_diag_base*U';
  Gamma = (Gamma+(Gamma'))/2.

  return Gamma;
end

function Thermal_fix_beta((Diag_H, U_H), beta)
 N_f   = convert(Int64, size(Diag_H,1)/2.);

 gamma = zeros(Complex{Float64}, size(Diag_H,1),size(Diag_H,1))
 for kiter=1:N_f
   e_k = Diag_H[kiter+N_f,kiter+N_f];
   gamma[kiter,kiter] = 1/(1+exp(2*beta*e_k));
   gamma[kiter+N_f,kiter+N_f] = 1/(1+exp(-2*beta*e_k))
 end

 gamma = U_H*gamma*U_H';

 return gamma;
end


function Thermal_fix_energy((Diag_H, U_H,), conserved_energy)
 N_f   = convert(Int64, size(Diag_H,1)/2.);
 a     = 0;
 b     = 100;
 beta = 0;


 temp_energy = 0;
 for iiter=1:1000
     beta  = (a+b)/2;
     temp_energy = 0
     for kiter=1:N_f
       e_k = Diag_H[kiter+N_f,kiter+N_f];
       temp_energy += e_k*(1/(1+exp(2*beta*e_k))-1/(1+exp(-2*beta*e_k)));
     end
     if (temp_energy > conserved_energy)
       a = beta;
       else
       b = beta;
     end
 end

 println("Thermal state with energy error of:", abs(temp_energy-conserved_energy)," and beta: ", beta);

 gamma = zeros(Complex{Float64}, size(Diag_H,1),size(Diag_H,1))
 for kiter=1:N_f
   e_k = Diag_H[kiter+N_f,kiter+N_f];
   gamma[kiter,kiter] = 1/(1+exp(2*beta*e_k));
   gamma[kiter+N_f,kiter+N_f] = 1/(1+exp(-2*beta*e_k))
 end

 gamma = U_H*gamma*U_H';

 return gamma, beta, abs(temp_energy-conserved_energy);
end


#Generate a random Hamiltonian with just nearest neighbour interactions
function Random_NNhamiltonian(N)
    ud  = rand(N-1).+im*rand(N-1);
    d   = rand(N).+im*rand(N);
    bd  = rand(N-1).+im*rand(N-1);
    A   = Tridiagonal(bd,d,ud);
    A   = (A+A')/2.;
    B   = Tridiagonal(bd, zeros(Complex{Float64}, N), -bd);
    H   = zeros(Complex{Float64}, 2*N, 2*N);
    H[(1:N),(1:N)]          = -conj(A);
    H[(1:N).+N,(1:N)]       = -conj(B);
    H[(1:N),(1:N).+N]       = B;
    H[(1:N).+N,(1:N).+N]    = A;

    return H;
end

function Energy(Γ,(D,U))
   N_f = convert(Int64, size(Γ,1)/2.);

   energy = 0;
   Γ = (Γ+Γ')/2.

   Γ_diag_base = real(U'*Γ*U);
   for iiter=1:(N_f)
       energy += Γ_diag_base[iiter,iiter]*D[iiter+N_f,iiter+N_f];
       energy += Γ_diag_base[iiter+N_f,iiter+N_f]*D[iiter,iiter];
   end

   return real(energy);
end

function Evolve(M,(D,U),t)
   N = div(size(M,1),2);

   M = (M+M')/2.;

   M_diag_base = U'*M*U;
   M_diag_base_evolv = exp(im*2*D*t)*M_diag_base*exp(-im*2*D*t);
   M_evolv = U*M_diag_base_evolv*(U');

   M_evolv = (M_evolv+M_evolv')/2.

   return M_evolv;
end

function Product(Γ1,Γ2)
  N = div(size(Γ1, 1),2);

  Ω = Build_Omega(N);
  γ1 = 2*Ω*Γ1*Ω'-I;
  γ2 = 2*Ω*Γ2*Ω'-I;
  γp = I-(I-γ2)*inv(I+γ1*γ2)*(I-γ1);
  return (Ω'*(1/2)*(γp+I)*Ω);
end

function Evolve_imag(Γ,D,U,t)
   N = div(size(Γ,1),2);

   Γ_diag_base = U'*Γ*U;
   Γβ = Thermal_fix_beta((D,I),t);
   Γ_diag_base_evolv = Product(Γ_diag_base,Γβ);
   Γ_diag_base_evolv = Product(Γβ,Γ_diag_base_evolv);
   Γ_evolv             = U*Γ_diag_base_evolv*U';

   return Γ_evolv;
end

function Build_hopping_hamiltonian(N,PBC=false);
  H = zeros(Float64, 2*N, 2*N);
  A = zeros(Float64, N, N);
  A[1:N,1:N] = 1/2*Tridiagonal(ones(Int64,N-1),zeros(Int64,N),ones(Int64,N-1));
  if PBC
    A[1,N] = 1/2.;
    A[N,1] = 1/2.;
  end
  H[(1:N),(1:N)]      = -A;
  H[(1:N).+N,(1:N).+N]  = A;

  return H;
end

function Build_Fourier_matrix(N)
  ω   = exp(-im*2*pi/N);
  W   = ones(Complex{Float64}, N, N);
  U_ω = zeros(Complex{Float64}, 2*N, 2*N);
  for i=1:(N-1)
    for j=1:(N-1)
      W[i,j] = ω^(i*j);
    end
  end
  W = 1/sqrt(N)*W;
  U_ω[(1:N),(1:N)]        = W;
  U_ω[(1:N).+N,(1:N).+N]  = conj.(W);

  return U_ω;
end

end
