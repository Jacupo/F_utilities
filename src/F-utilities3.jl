using PyPlot;

function Print_matrix(title, matrix)
 figure(title)
 pcolormesh(matrix)
 colorbar()
 ylim(size(matrix,1),0)
end

function Build_Omega(N_f)
#Build the matrix omega of dimension 2*N_f, that is for N_f fermions.
 Omega                                 = zeros(Complex{Float64}, 2*N_f, 2*N_f);
 Omega[1:N_f,1:N_f]                    = (1/(sqrt(2)))*eye(N_f);
 Omega[1:N_f,(1:N_f)+N_f]              = (1/(sqrt(2)))*eye(N_f);
 Omega[(1:N_f)+N_f,1:N_f]              = im*(1/(sqrt(2)))*eye(N_f);
 Omega[(1:N_f)+N_f,(1:N_f)+N_f]        = -im*(1/(sqrt(2)))*eye(N_f);
 return Omega
end

function Build_FxxTxp(N_f)
 FxxTxp = zeros(Int64, 2*N_f, 2*N_f);
 for iiter=1:N_f
   FxxTxp[2*iiter-1,iiter]    = 1;
   FxxTxp[2*iiter, iiter+N_f] = 1;
 end
 return FxxTxp
end

function Build_FxpTxx(N_f)
 FxpTxx = zeros(Int64, 2*N_f, 2*N_f)
 for iiter=1:N_f
   FxpTxx[iiter,2*iiter-1]    = 1;
   FxpTxx[iiter+N_f, 2*iiter] = 1;
 end

 return FxpTxx
end

function Diag_real_skew(M)

 N_f = convert(Int64, size(M,1)/2.)


 M = real((M-M')/2.)
 # M += diagm(rand(2*N_f)*eps()) #Se la diagonale è proprio zero ho problemi di
                                  #convergenza
 # random_M = 1*rand(2*N_f,2*N_f)*eps();#Per SETTORI MOMENTO questo bastava 10 e non 100
 # random_M = (random_M-random_M')/2.;
 # M += random_M;
 Schur_object   = schurfact(M);
 #Per non avere errori di convergenza
 # Schur_object   = schurfact(M);
 #Forse questo risolve dei problemi numerici
 Schur_ort_i    = Schur_object[:vectors];
 Schur_blocks_i = Schur_object[:Schur];



 Schur_adjust = zeros(Int64, 2*N_f, 2*N_f);
 for iiter=1:(convert(Int64,N_f))
   if (Schur_blocks_i[2*iiter,2*iiter-1] >= 0.)
     Schur_adjust[2*iiter-1, 2*iiter] = 1;
     Schur_adjust[2*iiter, 2*iiter-1] = 1;
   else
     Schur_adjust[2*iiter-1,2*iiter-1] = 1;
     Schur_adjust[2*iiter, 2*iiter]    = 1;
   end
 end

 # #AGGIUNTO! PER TEST
 # Schur_adjust = eye(2*N_f)
 # #

 M_temp   = Schur_adjust*Schur_blocks_i*Schur_adjust';
 O_temp   = (Schur_ort_i*Schur_adjust);


 not_sorted = true;    ##DEVE ESSERE LASCIATO A TRUE PER ORDINARE
 while not_sorted
     not_sorted = false;
     for jiter=2:(N_f)
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
     N_f = N_f-1;
 end

   M_f = M_temp;
   O_f = real(O_temp);

 return M_f, O_f
end

function Diag_ferm(M)
    N_f = convert(Int64, size(M,1)/2.)

    F_xptxx = Build_FxpTxx(N_f);

    Omega  = Build_Omega(N_f)
    M_temp = real(-im*Omega*M*Omega')
    # M_temp += 0.01*rand(2*N_f,2*N_f);
    M_temp = (M_temp-M_temp')/2.;
    #M_temp[1,1] += eps();
    M_temp, O = Diag_real_skew(M_temp)
    M_temp = F_xptxx*M_temp*(F_xptxx')
    M_temp = im*Omega'*M_temp*(Omega)

    M_f = M_temp;
    U_f = Omega'*O*(F_xptxx')*Omega;

    return real(M_f), U_f
end

function Diag_gamma(Γ)
 #Diagonalizzo cosi' invece che con eig cosi' ottengo gli autovettori a coppie
 #ordinate come voglio.
 #La generica unitaria tra Fermioni di Dirac e' fatta come
 #Ω*O*Ω', con O una generica ortogonale NxN specificata da N(N-1)/2 angoli
 Γ = (Γ+Γ')/2.;
 γ, U = Diag_ferm(Γ-0.5*eye(size(Γ,1)));

 return U'*Γ*U,U#real(γ+0.5*eye(size(Γ,1))),U
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

function Evolve_gamma(M,D,U,t)
   N_f = convert(Int64, size(M,1)/2.);

   #####AGGIUNTO
   M = (M+M')/2.
   # M[1,1] += eps();
   ####

   M_diag_base = U'*M*U;
   M_diag_base_evolv = expm(im*D*t)*M_diag_base*expm(-im*D*t);
   M_evolv = U*M_diag_base_evolv*(U');

   #####AGGIUNTO
   M_evolv = (M_evolv+M_evolv')/2.
   ####
   return M_evolv;
end

function Evolve_gamma_imag(M,D,U,t)
   N_f = convert(Int64, size(M,1)/2.);

   #####AGGIUNTO
   M = (M+M')/2.
   M[1,1] += eps();
   ####

   M_diag_base = U'*M*U;
   A           = 2*(D*M_diag_base-M_diag_base*D)
   U_t         = expm(A*t)
   M_diag_base_evolv = U_t*M_diag_base*U_t';
   M_evolv             = U*M_diag_base_evolv*U';

   #####AGGIUNTO
   M_evolv = (M_evolv+M_evolv')/2.
   ####

   return M_evolv;
end


function Energy_fermion(Mat,D,U)
   N_f = convert(Int64, size(Mat,1)/2.);

   energy = 0;

   M = deepcopy(Mat)
   #####AGGIUNTO
   M = (M+M')/2.
   M[1,1] += eps()
   ####

   M_diag_base = real(U'*M*U);
   for iiter=1:(N_f)
       energy += M_diag_base[iiter,iiter]*D[iiter+N_f,iiter+N_f];
       energy += M_diag_base[iiter+N_f,iiter+N_f]*D[iiter,iiter];
   end

   return real(energy);
end



function Build_A_NOPBC(Jx, Jy, lambdas)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1   ##SEMBRA CHE NON DEVO
   M_A       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_A[iiter, iiter]   = -2*lambdas[iiter];
       M_A[iiter, iiter-1] = -(Jx[iiter-1]+Jy[iiter-1]);
       M_A[iiter, iiter+1] = -(Jx[iiter]+Jy[iiter]);
   end
   M_A[1,1]         = -2*lambdas[1];#
   M_A[1,2]         = -(Jx[1]+Jy[1]);#

   # M_A[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#
   # M_A[dimension,1] = (-(-1)^(parity))*-(Jx[1]+Jy[1]);
   M_A[dimension, dimension-1] = -(Jx[dimension-1] + Jy[dimension-1]);#
   M_A[dimension,dimension]    = -2*lambdas[dimension];#

   return M_A
end


function Build_A(Jx, Jy, lambdas)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1   ##SEMBRA CHE NON DEVO
   M_A       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_A[iiter, iiter]   = -2*lambdas[iiter];
       M_A[iiter, iiter-1] = -(Jx[iiter-1]+Jy[iiter-1]);   #Si connettono sempre con il coeff. del più piccolo
       M_A[iiter, iiter+1] = -(Jx[iiter]+Jy[iiter]);       #Si connettono sempre con il coeff. del più piccolo
   end
   M_A[1,1]         = -2*lambdas[1];#
   M_A[1,2]         = -(Jx[1]+Jy[1]);#

   M_A[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#
   M_A[dimension,1] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#Questo è l'unico caso dove connetto con il coeff del più grande
   M_A[dimension, dimension-1] = -(Jx[dimension-1] + Jy[dimension-1]);   #
   M_A[dimension,dimension]    = -2*lambdas[dimension];#

   return M_A
end


function Build_A_APBC(Jx, Jy, lambdas)
   dimension = size(Jx, 1);
   parity = mod(dimension+1, 2); #!!ATTENZIONE HO CAMBIATO con +1   ##SEMBRA CHE NON DEVO
   M_A       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_A[iiter, iiter]   = -2*lambdas[iiter];
       M_A[iiter, iiter-1] = -(Jx[iiter-1]+Jy[iiter-1]);
       M_A[iiter, iiter+1] = -(Jx[iiter]+Jy[iiter]);
   end
   M_A[1,1]         = -2*lambdas[1];#
   M_A[1,2]         = -(Jx[1]+Jy[1]);#

   M_A[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);#
   M_A[dimension,1] = (-(-1)^(parity))*-(Jx[dimension]+Jy[dimension]);
   M_A[dimension, dimension-1] = -(Jx[dimension-1] + Jy[dimension-1]);#
   M_A[dimension,dimension]    = -2*lambdas[dimension];#

   return M_A
end


function Build_B_NOPBC(Jx, Jy)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1 ##SEMBRA CHE NON DEVO
   M_B       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_B[iiter, iiter-1] = -(Jx[iiter-1]-Jy[iiter-1]);
       M_B[iiter, iiter+1] = Jx[iiter]-Jy[iiter];
   end
   M_B[1,2]         = Jx[1]-Jy[1];

   # M_B[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]-Jy[dimension]);#
   # M_B[dimension,1] = (-(-1)^(parity))*(Jx[1]-Jy[1]);#

   M_B[dimension, dimension-1] = -(Jx[dimension-1]-Jy[dimension-1]);#

   return M_B
end


function Build_B(Jx, Jy)
   dimension = size(Jx, 1);
   parity = mod(dimension, 2); #!!ATTENZIONE HO CAMBIATO con +1 ##SEMBRA CHE NON DEVO
   M_B       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_B[iiter, iiter-1] = -(Jx[iiter-1]-Jy[iiter-1]);
       M_B[iiter, iiter+1] = Jx[iiter]-Jy[iiter];
   end
   M_B[1,2]         = Jx[1]-Jy[1];

   M_B[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]-Jy[dimension]);#
   M_B[dimension,1] = (-(-1)^(parity))*(Jx[dimension]-Jy[dimension]);#

   M_B[dimension, dimension-1] = -(Jx[dimension-1]-Jy[dimension-1]);#

   return M_B
end

function Build_B_APBC(Jx, Jy)
   dimension = size(Jx, 1);
   parity = mod(dimension+1, 2); #!!ATTENZIONE HO CAMBIATO con +1 ##SEMBRA CHE NON DEVO
   M_B       = zeros(Float64, dimension, dimension)
   for iiter=2:dimension-1
       M_B[iiter, iiter-1] = -(Jx[iiter-1]-Jy[iiter-1]);
       M_B[iiter, iiter+1] = Jx[iiter]-Jy[iiter];
   end
   M_B[1,2]         = Jx[1]-Jy[1];

   M_B[1,dimension] = (-(-1)^(parity))*-(Jx[dimension]-Jy[dimension]);#
   M_B[dimension,1] = (-(-1)^(parity))*(Jx[dimension]-Jy[dimension]);#

   M_B[dimension, dimension-1] = -(Jx[dimension-1]-Jy[dimension-1]);#

   return M_B
end

function TFI_Hamiltonian(N,Jx,Jy,lambda)
   A = 0.5*Build_A(Jx,Jy,lambda);
   B = 0.5*Build_B(Jx,Jy);

   H_TFI                          = zeros(Float64, 2*N, 2*N);
   H_TFI[1:N,1:N]             = -A;
   H_TFI[((1:N).+N),1:N]       = -B;
   H_TFI[1:N,(1:N).+N]       =  B;
   H_TFI[(1:N).+N,(1:N).+N] =  A;

   return H_TFI;
end



function TFI_Hamiltonian_NOPBC(N,Jx,Jy,lambda)
   A = 0.5*Build_A_NOPBC(Jx,Jy,lambda);
   B = 0.5*Build_B_NOPBC(Jx,Jy);

   H_TFI                          = zeros(Float64, 2*N, 2*N);
   H_TFI[1:N,1:N]             = -A;
   H_TFI[((1:N).+N),1:N]       = -B;
   H_TFI[1:N,(1:N).+N]       =  +B;
   H_TFI[(1:N).+N,(1:N).+N] =  A;

   return H_TFI;
end

function TFI_Hamiltonian_APBC(N,Jx,Jy,lambda)
   A = 0.5*Build_A_APBC(Jx,Jy,lambda);
   B = 0.5*Build_B_APBC(Jx,Jy);

   H_TFI                          = zeros(Float64, 2*N, 2*N);
   H_TFI[1:N,1:N]             = -A;
   H_TFI[((1:N).+N),1:N]       = -B;
   H_TFI[1:N,(1:N).+N]       =  B;
   H_TFI[(1:N).+N,(1:N).+N] =  A;

   return H_TFI;
end

function Gaussian_Hamiltonian(N,A,B)

   H_TFI                          = zeros(Complex{Float64}, 2*N, 2*N);
   H_TFI[1:N,1:N]             = -A;
   H_TFI[((1:N).+N),1:N]       = -B;
   H_TFI[1:N,(1:N).+N]       =  B;
   H_TFI[(1:N).+N,(1:N).+N] =  A;

   return H_TFI;
end


function GS_Gamma(D,U)
   N = convert(Int64, size(D,1)/2);

   Gamma_diag_base = zeros(Complex{Float64}, 2*N, 2*N);
   for iiter=1:N
       Gamma_diag_base[iiter+N, iiter+N] = 1;
   end
   Gamma = U*Gamma_diag_base*U';

   Gamma = (Gamma+(Gamma'))/2.
   return Gamma;
end


function Reduce_gamma(M, N_partition, first_index)
   N_f = convert(Int64, size(M,1)/2.)
   first_index = first_index-1;
   periodic_dimension = max((N_partition+first_index-N_f),0)
   dim_UL = N_partition-periodic_dimension;

   redgamma = zeros(Complex{Float64}, N_partition*2, N_partition*2);
   #Copy the upper left left part of the correlation matrix
   redgamma[1:dim_UL,1:dim_UL] = M[(1:dim_UL)+first_index,(1:dim_UL)+first_index];
   redgamma[(1:dim_UL)+N_partition,1:dim_UL] = M[(1:dim_UL)+N_f+first_index,(1:dim_UL)+first_index];
   redgamma[1:dim_UL,(1:dim_UL)+N_partition] = M[(1:dim_UL)+first_index,(1:dim_UL)+N_f+first_index];
   redgamma[(1:dim_UL)+N_partition,(1:dim_UL)+N_partition] = M[(1:dim_UL)+N_f+first_index,(1:dim_UL)+N_f+first_index];

   if (periodic_dimension>0)
     redgamma[(dim_UL+(1:periodic_dimension)),(dim_UL+(1:periodic_dimension))] = M[1:periodic_dimension,1:periodic_dimension];
     redgamma[1:dim_UL,(dim_UL+(1:periodic_dimension))] = M[(first_index+(1:dim_UL)),1:periodic_dimension];
     redgamma[(dim_UL+(1:periodic_dimension)),1:dim_UL] = M[1:periodic_dimension,(first_index+(1:dim_UL))];

     redgamma[(dim_UL+(1:periodic_dimension))+N_partition,(dim_UL+(1:periodic_dimension))] = M[(1:periodic_dimension)+N_f,1:periodic_dimension];
     redgamma[(1:dim_UL)+N_partition,(dim_UL+(1:periodic_dimension))] = M[(first_index+(1:dim_UL))+N_f,1:periodic_dimension];
     redgamma[(dim_UL+(1:periodic_dimension))+N_partition,1:dim_UL] = M[(1:periodic_dimension)+N_f,(first_index+(1:dim_UL))];

     redgamma[(dim_UL+(1:periodic_dimension)),(dim_UL+(1:periodic_dimension))+N_partition] = M[1:periodic_dimension,(1:periodic_dimension)+N_f];
     redgamma[1:dim_UL,(dim_UL+(1:periodic_dimension))+N_partition] = M[(first_index+(1:dim_UL)),(1:periodic_dimension)+N_f];
     redgamma[(dim_UL+(1:periodic_dimension)),(1:dim_UL)+N_partition] = M[1:periodic_dimension,(first_index+(1:dim_UL))+N_f];

     redgamma[(dim_UL+(1:periodic_dimension))+N_partition,(dim_UL+(1:periodic_dimension))+N_partition] = M[(1:periodic_dimension)+N_f,(1:periodic_dimension)+N_f];
     redgamma[(1:dim_UL)+N_partition,(dim_UL+(1:periodic_dimension))+N_partition] = M[(first_index+(1:dim_UL))+N_f,(1:periodic_dimension)+N_f];
     redgamma[(dim_UL+(1:periodic_dimension))+N_partition,(1:dim_UL)+N_partition] = M[(1:periodic_dimension)+N_f,(first_index+(1:dim_UL))+N_f];
   end


   return redgamma
end

function Inject_gamma(gamma, injection, first_index)
 dim_gamma     = convert(Int64, size(gamma, 1)/2.0);
 dim_injection = convert(Int64, size(injection, 1)/2.0);

 first_index = first_index-1;
 periodic_dimension = max((dim_injection+first_index-dim_gamma),0)
 dim_UL = dim_injection-periodic_dimension;

 #Injecto la parte Z nei 4 riquadri
 gamma[(1:dim_UL)+first_index,(1:dim_UL)+first_index]                      = injection[(1:dim_UL),(1:dim_UL)]
 gamma[(1:dim_UL)+first_index+dim_gamma,(1:dim_UL)+first_index]            = injection[(1:dim_UL)+dim_injection,(1:dim_UL)]
 gamma[(1:dim_UL)+first_index,(1:dim_UL)+first_index+dim_gamma]            = injection[(1:dim_UL),(1:dim_UL)+dim_injection]
 gamma[(1:dim_UL)+first_index+dim_gamma,(1:dim_UL)+first_index+dim_gamma]  = injection[(1:dim_UL)+dim_injection,(1:dim_UL)+dim_injection]


 if (periodic_dimension>0)
   #Injecto A,B,C  per ogni riquadro
   gamma[1:periodic_dimension,1:periodic_dimension] = injection[(dim_UL+(1:periodic_dimension)),(dim_UL+(1:periodic_dimension))];
   gamma[(first_index+(1:dim_UL)),1:periodic_dimension] = injection[1:dim_UL,(dim_UL+(1:periodic_dimension))];
   gamma[1:periodic_dimension,(first_index+(1:dim_UL))] = injection[(dim_UL+(1:periodic_dimension)),1:dim_UL];

   gamma[(1:periodic_dimension)+dim_gamma,1:periodic_dimension] = injection[(dim_UL+(1:periodic_dimension))+dim_injection,(dim_UL+(1:periodic_dimension))];
   gamma[(first_index+(1:dim_UL))+dim_gamma,1:periodic_dimension] = injection[(1:dim_UL)+dim_injection,(dim_UL+(1:periodic_dimension))];
   gamma[(1:periodic_dimension)+dim_gamma,(first_index+(1:dim_UL))] = injection[(dim_UL+(1:periodic_dimension))+dim_injection,1:dim_UL];

   gamma[1:periodic_dimension,(1:periodic_dimension)+dim_gamma] = injection[(dim_UL+(1:periodic_dimension)),(dim_UL+(1:periodic_dimension))+dim_injection];
   gamma[(first_index+(1:dim_UL)),(1:periodic_dimension)+dim_gamma] = injection[1:dim_UL,(dim_UL+(1:periodic_dimension))+dim_injection];
   gamma[1:periodic_dimension,(first_index+(1:dim_UL))+dim_gamma] = injection[(dim_UL+(1:periodic_dimension)),(1:dim_UL)+dim_injection];

   gamma[(1:periodic_dimension)+dim_gamma,(1:periodic_dimension)+dim_gamma] = injection[(dim_UL+(1:periodic_dimension))+dim_injection,(dim_UL+(1:periodic_dimension))+dim_injection];
   gamma[(first_index+(1:dim_UL))+dim_gamma,(1:periodic_dimension)+dim_gamma] = injection[(1:dim_UL)+dim_injection,(dim_UL+(1:periodic_dimension))+dim_injection];
   gamma[(1:periodic_dimension)+dim_gamma,(first_index+(1:dim_UL))+dim_gamma] = injection[(dim_UL+(1:periodic_dimension))+dim_injection,(1:dim_UL)+dim_injection];
 end

 return gamma
end

function Eigenvalues_of_rho(M)
   N = convert(Int64, size(M,1)/2);#H_ferm_f_Q, U_diag_f_Q

   evor = ones(Float64, 2^N)

   D,U = Diag_gamma((M+M')/2.)


   #eval = (1+sort(real(eigvals(M))))/2
   for iiter=1:2^N
       index = iiter-1;
       for jiter=1:N
           evor[iiter] = evor[iiter]*round((mod(index-1,2)*D[jiter,jiter]+(1-mod(index-1,2))*(D[jiter+N,jiter+N])),18)
           #evor[iiter] = evor[iiter]*round((mod(index-1,2)*eval[jiter]+(1-mod(index-1,2))*(1-eval[jiter])),15)
           index -= mod(index,2);
           index = index/2;
           # println("--> ", real(D[jiter,jiter]), " aa ", real(eval[jiter]))
           # println("--> ", real(D[jiter+N,jiter+N]), " aa ", real(eval[jiter+N]))
       end
   end

   return evor;
end

function VN_entropy(M)
   N = convert(Int64, size(M,1));

   D,U = eig((M+M')/2.);   #Se voglio fare con eig
   D = diagm(real(D));

   #Forzo la Hermitianità, se ho risultati strani controlla QUA
   #che M sia abbastanza vicina all'hermitiana prima.
   # M   = (M+M')/2.
   # D,U = Diag_gamma(M);
   S = 0;

   for iiter=1:N
       if (round(D[iiter,iiter],18)<-0.0000000000001)
        De,Ue = eig((M+M')/2.);
        for iiter=1:N
         println("DG: ", D[iiter,iiter]);
        end
        for iiter=1:N
         println("DE: ", De[iiter]);
        end
        save("/home/jacopo/Dropbox/ermatr.jld","M", M)
        error = string("Eigenvalue in VE not in [0,1]: ",round(D[iiter,iiter],18))
        throw(ArgumentError(error))
       end
       nu = abs(round(D[iiter,iiter],18))
       if (nu != 0 && nu != 1)
        #Invece di arrivare fino a N/2 nel ciclo e sommare nu e 1-nu, li passo tutti
        #perchè potrebbe essere che facendo una generica transformazione ortogonale non li abbia
        #ordinati in coppie, anche se Diag_gamma dovrebbe metterli in ordine
           S -= log(nu)*nu;
       end
   end

   return S;
end



function Project_diagonals(M4,off_diagonals)
 #Return a 4-blocks matrix, in wich in each block only
 #the first off_diagonals diagonal off diagonal are manteined
 #the rest is set to 0.
 #If off_diagonals=0 then it mantains only the diagonal of each block
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


function Is_phys(gamma)
 v = round.(eigvals(gamma),18);

 for iiter=1:size(v,1)
   if (real(v[iiter])>1 || real(v[iiter])<0)
     println("AV---> ", v[iiter])
     return false
   end
 end

 return true
end

function Physicalise(gamma)
 D,V = eig(gamma);

 for iiter=1:size(D,1)
   if (real(D[iiter])>1)
     D[iiter]=1
   elseif (real(D[iiter])<0)
     D[iiter]=0
   end
 end

 gamma = V*diagm(D)*V';

 return gamma
end



function Mutual_information(gamma, dim_a, start_a, dim_b, start_b)
 SA    = VN_entropy(Reduce_gamma(gamma, dim_a, start_a))
 SB    = VN_entropy(Reduce_gamma(gamma, dim_b, start_b))
 STot  = VN_entropy(gamma)

 return (SA+SB-STot)
end



##########Experimental Stuffs

function Fourier_matrix(N_f)
   F = zeros(Complex{Float64}, 2*N_f, 2*N_f)
   ND = (N_f-1)/2.;

   for iiter=1:(2*ND+1)
       index1=iiter-1-ND;
       iiter = convert(Int64, iiter);
       for jiter=1:(2*ND+1)
           index2=jiter-1-ND;
           jiter = convert(Int64, jiter);
           F[iiter,jiter]          = exp(-im*2*pi*index1*index2/(N_f))/sqrt(N_f)
           # F[iiter+N_f,jiter]      = exp(im*2*pi*index1*index2/(N_f))/sqrt(N_f)
           # F[iiter,jiter+N_f]      = exp(im*2*pi*index1*index2/(N_f))/sqrt(N_f)
           F[iiter+N_f,jiter+N_f]  = exp(im*2*pi*index1*index2/(N_f))/sqrt(N_f)
       end
   end

   return F;
end

function Bogoliubov_matrix(Jx,Jy,lambda,N_f)
   B = zeros(Complex{Float64}, 2*N_f, 2*N_f);

   theta = zeros(Float64, N_f);
   for kiter=1:N_f
       k = kiter-convert(Int64, (N_f-1)/2);
       theta[kiter] = 0.5*atan(-(((Jx-Jy)*sin(2*pi*k/N_f))/((Jx+Jy)*cos(2*pi*k/N_f)+lambda)));
   end

   for iiter=1:N_f
       B[iiter+N_f,2*N_f-iiter+1]  = -im*cos(theta[iiter]);
       B[iiter,iiter+N_f]          = sin(theta[iiter]);
       B[iiter+N_f,iiter]          = sin(theta[iiter]);
       B[iiter,N_f-iiter+1]        = im*cos(theta[iiter]);
   end

   return B;
end


function Build_effective_term_A(Jx, Jy, m)
   dimension = size(Jx, 1);
   M_A       = zeros(Float64, dimension, dimension)

   for iiter=1:dimension
       M_A[iiter, mod(iiter+m-1,dimension)+1]   = -(Jx[1]+Jy[1]);
       M_A[mod(iiter+m-1,dimension)+1, iiter] = -(Jx[1]+Jy[1]);
   end

   return M_A
end

function Build_effective_term_B(Jx, Jy, m)
   dimension = size(Jx, 1);
   M_B       = zeros(Float64, dimension, dimension)
   for iiter=1:dimension
       M_B[mod(iiter+m-1,dimension)+1, iiter] = -(Jx[1]-Jy[1]);
       M_B[iiter, mod(iiter+m-1,dimension)+1] = Jx[1]-Jy[1];
   end

   return M_B
end

function Effective_TFI_Hamiltonian(N,Jx,Jy,lambda, m)
   A = 0.5*Build_A(Jx,Jy,lambda)+0.5*Build_effective_term_A(Jx,Jy,m);
   B = 0.5*Build_B(Jx,Jy)+0.5*Build_effective_term_B(Jx,Jy,m);

   H_TFI                          = zeros(Float64, 2*N, 2*N);
   H_TFI[1:N_f,1:N_f]             = -A;
   H_TFI[(1:N_f)+N_f,1:N_f]       = -B;
   H_TFI[1:N_f,(1:N_f)+N_f]       =  B;
   H_TFI[(1:N_f)+N_f,(1:N_f)+N_f] =  A;

   return H_TFI;

end

function Build_thermal_gamma_fixed_beta(Diag_H, U_H, beta)
 N_f   = convert(Int64, size(Diag_H,1)/2.);

 gamma = zeros(Complex{Float64}, size(Diag_H,1),size(Diag_H,1))
 for kiter=1:N_f
   e_k = Diag_H[kiter+N_f,kiter+N_f];
   gamma[kiter,kiter] = 1/(1+exp(beta*e_k));
   gamma[kiter+N_f,kiter+N_f] = 1/(1+exp(-beta*e_k))
 end

 gamma = U_H*gamma*U_H';

 return gamma;
end


function Build_thermal_gamma(Diag_H, U_H, conserved_energy)
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
       temp_energy += e_k*(1/(1+exp(beta*e_k))-1/(1+exp(-beta*e_k)));
     end
     if (temp_energy > conserved_energy)
       a = beta;
       else
       b = beta;
     end
 end

 println("Thermal state with energy precision of:", abs(temp_energy-conserved_energy)," and beta: ", beta);

 gamma = zeros(Complex{Float64}, size(Diag_H,1),size(Diag_H,1))
 for kiter=1:N_f
   e_k = Diag_H[kiter+N_f,kiter+N_f];
   gamma[kiter,kiter] = 1/(1+exp(beta*e_k));
   gamma[kiter+N_f,kiter+N_f] = 1/(1+exp(-beta*e_k))
 end

 gamma = U_H*gamma*U_H';

 return gamma;
end

function E_contribute(U,i,k)
   N_f = convert(Int64, size(U,1)/2.)
   return real(0.5*(U[i,k]*conj(U[i,k])+U[i,k+N_f]*conj(U[i,k+N_f])+U[i+N_f,k]*conj(U[i+N_f,k])+U[i+N_f,k+N_f]*conj(U[i+N_f,k+N_f])));
end
#####

function exponentially_truncate(Γ,m)
 N_f = convert(Int64, size(Γ,1)/2);
 δs  = real(-log10(abs(Γ[1,m]))+log10(abs(Γ[1,(m-1)])));
 δd  = real(-log10(abs(Γ[N_f+1,m]))+log10(abs(Γ[N_f+1,(m-1)])));
 off_diagonals = convert(Int64,N_f/2);

 #Return a 4-blocks matrix, in wich in each block only
 #the first off_diagonals diagonal off diagonal are manteined
 #the rest is set to 0.
 #If off_diagonals=0 then it mantains only the diagonal of each block
 N_f = convert(Int64, size(Γ,1)/2.)
 M_finale = deepcopy(Γ);#zeros(Complex{Float64}, 2*N_f, 2*N_f)

 for iiter=1:N_f
   M_finale[iiter, iiter]          = Γ[iiter, iiter]
   M_finale[iiter+N_f, iiter]      = Γ[iiter+N_f, iiter]
   M_finale[iiter, iiter+N_f]      = Γ[iiter, iiter+N_f]
   M_finale[iiter+N_f, iiter+N_f]  = Γ[iiter+N_f, iiter+N_f]
   for jiter=(m-1):off_diagonals
     M_finale[iiter, mod(iiter+jiter-1,N_f)+1]          = min(1, round(10^(-δs*(jiter-(m-2))),12))*Γ[iiter, mod(iiter+jiter-2,N_f)+1]
     M_finale[iiter+N_f, mod(iiter+jiter-1,N_f)+1]      = min(1, round(10^(-δd*(jiter-(m-2))),12))*Γ[iiter+N_f, mod(iiter+jiter-2,N_f)+1]
     M_finale[iiter, mod(iiter+jiter-1,N_f)+1+N_f]      = min(1, round(10^(-δd*(jiter-(m-2))),12))*Γ[iiter, mod(iiter+jiter-2,N_f)+1+N_f]
     M_finale[iiter+N_f, mod(iiter+jiter-1,N_f)+1+N_f]  = min(1, round(10^(-δs*(jiter-(m-2))),12))*Γ[iiter+N_f, mod(iiter+jiter-2,N_f)+1+N_f]

     M_finale[iiter, mod(iiter-jiter-1,N_f)+1]          = min(1, round(10^(-δs*(jiter-(m-2))),12))*Γ[iiter, mod(iiter-jiter,N_f)+1]
     M_finale[iiter+N_f, mod(iiter-jiter-1,N_f)+1]      = min(1, round(10^(-δd*(jiter-(m-2))),12))*Γ[iiter+N_f, mod(iiter-jiter,N_f)+1]
     M_finale[iiter, mod(iiter-jiter-1,N_f)+1+N_f]      = min(1, round(10^(-δd*(jiter-(m-2))),12))*Γ[iiter, mod(iiter-jiter,N_f)+1+N_f]
     M_finale[iiter+N_f, mod(iiter-jiter-1,N_f)+1+N_f]  = min(1, round(10^(-δs*(jiter-(m-2))),12))*Γ[iiter+N_f, mod(iiter-jiter,N_f)+1+N_f]
   end
 end
 return M_finale

 # println("zz ", δs, "   ", round(10^(-δs*(1)),12))
 #
 # for iiter=1:N_f
 #   for jiter=3:convert(Int64, N_f/2)
 #     Γ[iiter,jiter] = round(10^(-δs*(jiter-2)),12)*Γ[iiter,jiter];
 #     Γ[iiter,N_f-jiter] = round(10^(-δs*(jiter-2)),12)*Γ[iiter,N_f-jiter];
 #     Γ[iiter+N_f,jiter+N_f] = round(10^(-δs*(jiter-2)),12)*Γ[iiter+N_f,jiter+N_f];
 #     Γ[iiter+N_f,2*N_f-jiter] = round(10^(-δs*(jiter-2)),12)*Γ[iiter+N_f,2*N_f-jiter];
 #
 #     # Γ[iiter,jiter+N_f] = round(10^(-δd*(jiter-2)),12)*Γ[iiter,jiter+N_f];
 #     # Γ[iiter,2*N_f-jiter] = round(10^(-δd*(jiter-2)),12)*Γ[iiter,2*N_f-jiter];
 #     # Γ[iiter+N_f,jiter] = round(10^(-δs*(jiter-2)),12)*Γ[iiter+N_f,jiter];
 #     # Γ[iiter+N_f,N_f-jiter] = round(10^(-δd*(jiter-2)),12)*Γ[iiter+N_f,N_f-jiter];
 #   end
 # end
 #
 # return Γ;
end


function exchange(U,column1,column2,element,site)
 element = element-1
 N_red = convert(Int64, size(U,1)/2.);

 ν1 = deepcopy(U[:,column1])
 c1 = (U[site+element*N_red,column1])
 ν2 = deepcopy(U[:,column2])
 c2 = (U[site+element*N_red,column2])
 norma = norm((1/c1)*ν1-(1/c2)*ν2);
 d1 = (1/c1)*(1/(norma));
 d2 = (1/c2)*(1/(norma));
 U[:,column1]   = d1*ν1-d2*ν2;
 U[:,column2]  = (conj(d2)*ν1+conj(d1)*ν2);

 return U
end

function zero_border(U, e)
 N_red = convert(Int64, size(U,1)/2)

 #Cancello bordo SX
 U = exchange(U,e,e+1,1,1);
 U = exchange(U,e+2,e+3,1,1);
 U = exchange(U,e,e+2,2,1);

 #AGGIUNGO cancellazione bordo dx
 U = exchange(U,e+1,e+3,1,1);
 U = exchange(U,e+1,e+2,2,1);
 U = exchange(U,e,e+1,1,N_red)
 U = exchange(U,e+3,e+4,1,1)
 U = exchange(U,e+2,e+3,2,1)
 U = exchange(U,e+1,e+2,1,N_red)
 U = exchange(U,e,e+1,2,N_red)

 ##Per le altre colonne col +N_red
 ######EX parziale senza ridefinizione malefica di e
 # U = exchange(U,e+N_red,e+1+N_red,1,1);
 # U = exchange(U,e+2+N_red,e+3+N_red,1,1);
 # U = exchange(U,e+N_red,e+2+N_red,2,1);
 ######
 e=e+N_red
 #Cancello bordo SX
 U = exchange(U,e,e+1,1,1);
 U = exchange(U,e+2,e+3,1,1);
 U = exchange(U,e,e+2,2,1);

 #AGGIUNGO cancellazione bordo dx
 U = exchange(U,e+1,e+3,1,1);
 U = exchange(U,e+1,e+2,2,1);
 U = exchange(U,e,e+1,1,N_red)
 U = exchange(U,e+3,e+4,1,1)
 U = exchange(U,e+2,e+3,2,1)
 U = exchange(U,e+1,e+2,1,N_red)
 U = exchange(U,e,e+1,2,N_red)

 return U
end

function decycle(M)
 dim = convert(Int64, size(M,1))

 for iiter=0:(dim-1)
   λ = M[1,iiter+1];
   for jiter=1:(dim-iiter)
     M[jiter,iiter+jiter] = M[jiter,iiter+jiter,]-λ
     # println("[",jiter,",",iiter+jiter,"]-[1,",iiter+1,"]");
   end
 end

 for iiter=1:(dim-1)
   λ = M[iiter+1,1]
   for jiter=1:(dim-iiter)
     M[jiter+iiter,jiter] = M[jiter+iiter,jiter]-λ
     # println(">>[",jiter+iiter,",",jiter,"]-[",iiter+1,",1]");
   end
 end

 return M
end

function Momentum_sector_H(L1,L2,γ,μ,k)
 #Ricordati che devi aggiungere il termine k=0 per l'hamiltoniana L1 dispari
 #Quaderno 5-18 p.7
 if (mod(L1,2)==1)
   println("_Momentum_sector_H_");
   println("ATTENZIONE L'HAMILTONIANA è FATTA PER I PARI");
   println("Per i dispari devi aggiungere un termine");
 end

 H = zeros(Complex{Float64},4*L2,4*L2);

 γ_k = γ*sin(((2*pi)/L1*k));
 ϵ_k = (μ-2*cos((2*pi)/L1*k));

 h_γ = zeros(Complex{Float64}, L2, L2);
 h_ϵ = zeros(Complex{Float64}, L2, L2);

 for iiter=1:(L2-1)
   h_γ[iiter,iiter]    = -im*2*γ_k;
   h_γ[iiter,iiter+1] = -γ;
   h_γ[iiter+1, iiter] = γ;

   h_ϵ[iiter,iiter]    = ϵ_k;
   h_ϵ[iiter,iiter+1]  = 1;
   h_ϵ[iiter+1,iiter]  = 1;
 end
 h_γ[L2,L2] = -im*2*γ_k;
 h_ϵ[L2,L2] = ϵ_k;

 A = zeros(Complex{Float64}, 2*L2, 2*L2);
 B = zeros(Complex{Float64}, 2*L2, 2*L2);

 A[1:L2,1:L2] = h_ϵ;
 A[(L2+1):(2*L2),(L2+1):(2*L2)] = h_ϵ;

 B[1:L2,(L2+1):(2*L2)] = h_γ;
 B[(L2+1):(2*L2),1:L2] = -transpose(h_γ);

 H[1:(2*L2),1:(2*L2)]                = -conj(A);
 H[(2*L2+1):(4*L2),(2*L2+1):(4*L2)]  = A;
 H[1:(2*L2),(2*L2+1):(4*L2)]         = B;
 H[(2*L2+1):(4*L2),1:(2*L2)]         = -conj(B);

 H = (1/2.)*(H/2. +L2*ϵ_k);

 return H;
end

function diagonalise_block(γi, starting_site, dimension)
 #Questa funzione prende una matrice γ e ritorna
 #la matrice γ_finale con il blocco grande "dimension" che parte dal
 #sito "starting_site" diagonalizzato.
 #Inoltre ritorna U_tot tale che
 #γ =  U_tot*γ_finale*U_tot'

 γ = deepcopy(γi)
 γ = (γ+γ')/2.
 dimension = convert(Int64, dimension)
 if ((starting_site+dimension-1)>size(γ,1)/2)
   println("_diagonalise_block_");
   println("ATTENZIONE IL BLOCCO SFORA IL BORDO");
   println("SS ", starting_site, " d ", dimension, " size(γ,1) ", size(γ,1));
 end
 γ_reduced    = Reduce_gamma(γ,dimension,starting_site);
 D_red, U_red = Diag_gamma(γ_reduced);
 ident        = convert(Array{Complex{Float64},2}, eye(size(γ,1)));
 U_tot        = Inject_gamma(ident, U_red, starting_site);
 γ_finale     = U_tot'*γ*U_tot;
 γ_finale     = (γ_finale+γ_finale')/2.



 return γ_finale,U_tot
end

function reduce_bond_dimension(Γ_in,L1,L2,m,μ)
 #Questa funzione prende un vettore di "N_momenti" Gamme e lo ritorna
 #come un vettore di "N_momenti" gamme delle quali il numero di
 #modi totale è "m".
 Γ_init = deepcopy(Γ_in)
 N_momenti = convert(Int64, size(Γ_init,1));
 #Metto in Γ_it la matrice già diagonalizzata nel primo blocchetto
 #così da riempire la tabella v ed avere le Γ di partenza
 Γ   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U_Γ   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 v   = zeros(Float64, 3, N_momenti);
 v[2,:] = 1;
 for iiter=1:N_momenti
     Γ[iiter,:,:], U_Γ[iiter,:,:]  = diagonalise_block(Γ_init[iiter,:,:],1, μ+1);
     v[1,iiter]                    = real(Γ[iiter,1,1]);
 end
 Tot_modi = L1*L2;
 while(Tot_modi>m)
    k = 0;
    check = true;
    perm = sortperm(v[1,:]);
    pos = 1;
    #Scelgo la piu piccola matrice da cui non ho già tolto tutti i modi
    while (check)
      #PRIMA C'ERA UN +2 come sul quaderno, controlla che vada bene toglierlo in
      #questo caso
        k = convert(Int64, perm[pos]);
        if ((convert(Int64,v[3,k]))==2*L2)
            pos +=1;
        else
            check = false;
        end
    end
    #Cancello le correlazioni e porto a zero l'entropia del modo
    v3k = convert(Int64, v[3,k])
     #println("k: ", k, "v1:", v[1,k], " Γ11: ", Γ[k,v3k+1,v3k+1], " 1-ν: ", Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)])
    Γ[k,v3k+1,:]                  = 0;
    Γ[k,:,v3k+1]                  = 0;
    Γ[k,v3k+1+(2*L2),:]             = 0;
    Γ[k,:, v3k+1+(2*L2)]            = 0;
    Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)]    = 1;

    v[3,k] += 1; #->In questa partizione ne ho tolto uno in più

    #Controllo che non sia già arrivato a diagonalizzare l'ultimo blocco
    if (v[2,k]<(2*L2-μ))
        v[2,k] += 1;
    end
    Γ[k,:,:], U  = diagonalise_block(Γ[k,:,:],convert(Int64, v[2,k]),(μ+1));
    U_Γ[k,:,:]   = U_Γ[k,:,:]*U;
    v3k = convert(Int64, v[3,k])
    v[1,k]       = real(Γ[k,v3k+1,v3k+1]); #Forse +2?
    Tot_modi -= 1;
 end
 for kiter=1:N_momenti
     Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 return Γ;
end


function trunker(Γ,m)
 #Questa tronca lasciando m autoval
 #m deve essere minore di N/2
 N_f = convert(Int64, size(Γ,1)/2.)
 contribute = 0;

 for i=(m+1):(N_f-(m+1))
     Γ, U_temp1 = diagonalise_block(Γ,1,i)
     Γ,U_temp2 = diagonalise_block(Γ ,i+1,N_f-i)
     ind_a = max(1,i-convert(Int64, N_f/2.)+1)

     contribute += Γ[ind_a, ind_a];

     Γ[ind_a,:]      = 0
     Γ[:,ind_a]      = 0
     Γ[ind_a+N_f,:]  = 0
     Γ[:, ind_a+N_f] = 0
     Γ[ind_a+N_f,ind_a+N_f]  = 1


     ind_b = i+max(1, N_f-i+1-convert(Int64, N_f/2.))
     Γ[ind_b,:]      = 0
     Γ[:,ind_b]      = 0
     Γ[ind_b+N_f,:]  = 0
     Γ[:, ind_b+N_f] = 0
     Γ[ind_b+N_f,ind_b+N_f]  = 1

     Γ = U_temp1*U_temp2*Γ*(U_temp2')*(U_temp1')
 end

 return Γ, contribute
end


function turbo_trunker(Γ,m)
 #Questa tronca lasciando m autoval
 #m deve essere minore di N/2
 N_f = convert(Int64, size(Γ,1)/2.)
 contribute = 0;

 i=(m+1)
 Γ, U_temp1 = diagonalise_block(Γ,1,i)
 Γ,U_temp2 = diagonalise_block(Γ ,i+1,N_f-i)
 ind_a = max(1,i-convert(Int64, N_f/2.)+1)
 contribute += Γ[ind_a, ind_a];
 Γ[ind_a,:]      = 0
 Γ[:,ind_a]      = 0
 Γ[ind_a+N_f,:]  = 0
 Γ[:, ind_a+N_f] = 0
 Γ[ind_a+N_f,ind_a+N_f]  = 1


 ind_b = i+max(1, N_f-i+1-convert(Int64, N_f/2.))
 Γ[ind_b,:]      = 0
 Γ[:,ind_b]      = 0
 Γ[ind_b+N_f,:]  = 0
 Γ[:, ind_b+N_f] = 0
 Γ[ind_b+N_f,ind_b+N_f]  = 1

 Γ = U_temp1*U_temp2*Γ*(U_temp2')*(U_temp1')

 i = (N_f-(m+1))
 Γ, U_temp1 = diagonalise_block(Γ,1,i)
 Γ,U_temp2 = diagonalise_block(Γ ,i+1,N_f-i)
 ind_a = max(1,i-convert(Int64, N_f/2.)+1)
 contribute += Γ[ind_a, ind_a];
 Γ[ind_a,:]      = 0
 Γ[:,ind_a]      = 0
 Γ[ind_a+N_f,:]  = 0
 Γ[:, ind_a+N_f] = 0
 Γ[ind_a+N_f,ind_a+N_f]  = 1


 ind_b = i+max(1, N_f-i+1-convert(Int64, N_f/2.))
 Γ[ind_b,:]      = 0
 Γ[:,ind_b]      = 0
 Γ[ind_b+N_f,:]  = 0
 Γ[:, ind_b+N_f] = 0
 Γ[ind_b+N_f,ind_b+N_f]  = 1

 Γ = U_temp1*U_temp2*Γ*(U_temp2')*(U_temp1')

 return Γ, contribute
end


function reduce_bond_dimension_2(Γ,L1,L2,m)
 N_momenti = convert(Int64, size(Γ,1));

 Γ_new   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U_Γ   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 v   = zeros(Float64, 2, N_momenti);

 for k=1:N_momenti
   Γ_new[k,:,:] = deepcopy(Γ[k,:,:])
   Γ_new[k,:,:], v[1,k] = turbo_trunker(Γ_new[k,:,:], convert(Int64, L2-v[2,k]-1))
 end

 Tot_conserved_mode = (L1*L2)/2.
 while (Tot_conserved_mode>m)
   #println("...",v[1,:],"...")
   perm = sortperm(v[1,:]);
   check = true;
   k = 1;
   pos = 1;
   while (check)
       k = convert(Int64, perm[pos]);
       if ((convert(Int64,v[2,k]))==(L2-1))
           pos +=1;
       else
           check = false;
       end
   end
   Γ[k,:,:] = Γ_new[k,:,:]
   v[2,k] += 1;
   Tot_conserved_mode -= 1;
   #if (Tot_conserved_mode>m)
     Γ_new[k,:,:] = deepcopy(Γ[k,:,:])
     Γ_new[k,:,:], v[1,k] = turbo_trunker(Γ_new[k,:,:], convert(Int64, L2-v[2,k]-1))
   #end

   #println("----> ", Tot_conserved_mode, " -- ", k, " --- ", v[1,:])
 end

 return Γ
end



function reduce_bond_dimension_3(Γ_in,L1,L2,m)
 #Questa funzione prende un vettore di "N_momenti" Gamme e lo ritorna
 #come un vettore di "N_momenti" gamme delle quali il numero di
 #modi totale è "m".
 Γ_init = deepcopy(Γ_in)
 N_momenti = convert(Int64, size(Γ_init,1));
 O = build_riordina_settore(L2)
 for iiter=1:N_momenti
   Γ_init[iiter,:,:] = O*Γ_init[iiter,:,:]*(O');
 end

 #Quanto fare grande il blocchetto in ogni gamma -2
 μ = convert(Int64, 2*(m/L2))
 #Metto in Γ_it la matrice già diagonalizzata nel primo blocchetto
 #così da riempire la tabella v ed avere le Γ di partenza
 Γ   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U_Γ   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 v   = zeros(Float64, 3, N_momenti);
 v[2,:] = 1;
 for iiter=1:N_momenti
     Γ[iiter,:,:], U_Γ[iiter,:,:]  = diagonalise_block(Γ_init[iiter,:,:],1, μ+2);
     v[1,iiter]                    = real(Γ[iiter,1,1])+real(Γ[iiter,2,2]);
 end
 Tot_modi = L1*L2;
 while(Tot_modi>m)
    k = 0;
    check = true;
    perm = sortperm(v[1,:]);
    pos = 1;
    #Scelgo la piu piccola matrice da cui non ho già tolto tutti i modi
    while (check)
      #PRIMA C'ERA UN +2 come sul quaderno, controlla che vada bene toglierlo in
      #questo caso
        k = convert(Int64, perm[pos]);
        if ((convert(Int64,v[3,k]))==2*L2)
            pos +=1;
        else
            check = false;
        end
    end
    #Cancello le correlazioni e porto a zero l'entropia del modo
    v3k = convert(Int64, v[3,k])
     #println("k: ", k, "v1:", v[1,k], " Γ11: ", Γ[k,v3k+1,v3k+1], " 1-ν: ", Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)])
    Γ[k,v3k+1,:]                  = 0;
    Γ[k,:,v3k+1]                  = 0;
    Γ[k,v3k+1+(2*L2),:]             = 0;
    Γ[k,:, v3k+1+(2*L2)]            = 0;
    Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)]    = 1;
    Γ[k,v3k+2,:]                  = 0;
    Γ[k,:,v3k+2]                  = 0;
    Γ[k,v3k+2+(2*L2),:]             = 0;
    Γ[k,:, v3k+2+(2*L2)]            = 0;
    Γ[k,v3k+2+(2*L2),v3k+2+(2*L2)]    = 1;

    v[3,k] += 2; #->In questa partizione ne ho tolto uno in più

    #rollo che non sia già arrivato a diagonalizzare l'ultimo blocco
    if (v[2,k]<(2*L2-μ-1))
        v[2,k] += 2;
    end
    Γ[k,:,:], U  = diagonalise_block(Γ[k,:,:],convert(Int64, v[2,k]),(μ+2));
    U_Γ[k,:,:]   = U_Γ[k,:,:]*U;
    v3k = convert(Int64, v[3,k])
    v[1,k]       = real(Γ[k,v3k+1,v3k+1])+real(Γ[k,v3k+2,v3k+2])
    Tot_modi -= 2;
 end
 for kiter=1:N_momenti
     Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 for iiter=1:N_momenti
   Γ[iiter,:,:] = (O')*Γ[iiter,:,:]*O;
 end

 return Γ;
end

function build_riordina_settore(L2)
 #Questa O mi trasforma da a†k1,a†k2,..,a†-k1,a†-k2,...,ak1,...,a-k1,....
 #all ordine a†k1,a†-k1,a†k2,a†-k2,.....ak1,a-k1,..
 O = zeros(Int64, 4*L2, 4*L2)

 for iiter=1:(2*L2)
     O[iiter,convert(Int64, ceil(iiter/2)+mod(iiter-1,2)*L2)] = 1;
 end
 for iiter=1:(2*L2)
     O[iiter+(2*L2),convert(Int64, ceil(iiter/2)+mod(iiter-1,2)*L2)+(2*L2)] = 1;
 end

 return O;
end


function build_riordina_settore_2(L2)
 #Questa O mi trasforma da a†k1,a†k2,..,a†-k1,a†-k2,...,ak1,...,a-k1,....
 #all ordine a†k1,a†-k1,a†k2,a†-k2,.....ak1,a-k1,..
 O = zeros(Int64, 4*L2, 4*L2)

 for iiter=1:(2*L2)
     O[iiter,convert(Int64, ceil(iiter/2)+mod(iiter-1,2)*L2)] = 1;
 end
 for iiter=1:(2*L2)
     O[iiter+(2*L2),convert(Int64, ceil(iiter/2)+mod(iiter-1,2)*L2)+(2*L2)] = 1;
 end

 return O;
end




function Complex_eye(N)
 return convert(Array{Complex{Float64},2}, eye(N));
end


function Purification(Γ_mixed)
 N = convert(Int64, size(Γ_mixed,1)/2.);

 Γ_pure_diag = zeros(Complex{Float64}, 4*N, 4*N);
 Γ_pure      = zeros(Complex{Float64}, 4*N, 4*N);
 U_pure      = Complex_eye(4*N);

 #Diagonalise the mixed Γ
 Γ_mixed_diag, U_mixed = Diag_ferm(Γ_mixed-0.5*eye(2*N));
 Γ_mixed_diag          = U_mixed'*Γ_mixed*U_mixed;

 #U_pure move Γ_pure_diag to the space basis
 U_pure = Inject_gamma(U_pure, U_mixed,1);
 U_pure = Inject_gamma(U_pure, U_mixed,N+1);

 Γ_pure_diag = Inject_gamma(Γ_pure_diag,Γ_mixed_diag,1);
 Γ_pure_diag = Inject_gamma(Γ_pure_diag,Γ_mixed_diag,N+1);
 for iiter=1:N
   ν   = real(Γ_mixed_diag[iiter,iiter]);
   Γ_pure_diag[iiter,iiter+3*N]    =   sqrt(ν*(1-ν));
   Γ_pure_diag[iiter+N,iiter+2*N]  =  -sqrt(ν*(1-ν));
   Γ_pure_diag[iiter+2*N,iiter+N]  =  -sqrt(ν*(1-ν));
   Γ_pure_diag[iiter+3*N,iiter]    =   sqrt(ν*(1-ν));
 end

 Γ_pure = U_pure*Γ_pure_diag*(U_pure');

 return Γ_pure;
end

function two_sites(Γ,s1,s2)
 N = convert(Int64,size(Γ,1)/2.);

 γ = zeros(Complex{Float64},4,4);
 γ[1,1] = Γ[s1,s1];
 γ[2,2] = Γ[s2,s2];
 γ[1,2] = Γ[s1,s2];
 γ[2,1] = Γ[s2,s1];

 γ[1,1+2] = Γ[s1,s1+N];
 γ[2,2+2] = Γ[s2,s2+N];
 γ[1,2+2] = Γ[s1,s2+N];
 γ[2,1+2] = Γ[s2,s1+N];

 γ[1+2,1] = Γ[s1+N,s1];
 γ[2+2,2] = Γ[s2+N,s2];
 γ[1+2,2] = Γ[s1+N,s2];
 γ[2+2,1] = Γ[s2+N,s1];

 γ[1+2,1+2] = Γ[s1+N,s1+N];
 γ[2+2,2+2] = Γ[s2+N,s2+N];
 γ[1+2,2+2] = Γ[s1+N,s2+N];
 γ[2+2,1+2] = Γ[s2+N,s1+N];

 return γ;
end


function negativita(gamma)
 D = eigvals(gamma);

 somma = 0;

 for iiter=1:size(D,1)
   if (real(D[iiter])<0)
     somma += abs(D[iiter]);
   end
 end

 return somma
end

function Build_ρ2x2(γ)
 I2  = eye(2);
 σp  = [0 1; 0 0];
 σm  = [0 0; 1 0];
 σz  = [1 0; 0 -1];
 βv  = zeros(Complex{Float64}, 4, 4, 4)
 βv[3,:,:] = -kron(σm,I2);
 βv[4,:,:] = -kron(σz,σm);
 βv[1,:,:] = -kron(σp,I2);
 βv[2,:,:] = -kron(σz,σp);

 D,U = Diag_gamma(γ);
 Hd   = zeros(Complex{Float64}, 4,4)
 for iiter=1:4
   Hd[iiter,iiter] = log((1-D[iiter,iiter])/D[iiter,iiter])
 end
 H = U*Hd*U';

 αv  = zeros(Complex{Float64}, 4, 4, 4)
 αh  = zeros(Complex{Float64}, 4, 4, 4)
 for iiter=1:4
   for jiter=1:4
     αv[iiter,:,:] += U[iiter,jiter]*βv[jiter,:,:]
   end
 end
 αh[1,:,:] = αv[3,:,:]
 αh[2,:,:] = αv[4,:,:]
 αh[3,:,:] = αv[1,:,:]
 αh[4,:,:] = αv[2,:,:]

 τv  = zeros(Complex{Float64}, 4, 4, 4)
 for iiter=1:4
   for jiter=1:4
     τv[iiter,:,:] += H[iiter,jiter]*αv[jiter,:,:]
   end
 end

 Me = zeros(Complex{Float64}, 4, 4)
 for iiter=1:4
     Me += αh[iiter,:,:]*τv[iiter,:,:]
 end
 Me = -0.5*Me;

 Z = trace(expm(Me))
 ρ = expm(Me)/Z

 return ρ
end

function Unitary_composite(λ)
   #From https://arxiv.org/pdf/1103.3408.pdf
   #Le lambde sono ordinate come λ_{m,n}=λ[(m-1)*d+n];
   #Il numero di parametri e' d^2
   #La unitaria ha parametri {0,0,0,0,0,0,...,0}
   d   = convert(Int64, sqrt(size(λ,1)));
   size_y  = convert(Int64, (d*(d-1))/2.);

   onb     = eye(d)                #Orthonormal basis
   ULHS    = convert(Array{Complex{Float64},2}, eye(d));
   URHS    = convert(Array{Complex{Float64},2}, eye(d));
   for miter=1:d
       for niter=(miter+1):d
           P_exponent  = expm(im*(onb[niter,:]*onb[niter,:]')*λ[(niter-1)*d+miter]);
           Y_matrix    = -im*(onb[miter,:]*onb[niter,:]')+im*(onb[niter,:]*onb[miter,:]');
           Y_exponent  = expm(im*Y_matrix*λ[(miter-1)*d+niter]);
           ULHS        = ULHS*P_exponent*Y_exponent;
       end
       URHS    = URHS*expm(im*(onb[miter,:]*onb[miter,:]')*λ[(miter-1)*d+miter]);
   end

   return (ULHS*URHS)
end



function Periodic_rotation(θ,p,n,N)
 #Costruisce una matrice NxN con una rotazione tra il sito p ed il sito p+n
 #n = numero di salti
 #p = posizione iniziale
 r_nxn      = eye(N);
 r_nxn[mod1(n+p,N),mod1(n+p,N)] = cos(θ);
 r_nxn[p,p] = cos(θ);
 r_nxn[mod1(n+p,N),p] = -sin(θ);
 r_nxn[p,mod1(n+p,N)] = sin(θ);

 return r_nxn
end


function Orthogonal_composite(ϕ);
 #Costruisce una matrice ortogonale parametrizzata dal vettore degli
 #angoli delle rotazioni ϕ. La matrice è NxN dove N è calcolato dal numero di
 #parametri che ho passato alla funzione

 #Da Wikipedia, una ortogonale NxN ha in genere N(N-1)/2 parametri liberi
 #θ è il vettore di questi parametri liberi

#Questa è semplicemente la formula inversa per trovare N a partire dal numero di
#di parametri size(θ,1)
N = convert(Int64,0.5*(sqrt(4*2*size(ϕ,1)+1)+1));
O = eye(N);

index = 0;
for p=1:(N-1)
 for n=1:(N-p)
   index+=1;
   O = O*Periodic_rotation(ϕ[index],p,n,N);
  end
 end

 return O;
end





















function Maximally_decohere_2_sites(Γ)
   U_tot = convert(Array{Complex{Float64},2}, eye(2*N_f));

   check = true
   iiter = 1;
   while (check)
       backup = Coherence(Reduce_gamma(Γ,2,1));
       U_semi_tot, Γ  = Minimise_coherence_4_TI(Γ)
       Δ = abs(Coherence(Reduce_gamma(Γ,2,1))-backup);
       if (Δ<10. ^(-14.) || iiter>100)
           check = false;
       end
       U_tot = U_semi_tot*U_tot;
       iiter += 1;
       if(mod(iiter,10)==0)
         println("----> substep: ", iiter, " Δ: ", Δ, "Tot Coh: ", Sum_2_site_coherences(Γ))
       end
   end

   return U_tot, (Γ+Γ')/2.
end

function Minimise_coherence_4_TI(Γ)
   n   = convert(Int64, size(Γ,1)/2.); #Questo deve essere un multiplo di 4

   Ω = Build_Omega(div(size(Γ,1),2))
   UA = convert(Array{Complex{Float64},2}, eye(2*n));
   #A
   res = optimize(λ -> f(λ,Reduce_gamma(Γ,4,4)), zeros(6));
   λ   = Optim.minimizer(res);
   U2   = Orthogonal_composite(λ);
   back = S1S2S3(Reduce_gamma(Γ,4,4))
   for iiter=1:(convert(Int64, n/4.))
       UA = Ω'*Inject_gamma(convert(Array{Complex{Float64},2}, eye(2*n)), U2, (iiter-1)*4+1)*Ω*UA;
   end
   Γ   = UA*Γ*UA';
   Γ   = (Γ+Γ')/2.;

   UB = convert(Array{Complex{Float64},2}, eye(2*n));
   #B
   res = optimize(λ -> f(λ,Reduce_gamma(Γ,4,2)), zeros(6));
   λ   = Optim.minimizer(res);
   U2   = Orthogonal_composite(λ);
   back = S1S2S3(Reduce_gamma(Γ,4,2))
   for iiter=1:(convert(Int64, n/4.))
       UB = Ω'*Inject_gamma(convert(Array{Complex{Float64},2}, eye(2*n)), U2, (iiter-1)*4+3)*Ω*UB;
   end
   Γ   = UB*Γ*UB';LinAlg.LAPACKException(10)
   Γ   = (Γ+Γ')/2.;

   UC = convert(Array{Complex{Float64},2}, eye(2*n));
   #C
   res = optimize(λ -> f(λ,Reduce_gamma(Γ,4,1)), zeros(6));
   λ   = Optim.minimizer(res);
   U2   = Orthogonal_composite(λ);
   back = S1S2S3(Reduce_gamma(Γ,4,1))
   for iiter=1:(convert(Int64, n/4.))
       UC = Ω'*Inject_gamma(convert(Array{Complex{Float64},2}, eye(2*n)), U2, (iiter-1)*4+2)*Ω*UC;
   end
   Γ = UC*Γ*UC';
   Γ   = (Γ+Γ')/2.;

   UD = convert(Array{Complex{Float64},2}, eye(2*n));
   # #D
   back = S1S2S3(Reduce_gamma(Γ,4,3))
   res = optimize(λ -> f(λ,Reduce_gamma(Γ,4,3)), zeros(6));
   λ   = Optim.minimizer(res);
   U2   = Orthogonal_composite(λ);
   for iiter=1:(convert(Int64, n/4.))
       UD = Ω'*Inject_gamma(convert(Array{Complex{Float64},2}, eye(2*n)), U2, iiter*4)*Ω*UD;
   end
   Γ   = UD*Γ*UD';
   Γ   = (Γ+Γ')/2.;

   return (UD*UC*UB*UA),(Γ+Γ')/2.;
end

function Coherence(γ)
   return (VN_entropy(Project_just_diagonal(γ))-VN_entropy(γ))
end

function Sum_2_site_coherences(Γ)
   n   = convert(Int64, size(Γ,1)/2.);

   sum_coh = 0;
   for iiter=1:n
       sum_coh += Coherence(Reduce_gamma(Γ,2,iiter))
   end

   return sum_coh;
end

function S1S2S3(γ)
   n = size(γ,1)

   S1 = Coherence(Reduce_gamma(γ,2,1))
   S2 = Coherence(Reduce_gamma(γ,2,2))
   S3 = Coherence(Reduce_gamma(γ,2,3))

   return (S1+S2+S3)
end


function f(ϕ,γ)
   Ω   = Build_Omega(4)
   U2  = Orthogonal_composite(ϕ);
   U   = convert(Array{Complex{Float64},2}, eye(8));
   U   = Inject_gamma(U,U2,2);
   #La generica unitaria che posso fare sulla symbol matrix è una generica ortogonale
   #tra due gamma, cioè posso parametrizzare la generica sui majorana ed aggiungere la
   #majoranizzazione e la demajoranizzazione
   U   = Ω'*U*Ω;
   γ   = U*γ*U';
   Coh = S1S2S3(γ);
   γ   = U'*γ*U;

   return Coh
end


function Project_just_diagonal(M_in)
 N = size(M_in,1)

 M_out = zeros(typeof(M_in[1,1]),N,N)
 for iiter=1:N
  M_out[iiter,iiter] = M_in[iiter,iiter];
 end

 return M_out;
end

function unentangle_2_sites(Γ,s,d)
 N = convert(Int64, size(Γ,1)/2.);

 #I due blocchi dovrebbero essere uguali (per una pseudo T.I. totale symbol matrix)
 #Qua faccio lo swipe da sinistra a destra, ma potrei fare anche due layer
 #shiftati se vedo che ho problemi con le correlazioni lunghe (non ho pensato se
 #c'entra qualcosa)
 γsx = Reduce_gamma(Γ,2,s);
 γdx = Reduce_gamma(Γ,2,mod1(s+4+d,N));

 dsx,usx = Diag_gamma(γsx);
 ddx,udx = Diag_gamma(γdx);
 #Queste due dovrebbero essere uguali se lo stato è TI, e con uno SWIPE da sx a dx
 #dovrebbero essere uguali anche ad ogni turno

 U = convert(Array{Complex{Float64},2}, eye(2*N));
 U = Inject_gamma(U,usx,s);
 U = Inject_gamma(U,udx,mod1(s+2,N));

 Γ = U'*Γ*U;


 v = zeros(Int64,2)
 v[1] = mod1(s,N)
 v[2] = mod1(s+1,N)
#Un blocco a destra della prima delle due matricette
 Γ[v,mod1(s+4+d,N)] = 0;
 Γ[v,mod1(s+5+d,N)] = 0;

#Un blocco sotto della prima delle due matricette
 Γ[mod1(s+4+d,N),v] = 0;
 Γ[mod1(s+5+d,N),v] = 0;

#Torno alla base di prima
 Γ = U*Γ*U';

return Γ;
end

function swipe_untentangle_2_sites(Γ,d)
 N = convert(Int64, size(Γ,1)/2.);

 for s=1:N
  Γ = unentangle_2_sites(Γ,s,d)
 end

 return Γ;
end

function Majo_Greplova(Λ)
 N = convert(Int64, size(Γ,1)/2.);
 Ω = Build_Omega(N);

 return im*(Ω*(2*Λ-1)*(Omega'))
end



function Exchange_fermions(Λ,f1,f2)
  N = convert(Int64, size(Λ,1)/2.);

  # Γ = deepcopy(Λ);
  # temp      = Γ[f1,:];
  # Γ[f1,:]   = Γ[f2,:];
  # Γ[f2,:]   = temp;
  # Print_matrix("1", abs(Γ)[1:N,1:N])
  #
  # temp      = Γ[f1+N,:];
  # Γ[f1+N,:] = Γ[f2+N,:];
  # Γ[f2+N,:] = temp;
  #
  # temp      = Γ[:,f1];
  # Γ[:,f1]   = Γ[:,f2];
  # Γ[:,f2]   = temp;
  # Print_matrix("3", abs(Γ)[1:N,1:N])
  #
  # temp      = Γ[:,f1+N];
  # Γ[:,f1+N] = Γ[:,f2+N];
  # Γ[:,f2+N] = temp;

  U = eye(2*N)
  U[f1,f1] = 0;
  U[f2,f2] = 0;
  U[f1,f2] = 1;
  U[f2,f1] = 1;
  U[f1+N,f1+N] = 0;
  U[f2+N,f2+N] = 0;
  U[f1+N,f2+N] = 1;
  U[f2+N,f1+N] = 1;

  return U,U*Λ*U';
 end

function Exchange_blocks(Γ,g,s1,s2)
 N = convert(Int64, size(Γ,1)/2.);
 U = eye(2*N)
 for Ϟ=0:(g-1)
  u,Γ = Exchange_fermions(Γ,mod1(s1+Ϟ,N),mod1(s2+Ϟ,N))
  U = u*U;
 end

 return U,Γ
end


# function Unentangle_blocks(Γ,g,d,s) #Questa faceva problemi per quench da - a +
#  N = convert(Int64, size(Γ,1)/2.);
#
#  #Bring Γ in the diagonal basis for the two blocks
#  γsx = Reduce_gamma(Γ,g,s);
#  γdx = Reduce_gamma(Γ,g,mod1(s+g+d,N));
#  dsx,usx = Diag_gamma(γsx);
#  ddx,udx = Diag_gamma(γdx);
#  U = convert(Array{Complex{Float64},2}, eye(2*N));
#  U = Inject_gamma(U,usx,s);
#  U = Inject_gamma(U,udx,mod1(s+g+d,N));
#  Γ = U'*Γ*U;
#  # Γ = U*Γ*U';
#
# #I have to set to 0 the desired elements for each one of the 4 blocks.
# #Print_matrix("Γ", abs(Γ));
# #Ora ho una matrice
#
# end


function Unentangle_blocks(Γ,g,d,s)
 N = convert(Int64, size(Γ,1)/2.);


 T,Γ = Exchange_blocks(Γ,g,mod1(s+g,N),mod1(s+g+d,N));
 γ = Reduce_gamma(Γ,g*2,s);
 #Bring Γ in the diagonal basis for the two blocks
 γsx = Reduce_gamma(γ,g,1);
 γdx = Reduce_gamma(γ,g,g+1);
 dsx,usx = Diag_gamma(γsx);
 ddx,udx = Diag_gamma(γdx);
 U = convert(Array{Complex{Float64},2}, eye(4*g));
 U = Inject_gamma(U,usx,1);
 U = Inject_gamma(U,udx,g+1);
 # Print_matrix("γ1", abs(γ));
 γ = U'*γ*U;


 #QUESTO TAGLIA TUTTO->EQUIVALENTE a mtrunc
 γ = Project_just_diagonal(γ);

 #QUESTO TAGLIA SOLO TR E BL, SOLO SPERIMENTALE, NON HA SENSO
 # γ[(1:2*g),(1:2*g)+2*g] = zeros(2*g,2*g);
 # γ[(1:2*g)+2*g,(1:2*g)] = zeros(2*g,2*g);


 γ = U*γ*U';


# Print_matrix("γ2", abs(γ)); #Con questi due print vedo che cancellare queste
#vuol dire semplicemente prendere la somma diretta delle due ridotte! Quindi
#usare questo metodo probabilmente equivale a cancellare tutte le correlazioni per
#una grandezza superiore a m dove m è uguale a d.
#Questo però succede se io  cancello tutto, protrebbe essere che la cosa che si può
#anche fare sia tagliare solo in TR and BL.
Γ = Inject_gamma(Γ,γ,s)
# Γ = Exchange_blocks(Γ,g,mod1(s+g,N),mod1(s+g+d,N));
Γ = T'*Γ*T;
return Γ;
end

function Sweep_Unentangle(Γ,g,d)
 N = convert(Int64, size(Γ,1)/2.);

 for s=1:(N)
  Γ = Unentangle_blocks(Γ,g,d,s)
 end

 return Γ;
end

function Mod_Diag(M,d,ϵ)
 N = convert(Int64, size(M,1))

 for iiter=1:N
  M[iiter,mod1(iiter+d,N)] = ϵ*M[iiter,mod1(iiter+d,N)];
 end

 return M;
end

function Mod_Diag_Symm(M,d,ϵ)
 Mod_Diag(M,d,ϵ);
 Mod_Diag(M,-d,ϵ);

 return M;
end

function Mod_Diag_4(M,d,ϵ)
 N = convert(Int64, div(size(M,1),2))

 M[(1:N),(1:N)]     = Mod_Diag(M[(1:N),(1:N)],d,ϵ);
 M[(1:N)+N,(1:N)]   = Mod_Diag(M[(1:N)+N,(1:N)],d,ϵ);
 M[(1:N),(1:N)+N]   = Mod_Diag(M[(1:N),(1:N)+N],d,ϵ);
 M[(1:N)+N,(1:N)+N] = Mod_Diag(M[(1:N)+N,(1:N)+N],d,ϵ);

 return M;
end


function Mod_Diag_Symm_4(M,d,ϵ)
 N = convert(Int64, div(size(M,1),2))

 M[(1:N),(1:N)]     = Mod_Diag_Symm(M[(1:N),(1:N)],d,ϵ);
 M[(1:N)+N,(1:N)]   = Mod_Diag_Symm(M[(1:N)+N,(1:N)],d,ϵ);
 M[(1:N),(1:N)+N]   = Mod_Diag_Symm(M[(1:N),(1:N)+N],d,ϵ);
 M[(1:N)+N,(1:N)+N] = Mod_Diag_Symm(M[(1:N)+N,(1:N)+N],d,ϵ);

 return M;
end

function Diagonal_elements(M)
 N = convert(Int64, size(M,1));

 v = zeros(Complex{Float64},N)
 for iiter=1:N
  v[iiter] = M[iiter,iiter]
 end

 return v;
end









function reduce_bond_dimension_3(Γ_in,L1,L2,m)
 #Questa funzione prende un vettore di "N_momenti" Gamme e lo ritorna
 #come un vettore di "N_momenti" gamme delle quali il numero di
 #modi totale è "m".
 Γ_init = deepcopy(Γ_in)
 N_momenti = convert(Int64, size(Γ_init,1));
 O = build_riordina_settore(L2)
 for iiter=1:N_momenti
   Γ_init[iiter,:,:] = O*Γ_init[iiter,:,:]*(O');
 end

 #Quanto fare grande il blocchetto in ogni gamma -2
 μ = convert(Int64, 2*(m))
 #Metto in Γ_it la matrice già diagonalizzata nel primo blocchetto
 #così da riempire la tabella v ed avere le Γ di partenza
 Γ   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U_Γ   = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 v   = zeros(Float64, 3, N_momenti);
 v[2,:] = 1;
 for iiter=1:N_momenti
     Γ[iiter,:,:], U_Γ[iiter,:,:]  = diagonalise_block(Γ_init[iiter,:,:],1, μ+2);
     v[1,iiter]                    = real(Γ[iiter,1,1])+real(Γ[iiter,2,2]);
 end
 Tot_modi = L1*L2;
 while(Tot_modi>m)
    k = 0;
    check = true;
    perm = sortperm(v[1,:]);
    pos = 1;
    #Scelgo la piu piccola matrice da cui non ho già tolto tutti i modi
    while (check)
      #PRIMA C'ERA UN +2 come sul quaderno, controlla che vada bene toglierlo in
      #questo caso
        k = convert(Int64, perm[pos]);
        if ((convert(Int64,v[3,k]))==2*L2)
            pos +=1;
        else
            check = false;
        end
    end
    #Cancello le correlazioni e porto a zero l'entropia del modo
    v3k = convert(Int64, v[3,k])
     #println("k: ", k, "v1:", v[1,k], " Γ11: ", Γ[k,v3k+1,v3k+1], " 1-ν: ", Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)])
    Γ[k,v3k+1,:]                  = 0;
    Γ[k,:,v3k+1]                  = 0;
    Γ[k,v3k+1+(2*L2),:]             = 0;
    Γ[k,:, v3k+1+(2*L2)]            = 0;
    Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)]    = 1;
    Γ[k,v3k+2,:]                  = 0;
    Γ[k,:,v3k+2]                  = 0;
    Γ[k,v3k+2+(2*L2),:]             = 0;
    Γ[k,:, v3k+2+(2*L2)]            = 0;
    Γ[k,v3k+2+(2*L2),v3k+2+(2*L2)]    = 1;

    v[3,k] += 2; #->In questa partizione ne ho tolto uno in più

    #Controllo che non sia già arrivato a diagonalizzare l'ultimo blocco
    if (v[2,k]<(2*L2-μ-1))
        v[2,k] += 2;
    end
    Γ[k,:,:], U  = diagonalise_block(Γ[k,:,:],convert(Int64, v[2,k]),(μ+2));
    U_Γ[k,:,:]   = U_Γ[k,:,:]*U;
    v3k = convert(Int64, v[3,k])
    v[1,k]       = real(Γ[k,v3k+1,v3k+1])+real(Γ[k,v3k+2,v3k+2])
    Tot_modi -= 2;
 end
 for kiter=1:N_momenti
     Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 for iiter=1:N_momenti
   Γ[iiter,:,:] = (O')*Γ[iiter,:,:]*O;
 end

 return Γ;
end

function RBD(Γ_in,L1,L2,m)
 #If the bond dimension is maximal do nothing
 if (m>=L2)
  println("L2=", L2, " m=", m, " impossible to RBD")
  return Γ_in;
 end

 Γ = deepcopy(Γ_in)
 N_momenti = convert(Int64, size(Γ,1));
 O = build_riordina_settore(L2)
 for iiter=1:N_momenti
   Γ[iiter,:,:] = O*Γ[iiter,:,:]*(O');
 end

 U_Γ = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 for kiter=1:N_momenti
  U_Γ[kiter,:,:] = eye(4*L2);
 end
 v   = zeros(Float64, 3, N_momenti);


 for siter=1:(L2-m)
  for iiter=1:N_momenti
      Γ[iiter,:,:], U[iiter,:,:] = diagonalise_block(Γ[iiter,:,:],convert(Int64,v[2,iiter]+1), convert(Int64,2*(m+siter)-v[2,iiter]));
      U_Γ[iiter,:,:]             = U_Γ[iiter,:,:]*U[iiter,:,:];
      v3k                        = convert(Int64, v[2,iiter])
      v[1,iiter]                 = real(Γ[iiter,v3k+1,v3k+1])+real(Γ[iiter,v3k+2,v3k+2])
  end
  for miter=1:div(L1,2)
     k = 0;
     check = true;
     perm = sortperm(v[1,:]);
     pos = 1;
     #Scelgo la piu piccola matrice da cui non ho già tolto tutti i modi
     while (check)
       #PRIMA C'ERA UN +2 come sul quaderno, controlla che vada bene toglierlo in
       #questo caso
         k = convert(Int64, perm[pos]);
         if ((convert(Int64,v[2,k]))>=2*(m+siter))
             pos +=1;
         else
             check = false;
         end
     end
     #Cancello le correlazioni e porto a zero l'entropia del modo
     v3k = convert(Int64, v[2,k])
      #println("k: ", k, "v1:", v[1,k], " Γ11: ", Γ[k,v3k+1,v3k+1], " 1-ν: ", Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)])
     Γ[k,v3k+1,:]                   = 0;
     Γ[k,:,v3k+1]                   = 0;
     Γ[k,v3k+1+(2*L2),:]            = 0;
     Γ[k,:, v3k+1+(2*L2)]           = 0;
     Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)] = 1;
     Γ[k,v3k+2,:]                   = 0;
     Γ[k,:,v3k+2]                   = 0;
     Γ[k,v3k+2+(2*L2),:]            = 0;
     Γ[k,:, v3k+2+(2*L2)]           = 0;
     Γ[k,v3k+2+(2*L2),v3k+2+(2*L2)] = 1;

     v[2,k] += 2; #->In questa colonna ne ho tolti 2


     v3k          = convert(Int64, v[2,k])
     v[1,k]       = real(Γ[k,v3k+1,v3k+1])+real(Γ[k,v3k+2,v3k+2])
  end
 end

 # for k=1:N_momenti
 #   println("K-->",k)
 #   for iiter=1:(2*L2)
 #     if (abs(round(real(Γ[k,iiter,iiter]),13))!=0)
 #       print("__",iiter,"__")
 #     end
 #   end
 #   println(" ");
 # end

 for kiter=1:N_momenti
  Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 for iiter=1:N_momenti
  Γ[iiter,:,:] = (O')*Γ[iiter,:,:]*O;
 end

 return Γ;
end



function Half_RBD(Γ_in,L1,L2,m)
 #If the bond dimension is maximal do nothing
 if (m>=div(L2,2))
  println("L2=", L2, " m=", m, " impossible to RBD")
  return Γ_in;
 end

 Γ = deepcopy(Γ_in)
 N_momenti = convert(Int64, size(Γ,1));
 O = build_riordina_settore(L2)
 for iiter=1:N_momenti
   Γ[iiter,:,:] = O*Γ[iiter,:,:]*(O');
 end

 U_Γ = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 for kiter=1:N_momenti
  U_Γ[kiter,:,:] = eye(4*L2);
 end
 v   = zeros(Float64, 3, N_momenti);



 for iiter=1:N_momenti
     Γ[iiter,:,:], U[iiter,:,:] = diagonalise_block(Γ[iiter,:,:],1, div(L2,2));
     U_Γ[iiter,:,:]             = U_Γ[iiter,:,:]*U[iiter,:,:];
     v3k                        = convert(Int64, v[2,iiter])
     v[1,iiter]                 = real(Γ[iiter,v3k+1,v3k+1])+real(Γ[iiter,v3k+2,v3k+2])
 end
 for miter=1:div((div(L2,2)-m)*L1,2)   #Divido per due perchè ne tolgo due per volta
    k = 0;
    check = true;
    perm = sortperm(v[1,:]);
    pos = 1;
    #Scelgo la piu piccola matrice da cui non ho già tolto tutti i modi
    while (check)
      #PRIMA C'ERA UN +2 come sul quaderno, controlla che vada bene toglierlo in
      #questo caso
        k = convert(Int64, perm[pos]);
        if ((convert(Int64,v[2,k]))>=L2)
            pos +=1;
        else
            check = false;
        end
    end
    #Cancello le correlazioni e porto a zero l'entropia del modo
    v3k = convert(Int64, v[2,k])
     #println("k: ", k, "v1:", v[1,k], " Γ11: ", Γ[k,v3k+1,v3k+1], " 1-ν: ", Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)])
    Γ[k,v3k+1,:]                   = 0;
    Γ[k,:,v3k+1]                   = 0;
    Γ[k,v3k+1+(2*L2),:]            = 0;
    Γ[k,:, v3k+1+(2*L2)]           = 0;
    Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)] = 1;
    Γ[k,v3k+2,:]                   = 0;
    Γ[k,:,v3k+2]                   = 0;
    Γ[k,v3k+2+(2*L2),:]            = 0;
    Γ[k,:, v3k+2+(2*L2)]           = 0;
    Γ[k,v3k+2+(2*L2),v3k+2+(2*L2)] = 1;

    v[2,k] += 2; #->In questa colonna ne ho tolti 2


    v3k          = convert(Int64, v[2,k])
    v[1,k]       = real(Γ[k,v3k+1,v3k+1])+real(Γ[k,v3k+2,v3k+2])
 end

 for kiter=1:N_momenti
  Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 for iiter=1:N_momenti
  Γ[iiter,:,:] = (O')*Γ[iiter,:,:]*O;
 end

 return Γ;
end








function RBD_uno_per_volta(Γ_in,L1,L2,mdL)
 #Qua con m intendo m/L1
 #If the bond dimension is maximal do nothing
 if (mdL>=2*L2)
  println("L2=", L2, " m/L=", mdL, " impossible to RBD")
  return Γ_in;
 end

 Γ = deepcopy(Γ_in)
 N_momenti = convert(Int64, size(Γ,1));
 O = build_riordina_settore(L2)
 for iiter=1:N_momenti
   Γ[iiter,:,:] = O*Γ[iiter,:,:]*(O');
 end

 U_Γ = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 for kiter=1:N_momenti
  U_Γ[kiter,:,:] = eye(4*L2);
 end
 v   = ones(Float64, 2, N_momenti);


 for siter=1:(L2-mdL)
  for iiter=1:N_momenti
      Γ[iiter,:,:], U[iiter,:,:] = diagonalise_block(Γ[iiter,:,:],1, 2*(mdL+siter));
      U_Γ[iiter,:,:]             = U_Γ[iiter,:,:]*U[iiter,:,:];
      v2k                        = convert(Int64, v[2,iiter])
      v[1,iiter]                 = real(Γ[iiter,v2k,v2k])
  end
  for miter=1:(2*N_momenti)
     k = 0;
     check = true;
     perm = sortperm(v[1,:]);
     pos = 1;
     #Scelgo la piu piccola matrice da cui non ho già tolto tutti i modi
     while (check)
       #PRIMA C'ERA UN +2 come sul quaderno, controlla che vada bene toglierlo in
       #questo caso
         k = convert(Int64, perm[pos]);
         if ((convert(Int64,v[2,k]))>2*(mdL+siter))
             pos +=1;
         else
             check = false;
         end
     end
     #Cancello le correlazioni e porto a zero l'entropia del modo
     v2k = convert(Int64, v[2,k])
      #println("k: ", k, "v1:", v[1,k], " Γ11: ", Γ[k,v3k+1,v3k+1], " 1-ν: ", Γ[k,v3k+1+(2*L2),v3k+1+(2*L2)])
     Γ[k,v2k,:]                   = 0;
     Γ[k,:,v2k]                   = 0;
     Γ[k,v2k+(2*L2),:]            = 0;
     Γ[k,:,v2k+(2*L2)]            = 0;
     Γ[k,v2k+(2*L2),v2k+(2*L2)]   = 1;
     # Γ[k,v3k+2,:]                   = 0;
     # Γ[k,:,v3k+2]                   = 0;
     # Γ[k,v3k+2+(2*L2),:]            = 0;
     # Γ[k,:, v3k+2+(2*L2)]           = 0;
     # Γ[k,v3k+2+(2*L2),v3k+2+(2*L2)] = 1;

     v[2,k] += 1; #->In questa colonna ne ho tolti 2


     v2k          = convert(Int64, v[2,k])
     v[1,k]       = real(Γ[k,v2k,v2k])
  end
 end

 # for k=1:N_momenti
 #   println("K-->",k)
 #   for iiter=1:(2*L2)
 #     if (abs(round(real(Γ[k,iiter,iiter]),13))!=0)
 #       print("__",iiter,"__")
 #     end
 #   end
 #   println(" ");
 # end

 for kiter=1:N_momenti
  Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 for iiter=1:N_momenti
  Γ[iiter,:,:] = (O')*Γ[iiter,:,:]*O;
 end



 return Γ;
end





function RBD_nuovo(Γ_in,L1,L2,mdL)
 #Qua con m intendo m/L1
 #If the bond dimension is maximal do nothing
 if (mdL>=L2)
  println("L2=", L2, " m/L=", mdL, " impossible to RBD")
  return Γ_in;
 end

 Γ = deepcopy(Γ_in)
 N_momenti = convert(Int64, size(Γ,1));
 O = build_riordina_settore(L2)
 for iiter=1:N_momenti
   Γ[iiter,:,:] = O*Γ[iiter,:,:]*(O');
 end

 U_Γ = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 U = zeros(Complex{Float64}, N_momenti, 4*L2, 4*L2);
 for kiter=1:N_momenti
  U_Γ[kiter,:,:] = eye(4*L2);
 end
 v   = ones(Float64, 2, N_momenti);

##################
for L2iter=1:(L2-mdL)
  for k=1:N_momenti
      Γ[k,:,:], U[k,:,:]  = diagonalise_block(Γ[k,:,:],convert(Int64,v[2,k]), 2*(mdL+L2iter)-v[2,k]+1);
      U_Γ[k,:,:]          = U_Γ[k,:,:]*U[k,:,:];
      TBz                 = convert(Int64, v[2,k])
      v[1,k]              = real(Γ[k,TBz,TBz])
  end
  for L1iter=1:L1
    k = sortperm(v[1,:])
    i = 0;
    for j=1:N_momenti
      i=j;
      if (v[2,k[i]]<=(2*(mdL+L2iter)))
        break;
      end
    end
    Γ[k[i],Int(v[2,k[i]]),:] = 0;
    Γ[k[i],:,Int(v[2,k[i]])] = 0;
    Γ[k[i],Int(v[2,k[i]])+2*L2,:] = 0;
    Γ[k[i],:,Int(v[2,k[i]])+2*L2] = 0;
    Γ[k[i],Int(v[2,k[i]])+2*L2,Int(v[2,k[i]])+2*L2] = 1;

    v[2,k[i]] +=1;
    TBz                 = convert(Int64, v[2,k[i]])
    v[1,k[i]]           = real(Γ[k[i],TBz,TBz])
  end
end
##################

 for kiter=1:N_momenti
  Γ[kiter,:,:] = U_Γ[kiter,:,:]*Γ[kiter,:,:]*(U_Γ[kiter,:,:]');
 end

 for iiter=1:N_momenti
  Γ[iiter,:,:] = (O')*Γ[iiter,:,:]*O;
 end

 return Γ;
end























function Build_hopping_H(N)
    A = zeros(Float64,N,N)
    B = zeros(Float64,N,N)

    for i=1:(N-1)
        A[i,i+1]=1;
    end
    A[1,N] = 1;
    A = (A+A')/2.;

    H = [-A B; -B A];

    return H;
end

function Build_fourier_matrix(N)
    U   = zeros(Complex{Float64},N,N);
    FU  = zeros(Complex{Float64},2*N,2*N);

    for k=1:N
        for x=1:N
            U[k,x]= exp(-im*2*pi*(k-(N-1)/2.)*(x-(N-1)/2.)/N);
        end
    end

    FU = 1/sqrt(N)*[U 0*eye(N); 0*eye(N) (conj(U))]
    return FU;
end


function tilde_det(M)
    return sqrt(abs(det((M))));
end

function Fermionic_Negativity(Λ,N_A)
    # Compute the Fermionic Negativity (Eisler-Zimboras) of the system with
    # first partition from 1 to N_A
    dim_Λ   = size(Λ,1);
    N       = div(dim_Λ,2);
    N_B     = N-N_A;
    F_xxtxp = Build_FxxTxp(N);
    Ω       = Build_Omega(N);
    I       = eye(dim_Λ);

    #### Build the Majorana correlation matrix M ####
    γ = Λ-0.5*I;
    γ = real(-im*Ω*γ*Ω');
    # γ = (γ-γ')/2.;
    γ = F_xxtxp*γ*F_xxtxp';
    ####    ####

    Tp = [eye(2*N_A) zeros(Int64, 2*N_A, 2*N_B); zeros(Int64, 2*N_B, 2*N_A) im*eye(2*N_B)];
    Tm = [eye(2*N_A) zeros(Int64, 2*N_A, 2*N_B); zeros(Int64, 2*N_B, 2*N_A) -im*eye(2*N_B)];

    γp = Tp*γ*Tp;
    γm = Tm*γ*Tm;

    γx = im*(I-(I+im*γm)*inv(I-γp*γm)*(I+im*γp))

    ZZZ = tilde_det((I-γ*γ)/2.);
    TR  = tilde_det(sqrt((I+im*γx)/2.)+sqrt((I-im*γx)/2.));

    return log(TR*sqrt(ZZZ));

end


function Reduce_double_gamma(Λ,s1,d1,s2,d2)
  dim_Λ     = size(Λ,1);
  N         = div(dim_Λ,2);
  F_xxtxp   = Build_FxxTxp(N);
  Ω         = Build_Omega(N);
  I         = eye(dim_Λ);
  γ_reduced = zeros(Complex{Float64},2*(d1+d2),2*(d1+d2))


  γ = F_xxtxp*Λ*F_xxtxp';


  γ_reduced[1:(2*d1),1:(2*d1)] = γ[((1:(2*d1))+2*s1-2),((1:(2*d1))+2*s1-2)];
  γ_reduced[((1:(2*d2))+2*d1),((1:(2*d2))+2*d1)] = γ[((1:(2*d2))+2*s2-2),((1:(2*d2))+2*s2-2)];
  γ_reduced[1:(2*d1),((1:(2*d2))+2*d1)] = γ[((1:(2*d1))+2*s1-2),((1:(2*d2))+2*s2-2)];
  γ_reduced[((1:(2*d2))+2*d1),1:(2*d1)] = γ[((1:(2*d2))+2*s2-2),((1:(2*d1))+2*s1-2)];

  F_xxtxp_reduced = Build_FxxTxp(d1+d2);
  γ_reduced       = F_xxtxp_reduced'*γ_reduced*F_xxtxp_reduced;

  return γ_reduced;

end


function Modes_spatial_contribution(Λ)
    #Return a vectore contributions[s,e] of the weight of the eigenmode
    #to the spatial mode s or viceversa.
    dim_Λ   = size(Λ,1);
    N       = div(dim_Λ,2);

    contributions = zeros(Float64,N,N);
    D,U = Diag_gamma(Λ);
    valy = real(diag(D[1:N,1:N]));
    for s=1:N
        for e=1:N
            contributions[s,e] += real(U[s,e]*conj(U[s,e])+U[s,e+N]*conj(U[s,e+N]));
            contributions[s,e] += real(U[s+N,e]*conj(U[s+N,e])+U[s+N,e+N]*conj(U[s+N,e+N]));
        end
    end
    contributions = contributions;

    return contributions, valy;
end


function U_spatial_contribution(U)
    #Return a vectore contributions[s,e] of the weight of the eigenmode
    #to the spatial mode s or viceversa.
    dim_Λ   = size(U,1);
    N       = div(dim_Λ,2);

    contributions = zeros(Float64,N,N);
    for s=1:N
        for e=1:N
            contributions[s,e] += real(U[s,e]*conj(U[s,e])+U[s,e+N]*conj(U[s,e+N]));
            contributions[s,e] += real(U[s+N,e]*conj(U[s,e])+U[s+N,e+N]*conj(U[s,e+N]));
        end
    end
    contributions = contributions;

    return contributions;
end


function Mutual_information(Λ,N_A)
    # Compute the Mutual Information of the system with
    # first partition from 1 to N_A
    dim_Λ   = size(Λ,1);
    N       = div(dim_Λ,2);
    N_B     = N-N_A;

    vneA  = VN_entropy(Reduce_gamma(Λ,N_A,1));
    vneB  = VN_entropy(Reduce_gamma(Λ,N_B,N_A+1));
    vneT  = VN_entropy(Λ);

    return vneA+vneB-vneT;
end


function myhist(data, min, max, nbins)
  #plot(bin,out)

  N = length(data)             # How many elements in the input vector 'data' ?
  delta = (max-min)/nbins      # Bin size is inferred here from the maximal, minimal, and bin number
  out = zeros(nbins)           # Let's initialize the output data structures for the bin count
  bin = zeros(nbins)           # and for the bin centres...

  start = min                  # Left edge
  for k=1:nbins
    stop   = start + delta   # Right edge
    out[k] = length(find((start.<=data.<stop))) # Count how many elements are between left and right
    bin[k] = start + delta/2. # Centre of the bin
    start  = stop            # New left edge
   end
   return out, bin
  end

  function Ent_cont(Λ)
    dim_Λ   = size(Λ,1);
    N       = div(dim_Λ,2);

    F_xxtxp = Build_FxxTxp(N);
    Ω       = Build_Omega(N);
    I       = eye(dim_Λ);

    #### Build the Majorana correlation matrix M ####
    γ = Λ-0.5*I;
    γ = real(-im*Ω*γ*Ω');
    γ = F_xxtxp*γ*F_xxtxp';
    γ = (γ-γ')/2.;
    #########

    M,O = Diag_real_skew(γ);

    p = zeros(Float64,N,N);

    for i=1:N
      for k=1:N
        p[i,k] = 0.5*((O[2*i-1,2*k-1]*conj(O[2*i-1,2*k-1]))+(O[2*i,2*k-1]*conj(O[2*i,2*k-1]))
        +(O[2*i-1,2*k]*conj(O[2*i-1,2*k]))+(O[2*i,2*k]*conj(O[2*i,2*k])))
      end
    end

    #### start debug ####
    # println("--> Devono esserci solo 1:")
    # for i=1:N
    #   print("- ",sum(p[i,:]));
    # end
    #### end debug ####

    Saffi = zeros(Float64, N);
    for i=1:N
      for k=1:N
        ν = round(0.5+M[2*k-1,2*k],15);
        if (ν<0.)
          # println("ν < ZERO!!!! ", ν)
          ν = 0;
        end
        if (ν>1.)
          # println("ν > UNO!!!!", ν)
          ν = 1;
        end
        if (ν!=0.0 && ν!=1.0)
          Saffi[i] -= p[i,k]*(ν*log2(ν)+(1-ν)*log2(1-ν));
        end
      end
    end

    return Saffi;
  end




  function Ent_cont_distr(Λ)
    dim_Λ   = size(Λ,1);
    N       = div(dim_Λ,2);

    F_xxtxp = Build_FxxTxp(N);
    Ω       = Build_Omega(N);
    I       = eye(dim_Λ);

    #### Build the Majorana correlation matrix M ####
    γ = Λ-0.5*I;
    γ = real(-im*Ω*γ*Ω');
    γ = F_xxtxp*γ*F_xxtxp';
    γ = (γ-γ')/2.;
    #########

    M,O = Diag_real_skew(γ);

    p = zeros(Float64,N,N);

    for i=1:N
      for k=1:N
        p[i,k] = 0.5*((O[2*i-1,2*k-1]*conj(O[2*i-1,2*k-1]))+(O[2*i,2*k-1]*conj(O[2*i,2*k-1]))
        +(O[2*i-1,2*k]*conj(O[2*i-1,2*k]))+(O[2*i,2*k]*conj(O[2*i,2*k])))
      end
    end


    return p;
  end






  function Border_Contour(Λ,ns)
      #Ritorna la somma del contour sui primi ns siti a destra ed a sinistra
      cont = Ent_cont(Λ);
      return (sum(cont[1:ns])+sum(cont[(end-ns+1):end]))
  end


  function fg(Λ,NA,NB,ns,λ)
      #Le trasformazioni fermioniche che posso fare sono
      #facili da caratterizzare sui Majorana, qua sono delle
      #ortogonali
      N   = div(size(Λ,1),2);
      Oab = Orthogonal_composite(λ);
      O = eye(2*N)
      O[(2*(NA-ns+1)-1):(2*(NA+ns)),(2*(NA-ns+1)-1):(2*(NA+ns))] = Oab;
      O[(2*((NA+NB)-ns+1)-1):(2*((NA+NB)+ns)),(2*((NA+NB)-ns+1)-1):(2*((NA+NB)+ns))] = Oab;

      F_xxtxp = Build_FxxTxp(N);
      Ω       = Build_Omega(N);
      I       = eye(2*N);

      #### Build the Majorana correlation matrix γ ####
      γ = deepcopy(Λ-0.5*I);
      γ = real(-im*Ω*γ*Ω');
      # γ = (γ-γ')/2.;
      γ = F_xxtxp*γ*F_xxtxp';
      ####    ####

      γ = O*γ*transpose(O);

      #### Build the symbol matrix Λ ####
      Λ = F_xxtxp'*γ*F_xxtxp;
      Λ = im*(Ω'*Λ*Ω);
      Λ += 0.5*I;
      ####    ####


      return Border_Contour(Reduce_gamma(Λ,NB,NA+1),ns);
  end

  function Transformed_Λ(Λ,NA,NB,ns,λ)
      N   = div(size(Λ,1),2);
      Oab = Orthogonal_composite(λ);
      O = eye(2*N)
      O[(2*(NA-ns+1)-1):(2*(NA+ns)),(2*(NA-ns+1)-1):(2*(NA+ns))] = Oab;
      O[(2*((NA+NB)-ns+1)-1):(2*((NA+NB)+ns)),(2*((NA+NB)-ns+1)-1):(2*((NA+NB)+ns))] = Oab;

      F_xxtxp = Build_FxxTxp(N);
      Ω       = Build_Omega(N);
      I       = eye(2*N);

      #### Build the Majorana correlation matrix γ ####
      γ = deepcopy(Λ-0.5*I);
      γ = real(-im*Ω*γ*Ω');
      # γ = (γ-γ')/2.;
      γ = F_xxtxp*γ*F_xxtxp';
      ####    ####

      γ = O*γ*transpose(O);

      #### Build the symbol matrix Λ ####
      Λ = F_xxtxp'*γ*F_xxtxp;
      Λ = im*(Ω'*Λ*Ω);
      Λ += 0.5*I;
      ####    ####

      return Λ;
  end


  function Back_from_Transformed_Λ(Λ,NA,NB,ns,λ)
      N   = div(size(Λ,1),2);
      Oab = Orthogonal_composite(λ);
      O = eye(2*N)
      O[(2*(NA-ns+1)-1):(2*(NA+ns)),(2*(NA-ns+1)-1):(2*(NA+ns))] = Oab;
      O[(2*((NA+NB)-ns+1)-1):(2*((NA+NB)+ns)),(2*((NA+NB)-ns+1)-1):(2*((NA+NB)+ns))] = Oab;

      F_xxtxp = Build_FxxTxp(N);
      Ω       = Build_Omega(N);
      I       = eye(2*N);

      #### Build the Majorana correlation matrix γ ####
      γ = deepcopy(Λ-0.5*I);
      γ = real(-im*Ω*γ*Ω');
      # γ = (γ-γ')/2.;
      γ = F_xxtxp*γ*F_xxtxp';
      ####    ####

      γ = transpose(O)*γ*O;

      #### Build the symbol matrix Λ ####
      Λ = F_xxtxp'*γ*F_xxtxp;
      Λ = im*(Ω'*Λ*Ω);
      Λ += 0.5*I;
      ####    ####

      return Λ;
  end














  function ASYM_Border_Contour(Λ,ns,😈)
      #Ritorna la somma del contour sui primi ns siti a destra ed a sinistra
      cont = Ent_cont(Λ);
      if (😈==0)
        return sum(cont[1:ns])
      end
      if (😈==1)
        return sum(cont[(end-ns+1):end])
      end
  end


  function ASYM_fg(Λ,NA,NB,ns,λ,😈)
      #Le trasformazioni fermioniche che posso fare sono
      #facili da caratterizzare sui Majorana, qua sono delle
      #ortogonali
      Λ = ASYM_Transformed_Λ(Λ,NA,NB,ns,λ,😈)

      return ASYM_Border_Contour(Reduce_gamma(Λ,NB,NA+1),ns,😈);
  end

  function ASYM_Transformed_Λ(Λ,NA,NB,ns,λ,😈)
    N   = div(size(Λ,1),2);
    Oab = Orthogonal_composite(λ);
    O = eye(2*N)
    O[(2*(NA+(😈*NB)-ns+1)-1):(2*(NA+(😈*NB)+ns)),(2*(NA+(😈*NB)-ns+1)-1):(2*(NA+(😈*NB)+ns))] = Oab;

    F_xxtxp = Build_FxxTxp(N);
    Ω       = Build_Omega(N);
    I       = eye(2*N);

    #### Build the Majorana correlation matrix γ ####
    γ = deepcopy(Λ-0.5*I);
    γ = real(-im*Ω*γ*Ω');
    # γ = (γ-γ')/2.;
    γ = F_xxtxp*γ*F_xxtxp';
    ####    ####

    γ = O*γ*transpose(O);

    #### Build the symbol matrix Λ ####
    Λ = F_xxtxp'*γ*F_xxtxp;
    Λ = im*(Ω'*Λ*Ω);
    Λ += 0.5*I;
    ####    ####

      return Λ;
  end


  function ASYM_Back_from_Transformed_Λ(Λ,NA,NB,ns,λ,😈)
      N   = div(size(Λ,1),2);
      Oab = Orthogonal_composite(λ);
      O = eye(2*N)
      O[(2*(NA+(😈*NB)-ns+1)-1):(2*(NA+(😈*NB)+ns)),(2*(NA+(😈*NB)-ns+1)-1):(2*(NA+(😈*NB)+ns))] = Oab;

      F_xxtxp = Build_FxxTxp(N);
      Ω       = Build_Omega(N);
      I       = eye(2*N);

      #### Build the Majorana correlation matrix γ ####
      γ = deepcopy(Λ-0.5*I);
      γ = real(-im*Ω*γ*Ω');
      # γ = (γ-γ')/2.;
      γ = F_xxtxp*γ*F_xxtxp';
      ####    ####

      γ = transpose(O)*γ*O;

      #### Build the symbol matrix Λ ####
      Λ = F_xxtxp'*γ*F_xxtxp;
      Λ = im*(Ω'*Λ*Ω);
      Λ += 0.5*I;
      ####    ####

      return Λ;
  end

function Envelope(v)
  d = size(v,1);

  envelope_top    = zeros(Float64, d);
  envelope_bottom = zeros(Float64, d);
  width           = zeros(Float64, d);
  for i=1:d
      envelope_top[i]     = maximum(v[i:end]);
      envelope_bottom[i]  = minimum(v[i:end]);
      width[i]            = abs(envelope_top[i]-envelope_bottom[i]);
  end

  return envelope_top,envelope_bottom,width;
end

@. model_lin(x,p)     = p[1]+p[2]*x;
function jacobian_model_lin(x,p)
    J = Array{Float64}(undef, length(x),length(p))
    J[:,1] = 1;    #dmodel/dp[1]
    J[:,2] = x;  #dmodel/dp[2]
    J
end


function Piattezza(v)

        n = size(v,1);
        media = mean(v);

        p = 0;
        for i=(1):(n)
                p += abs(v[i]-media);
        end
        p = p/(n);

        return p;
end

function Max_oscillation(v)

        n = size(v,1);
        media = mean(v);

        p = abs(v[1]-media);
        for i=(2):(n)
          if abs(v[i]-media)>p;
                p = abs(v[i]-media);
              end
        end

        return p;
end


@. model_exp(x,p)     = p[1]+e^(-p[2]*x);
function jacobian_model_exp(x,p)
    J = Array{Float64}(undef, length(x),length(p))
    J[:,1] = 1;    #dmodel/dp[1]
    J[:,2] = -p[2]*e.^(-p[2]*x);  #dmodel/dp[2]
    J
end




function xc(Exact_E, c)
    N_points   = size(Exact_E,1)+1;
    renorm     =   zeros(Float64, N_points-1);
    for n=1:(N_points-1)
        renorm[n] = ((log(abs(Exact_E[n]+c)))-minimum(log.(abs.(Exact_E[:]+c))))/(maximum(log.(abs.(Exact_E[:]+c)))-minimum(log.(abs.(Exact_E[:]+c))));
    end
    der_seconde = zeros(Float64, N_points-3);
    for n=2:(N_points-2)
        der_seconde[n-1] = (renorm[n+1]-2*renorm[n]+renorm[n-1])/2.;
    end
    mean_rapp_incr = mean(der_seconde);
    Δ_rap_incr      = zeros(Float64, N_points-3);
    for n=2:(N_points-2)
        Δ_rap_incr[n-1] = ((der_seconde[n-1]))^2;
    end
    sum_rap_incr = sum(Δ_rap_incr);

    return sum_rap_incr
    # fitt = curve_fit(model_lin, jacobian_model_lin, 2:N_points, log.(abs.(Exact_E+c)), [1.1,-1.1]);
    # return stderror(fitt)[2];
end

function Find_decay_constant(v,a,b)
  #v é un vettore che mi aspetto decada come exp[-αx]+c
  #questa funzione trova c nell'intervallo [a,b]
  #!!Non so quanto sia precisa e tante volte è meglio restringere agli ultimi punti  di v
  C_stimato = optimize(c->xc(v,c),a,b);
  C = Optim.minimizer(C_stimato)

  return C;
end


function Z_function(H_D,β);
  N = div(size(H_D,1),2);

  Z = 1;
  for i=1:N
    Z = Z*(1+exp(-β*H_D[i,i]));
  end

  return Z;
end



function Biggest_eig_rho(M,n_gaps)
    #La grandezza del vettore restituito scala come 2^(n_gaps-1);
    N = convert(Int64, size(M,1)/2);

    D,U = Diag_gamma(M)

    d_rev = sort(real(diag(D)),rev=true); #Il primo è il più grande, I primi N sono maggiori di 0.5
    # d_norev = sort(real(diag(D)));        #Questo è come 1-d_rev

    initial_coeff = 1;
    for i=1:(N-n_gaps)
        initial_coeff = initial_coeff*d_rev[i]
    end
    evor = initial_coeff*ones(Float64,convert(Int64,2^(n_gaps-1))+1);

    for i=0:convert(Int64,2^(n_gaps-1))
        bin_string = i;
        for j=0:(n_gaps-1)
            evor[i+1] = evor[i+1]*abs(mod(bin_string,2)-d_rev[N-j]);
            bin_string = div(bin_string,2);
        end
    end

    return sort(evor,rev=true);
end


function Biggest_eig_H(M,n_eig)
    #La grandezza del vettore restituito scala come 2^(n_gaps-1);
    N = convert(Int64, size(M,1)/2);

    D,U = Diag_ferm(M)

    # d_rev = sort(real(diag(D)),rev=true); #Il primo è il più grande, sono positivi
    d_norev = sort(real(diag(D)));          #Il primo è il più piccolo, sono negativi

    # println(sum(d_norev[1:N]))

    initial_coeff = 0;
    for i=1:(N-n_eig)
        initial_coeff += d_norev[i]
    end
    evor = initial_coeff*ones(Float64,convert(Int64,2^(n_eig-1))+1);

    for i=0:convert(Int64,2^(n_eig-1))
        bin_string = i;
        for j=0:(n_eig-1)
            evor[i+1] += (((-1)^(mod(bin_string,2)))*d_norev[N-j]);
            bin_string = div(bin_string,2);
        end
    end

    return sort(evor,rev=true);
end



function Number_fermion(Mat,D,U)
   N_f = convert(Int64, size(Mat,1)/2.);

   number_fermion = 0;

   M = deepcopy(Mat)
   #####AGGIUNTO
   M = (M+M')/2.
   M[1,1] += eps()
   ####

   M_diag_base = real(U'*M*U);
   for iiter=1:(N_f)
       number_fermion += M_diag_base[iiter,iiter];
   end

   return real(number_fermion);
end


function Majoranise(M)
  N   = div(size(M,1),2);
  F_xxtxp = Build_FxxTxp(N);
  Ω       = Build_Omega(N);
  I       = eye(2*N);

  #### Build the Majorana correlation matrix γ ####
  γ = deepcopy(M-0.5*I);
  γ = real(-im*Ω*γ*Ω');
  # γ = (γ-γ')/2.;
  γ = F_xxtxp*γ*F_xxtxp';
  γ = (γ-γ')/2.;
  #########

  return γ
end


function DeMajoranise(M)
  N   = div(size(M,1),2);
  F_xxtxp = Build_FxxTxp(N);
  Ω       = Build_Omega(N);
  I       = eye(2*N);

  Λ = F_xxtxp'*M*F_xxtxp;
  Λ = im*(Ω'*Λ*Ω);
  Λ += 0.5*I;


  return Λ;
end


function Parity(M)
  #I multiply by two the matrix to rinormalise to 1 or -1
  #If M has dimension 2nx2n then, for a scalar λ, Pf(λA)=λ^n*Pf(A)
  #Suppose for example that the M is the Majorana correlation matrix
  #for a pure state, thus is the direct sum of {{0,0.5},{-0.5,0}}.
  #The Pfaffian would be 0.5^n, multiplying M by 2, the Pfaffian become 1.
  return (Pfaffian(Majoranise(2*M)));
end

function Parity2(M)
  N   = div(size(M,1),2);

  x = Complex(0);
  for i=1:N
    x += (complex(-1))^M[i,i];
  end

  return x;
end

function Parity3(M)
  N   = div(size(M,1),2);

  x = Complex(0);
  for i=1:N
    x += 1-2*M[i,i];
  end

  return x;
end


function GDE(g_ferm_i,U_diag_f_Q)
  return (U_diag_f_Q*Project_diagonals(U_diag_f_Q'*g_ferm_i*U_diag_f_Q,0)*U_diag_f_Q');;
end



function VN_entropy_basse_energie(M,χ)
   N = convert(Int64, div(size(M,1),2));

   D,U = Diag_gamma(M);   #Se voglio fare con eig
   D = (real(D));

   #Forzo la Hermitianità, se ho risultati strani controlla QUA
   #che M sia abbastanza vicina all'hermitiana prima.
   # M   = (M+M')/2.
   # D,U = Diag_gamma(M);
   S = 0;



   for chiiter=1:χ
     iiter = N-chiiter+1;
       if (round(D[iiter,iiter],18)<-0.0000000000001)
        De,Ue = eig((M+M')/2.);
        for iiter=1:N
         println("DG: ", D[iiter,iiter]);
        end
        for iiter=1:N
         println("DE: ", De[iiter]);
        end
        save("/home/jacopo/Dropbox/ermatr.jld","M", M)
        error = string("Eigenvalue in VE not in [0,1]: ",round(D[iiter,iiter],18))
        throw(ArgumentError(error))
       end
       nu = abs(round(D[iiter,iiter],18))
       if (nu != 0 && nu != 1)
        #Invece di arrivare fino a N/2 nel ciclo e sommare nu e 1-nu, li passo tutti
        #perchè potrebbe essere che facendo una generica transformazione ortogonale non li abbia
        #ordinati in coppie, anche se Diag_gamma dovrebbe metterli in ordine
           S -= log(nu)*nu+log(1-nu)*(1-nu);
       end
   end

   return S;
end



function Ent_cont_basse_energie(Λ,χ)
  dim_Λ   = size(Λ,1);
  N       = div(dim_Λ,2);

  F_xxtxp = Build_FxxTxp(N);
  Ω       = Build_Omega(N);
  I       = eye(dim_Λ);

  #### Build the Majorana correlation matrix M ####
  γ = Λ-0.5*I;
  γ = real(-im*Ω*γ*Ω');
  γ = F_xxtxp*γ*F_xxtxp';
  γ = (γ-γ')/2.;
  #########

  M,O = Diag_real_skew(γ);

  p = zeros(Float64,N,N);

  for i=1:N
    for k=1:N
      p[i,k] = 0.5*((O[2*i-1,2*k-1]*conj(O[2*i-1,2*k-1]))+(O[2*i,2*k-1]*conj(O[2*i,2*k-1]))
      +(O[2*i-1,2*k]*conj(O[2*i-1,2*k]))+(O[2*i,2*k]*conj(O[2*i,2*k])))
    end
  end

  #### start debug ####
  # println("--> Devono esserci solo 1:")
  # for i=1:N
  #   print("- ",sum(p[i,:]));
  # end
  #### end debug ####

  Saffi = zeros(Float64, N);
  for i=1:N
    for s=1:χ
      k = N-s+1;
      ν = round(0.5+M[2*k-1,2*k],15);
      if (ν<0.)
        # println("ν < ZERO!!!! ", ν)
        ν = 0;
      end
      if (ν>1.)
        # println("ν > UNO!!!!", ν)
        ν = 1;
      end
      if (ν!=0.0 && ν!=1.0)
        Saffi[i] -= p[i,k]*(ν*log2(ν)+(1-ν)*log2(1-ν));
      end
    end
  end

  return Saffi;
end


function VN_entropy_singolo_modo(M,χ)
   N = convert(Int64, div(size(M,1),2));

   D,U = Diag_gamma(M);   #Se voglio fare con eig
   D = (real(D));

   #Forzo la Hermitianità, se ho risultati strani controlla QUA
   #che M sia abbastanza vicina all'hermitiana prima.
   # M   = (M+M')/2.
   # D,U = Diag_gamma(M);
   S = 0;



    iiter = N-χ+1;
     if (round(D[iiter,iiter],18)<-0.0000000000001)
      De,Ue = eig((M+M')/2.);
      for iiter=1:N
       println("DG: ", D[iiter,iiter]);
      end
      for iiter=1:N
       println("DE: ", De[iiter]);
      end
      save("/home/jacopo/Dropbox/ermatr.jld","M", M)
      error = string("Eigenvalue in VE not in [0,1]: ",round(D[iiter,iiter],18))
      throw(ArgumentError(error))
     end
     nu = abs(round(D[iiter,iiter],18))
     if (nu != 0 && nu != 1)
      #Invece di arrivare fino a N/2 nel ciclo e sommare nu e 1-nu, li passo tutti
      #perchè potrebbe essere che facendo una generica transformazione ortogonale non li abbia
      #ordinati in coppie, anche se Diag_gamma dovrebbe metterli in ordine
         S -= log(nu)*nu+log(1-nu)*(1-nu);
     end

   return S;
end


function Ent_cont_singolo_modo(Λ,χ)
  dim_Λ   = size(Λ,1);
  N       = div(dim_Λ,2);

  F_xxtxp = Build_FxxTxp(N);
  Ω       = Build_Omega(N);
  I       = eye(dim_Λ);

  #### Build the Majorana correlation matrix M ####
  γ = Λ-0.5*I;
  γ = real(-im*Ω*γ*Ω');
  γ = F_xxtxp*γ*F_xxtxp';
  γ = (γ-γ')/2.;
  #########

  M,O = Diag_real_skew(γ);

  p = zeros(Float64,N,N);

  for i=1:N
    for k=1:N
      p[i,k] = 0.5*((O[2*i-1,2*k-1]*conj(O[2*i-1,2*k-1]))+(O[2*i,2*k-1]*conj(O[2*i,2*k-1]))
      +(O[2*i-1,2*k]*conj(O[2*i-1,2*k]))+(O[2*i,2*k]*conj(O[2*i,2*k])))
    end
  end

  #### start debug ####
  # println("--> Devono esserci solo 1:")
  # for i=1:N
  #   print("- ",sum(p[i,:]));
  # end
  #### end debug ####

  Saffi = zeros(Float64, N);
  for i=1:N
    k = N-χ+1;
    ν = round(0.5+M[2*k-1,2*k],15);
    if (ν<0.)
      # println("ν < ZERO!!!! ", ν)
      ν = 0;
    end
    if (ν>1.)
      # println("ν > UNO!!!!", ν)
      ν = 1;
    end
    if (ν!=0.0 && ν!=1.0)
      Saffi[i] -= p[i,k]*(ν*log2(ν)+(1-ν)*log2(1-ν));
    end
  end

  return Saffi;
end







function Diag_real_skew_eig(M)
 N_f = convert(Int64, size(M,1)/2.)

  DD,UU = eig(M);
  vec_sort = sortperm(abs.(DD), rev=true);
  mat_sort = zeros(Complex{Float64}, 2*N_f, 2*N_f);
  for i=1:2*N_f
      mat_sort[i, vec_sort[i]] = 1.0;
  end
  #questa li mette in ordine ma non sono sicuro siano sempre +-,+-,+-. Devo fare un check o assicurarmi che cmq questo algoritmo faccia effettivamente così
  block_antidiag = zeros(Complex{Float64}, 2*N_f, 2*N_f);
  for i=1:N_f
      block_antidiag[2*i-1,2*i-1] = -(1/sqrt(2))*im;
      block_antidiag[2*i-1,2*i]   = *(1/sqrt(2))*im;
      block_antidiag[2*i,2*i-1]   = (1/sqrt(2))*1;
      block_antidiag[2*i,2*i]     = (1/sqrt(2))*1;
  end
  UOU = im*block_antidiag*mat_sort*UU';

  return real.(im*block_antidiag*mat_sort*diagm(DD)*(im*block_antidiag*mat_sort)'), UOU';
end

function Diag_ferm_eig(M)
    N_f = convert(Int64, size(M,1)/2.)

    F_xptxx = Build_FxpTxx(N_f);

    Omega  = Build_Omega(N_f)
    M_temp = real(-im*Omega*M*Omega')
    # M_temp += eps()*rand(2*N_f,2*N_f);
    M_temp = (M_temp-M_temp')/2.;
    # M_temp[1,1] += eps();
    M_temp, O = Diag_real_skew_eig(M_temp)
    M_temp = F_xptxx*M_temp*(F_xptxx')
    M_temp = im*Omega'*M_temp*(Omega)

    M_f = M_temp;
    U_f = Omega'*O'*(F_xptxx')*Omega;

    return real(M_f), U_f
end
