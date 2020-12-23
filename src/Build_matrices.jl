   function Circulant(cv)
      cv = reverse(cv)
      return Toeplitz(vcat(cv,cv[1:(end-1)]));
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

  function Build_hopping_hamiltonian(N;PBC=false);
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


  function Build_A_TFI(N, θ, PBC)
     M_A      = LinearAlgebra.diagm(-1 => ones(Float64,N-1), 0 => 2*cot(θ)*ones(Float64,N), 1 => ones(Float64,N-1));
     M_A[1,N] = PBC;
     M_A[N,1] = PBC;

     return -1/2. .*M_A;
  end

  function Build_B_TFI(N, PBC)
     M_B     = LinearAlgebra.diagm(-1 => ones(Float64,N-1), 0 => zeros(Float64,N), 1 => -ones(Float64,N-1));
     M_B[1,N] = PBC;
     M_B[N,1] = -PBC;

     return -1/2. .*M_B;
  end

  function TFI_Hamiltonian(N, θ; PBC=+1)
     A = Build_A_TFI(N, θ, PBC);
     B = Build_B_TFI(N, PBC);

     H_TFI                    = zeros(Float64, 2*N, 2*N);
     H_TFI[1:N,1:N]           = -A;
     H_TFI[((1:N).+N),1:N]    = -B;
     H_TFI[1:N,(1:N).+N]      =  B;
     H_TFI[(1:N).+N,(1:N).+N] =  A;

     return H_TFI;
  end

  function Build_A_TFI_FIX(N, θ, PBC)
     M_A      = LinearAlgebra.diagm(-1 => ones(Float64,N-1), 0 => 2*cot(θ)*ones(Float64,N), 1 => ones(Float64,N-1));
     M_A[1,N] = PBC;
     M_A[N,1] = PBC;
     M_A[1,1] = 0;
     M_A[1,2] = 0;
     M_A[2,1] = 0;

     return -1/2. .*M_A;
  end
  function Build_B_TFI_FIX(N, PBC)
     M_B     = LinearAlgebra.diagm(-1 => ones(Float64,N-1), 0 => zeros(Float64,N), 1 => -ones(Float64,N-1));
     M_B[1,N] = 0;
     M_B[N,1] = -0;
     M_B[1,2] = 0;
     M_B[2,1] = -0;

     return -1/2. .*M_B;
  end

  function TFI_Hamiltonian_FIX(N, θ; PBC=+1)
     A = Build_A_TFI_FIX(N, θ, PBC);
     B = Build_B_TFI_FIX(N, PBC);

     H_TFI                    = zeros(Float64, 2*N, 2*N);
     H_TFI[1:N,1:N]           = -A;
     H_TFI[((1:N).+N),1:N]    = -B;
     H_TFI[1:N,(1:N).+N]      =  B;
     H_TFI[(1:N).+N,(1:N).+N] =  A;

     return H_TFI;
  end







  function Build_A_JXJY(Jx, Jy, lambdas)
     dimension = size(Jx, 1);
     M_A       = zeros(Float64, dimension, dimension)
     for iiter=2:dimension-1
         M_A[iiter, iiter]   = -2*lambdas[iiter];
         M_A[iiter, iiter-1] = -(Jx[iiter-1]+Jy[iiter-1]);
         M_A[iiter, iiter+1] = -(Jx[iiter]+Jy[iiter]);
     end
     M_A[1,1]         = -2*lambdas[1];#
     M_A[1,2]         = -(Jx[1]+Jy[1]);#

     M_A[1,dimension] = -(Jx[dimension]+Jy[dimension]);#
     M_A[dimension,1] = -(Jx[dimension]+Jy[dimension]);
     M_A[dimension, dimension-1] = -(Jx[dimension-1] + Jy[dimension-1]);   #
     M_A[dimension,dimension]    = -2*lambdas[dimension];#

     return M_A
  end



  function Build_B_JXJY(Jx, Jy)
     dimension = size(Jx, 1);
     M_B       = zeros(Float64, dimension, dimension)
     for iiter=2:dimension-1
         M_B[iiter, iiter-1] = -(Jx[iiter-1]-Jy[iiter-1]);
         M_B[iiter, iiter+1] = Jx[iiter]-Jy[iiter];
     end
     M_B[1,2]         = Jx[1]-Jy[1];

     M_B[1,dimension] = -(Jx[dimension]-Jy[dimension]);#
     M_B[dimension,1] = (Jx[dimension]-Jy[dimension]);#

     M_B[dimension, dimension-1] = -(Jx[dimension-1]-Jy[dimension-1]);

     return M_B
  end


  function JXJY_Hamiltonian(N,Jx,Jy,lambda)
     A = 0.5*Build_A_JXJY(Jx,Jy,lambda);
     B = 0.5*Build_B_JXJY(Jx,Jy);

     H_JXJY                          = zeros(Float64, 2*N, 2*N);
     H_JXJY[1:N,1:N]             = -A;
     H_JXJY[((1:N).+N),1:N]       = -B;
     H_JXJY[1:N,(1:N).+N]       =  B;
     H_JXJY[(1:N).+N,(1:N).+N] =  A;

     return H_JXJY;
  end
