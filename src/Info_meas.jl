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
