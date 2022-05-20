#### Preferred indexing:
# -----------------------------
# | 1    | 2    | 3    | 4    |
# -----------------------------
# | 5    | 6    | 7    | 8    |
# -----------------------------
# | 9    | 10   | 11   | 12   |
# -----------------------------
# | 13   | 14   | 15   | 16   |
# -----------------------------

function T1D(Lx,Ly,r)
  return (r[2]-1)*Lx+r[1]
end

function T2D(Lx,Ly,N)
  return (rem(N-1,Lx)+1,div(N,Lx+1)+1)
end


"""
  Reduce_gamma2D(M, Lx,Ly, lx,ly, x,y)

  If `M` is a Dirac correlation matrix of a 2D system of dimensions '(Lx,Ly)' , returns the reduced density matrix of the system
   of dimension '(lx,ly)' starting with the top-left corner in '(x,y)'.
"""
function Reduce_gamma2D(M, Lx,Ly, lx,ly, r)
  #Check if the partition fits inside the system
  if !((r[1]+lx)<=Lx&&(r[2]+ly)<=Ly)
    error("In Reduce_gamma2D the selected partition exceeds the boundaries of the system")
  end

  redgammaTL = zeros(Complex{Float64}, lx*ly, lx*ly);
  redgammaTR = zeros(Complex{Float64}, lx*ly, lx*ly);
  redgammaBL = zeros(Complex{Float64}, lx*ly, lx*ly);
  redgammaBR = zeros(Complex{Float64}, lx*ly, lx*ly);
  MTL=M[(1:(Lx*Ly)),(1:(Lx*Ly))];
  MTR=M[(1:(Lx*Ly)).+(Lx*Ly),(1:(Lx*Ly))];
  MBL=M[(1:(Lx*Ly)),(1:(Lx*Ly)).+(Lx*Ly)];
  MBR=M[(1:(Lx*Ly)).+(Lx*Ly),(1:(Lx*Ly)).+(Lx*Ly)];
  for i=0:(ly-1)
    for j=0:(ly-1)
      redgammaTL[(i*lx+1):((i+1)*lx),(j*lx+1):((j+1)*lx)]=MTL[T1D(Lx,Ly,(r+[0,i])):(T1D(Lx,Ly,(r+[0,i]))+lx-1),T1D(Lx,Ly,(r+[0,j])):(T1D(Lx,Ly,(r+[0,j]))+lx-1)]
      redgammaTR[(i*lx+1):((i+1)*lx),(j*lx+1):((j+1)*lx)]=MTR[T1D(Lx,Ly,(r+[0,i])):(T1D(Lx,Ly,(r+[0,i]))+lx-1),T1D(Lx,Ly,(r+[0,j])):(T1D(Lx,Ly,(r+[0,j]))+lx-1)]
      redgammaBL[(i*lx+1):((i+1)*lx),(j*lx+1):((j+1)*lx)]=MBL[T1D(Lx,Ly,(r+[0,i])):(T1D(Lx,Ly,(r+[0,i]))+lx-1),T1D(Lx,Ly,(r+[0,j])):(T1D(Lx,Ly,(r+[0,j]))+lx-1)]
      redgammaBR[(i*lx+1):((i+1)*lx),(j*lx+1):((j+1)*lx)]=MBR[T1D(Lx,Ly,(r+[0,i])):(T1D(Lx,Ly,(r+[0,i]))+lx-1),T1D(Lx,Ly,(r+[0,j])):(T1D(Lx,Ly,(r+[0,j]))+lx-1)]
    end
  end

  return [redgammaTL redgammaTR; redgammaBL redgammaBR]
end



"""
  Biggest_eig_rho(M, n)
  Return the first 2^n-1 eigenvalues of the density matrix rho described by the correlation matrix 'M'.
  In the SciPost notation, it computes the eigenvalues vec{x}=Binary(0),Binary(1),dots,Binary(2^n-1).
  The first 'n-1' eigenvalues are assured to be the first 'n-1' biggest eigenvalues of rho.

"""
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
