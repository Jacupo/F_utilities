function Diag_real_skew(M, rand_perturbation::Int64=0)
  N = div(size(M, 1), 2)

  # Random perturbation before forcing skew symmetrisation
  if (rand_perturbation != 0)
    if (rand_perturbation == 1)
        random_M = rand(2N, 2N) * eps()
        random_M = (random_M - random_M') / 2.0
        M += random_M
      end
      if (rand_perturbation == 4)
        random_M = zeros(Complex{Float64}, 2N, 2N)
        random_M[1, N+2] = eps()
        random_M[2, N+1] = -random_M[1, N+2]
        random_M = (random_M - random_M') / 2.0
        M += random_M
      end
      if (rand_perturbation == 5)
        random_M = zeros(Complex{Float64}, 2N, 2N)
        r = eps()
        M[1, 2] += r
        M[2, 1] -= r
      end
  end

  # M = real((M-M')/2.); #Force skew-symmetry
  # #Random pertubation after the skew symmetrization
  # if (rand_perturbation != 0)
  #   if (rand_perturbation == 2)    #Perturb the diagonal elements (loose perfect skew-symmetry)
  #     M += diagm(rand(2*N)*eps())
  #   end
  #   if (rand_perturbation == 3)  #Perturb the whole matrix (loose perfect skew-symmetry)
  #     random_M = 1*rand(2*N,2*N)*eps();
  #     random_M = (random_M-random_M')/2.;
  #     M += random_M;
  #   end
  # end

  Schur_object = LinearAlgebra.schur(M)

  Schur_ort_i = Schur_object.vectors
  Schur_blocks_i = Schur_object.Schur

  # Reorder so that all 2x2 blocks have positive value in top right, accounting for 1x1 blocks.
  swap_perm = Vector{Int64}(undef, 0)
  sizehint!(swap_perm, 2N)
  iiter = 1
  while iiter < 2N
      if abs(Schur_blocks_i[iiter, iiter]) >= abs(Schur_blocks_i[iiter+1, iiter]) # alternative: s_b_i[i+1, i] == 0
          # We have a 1x1 block - move it to the front. These always come in pairs, but not always sequentially.
          pushfirst!(swap_perm, iiter)
          iiter += 1
      elseif (Schur_blocks_i[iiter + 1, iiter] >= 0.0)
          # Flipped 2x2 block
          push!(swap_perm, iiter + 1, iiter)
          iiter += 2
      else
          # Unflipped 2x2 block
          push!(swap_perm, iiter, iiter + 1)
          iiter += 2
      end
  end
  if iiter == 2N  # catch the final 1x1 block if it exists
    pushfirst!(swap_perm, 2N)
  end
  
  M_temp = Schur_blocks_i[swap_perm, swap_perm]
  O_temp = Schur_ort_i[:, swap_perm]

  # Sort the blocks, λ_1>=λ_2>=...>=λ_N with λ_1 the coefficient in the upper left block
  psort = sortperm(diag(M_temp, 1)[begin:2:end], rev=true)
  full_psort = zeros(Int64, 2N)
  full_psort[begin:2:end] .= 2 .* psort .- 1
  full_psort[begin + 1:2:end] .= 2 .* psort

  M_f = M_temp[full_psort, full_psort]
  O_f = real(O_temp[:, full_psort])

  return M_f, O_f
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

  return real.(U_f'*M*U_f), U_f;#real(M_f), U_f;
end

function Diag_gamma(Γ,rand_perturbation::Int64=0)
Γ = (Γ+Γ')/2.;
γ, U = Diag_h(Γ-0.5*I,rand_perturbation);

return U'*Γ*U,U;#real(γ+0.5*eye(size(Γ,1))),U
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
    if real(D[index,index])<0
      Gamma_diag_base[index+N,index+N] = 1;
    else
      Gamma_diag_base[index,index] = 1;
    end
    # if real(D[index+N,index+N])<=0
    #   Gamma_diag_base[index,index] = 1;
    # end
  end
  Gamma = U*Gamma_diag_base*U';
  Gamma = (Gamma+(Gamma'))/2.

  return Gamma;
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
