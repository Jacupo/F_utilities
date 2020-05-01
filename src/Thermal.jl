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
