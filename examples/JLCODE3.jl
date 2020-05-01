using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;


N =127;

H   = Fu.Build_hopping_hamiltonian(N,true);

U_omega = Fu.Build_Fourier_matrix(N);
D_omega = U_omega'*H*U_omega;
D,U = Fu.Diag_h(H);


figure("Energies")
plot(diag(real.(D_omega))[(N+1):(2*N)],label="Method Fourier");
plot(real.(diag(D))[(N+1):(2*N)], label="Method Diag_h");
xlabel(L"$k$");
ylabel(L"$\epsilon_k$");
legend();

Gamma_omega = Fu.GS_gamma(D_omega,U_omega);
Gamma   = Fu.GS_gamma(D,U);
println("")
println("Energy GS Method Fourier:      ",Fu.Energy(Gamma_omega,(D_omega,U_omega)))
println("En GS Method Diag_h:           ",Fu.Energy(Gamma,(D,U)))

Fu.Print_complex_matrix("Differenza Gamma", Gamma-Gamma_omega)
