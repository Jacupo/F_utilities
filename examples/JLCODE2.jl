using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;

N = 64;
H = Fu.Random_NNhamiltonian(N);
H_D, U_D = Fu.Diag_h(H);

Gamma = Fu.GS_Gamma(H_D, U_D);
println("The energy of the ground state is: ", Fu.Energy_fermions(Gamma,H_D, U_D));

N_A = 32;
#I consider the reduced state over the sites 17,2,...,48
Gamma_A = Fu.Reduce_gamma(Gamma,N_A,17);
#I compute the entangement entropy
S_A = Fu.VN_entropy(Gamma_A);
#I compute the contour of partition A
c_A = Fu.Contour(Gamma_A);

lbl_title   = string(L"$S(A) = $", S_A);
lbl_legend  = string(L"$\sum_{i=1}^{N_A} c_A(i) = $", sum(c_A));
figure("Contour of A")
title(lbl_title)
plot(1:N_A, c_A, marker="o", label=lbl_legend);
xlabel("i")
ylabel(L"$c_A(i)$")
legend();
