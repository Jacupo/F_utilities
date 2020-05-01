using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;

N = 64;

#Generate and diagonalise the hamiltonian
H = Fu.Random_NNhamiltonian(N)
H_D, U_D = Fu.Diag_h(H)

#Print the energy modes epsilon_k
figure("Energies")
plot(1:N,diag(H_D)[1:N])
xlabel(L"$k$")
ylabel(L"$\epsilon_k$")
