using F_utilities;
using LinearAlgebra;
using Test;

@testset "Hopping model" begin
    N = 128;
    H   = Build_hopping_hamiltonian(N,true);
    U_ω = Build_Fourier_matrix(N);
    D_Fourier = U_ω'*H*U_ω;
    D,U = Diag_h(H);
    @testset "Energies" begin
        @test all(sort(diag(real.(D_Fourier))).-sort(real.(diag(D))) .< 10^(-12));
    end
end
