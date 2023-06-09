using F_utilities;
using LinearAlgebra;
using Test;

@testset "TFIM zero mode" begin
    N = 200
    J = 1
    Γ = 1e-3
    J_coupling = rand(Float64, N) * J
    J_coupling[end] = 0
    Γ_coupling = rand(Float64, N) * Γ
    offdiagonal = [reduce(vcat, [-Γ_coupling[idx], J_coupling[idx]] for idx in 1:N-1); -Γ_coupling[end]]
    h_free = Tridiagonal(-offdiagonal, zeros(2N), offdiagonal)
    h_d, O_free = F_utilities.Diag_real_skew(h_free, 0)
    free_modes = diag(h_d, 1)[begin:2:end]  # Biggest to smallest convention
    @test all(abs.(h_free - (O_free * h_d * O_free')) .<= 1e-12)
    @test minimum(free_modes) < 1e-12
    @test count(free_modes .< 1e-12) == 1
end

@testset "Hopping model" begin
    N = 128;
    H   = F_utilities.Build_hopping_hamiltonian(N,PBC=true);
    U_ω = F_utilities.Build_Fourier_matrix(N);
    D_Fourier = U_ω'*H*U_ω;
    D,U = F_utilities.Diag_h(H);
    @testset "Energies" begin
        @test all(sort(diag(real.(D_Fourier))).-sort(real.(diag(D))) .< 10^(-12));
    end
end