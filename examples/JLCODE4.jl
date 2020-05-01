using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;


N=50;
N_steps = 100;
delta_steps = 0.1;

#Build the circulant vector for the adaa part of the Gamma with exponential decaying correlations
adaa = zeros(Complex{Float64},N);#
for i=1:div(N,2)
    adaa[i] = exp(-i*0.15)*(rand()+im*rand())
end
adaa[((div(N,2))+1):N]= reverse(adaa[1:div(N,2)]);
#Build the translational invariant adaa part of the Gamma
Gamma_adaa = Fu.Circulant(adaa);
Gamma_adaa = (Gamma_adaa+Gamma_adaa')/2.

#Build the circulant vector for the aa part of the Gamma
aa = zeros(Complex{Float64},N);
aa[2] = rand()+im*rand();
aa[3] = rand()+im*rand();
#Build the translational invariant aa part of the Gamma
Gamma_aa = Fu.Circulant(aa)
Gamma_aa = (Gamma_aa-transpose(Gamma_aa))/2.;

#Build the translational invariant Gamma
Gamma= zeros(Complex{Float64}, 2N,2N);
Gamma[(1:N),(1:N)]          = Gamma_adaa;
Gamma[(1:N).+N,(1:N).+N]    = (I-Gamma_adaa)';
Gamma[(1:N).+N,(1:N)]       = Gamma_aa;
Gamma[(1:N),(1:N).+N]       = -conj(Gamma_aa);
Fu.Print_complex_matrix("Gamma",Gamma)

H   = Fu.Build_hopping_hamiltonian(N,true);
D,U = Diag_h(H);


Gamma_evolved   = copy(Gamma);
adaa        = zeros(Complex{Float64}, N_steps)
aa          = similar(adaa);

#Start the time evolution cycle
#at each cycle it saves the value of two correalotors
adaa[1]     = Gamma_evolved[1,2];
aa[1]       = Gamma_evolved[N+1,2];
for t=2:N_steps
    global Gamma_evolved = Fu.Evolve(Gamma_evolved,(D,U),delta_steps);
    adaa[t]     = Gamma_evolved[1,2];
    aa[t]       = Gamma_evolved[N+1,2];
end

figure("Evolutions")
plot(real.(adaa), color="black", label=L"$\mathfrak{R}(\langle a_1^{\dagger}a_2 \rangle)(t)$");
plot(imag.(adaa), color="black",linestyle="--", label=L"$\mathfrak{I}(\langle a_1^{\dagger} a_2 \rangle)(t)$");
plot(real.(aa), color="purple", label=L"$\mathfrak{R}(\langle a_1a_2 \rangle|(t))$");
plot(imag.(aa), color="purple", linestyle="--", label=L"$\mathfrak{I}(\langle a_1 a_2 \rangle)(t)$");
legend()
xlabel("timestep")
