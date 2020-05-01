using F_utilities;
using PyPlot;

const Fu = F_utilities;

theta   = pi/8;
Delta_E = zeros(Float64, 47)

for N=4:50
        H_APBC          = Fu.TFI_Hamiltonian(N, theta; PBC=-1);
        HD_APBC, U_APBC = Fu.Diag_h(H_APBC);
        H_PBC           = Fu.TFI_Hamiltonian(N, theta; PBC=+1);
        HD_PBC, U_PBC   = Fu.Diag_h(H_PBC);

        E_GS_APBC   = Fu.Energy(Fu.GS_gamma(HD_APBC,U_APBC),(HD_APBC,U_APBC));
        E_GS_PBC    = Fu.Energy(Fu.GS_gamma(HD_PBC,U_PBC),(HD_PBC,U_PBC));

        global Delta_E[N-3]= abs(E_GS_APBC-E_GS_PBC);
end

figure("|E_GS(GF=+1)-E_GS(GF=-1)|")
plot(4:50, log10.(abs.(Delta_E)));
xlabel(L"$N$");
ylabel(L"$|E_{GS}(g_F=+1,N)-E_{GS}(g_F=-1,N)|$");
tight_layout();
