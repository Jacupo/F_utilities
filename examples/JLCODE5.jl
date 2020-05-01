using F_utilities;
using PyPlot;
using LinearAlgebra;

const Fu = F_utilities;

N   = 10;
theta   = pi/8;

H_APBC          = Fu.TFI_Hamiltonian(N, theta; PBC=-1);
HD_APBC, U_APBC = Fu.Diag_h(H_APBC);
NE_APBC         = diag(HD_APBC);

H_PBC           = Fu.TFI_Hamiltonian(N, theta; PBC=+1);
HD_PBC, U_PBC   = Fu.Diag_h(H_PBC);
NE_PBC          = diag(HD_PBC);

H_OBC           = Fu.TFI_Hamiltonian(N, theta; PBC=0);
HD_OBC, U_OBC   = Fu.Diag_h(H_OBC,2);
NE_OBC          = diag(HD_OBC);

AE_APBC = zeros(Float64, 2*N);
AE_PBC  = similar(AE_APBC);
AE_OBC  = similar(AE_APBC);

#The solutions of equation sin((N+1)*phi)/sin(N*phi)=-1/cot(theta)
phi = [0.293377974249272,0.586547314382234,
        0.879273168816649,1.17126278605144,
        1.46212217642804,1.75129510871389,
        2.03798675412035,2.32111594487769,
        2.59947341172037,2.87247738375037]
for n=1:N
    AE_APBC[n]  = sqrt(1+cot(theta)^2+2*cot(theta)*cos(2*π*((1-N)/2+n-1)/N));
    AE_PBC[n]   = sqrt(1+cot(theta)^2+2*cot(theta)*cos(2*π*((-N)/2+n-1)/N));
    AE_OBC[n]   = sqrt(1+cot(theta)^2+2*cot(theta)*cos(phi[n]));
end

AE_APBC[(N+1):(2*N)]    = -AE_APBC[1:N];
AE_PBC[(N+1):(2*N)]     = -AE_PBC[1:N];
AE_OBC[(N+1):(2*N)]     = -AE_OBC[1:N];


fig = plt.figure("Comparison Analitical and Numerical Results",figsize=(10, 6), dpi=80)
plt.subplots_adjust(wspace=0, hspace=0)
ax1 = plt.subplot2grid((21,10), (0,0), colspan=10, rowspan=7);
ax1.set_title("Comparison Analitical and Numerical ϵ_p/2")
ax1.plot(sort(AE_APBC), color="black",  marker="o",
        markersize=10, mfc="none" ,   label="Analytical APBC");
ax1.plot(sort(NE_APBC), color="red", marker="+",
        markersize=10,  linestyle="None", label="Numerical APBC");
ax1.xaxis.set_ticklabels([])
legend();
ax2 = plt.subplot2grid((21,10), (7,0), colspan=10, rowspan=7);
ax2.plot(sort(NE_PBC),  color="purple",  marker="o",
        markersize=10, mfc="none" , label="Analytical PBC" );
ax2.plot(sort(AE_PBC), color="green", marker="+",
        markersize=10,  linestyle="None", label="Numerical PBC");
ax2.xaxis.set_ticklabels([])
legend();
ax3 = plt.subplot2grid((21,10), (14,0), colspan=10, rowspan=7);
ax3.plot(sort(NE_OBC),   color="blue",  marker="o",
        markersize=10, mfc="none" , label="Analytical OBC" );
ax3.plot(sort(AE_OBC),  color="orange", marker="+",
        markersize=10,  linestyle="None", label="Numerical OBC");
ax3.set_xlabel("p")
legend();
tight_layout();

GS_APBC     = Fu.GS_gamma(HD_APBC,U_APBC);
GS_PBC      = Fu.GS_gamma(HD_PBC,U_PBC);
E_GS_APBC   = Fu.Energy(GS_APBC,(HD_APBC,U_APBC));
E_GS_PBC    = Fu.Energy(GS_PBC,(HD_PBC,U_PBC));

println("Ground State Energies");
println("G_F=-1 :       ", E_GS_APBC);
println("G_F=+1 :       ", E_GS_PBC);
