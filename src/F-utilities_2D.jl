function site_1d(i,j,M,N)
    return convert(Int64, (i-1)*N+j)
end

function Map_to_2d(index,M,N)
    return convert(Int64,ceil(index/N)),convert(Int64,index-(N*(ceil(index/N)-1)))
end

function right_site_1d(i,j,M,N)
    return (i-1)*N+mod(j,N)+1
end

function bottom_site_1d(i,j,M,N)
    return mod(i,M)*N+j
end

function Build_2DA(M,N,J_x,J_y,lambda)
    A   = zeros(Float64, M*N, M*N);
    J_a = J_x+J_y;

    for iiter=1:M
        for jiter=1:N
            current_site = site_1d(iiter,jiter,M,N);
            right_site   = right_site_1d(iiter,jiter,M,N);
            bottom_site  = bottom_site_1d(iiter,jiter,M,N);
            if (right_site != current_site)
                A[current_site, right_site]     = -J_a;
                A[right_site, current_site]     = -J_a;
            end
            if (bottom_site != current_site)
                A[current_site, bottom_site]    = -J_a;
                A[bottom_site, current_site]    = -J_a;
            end
            A[current_site, current_site]   = -2*lambda;
        end
    end

    return A;
end

function Build_2DB(M,N,J_x,J_y)
    B   = zeros(Float64, M*N,M*N);
    J_b = J_x-J_y;

    for iiter=1:M
        for jiter=1:N
            current_site = site_1d(iiter,jiter,M,N);
            right_site   = right_site_1d(iiter,jiter,M,N);
            bottom_site  = bottom_site_1d(iiter,jiter,M,N);
            if (right_site != current_site)
                B[current_site,right_site]      = J_b;
                B[right_site, current_site]     = -J_b;
            end
            if (bottom_site != current_site)
                B[current_site, bottom_site]    = J_b;
                B[bottom_site, current_site]    = -J_b;
            end
        end
    end

    return B;
end

function Hamiltonian_2D(M,N,J_a,J_b,lambda)
    A = 0.5*Build_2DA(M,N,J_a,J_b,lambda);
    B = 0.5*Build_2DB(M,N,J_a,J_b);

    H                          = zeros(Float64, 2*M*N, 2*M*N);
    H[1:N_f,1:N_f]             = -A;
    H[(1:N_f)+N_f,1:N_f]       = -B;
    H[1:N_f,(1:N_f)+N_f]       =  B;
    H[(1:N_f)+N_f,(1:N_f)+N_f] =  A;

    return H;
end
