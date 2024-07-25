# Dieno Diba 2022
# FITE2DMT
# Calculate the field at side boundaries for TE Mode
# 1D analytical solution

function CalculateFieldAtSideBoundaryTE(rho::Any, T::Number, z1::Vector, Utop::Number)
    mu = 4*pi*1e-07
    freq = 1/T
    omega = 2*pi*freq

    rhoki = rho[:,1]
    rhoka = rho[:,end]

    Uleft = zeros(ComplexF64,length(z1))
    Uright = zeros(ComplexF64,length(z1))

    z2 = z1 .- z1[1]        # z2 starts from 0 

    for idk = 1:2
        if idk == 1
            rhol = rhoki
        else
            rhol = rhoka
        end
        
        N = 1
        for idr = 1:length(rhol)-1
            if rhol[idr] != rhol[idr+1]
                N += 1      # number of layer
            end
        end
        
        if N > 1
            rho1d = zeros(Float64,N)
            tmp = 1
            for idr = 1:length(rhol)
                if idr == 1 || rhol[idr] != rhol[idr-1]
                    rho1d[tmp] = rhol[idr]     # resistivity of each layer
                    tmp += 1
                end
            end
            h = zeros(Float64,N-1)
            tmp = 1
            za = z2[1]
            for idz = 1:length(z2)
                if tmp < N
                    if rhol[idz] != rhol[idz+1]
                        h[tmp] = z2[idz+1] - za     # thickness of each layer
                        tmp += 1
                        za = z2[idz+1]
                    end
                end
            end

            Z = zeros(ComplexF64,N)
            Z0 = zeros(ComplexF64,N-1) 
            R = zeros(ComplexF64,N-1)
            k = zeros(ComplexF64,N-1)
            Z[N] = sqrt(im*omega*mu*rho1d[N])
            for m = N-1:-1:1        # looping layer
                Z0[m] = sqrt(im*omega*mu*rho1d[m])
                R[m] = (Z0[m] - Z[m+1]) / (Z0[m] + Z[m+1])
                k[m] = Z0[m] / rho1d[m]
                Z[m] = Z0[m] * ((1 - (R[m]*exp(-2*k[m]*h[m]))) / (1 + (R[m]*exp(-2*k[m]*h[m]))))
            end

            # E(z) = A*exp(ikz) + B*exp(-ikz)
            A = zeros(ComplexF64,N)
            B = zeros(ComplexF64,N)

            A[1] = 1.0 + 0.0*im
            A[N] = 0.0 + 0.0*im
            
            dcy = zeros(ComplexF64,N)
            for idx = 1:N
                dcy[idx] = sqrt(-im*omega*mu/rho1d[idx])       # decay factor (k)
            end
            Cof = zeros(ComplexF64,N-1)     # Cof = A/B
            dep = 0.0
            for idx = 1:N-1
                dep = dep + h[idx]
                Cof[idx] = ((Z[idx+1] - omega*mu/dcy[idx])/(Z[idx+1] + omega*mu/dcy[idx]))*exp(-2*im*dcy[idx]*dep)
            end
            dep = 0.0
            for idx = 1:N-1
                dep += h[idx]
                B[idx] = A[idx]/Cof[idx]
                F = A[idx]*exp(im*dcy[idx]*dep) + B[idx]*exp(-im*dcy[idx]*dep)     # F is E
                if idx < N-1
                    A[idx+1] = F / (exp(im*dcy[idx+1]*dep) + (Cof[idx+1]^-1)*exp(-im*dcy[idx+1]*dep))
                end
            end
            B[N] = (A[N-1]*exp(im*dcy[N-1]*dep) + B[N-1]*exp(-im*dcy[N-1]*dep))/exp(-im*dcy[N]*dep)
            
            dpth = zeros(Float64,N)
            for idh = 1:length(h)
                dpth[idh+1] = dpth[idh] + h[idh]
            end
            dpth = filter!(!iszero,dpth)
            
            tmp = 1
            U = zeros(ComplexF64,length(z2))
            for idz = 1:length(z2)
                if tmp < N && Cof[tmp] != 0
                    U[idz] = A[tmp]*exp(im*dcy[tmp]*z2[idz]) + B[tmp]*exp(-im*dcy[tmp]*z2[idz])      
                    if z2[idz] >= dpth[tmp]
                        tmp += 1
                    end
                else
                    U[idz] = U[idz-1]*exp(-im*dcy[tmp]*(z2[idz]-z2[idz-1]))
                end
            end
            rscU = U / U[1] * Utop
        end
        
        # NaN may result from a division by complex inf
        # Anticipate NaN value
        for idu = 1:length(U)       
            if isnan(rscU[idu]) == 1
                rscU[idu] = 0
            end
        end
        
        if idk == 1
            Uleft = rscU
        else
            Uright = rscU
        end
    end

    return Uleft,Uright
end

