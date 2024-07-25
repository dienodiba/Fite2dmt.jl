# Dieno Diba 2022
# FITE2DMT
# Calculate the fields at side boundary for TM mode
# 1D analytical solution

function CalculateFieldAtSideBoundaryTM(rho::Any, T::Number, z1::Vector, Xtop::Number)
    mu = 4*pi*1e-07
    freq = 1/T
    omega = 2*pi*freq

    rhoki = rho[:,1]
    rhoka = rho[:,end]

    Xleft = zeros(ComplexF64,length(z1))
    Xright = zeros(ComplexF64,length(z1))

    for idk = 1:2
        if idk == 1
            rhol = rhoki
        else
            rhol = rhoka
        end

        # Count the number of layer
        N = 1
        for idr = 1:length(rhol)-1
            if rhol[idr] != rhol[idr+1]
                N += 1
            end
        end

        
        if N == 1   # If the resistivity in the boundary is homogeneous
            rhohom = rhol[1]
            X = zeros(ComplexF64,length(z1))
            X[1] = Xtop
            for idz = 1:length(z1)-1
                dcy = sqrt(-im*omega*mu/rhohom)
                X[idz+1] = X[1]*exp(-im*dcy*z1[idz+1])
            end
            rscX = X
        else
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
            za = 0
            for idz = 1:length(z1)
                if tmp < N
                    if rhol[idz] != rhol[idz+1]
                        h[tmp] = z1[idz+1] - za    # thickness of each layer 
                        tmp += 1
                        za = z1[idz+1]
                    end
                end
            end

            # 1D forward calculation to get impedance (Z)
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

            # X(z) = A*exp(ikz) + B*exp(-ikz)
            A = zeros(ComplexF64,N)
            B = zeros(ComplexF64,N)

            A[1] = 1.0+0.0*im
            A[N] = 0.0+0.0*im
            
            dcy = zeros(ComplexF64,N)
            for idx = 1:N
                dcy[idx] = sqrt(-im*omega*mu/rho1d[idx])   # decay factor (k)
            end
            Cof = zeros(ComplexF64,N-1)     # Cof = A/B
            dep = 0
            for idx = 1:N-1
                dep = dep + h[idx]
                Cof[idx] = ((Z[idx+1] - im*rho1d[idx]*dcy[idx])/(-im*rho1d[idx]*dcy[idx] - Z[idx+1]))*exp(-2*im*dcy[idx]*dep)
            end
            
            dep = 0
            for idx = 1:N-1
                dep = dep + h[idx]
                B[idx] = A[idx]/Cof[idx]
                F = A[idx]*exp(im*dcy[idx]*dep) + B[idx]*exp(-im*dcy[idx]*dep)
                if idx < N-1
                    A[idx+1] = F / (exp(im*dcy[idx+1]*dep) + (Cof[idx+1]^-1)*exp(-im*dcy[idx+1]*dep))
                end
            end
            B[N] = (A[N-1]*exp(im*dcy[N-1]*dep) + B[N-1]*exp(-im*dcy[N-1]*dep))/exp(-im*dcy[N]*dep)

            dpth = zeros(Float64,N)
            for idh = 1:length(h)
                dpth[idh+1] = dpth[idh] + h[idh]
            end
            dpth = filter(!iszero,dpth)

            tmp = 1
            X = zeros(ComplexF64,length(z1))
            for idz = 1:length(z1)
                if tmp < N && Cof[tmp] != 0
                    X[idz] = A[tmp]*exp(im*dcy[tmp]*z1[idz]) + B[tmp]*exp(-im*dcy[tmp]*z1[idz])
                    if z1[idz] >= dpth[tmp]
                        tmp += 1
                    end
                else
                    if idz == 1
                        X[idz] = exp(-im*dcy[tmp]*z1[idz]);
                    else
                        X[idz] = X[idz-1]*exp(-im*dcy[tmp]*(z1[idz]-z1[idz-1]))
                    end
                end
            end
            rscX = X / X[1] * Xtop
        end

        # NaN may result from a division by complex inf
        # Replace NaN with 0
        for idx = 1:length(X)
            if isnan(rscX[idx]) == 1
                rscX[idx] = 0
            end
        end
        
        if idk == 1
            Xleft = rscX
        else
            Xright = rscX
        end
    end

    return Xleft, Xright
end
