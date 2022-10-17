% Dieno Diba 2022
% 1D analytical solution for 2D boundary condition: TM mode

function [Xleft, Xright] = BoundaryTM (rho, T, z1, Xtop)

mu = 4*pi*10^(-7);
freq = 1/T;
n = length(freq);
omega = 2*pi*freq;

rhoki = rho(:,1);
rhoka = rho(:,end);

for idk = 1:2
    if idk == 1
        rhol = rhoki;
    else
        rhol = rhoka;
    end
    
    %Count the number of layer
    N = 1;
    for idr = 1:length(rhol)-1
        if rhol(idr) ~= rhol(idr+1)
            N = N + 1;
        end
    end
    
    %The case when resistivity in the boundary is non-uniform
    if N > 1
        rho1d = zeros(N,1);
        tmp = 1;
        for idr = 1:length(rhol)
            if idr == 1 || rhol(idr) ~= rhol(idr-1)
                rho1d(tmp) = rhol(idr);     %resistivity of each layer
                tmp = tmp + 1;
            end
        end
        h = zeros(N-1,1);
        tmp = 1;
        za = 0;
        for idz = 1:length(z1)
            if tmp < N
                if rhol(idz) ~= rhol(idz+1)
                    h(tmp) = z1(idz+1) - za;    %layer thickness 
                    tmp = tmp + 1;
                    za = z1(idz+1);
                end
            end
        end

        rhoa = zeros(n,1);
        pha = zeros(n,1);
        z = zeros(n,1);
        for x = 1:n
            Z = zeros(N,1);
            Z(N) = sqrt(1i*omega(x)*mu*rho1d(N));
            for m = N-1:-1:1
                Z0(m) = sqrt(1i*omega(x)*mu*rho1d(m));
                R(m) = (Z0(m) - Z(m+1)) / (Z0(m) + Z(m+1)); 
                k(m) = Z0(m) / rho1d(m);
                Z(m) = Z0(m) * ((1 - (R(m)*exp(-2*k(m)*h(m)))) / (1 + (R(m)*exp(-2*k(m)*h(m)))));
            end
            rhoa(x) = (1/(omega(x)*mu)) * (abs(Z(1)))^2;
            pha(x) = (180/pi) * atan2(imag(Z(1)),real(Z(1)));
            z(x) = Z(1);    %surface impedance
        end

        %X(z) = A*exp(ikz) + B*exp(-ikz)
        A = zeros(N,1);
        B = zeros(N,1);
        A(1) = 1e-200;
        A(N) = 0;
        
        dcy = zeros(N,1);
        for idx = 1:N
            dcy(idx) = sqrt(-1i*omega*mu*(1/rho1d(idx)));   %(k)
        end
        Cof = zeros(N-1,1);     %A/B
        dep = 0;
        for idx = 1:N-1
            dep = dep + h(idx);
            Cof(idx) = ((Z(idx+1) - 1i*rho1d(idx)*dcy(idx))/(-1i*rho1d(idx)*dcy(idx) - Z(idx+1)))*exp(-2i*dcy(idx)*dep);
        end
        dep = 0;
        for idx = 1:N-1
            dep = dep + h(idx);
            B(idx) = A(idx)/Cof(idx);
            F = A(idx)*exp(1i*dcy(idx)*dep) + B(idx)*exp(-1i*dcy(idx)*dep);
            if idx < N-1
                A(idx+1) = F / (exp(1i*dcy(idx+1)*dep) + (Cof(idx+1)^-1)*exp(-1i*dcy(idx+1)*dep));
            end
        end
        B(N) = (A(N-1)*exp(1i*dcy(N-1)*dep) + B(N-1)*exp(-1i*dcy(N-1)*dep))/exp(-1i*dcy(N)*dep);

        dpth = zeros(N,1);
        for idh = 1:length(h)
            dpth(idh + 1) = dpth(idh) + h(idh);
        end
        dpth = nonzeros(dpth);

        tmp = 1;
        X = zeros(length(z1),1);
        for idz = 1:length(z1)
            if tmp < N && Cof(tmp) ~= 0
                X(idz) = A(tmp)*exp(1i*dcy(tmp)*z1(idz)) + B(tmp)*exp(-1i*dcy(tmp)*z1(idz));
                if z1(idz) >= dpth(tmp)
                    tmp = tmp + 1;
                end
            else
                if idz == 1
                    X(idz) = exp(-1i*dcy(tmp)*z1(idz));
                else
                    X(idz) = X(idz-1)*exp(-1i*dcy(tmp)*(z1(idz)-z1(idz-1)));
                end
            end
        end
        rscX = X / X(1) * Xtop;
    end
    %The case when resistivity in the boundary is uniform
    if N == 1       
        rhohom = rhol(1);
        X = zeros(length(z1),1);
        X(1) = Xtop;
        for idz = 1:length(z1)-1
            dcy = sqrt(-1i*omega*mu*(1/rhohom));
            X(idz+1) = X(1)*exp(-1i*dcy*z1(idz+1));
        end
        rscX = X;
    end
    
    %Replace NaN with 0
    %NaN resulted from a division by complex inf
    for idx = 1:length(X)
        if isnan(rscX(idx)) == 1
            rscX(idx) = 0;
        end
    end
    
    if idk == 1
        Xleft = rscX;
    else
        Xright = rscX;
    end
end

end