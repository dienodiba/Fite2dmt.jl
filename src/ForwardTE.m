% Dieno Diba 2022
% 2D finite element magnetotelluric modelling: TE Mode

function [rxy,pxy,iyy,tzy] = ForwardTE(el2no,no2yz,rho,period,sta,topo)

%km to m
no2yz = 1e+3 * no2yz;
topo = 1e+3 * topo;
sta = 1e+3 * sta;

nel = length(el2no);
nno = length(no2yz);

% Minimum edge length over all the elements
for ide = 1:nel
    edge1 = sqrt((no2yz(el2no(ide,3),1) - no2yz(el2no(ide,2),1))^2 + (no2yz(el2no(ide,3),2) - no2yz(el2no(ide,2),2))^2);
    edge2 = sqrt((no2yz(el2no(ide,1),1) - no2yz(el2no(ide,3),1))^2 + (no2yz(el2no(ide,1),2) - no2yz(el2no(ide,3),2))^2);
    edge3 = sqrt((no2yz(el2no(ide,2),1) - no2yz(el2no(ide,1),1))^2 + (no2yz(el2no(ide,2),2) - no2yz(el2no(ide,1),2))^2);
    tmp = min([edge1,edge2,edge3]);
    if ide == 1
        minedge = tmp;
    else
        if tmp < minedge
            minedge = tmp;
        end
    end
end

omega = 2*pi/period;
mu0 = 4*pi*1e-7;

% U: Electric field vector
U = zeros(nno,1);

% Boundary nodes and elements
ymin = min(no2yz(:,1));
ymax = max(no2yz(:,1));
zmin = min(no2yz(:,2));
zmax = max(no2yz(:,2));
zlef = [];
zrig = [];
ytop = [];
ybot = [];
nole = [];
nori = [];
notp = [];
nobt = [];
ysur = [];
nsur = [];
for idno = 1:nno
    if no2yz(idno,1) == ymin
        zlef = [zlef;no2yz(idno,2)];
        nole = [nole;idno];
    elseif no2yz(idno,1) == ymax
        zrig = [zrig;no2yz(idno,2)];
        nori = [nori;idno];
    elseif no2yz(idno,2) == zmin
        ytop = [ytop;no2yz(idno,1)];
        notp = [notp;idno];
    elseif no2yz(idno,2) == zmax
        ybot = [ybot;no2yz(idno,1)];
        nobt = [nobt;idno];
    end
    if abs(no2yz(idno,2) - interp1(topo(:,1),topo(:,2),no2yz(idno,1))) < minedge/2
        ysur = [ysur;no2yz(idno,1)];
        nsur = [nsur;idno];
    end
end
nole = sortrows([nole zlef],2);
nole = nole(:,1);
zlef = sort(zlef);
nori = sortrows([nori zrig],2);
nori = nori(:,1);
zrig = sort(zrig);
notp = sortrows([notp ytop],2);
notp = notp(:,1);
ytop = sort(ytop);
nobt = sortrows([nobt ybot],2);
nobt = nobt(:,1);
ybot = sort(ybot);
elle = [];
elri = [];
nsur = sortrows([nsur ysur],2);
nsur = nsur(:,1);
ysur = sort(ysur);
for ide = 1:nel
    if sum(ismember(el2no(ide,:),nole)) == 2
        tmp = ismember(el2no(ide,:),nole);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz(el2no(ide,pt1),2) no2yz(el2no(ide,pt2),2)]);        
        elle = [elle;[ide tmd]];
    end
    if sum(ismember(el2no(ide,:),nori)) == 2
        tmp = ismember(el2no(ide,:),nori);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz(el2no(ide,pt1),2) no2yz(el2no(ide,pt2),2)]);        
        elri = [elri;[ide tmd]];
    end
end
elle = sortrows(elle,2);
role = rho(elle(:,1));
elri = sortrows(elri,2);
rori = rho(elri(:,1));
% - Top boundary
Utop = 5;
U(notp) = Utop * ones(length(notp),1);
% - Left boundary
U(nole) = BoundaryTE(role,period,zlef,Utop);
% - Right boundary
U(nori) = BoundaryTM(rori,period,zrig,Utop);
% - Bottom boundary
U(nobt) = linspace(U(nole(end)),U(nori(end)),length(nobt));

% Assemble coefficient matrix De
% Compute element matrix and add the contribution to the global matrix
De = zeros(nno);
for ide = 1:nel
    % Nodes and their coordinates for the element ide
    no_ide = el2no(ide,:)';
    yz_ide = no2yz(no_ide,:)';
    % Edges
    s1_ide = yz_ide(:,3) - yz_ide(:,2);
    s2_ide = yz_ide(:,1) - yz_ide(:,3);
    s3_ide = yz_ide(:,2) - yz_ide(:,1);
    % Area of the element
    Ar_ide = abs(0.5 * (s2_ide(1)*s3_ide(2) - s2_ide(2)*s3_ide(1)));
    % Gradient of test functions
    grad_phi1 = [-s1_ide(2);s1_ide(1)]/(2 * Ar_ide);
    grad_phi2 = [-s2_ide(2);s2_ide(1)]/(2 * Ar_ide);
    grad_phi3 = [-s3_ide(2);s3_ide(1)]/(2 * Ar_ide);
    grad_phi = [grad_phi1 grad_phi2 grad_phi3];
    % Compute all the integrals for this particular element
    De_ide = zeros(3);
    for i = 1:3
        for j = 1:3
            if i == j
                De_ide(i,j) = (-(grad_phi(:,i)'*grad_phi(:,j)) - 1i*omega*mu0/rho(ide)*1/3) * Ar_ide;
            else
                De_ide(i,j) = -(grad_phi(:,i)'*grad_phi(:,j)) * Ar_ide;
            end
        end
    end
    De(no_ide,no_ide) = De(no_ide,no_ide) + De_ide;
end

% Known boundary nodes (noex), all nodes (noal), inner unknown nodes (noin)
noex = ([nole;nori;notp;nobt]);
noal = 1:nno;
noin = setdiff(noal,noex);

% Solve the linear system
U(noin) = De(noin,noin)\(-De(noin,noex)*U(noex));

% Compute transfer functions at the stations
nst = size(sta,1);
rxy = zeros(nst,1);
pxy = zeros(nst,1);
iyy = zeros(nst,1);
tzy = zeros(nst,1);
[~,nl1] = ismember([ymin 0],no2yz,'rows');
[~,nl2] = ismember([ymin zlef(find(zlef == 0)+1)],no2yz,'rows');
% Y at y=min(y),z=0 [for COMMEMI models]
Y0 = -1/1i/omega/mu0 * (U(nl2) - U(nl1))/no2yz(nl2,2);
for ids = 1:nst
    itp1 = find(abs(ysur-sta(ids,1)) == min(abs(ysur-sta(ids,1))));
    itp1 = itp1(1);
    if ysur(itp1) > sta(ids,1)
        itp2 = itp1 - 1;
    else
        itp2 = itp1 + 1;
    end
    tmp = find(no2yz(:,1) == ysur(itp1));
    ns1 = tmp(abs(no2yz(tmp,2)-sta(ids,2)) == min(abs(no2yz(tmp,2)-sta(ids,2))));
    tmp = find(no2yz(:,1) == ysur(itp2));
    ns2 = tmp(abs(no2yz(tmp,2)-sta(ids,2)) == min(abs(no2yz(tmp,2)-sta(ids,2))));
    for ide = 1:nel
        if sum(ismember([ns1 ns2],el2no(ide,:))) == 2
            tmpa = setdiff(el2no(ide,:),[ns1 ns2]);
            if no2yz(tmpa,2) > interp1(topo(:,1),topo(:,2),no2yz(tmpa,1))
                nbsa = tmpa;
            end
        end
    end
    a = sqrt((no2yz(ns1,1)-no2yz(nbsa,1))^2 + (no2yz(ns1,2)-no2yz(nbsa,2))^2);
    b = sqrt((no2yz(ns2,1)-no2yz(nbsa,1))^2 + (no2yz(ns2,2)-no2yz(nbsa,2))^2);
    c = sqrt((no2yz(ns1,1)-no2yz(ns2,1))^2 + (no2yz(ns1,2)-no2yz(ns2,2))^2);
    C = acos((a^2 + b^2 - c^2)/(2 * a * b));
    h = a*b/c * sin(C);
    d = a * cos(asin(h/a));
    if d > c
        disp(strcat("Obtuse scalene trianglular element under station ",num2str(ids)))
        disp("The response functions might be innacurate")
        Up = mean([U(ns1) U(ns2)]);
    else
        Up = interp1([0 c],[U(ns1) U(ns2)],d);
    end
    Y = -1/1i/omega/mu0 * (U(nbsa) - Up)/h;
    % Impedance, apparent resistivity, phase
    Zxy = mean([U(ns1) U(ns2)])/Y;
    rxy(ids) = 1/omega/mu0 * abs(Zxy).^2;
    pxy(ids) = 180/pi * atan2(imag(Zxy),real(Zxy));
    % Tzy
    if no2yz(ns1,1) > no2yz(ns2,1)
        Z = 1/1i/omega/mu0 * (U(ns1) - U(ns2))/c;
    else
        Z = 1/1i/omega/mu0 * (U(ns2) - U(ns1))/c;
    end
    tzy(ids) = Z/Y;
    % COMMEMI Iyy
    iyy(ids) = Y/Y0;
end

end