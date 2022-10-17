% Dieno Diba 2022
% 2D finite element magnetotelluric modelling: TM Mode

function [ryx,pyx] = ForwardTM(el2no,no2yz,rho,period,sta,topo)

%km to m
no2yz = 1e+3*no2yz; 
topo = 1e+3*topo;   
sta = 1e+3*sta;     

% Remove elements and nodes in the air half-space
rmel = [];
for ide = 1:length(el2no)
    tmp = no2yz(el2no(ide,:),:);
    if mean(tmp(:,2)) < interp1(topo(:,1),topo(:,2),mean(tmp(:,1)))
        rmel = [rmel;ide];
    end
end
syel = setdiff(1:length(el2no),rmel);
el2no_h = el2no(syel,:);
rho(rmel,:) = [];
nel_h = length(el2no_h);
rmno = [];
tmp = 0;
for idn = 1:length(no2yz)
    if ismember(idn,unique(el2no_h(:))) == 0
        rmno = [rmno;idn];
    else
        tmp = tmp + 1;
        trft(tmp,:) = [idn tmp];
    end
end
no2yz_h = no2yz(setdiff(1:length(no2yz),rmno),:);
nno_h = length(no2yz_h);
for ide = 1:nel_h
    for idv = 1:3
        el2no_h(ide,idv) = trft(trft(:,1) == el2no_h(ide,idv),2);
    end
    edge1 = sqrt((no2yz_h(el2no_h(ide,3),1) - no2yz_h(el2no_h(ide,2),1))^2 + (no2yz_h(el2no_h(ide,3),2) - no2yz_h(el2no_h(ide,2),2))^2);
    edge2 = sqrt((no2yz_h(el2no_h(ide,1),1) - no2yz_h(el2no_h(ide,3),1))^2 + (no2yz_h(el2no_h(ide,1),2) - no2yz_h(el2no_h(ide,3),2))^2);
    edge3 = sqrt((no2yz_h(el2no_h(ide,2),1) - no2yz_h(el2no_h(ide,1),1))^2 + (no2yz_h(el2no_h(ide,2),2) - no2yz_h(el2no_h(ide,1),2))^2);
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

% X: Magnetic field vector
X = zeros(nno_h,1);

% Boundary nodes and elements
ymin = min(no2yz_h(:,1));
ymax = max(no2yz_h(:,1));
zmax = max(no2yz_h(:,2));
zlef = [];
zrig = [];
ysur = [];
ybot = [];
nole = [];
nori = [];
nsur = [];
nobt = [];
for idno = 1:nno_h
    if no2yz_h(idno,1) == ymin
        zlef = [zlef;no2yz_h(idno,2)];
        nole = [nole;idno];
    elseif no2yz_h(idno,1) == ymax
        zrig = [zrig;no2yz_h(idno,2)];
        nori = [nori;idno];
    elseif no2yz_h(idno,2) == zmax
        ybot = [ybot;no2yz_h(idno,1)];
        nobt = [nobt;idno];
    elseif abs(no2yz_h(idno,2) - interp1(topo(:,1),topo(:,2),no2yz_h(idno,1))) < minedge/2
        ysur = [ysur;no2yz_h(idno,1)];
        nsur = [nsur;idno];
    end
end
nole = sortrows([nole zlef],2);
nole = nole(:,1);
zlef = sort(zlef);
nori = sortrows([nori zrig],2);
nori = nori(:,1);
zrig = sort(zrig);
nsur = sortrows([nsur ysur],2);
nsur = nsur(:,1);
ysur = sort(ysur);
nobt = sortrows([nobt ybot],2);
nobt = nobt(:,1);
ybot = sort(ybot);
elle = [];
elri = [];
for ide = 1:nel_h
    if sum(ismember(el2no_h(ide,:),nole)) == 2
        tmp = ismember(el2no_h(ide,:),nole);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz_h(el2no_h(ide,pt1),2) no2yz_h(el2no_h(ide,pt2),2)]);        
        elle = [elle;[ide tmd]];
    end
    if sum(ismember(el2no_h(ide,:),nori)) == 2
        tmp = ismember(el2no_h(ide,:),nori);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz_h(el2no_h(ide,pt1),2) no2yz_h(el2no_h(ide,pt2),2)]);        
        elri = [elri;[ide tmd]];
    end
end
elle = sortrows(elle,2);
role = rho(elle(:,1));
elri = sortrows(elri,2);
rori = rho(elri(:,1));
% - Top boundary
Xtop = 5;
X(nsur) = Xtop * ones(length(nsur),1);
% - Left boundary
X(nole) = BoundaryTM(role,period,zlef,Xtop);
% - Right boundary
X(nori) = BoundaryTM(rori,period,zrig,Xtop);
% - Bottom boundary
X(nobt) = linspace(X(nole(end)),X(nori(end)),length(nobt));

% Assemble the coefficient matrix Dh
% Compute each element matrix and add the contribution to the global matrix
Dh = zeros(nno_h);
for ide = 1:nel_h
    % Get nodes and their coordinates for the element ide
    no_ide = el2no_h(ide,:)';
    yz_ide = no2yz_h(no_ide,:)';
    % Edges
    s1 = yz_ide(:,3) - yz_ide(:,2);
    s2 = yz_ide(:,1) - yz_ide(:,3);
    s3 = yz_ide(:,2) - yz_ide(:,1);
    % Area of the triangle
    Ar_ide = abs(0.5 * (s2(1)*s3(2) - s2(2)*s3(1)));
    % Gradient of test functions
    grad_phi1 = [-s1(2);s1(1)]/(2 * Ar_ide);
    grad_phi2 = [-s2(2);s2(1)]/(2 * Ar_ide);
    grad_phi3 = [-s3(2);s3(1)]/(2 * Ar_ide);
    grad_phi = [grad_phi1 grad_phi2 grad_phi3];
    % Compute all the integrals for this particular element
    Dh_ide = zeros(3);
    for i = 1:3
        for j = 1:3
            if i == j
                Dh_ide(i,j) = (-rho(ide)*(grad_phi(:,i)'*grad_phi(:,j)) - 1i*omega*mu0*1/3) * Ar_ide;
            else
                Dh_ide(i,j) = (-rho(ide)*(grad_phi(:,i)'*grad_phi(:,j))) * Ar_ide;
            end
        end
    end

    Dh(no_ide,no_ide) = Dh(no_ide,no_ide) + Dh_ide;
end

% Known boundary nodes (noex), all nodes (noal), inner unknown nodes (noin)
noex = ([nole;nori;nsur;nobt]);
noal = 1:nno_h;
noin = setdiff(noal,noex);

% Solve the linear system
X(noin) = Dh(noin,noin)\(-Dh(noin,noex)*X(noex));

% Compute transfer functions at the stations
nst = size(sta,1);
ryx = zeros(nst,1);
pyx = zeros(nst,1);
for ids = 1:nst
    itp1 = find(abs(ysur-sta(ids,1)) == min(abs(ysur-sta(ids,1))));
    itp1 = itp1(1);
    if ysur(itp1) > sta(ids,1)
        itp2 = itp1 - 1;
    else
        itp2 = itp1 + 1;
    end
    tmp = find(no2yz_h(:,1) == ysur(itp1));
    ns1 = tmp(abs(no2yz_h(tmp,2)-sta(ids,2)) == min(abs(no2yz_h(tmp,2)-sta(ids,2))));
    tmp = find(no2yz_h(:,1) == ysur(itp2));
    ns2 = tmp(abs(no2yz_h(tmp,2)-sta(ids,2)) == min(abs(no2yz_h(tmp,2)-sta(ids,2))));
    for ide = 1:nel_h
        if sum(ismember([ns1 ns2],el2no_h(ide,:))) == 2
            els = ide;
            ns3 = setdiff(el2no_h(els,:),[ns1 ns2]);
            break;
        end
    end
    a = sqrt((no2yz_h(ns1,1)-no2yz_h(ns3,1))^2 + (no2yz_h(ns1,2)-no2yz_h(ns3,2))^2);
    b = sqrt((no2yz_h(ns2,1)-no2yz_h(ns3,1))^2 + (no2yz_h(ns2,2)-no2yz_h(ns3,2))^2);
    c = sqrt((no2yz_h(ns1,1)-no2yz_h(ns2,1))^2 + (no2yz_h(ns1,2)-no2yz_h(ns2,2))^2);
    C = acos((a^2 + b^2 - c^2)/(2 * a * b));
    h = a*b/c * sin(C);
    % Impedance, apparent resistivity, phase
    Zyx = (X(ns3)-Xtop)/h*rho(els)/Xtop;
    ryx(ids) = 1/omega/mu0 * abs(Zyx).^2;
    pyx(ids) = 180/pi * atan2(imag(Zyx),real(Zyx));
end

end
