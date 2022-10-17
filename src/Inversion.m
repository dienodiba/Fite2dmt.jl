% Dieno Diba 2022
% 2D finite element magnetotelluric inversion
% Inversion

function Inversion(filem0,filedat,filetopo,filestg)

% Import initial model
fid = fopen(filem0);
tmp = fgetl(fid);
nel = str2double(tmp(6:end));
tmp = fgetl(fid);
nno = str2double(tmp(6:end));
tmp = fgetl(fid);
tmp = fgetl(fid);
el2no = zeros(nel,3);
for ide = 1:nel
    tmp = fgetl(fid);
    el2no(ide,1) = str2double(tmp(1:12));
    el2no(ide,2) = str2double(tmp(13:24));
    el2no(ide,3) = str2double(tmp(25:36));
end
tmp = fgetl(fid);
tmp = fgetl(fid);
no2yz = zeros(nno,2);
for idn = 1:nno
    tmp = fgetl(fid);
    no2yz(idn,1) = 1e+3 * str2double(tmp(1:12));        %km to m
    no2yz(idn,2) = 1e+3 * str2double(tmp(13:24));
end
tmp = fgetl(fid);
tmp = fgetl(fid);
rho0 = zeros(nel,1);
for ide = 1:nel
    tmp = fgetl(fid);
    rho0(ide) = str2double(tmp(1:12));
end
fclose(fid);
% Import dataset
fid = fopen(filedat);
tmp = fgetl(fid);
nst = str2double(tmp(13:end));
tmp = fgetl(fid);
nper = str2double(tmp(13:end));
tmp = fgetl(fid);   % yref not needed
tmp = fgetl(fid);
sta = zeros(nst,2);
period = zeros(nper,1);
rhe_o = zeros(nper,nst);
ere_o = zeros(nper,nst);
phe_o = zeros(nper,nst);
epe_o = zeros(nper,nst);
rhh_o = zeros(nper,nst);
erh_o = zeros(nper,nst);
phh_o = zeros(nper,nst);
eph_o = zeros(nper,nst);
for ids = 1:nst
    tmp = fgetl(fid);
    tmp = fgetl(fid);
    sta(ids,1) = 1e+3 * str2double(tmp(8:end));
    tmp = fgetl(fid);
    sta(ids,2) = 1e+3 * str2double(tmp(8:end));
    tmp = fgetl(fid);
    for idp = 1:nper
        tmp = fgetl(fid);
        period(idp) = str2double(tmp(1:12));
        rhe_o(idp,ids) = log10(str2double(tmp(13:24)));
        phe_o(idp,ids) = pi/180*str2double(tmp(37:48));
        if rhe_o(idp,ids) ~= log10(99999)
            ere_o(idp,ids) = str2double(tmp(25:36))/str2double(tmp(13:24))/log(10);
            epe_o(idp,ids) = pi/180*str2double(tmp(49:60));
        else
            ere_o(idp,ids) = 99999;
            epe_o(idp,ids) = 99999;
        end
        rhh_o(idp,ids) = log10(str2double(tmp(61:72)));
        phh_o(idp,ids) = pi/180*str2double(tmp(85:96));
        if rhh_o(idp,ids) ~= log10(99999)
            erh_o(idp,ids) = str2double(tmp(73:84))/str2double(tmp(61:72))/log(10);
            eph_o(idp,ids) = pi/180*str2double(tmp(97:108));
        else
            erh_o(idp,ids) = 99999;
            eph_o(idp,ids) = 99999;
        end
    end
    tmp = fgetl(fid);
end
fclose(fid);
% Import setting
fid = fopen(filestg);
while feof(fid) == 0
    tmp = fgetl(fid);
    sheg = split(tmp);
    if strcmp(sheg(1),'LAMBDA') == 1
        tmp = fgetl(fid);
        lambda = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'DAMPING') == 1
        tmp = fgetl(fid);
        damping = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'MAXITER') == 1
        tmp = fgetl(fid);
        max_iter = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'USE_RTE') == 1
        tmp = fgetl(fid);
        wg_rte = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'USE_PTE') == 1
        tmp = fgetl(fid);
        wg_pte = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'USE_RTM') == 1
        tmp = fgetl(fid);
        wg_rtm = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'USE_PTM') == 1
        tmp = fgetl(fid);
        wg_ptm = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'ERF_ZTE') == 1
        tmp = fgetl(fid);
        ef_zte = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'ERF_ZTM') == 1
        tmp = fgetl(fid);
        ef_ztm = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'FIX_N') == 1
        tmp = fgetl(fid);
        fix_n = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'FIX_RES') == 1
        tmp = fgetl(fid);
        fix_el = str2double(tmp(1:end));
    end
    if strcmp(sheg(1),'OUTFILE') == 1
        tmp = fgetl(fid);
        outroot = tmp(1:end);
    end
end
fclose(fid);
rai = fix_el(1);
% Import topography
topo = importdata(filetopo);
topo = 1e+3 * topo.data;

% Data: 
% - 0.5*log10(rhoa)
% - phase/log(10)
N = nst * nper * 4;     % rte pte rtm ztm
Nh = N / 2;

% Construct vector of observed data considering error floors
dobs = zeros(N,1);
eobs = zeros(N,1);
idd = 0;
for idm = 1:2
    for idp = 1:nper
        for ids = 1:nst
            if idm == 1
                for idj = 1:2
                    idd = idd + 1;
                    if idj == 1
                        dobs(idd) = 0.5 * rhe_o(idp,ids);
                        eobs(idd) = ere_o(idp,ids);
                        if eobs(idd) < 2*ef_zte/log(10)
                            eobs(idd) = 2*ef_zte/log(10);
                        end
                    end
                    if idj == 2
                        dobs(idd) = phe_o(idp,ids) / log(10);
                        eobs(idd) = epe_o(idp,ids);
                        if eobs(idd) < asin(ef_zte)
                            eobs(idd) = asin(ef_zte);
                        end
                    end
                end
            else
                for idj = 1:2
                    idd = idd + 1;
                    if idj == 1
                        dobs(idd) = 0.5 * rhh_o(idp,ids);
                        eobs(idd) = erh_o(idp,ids);
                        if eobs(idd) < 2*ef_ztm/log(10)
                            eobs(idd) = 2*ef_ztm/log(10);
                        end
                    end
                    if idj == 2
                        dobs(idd) = phh_o(idp,ids) / log(10);
                        eobs(idd) = eph_o(idp,ids);
                        if eobs(idd) < asin(ef_ztm)
                            eobs(idd) = asin(ef_ztm);
                        end
                    end
                end
            end
        end
    end
end

% Diagonal matrix of error inverse for weighting
W = zeros(N);
for idN = 1:N
    if eobs(idN) ~= 99999
        W(idN,idN) = (eobs(idN))^-2;
    end
end

% Count the number of each data component
T = zeros(N,1);
idd1 = 0;
idd4 = 0;
idN = 0;
for idm = 1:2
    for idp = 1:nper
        for ids = 1:nst
            if idm == 1
                for idj = 1:2
                    idN = idN + 1;
                    if idj == 1
                        T(idN) = wg_rte;
                        if eobs(idN) ~= 99999
                            idd1 = idd1 + wg_rte;
                        end
                    else
                        T(idN) = wg_pte;
                        if eobs(idN) ~= 99999
                            idd1 = idd1 + wg_pte;
                        end
                    end
                end
            else
                for idj = 1:2
                    idN = idN + 1;
                    if idj == 1
                        T(idN) = wg_rtm;
                        if eobs(idN) ~= 99999
                            idd4 = idd4 + wg_rtm;
                        end
                    else
                        T(idN) = wg_pte;
                        if eobs(idN) ~= 99999
                            idd4 = idd4 + wg_ptm;
                        end
                    end
                end
            end
        end
    end
end

% Number of parameter M is the number of free resistivity elements
M = 0;
m0 = zeros(nel,1);
mbin = zeros(nel,1);
m2ee = zeros(nel,1);
Area = zeros(nel,1);
for ide = 1:nel
    if ismember(rho0(ide),fix_el) == 0
        M = M + 1;
        mbin(ide) = 1;      % 1 for free element, 0 for fixed element
        m0(ide) = log10(rho0(ide));
        m2ee(M) = ide;
    end
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
    no_ide = el2no(ide,:)';
    yz_ide = no2yz(no_ide,:)';
    % Edges
    s1_ide = yz_ide(:,3) - yz_ide(:,2);
    s2_ide = yz_ide(:,1) - yz_ide(:,3);
    s3_ide = yz_ide(:,2) - yz_ide(:,1);
    % Area of the element
    Area(ide)= abs(0.5 * (s2_ide(1)*s3_ide(2) - s2_ide(2)*s3_ide(1)));
end
m0 = nonzeros(m0);
m2ee = nonzeros(m2ee);

% Assign the first model
m = m0;
rho = rho0;

% Roughness: penalizing elements sharing an edge, index only, no weighting
L = zeros(nel);
for ide = 1:nel
    L(ide,ide) = -1;
    tmp = el2no(ide,:);
    sheg = zeros(3,1);
    count = 0;
    for ideg = 1:3
        idv1 = ideg;
        if ideg < 3
            idv2 = ideg + 1;
        else            
            idv2 = 1;
        end
        va1 = find(el2no(:,1) == tmp(idv1));
        vb1 = find(el2no(:,2) == tmp(idv1));
        vc1 = find(el2no(:,3) == tmp(idv1));
        vu1 = unique([va1;vb1;vc1]);
        va2 = find(el2no(:,1) == tmp(idv2));
        vb2 = find(el2no(:,2) == tmp(idv2));
        vc2 = find(el2no(:,3) == tmp(idv2));
        vu2 = unique([va2;vb2;vc2]);
        if isempty(setdiff(intersect(vu1,vu2),ide)) == 0
            sheg(ideg) = setdiff(intersect(vu1,vu2),ide);
            if mbin(sheg(ideg)) == 1
                count = count + 1;
            end
        end
    end
    L(ide,nonzeros(sheg)) = 1/count;
end
L = L(mbin == 1,mbin == 1);
K = sparse(L') * sparse(L);

Phid1 = zeros(max_iter,1);
Phid4 = zeros(max_iter,1);
RMSd1 = zeros(max_iter,1);
RMSd4 = zeros(max_iter,1);
Phi = zeros(max_iter,1);
Psi = zeros(max_iter,1);
Stab = zeros(max_iter,1);
RMS = zeros(max_iter,1);

mpre = zeros(M,1);
dpre = zeros(N,1);
frho = zeros(nel,1);
miscom = 0;
epoch = 0;

ymin = min(no2yz(:,1));
ymax = max(no2yz(:,1));
zmin = min(no2yz(:,2));
zmax = max(no2yz(:,2));

% Boundary nodes and elements for TE mode computation
zlef_e = []; zrig_e = [];
ytop_e = []; ybot_e = [];
nole_e = []; nori_e = [];
notp_e = []; nobt_e = [];
ysur_e = []; nsur_e = [];
for idno = 1:nno
    if no2yz(idno,1) == ymin
        zlef_e = [zlef_e;no2yz(idno,2)];
        nole_e = [nole_e;idno];
    elseif no2yz(idno,1) == ymax
        zrig_e = [zrig_e;no2yz(idno,2)];
        nori_e = [nori_e;idno];
    elseif no2yz(idno,2) == zmin
        ytop_e = [ytop_e;no2yz(idno,1)];
        notp_e = [notp_e;idno];
    elseif no2yz(idno,2) == zmax
        ybot_e = [ybot_e;no2yz(idno,1)];
        nobt_e = [nobt_e;idno];
    end
    if abs(no2yz(idno,2) - interp1(topo(:,1),topo(:,2),no2yz(idno,1))) < minedge/2
        ysur_e = [ysur_e;no2yz(idno,1)];
        nsur_e = [nsur_e;idno];
    end
end
nole_e = sortrows([nole_e zlef_e],2);
nole_e = nole_e(:,1);
zlef_e = sort(zlef_e);
nori_e = sortrows([nori_e zrig_e],2);
nori_e = nori_e(:,1);
zrig_e = sort(zrig_e);
notp_e = sortrows([notp_e ytop_e],2);
notp_e = notp_e(:,1);
ytop_e = sort(ytop_e);
nobt_e = sortrows([nobt_e ybot_e],2);
nobt_e = nobt_e(:,1);
ybot_e = sort(ybot_e);
elle_e = [];
elri_e = [];
nsur_e = sortrows([nsur_e ysur_e],2);
nsur_e = nsur_e(:,1);
ysur_e = sort(ysur_e);
for ide = 1:nel
    if sum(ismember(el2no(ide,:),nole_e)) == 2
        tmp = ismember(el2no(ide,:),nole_e);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz(el2no(ide,pt1),2) no2yz(el2no(ide,pt2),2)]);        
        elle_e = [elle_e;[ide tmd]];
    end
    if sum(ismember(el2no(ide,:),nori_e)) == 2
        tmp = ismember(el2no(ide,:),nori_e);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz(el2no(ide,pt1),2) no2yz(el2no(ide,pt2),2)]);        
        elri_e = [elri_e;[ide tmd]];
    end
end
elle_e = sortrows(elle_e,2);
elri_e = sortrows(elri_e,2);

% For TM, remove elements and nodes in the air half-space
rmel = [];
eh2ee = zeros(nel,1);
count = 0;
for ide = 1:nel
    tmp = no2yz(el2no(ide,:),:);
    if mean(tmp(:,2)) < interp1(topo(:,1),topo(:,2),mean(tmp(:,1)))
        rmel = [rmel;ide];        
    else
        count = count + 1;
        eh2ee(count) = ide;
    end
end
eh2ee = nonzeros(eh2ee);
syel = setdiff(1:nel,rmel);
el2no_h = el2no(syel,:);
rho_h = rho(syel,:);
nel_h = length(el2no_h);
rmno = [];
tmp = 0;
nh2ne = zeros(nno);
for idn = 1:nno
    if ismember(idn,unique(el2no_h(:))) == 0
        rmno = [rmno;idn];
    else
        tmp = tmp + 1;
        nh2ne(tmp) = idn;
    end
end
nh2ne = nonzeros(nh2ne);
no2yz_h = no2yz(setdiff(1:nno,rmno),:);
nno_h = length(no2yz_h);
for ide = 1:nel_h
    for idv = 1:3
        el2no_h(ide,idv) = find(nh2ne == el2no_h(ide,idv));
    end
end
% Boundary nodes and elements for TM mode computation
zlef_h = [];
zrig_h = [];
ysur_h = [];
ybot_h = [];
nole_h = [];
nori_h = [];
nsur_h = [];
nobt_h = [];
for idno = 1:nno_h
    if no2yz_h(idno,1) == ymin
        zlef_h = [zlef_h;no2yz_h(idno,2)];
        nole_h = [nole_h;idno];
    elseif no2yz_h(idno,1) == ymax
        zrig_h = [zrig_h;no2yz_h(idno,2)];
        nori_h = [nori_h;idno];
    elseif no2yz_h(idno,2) == zmax
        ybot_h = [ybot_h;no2yz_h(idno,1)];
        nobt_h = [nobt_h;idno];
    elseif abs(no2yz_h(idno,2) - interp1(topo(:,1),topo(:,2),no2yz_h(idno,1))) < minedge/2
        ysur_h = [ysur_h;no2yz_h(idno,1)];
        nsur_h = [nsur_h;idno];
    end
end
nole_h = sortrows([nole_h zlef_h],2);
nole_h = nole_h(:,1);
zlef_h = sort(zlef_h);
nori_h = sortrows([nori_h zrig_h],2);
nori_h = nori_h(:,1);
zrig_h = sort(zrig_h);
nsur_h = sortrows([nsur_h ysur_h],2);
nsur_h = nsur_h(:,1);
ysur_h = sort(ysur_h);
nobt_h = sortrows([nobt_h ybot_h],2);
nobt_h = nobt_h(:,1);
ybot_h = sort(ybot_h);
elle_h = [];
elri_h = [];
for ide = 1:nel_h
    if sum(ismember(el2no_h(ide,:),nole_h)) == 2
        tmp = ismember(el2no_h(ide,:),nole_h);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz_h(el2no_h(ide,pt1),2) no2yz_h(el2no_h(ide,pt2),2)]);        
        elle_h = [elle_h;[ide tmd]];
    end
    if sum(ismember(el2no_h(ide,:),nori_h)) == 2
        tmp = ismember(el2no_h(ide,:),nori_h);
        pt1 = find(tmp,1,'first');
        pt2 = find(tmp,1,'last');
        tmd = mean([no2yz_h(el2no_h(ide,pt1),2) no2yz_h(el2no_h(ide,pt2),2)]);        
        elri_h = [elri_h;[ide tmd]];
    end
end
elle_h = sortrows(elle_h,2);
elri_h = sortrows(elri_h,2);
m2eh = zeros(nel_h,1);
tmp = 0;
for ide = 1:nel_h
    if ismember(rho_h(ide),fix_el) == 0
        tmp = tmp + 1;
        m2eh(tmp) = ide;
    end
end
m2eh = nonzeros(m2eh);

% Gauss Newton Iteration
mu0 = 4*pi*1e-7;
tStart = cputime;
for itr = 1:max_iter
    epoch = epoch + 1;
    disp(strcat("Iteration ",num2str(itr)))    
    
    % ========= TE mode =========
    dcale1 = zeros(nper,nst);
    Aheg1 = zeros(nper*nst,M);
    for idp = 1:nper
        omega = 2*pi/period(idp);
        % Define ve: field solution vector
        ve = zeros(nno,1);
        % Set the values at the boundary
        % - Top boundary
        ve_top = 1;
        ve(notp_e) = ve_top * ones(length(notp_e),1);
        % - Left boundary
        role_e = rho(elle_e(:,1));
        ve(nole_e) = BoundaryTE(role_e,period(idp),zlef_e,ve_top);
        % - Right boundary
        rori_e = rho(elri_e(:,1));
        ve(nori_e) = BoundaryTE(rori_e,period(idp),zrig_e,ve_top);
        % - Bottom boundary
        ve(nobt_e) = linspace(ve(nole_e(end)),ve(nori_e(end)),length(nobt_e));

        % Assemble coefficient matrix De. It is sparse, symmeteric, banded
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
            % Gradient of test functions
            grad_phi1 = [-s1_ide(2);s1_ide(1)]/(2 * Area(ide));
            grad_phi2 = [-s2_ide(2);s2_ide(1)]/(2 * Area(ide));
            grad_phi3 = [-s3_ide(2);s3_ide(1)]/(2 * Area(ide));
            grad_phi = [grad_phi1 grad_phi2 grad_phi3];
            % Compute all the integrals for this particular element
            De_ide = zeros(3);
            for i = 1:3
                for j = 1:3
                    if i == j
                        De_ide(i,j) = (-(grad_phi(:,i)'*grad_phi(:,j)) - 1i*omega*mu0/rho(ide)*1/3) * Area(ide);
                    else
                        De_ide(i,j) = -(grad_phi(:,i)'*grad_phi(:,j)) * Area(ide);
                    end
                end
            end
            De(no_ide,no_ide) = De(no_ide,no_ide) + De_ide;
        end

        % Known boundary nodes (noex), all nodes (noal), inner unknown nodes (noin)
        noex_e = ([nole_e;nori_e;notp_e;nobt_e]);
        noal_e = 1:nno;
        noin_e = setdiff(noal_e,noex_e);

        % Solve the linear system
        dcDe = decomposition(De(noin_e,noin_e),'banded');
        ve(noin_e) = dcDe\(-De(noin_e,noex_e)*ve(noex_e));
        
        % ae1, be1: coefficients to interpolate and/or differentiate fields
        % ce1: required for pseudo-forward computation
        ce1 = complex(zeros(nno,nst));
        for ids = 1:nst
            itp1 = find(abs(ysur_e-sta(ids,1)) == min(abs(ysur_e-sta(ids,1))));
            itp1 = itp1(1);
            if ysur_e(itp1) > sta(ids,1)
                itp2 = itp1 - 1;
            else
                itp2 = itp1 + 1;
            end
            tmp = find(no2yz(:,1) == ysur_e(itp1));
            ns1 = tmp(abs(no2yz(tmp,2)-sta(ids,2)) == min(abs(no2yz(tmp,2)-sta(ids,2))));
            tmp = find(no2yz(:,1) == ysur_e(itp2));
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
            % rho-app and phase
            ae1 = zeros(3,1);
            be1 = zeros(3,1);
            ae1(1) = 0.5;
            ae1(2) = 0.5;
            if d > c    % obtuse triangle beneath the station ids
                be1(1) = 0.5;
                be1(2) = 0.5;
            else        % acute triangle beneath the station ids
                be1(1) = -1 * -1/1i/omega/mu0/h * (1 - d/c);
                be1(2) = -1 * -1/1i/omega/mu0/h * d/c;
            end
            be1(3) = -1/1i/omega/mu0/h;
            ce1([ns1 ns2 nbsa],ids) = ae1/(ae1.'*ve([ns1 ns2 nbsa])) - be1/(be1.'*ve([ns1 ns2 nbsa]));
            dcale1(idp,ids) = log10(sqrt(1/omega/mu0)*(ae1.'*ve([ns1 ns2 nbsa]))/(be1.'*ve([ns1 ns2 nbsa])));
        end
        
        % Pseudo-forward computation
        ue1 = complex(zeros(nno,nst));
        for ids = 1:nst
            ue1(noin_e,ids) = dcDe\ce1(noin_e,ids);
        end
        
        % Assemble Jacobian
        Ahe1 = complex(zeros(nst,M));
        for idM = 1:M
            kiw = el2no(m2ee(idM),:);
            pDe = 1i*omega*mu0*Area(m2ee(idM))/3/rho(m2ee(idM))^2;
            Ahe1(:,idM) = -ue1(kiw,:).'*(pDe*ve(kiw)) * rho(m2ee(idM));
        end
        Aheg1((idp-1)*nst+1:idp*nst,:) = Ahe1;
    end
    
    Aheg = [];
    for idG = 1:nper*nst
        Aheg = [Aheg;Aheg1(idG,:)];
    end
    
    dcaleg = zeros(nper*nst*2,1);
    idd = 0;
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:2
                idd = idd + 1;
                if idj == 1
                    dcaleg(idd) = real(dcale1(idp,ids));
                end
                if idj == 2
                    dcaleg(idd) = imag(dcale1(idp,ids));
                end
            end
        end
    end
    
    % ========= TM mode =========
    rho_h = rho(syel,:);
    role_h = rho_h(elle_h(:,1));
    rori_h = rho_h(elri_h(:,1));    
    dcalh = zeros(nper,nst);
    Ahhg = zeros(nper*nst,M);
    for idp = 1:nper
        omega = 2*pi/period(idp);
        % Define X: field solution vector
        vh = zeros(nno_h,1);
        % Set the values at the boundary
        % Top boundary
        vh_top = 1;
        vh(nsur_h) = vh_top * ones(length(nsur_h),1);
        % Left boundary
        vh(nole_h) = BoundaryTM(role_h,period(idp),zlef_h,vh_top);
        % Right boundary
        vh(nori_h) = BoundaryTM(rori_h,period(idp),zrig_h,vh_top);
        % Bottom boundary
        vh(nobt_h) = linspace(vh(nole_h(end)),vh(nori_h(end)),length(nobt_h));

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
                        Dh_ide(i,j) = (-rho_h(ide)*(grad_phi(:,i)'*grad_phi(:,j)) - 1i*omega*mu0*1/3) * Ar_ide;
                    else
                        Dh_ide(i,j) = (-rho_h(ide)*(grad_phi(:,i)'*grad_phi(:,j))) * Ar_ide;
                    end
                end
            end
            Dh(no_ide,no_ide) = Dh(no_ide,no_ide) + Dh_ide;
        end

        % Known boundary nodes (noex), all nodes (noal), inner unknown nodes (noin)
        noex_h = ([nole_h;nori_h;nsur_h;nobt_h]);
        noal_h = 1:nno_h;
        noin_h = setdiff(noal_h,noex_h);

        % Solve the linear system
        dcDh = decomposition(Dh(noin_h,noin_h),'banded');
        vh(noin_h) = dcDh\(-Dh(noin_h,noex_h)*vh(noex_h));
        
        % ah, bh: coefficients to interpolate and/or differentiate fields
        % ch: required for pseudo-forward computation
        ch = complex(zeros(nno_h,1));
        elsur = zeros(nst,1);
        for ids = 1:nst
            itp1 = find(abs(ysur_h-sta(ids,1)) == min(abs(ysur_h-sta(ids,1))));
            itp1 = itp1(1);
            if ysur_h(itp1) > sta(ids,1)
                itp2 = itp1 - 1;
            else
                itp2 = itp1 + 1;
            end
            tmp = find(no2yz_h(:,1) == ysur_h(itp1));
            ns1 = tmp(abs(no2yz_h(tmp,2)-sta(ids,2)) == min(abs(no2yz_h(tmp,2)-sta(ids,2))));
            tmp = find(no2yz_h(:,1) == ysur_h(itp2));
            ns2 = tmp(abs(no2yz_h(tmp,2)-sta(ids,2)) == min(abs(no2yz_h(tmp,2)-sta(ids,2))));
            for ide = 1:nel_h
                if sum(ismember([ns1 ns2],el2no_h(ide,:))) == 2
                    els = ide;
                    elsur(ids) = els;
                    nbs = setdiff(el2no_h(els,:),[ns1 ns2]);
                    break;
                end
            end
            a = sqrt((no2yz_h(ns1,1)-no2yz_h(nbs,1))^2 + (no2yz_h(ns1,2)-no2yz_h(nbs,2))^2);
            b = sqrt((no2yz_h(ns2,1)-no2yz_h(nbs,1))^2 + (no2yz_h(ns2,2)-no2yz_h(nbs,2))^2);
            c = sqrt((no2yz_h(ns1,1)-no2yz_h(ns2,1))^2 + (no2yz_h(ns1,2)-no2yz_h(ns2,2))^2);
            C = acos((a^2 + b^2 - c^2)/(2 * a * b));
            h = a*b/c * sin(C);
            ah = zeros(3,1);
            bh = zeros(3,1);
            ah(1) = rho_h(els)/h * -0.5;
            ah(2) = rho_h(els)/h * -0.5;
            ah(3) = rho_h(els)/h;
            bh(1) = 0.5;
            bh(2) = 0.5;
            dcalh(idp,ids) = log10(sqrt(1/omega/mu0) * (ah.'*vh([ns1 ns2 nbs])) / (bh.'*vh([ns1 ns2 nbs])));
            ch([ns1 ns2 nbs],ids) = ah/(ah.'*vh([ns1 ns2 nbs])) - bh/(bh.'*vh([ns1 ns2 nbs]));
        end
        
        %Pseduo-forward computation
        uh = complex(nno_h,nst);
        for ids = 1:nst
            uh(noin_h,ids) = dcDh\ch(noin_h,ids);
        end

        % Assemble Jacobian
        Ahh = complex(zeros(nst,M));
        for idM = 1:M
            kiw = el2no_h(m2eh(idM),:);
            yz = no2yz(el2no(m2ee(idM),:),:)';
            % Edges
            s1 = yz(:,3) - yz(:,2);
            s2 = yz(:,1) - yz(:,3);
            s3 = yz(:,2) - yz(:,1);
            % Area of the triangle
            Ar = abs(0.5 * (s2(1)*s3(2) - s2(2)*s3(1)));
            % Gradient of test functions
            grad_phi1 = [-s1(2);s1(1)]/(2 * Ar);
            grad_phi2 = [-s2(2);s2(1)]/(2 * Ar);
            grad_phi3 = [-s3(2);s3(1)]/(2 * Ar);
            grad_phi = [grad_phi1 grad_phi2 grad_phi3];
            pDh = zeros(3);            
            for i = 1:3
                for j = 1:3
                    pDh(i,j) = -1 * (grad_phi(:,i)'*grad_phi(:,j)) * Ar;
                end
            end            
            Ahh(:,idM) = -uh(kiw,:).'*(pDh*vh(kiw)) * rho_h(m2eh(idM));
        end
        % Delta Jacobian for elements under the stations
        % Small perturbation
        ptb = 1e-3;
        rptb = rho;
        rptb(eh2ee(elsur)) = rptb(eh2ee(elsur)) + ptb;
        [ryx_ptb,pyx_ptb] = ForwardTM(el2no,no2yz*1e-3,rptb,period(idp),sta*1e-3,topo*1e-3);
        d_ptb = 0.5*log10(ryx_ptb) + 1i* (pyx_ptb*pi/180/log(10));
        for ids = 1:nst
            Ahh(ids,m2eh == elsur(ids)) = (d_ptb(ids) - dcalh(idp,ids))/(log10(rptb(eh2ee(elsur(ids)))) - log10(rho(eh2ee(elsur(ids)))));
        end
        Ahhg((idp-1)*nst+1:idp*nst,:) = Ahh;
    end
    
    dcalhg = zeros(nst*nper*2,1);
    idd = 0;
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:2
                idd = idd + 1;
                if idj == 1
                    dcalhg(idd) = real(dcalh(idp,ids));
                else
                    dcalhg(idd) = imag(dcalh(idp,ids));
                end
            end
        end
    end
    
    %Merging the calculated data from TE mode and TM mode
    dcal = [dcaleg;dcalhg];
    %Merging the complex Jacobian
    Ah = [Aheg;Ahhg];
    %Convert Complex Jacobian into Standard Real Jacobian
    A = zeros(N,M);
    for idM = 1:M
        idN = 0;
        for idNh = 1:Nh
            for idj = 1:2
                idN = idN + 1;
                if idj == 1
                    A(idN,idM) = real(Ah(idNh,idM));
                else
                    A(idN,idM) = imag(Ah(idNh,idM));
                end
            end
        end
    end
    
    % Calculate the objective function
    idN = 0;
    for idm = 1:2
        for idp = 1:nper
            for ids = 1:nst
                if idm == 1
                    for idj = 1:2
                        idN = idN + 1;
                        if idj == 1 || idj == 2
                            Phid1(itr) = Phid1(itr) + (dcal(idN)-dobs(idN))^2 * W(idN,idN) * T(idN);
                        end
                    end
                else
                    for idj = 1:2
                        idN = idN + 1;
                        if idj == 1 || idj == 2
                            Phid4(itr) = Phid4(itr) + (dcal(idN)-dobs(idN))^2 * W(idN,idN) * T(idN);
                        end
                    end
                end
            end
        end
    end
    %Total data misfit
    Phi(itr) = Phid1(itr) + Phid4(itr);
    %Model roughness
    Stab(itr) = lambda * m.'*K*m;
    %Total objective function
    Psi(itr) = Phi(itr) + Stab(itr);       
    % RMS combined
    RMS(itr) = sqrt((Phid1(itr)+Phid4(itr))/(idd1+idd4));
    % RMS of each dataset
    RMSd1(itr) = sqrt(Phid1(itr)/idd1);
    RMSd4(itr) = sqrt(Phid4(itr)/idd4);
    
    % Termination criteria:
    % 1. Objective function increases (inversion diverges)
    if itr > 1 && Psi(itr) - Psi(itr-1) >= 0
        disp("Objective function is no longer decreasing")
        disp("The inversion is terminated here")
        break
    % 2. Objective function decreases insignificantly (converges)
    elseif itr > 1 && Psi(itr-1) - Psi(itr) < 0.01*Psi(itr)
        dpre = dcal;
        frho = rho;
        mpre = m;
        miscom = Phid1(itr) + Phid4(itr);        
        disp("The inversion converged to a solution")
        disp("Finished")
        break
    end
    
    % If the process continues, save data, model, and misfit
    dpre = dcal;
    frho = rho;
    mpre = m;
    miscom = Phid1(itr) + Phid4(itr);
    
    %Update the model when itr < max_iter
    if itr < max_iter
        % Assemble Gradient and Hessian
        Grd = -2*(A'*(T.*sparse(W))*(dobs-dcal)) + 2*lambda*K*m;
        Hes = 2*(A'*(T.*sparse(W)*A)) + 2*lambda*K;
        %Levenberg-Marquardt damping
        epl = damping*Psi(itr); 
        %Damped Hessian
        Hes_dmp = (Hes + epl*speye(M));      
        dcdmpH = decomposition(Hes_dmp,'lu');
        m = m  - dcdmpH\Grd;
        for idM = 1:M
            rho(m2ee(idM)) = 10^m(idM);
        end
    else
        % when itr = max_iter, not necessary to update the model
        break
    end
end
tEnd = cputime - tStart;

% Export inversion report
fid = fopen(strcat(outroot,'_log.txt'),'w');
fprintf(fid,"Inversion Report\n");
fprintf(fid,strcat(string(datetime),"\n\n"));
fprintf(fid,"===============================\n");
fprintf(fid,"Input:\n");
fprintf(fid,"%s\n",filedat);
fprintf(fid,"%s\n",filem0);
fprintf(fid,"%s\n",filestg);
fprintf(fid,"%s\n",filetopo);
fprintf(fid,"===============================\n");
fprintf(fid,"Final combined RMS = %4.2f\n",RMS(find(RMS,1,'last')));
fprintf(fid,"Final RMS_TE = %4.2f\n",RMSd1(find(RMSd1,1,'last')));
fprintf(fid,"Final RMS_TM = %4.2f\n",RMSd4(find(RMSd4,1,'last')));
fprintf(fid,"Final data misfit = %3.3e\n",miscom);
fprintf(fid,"Final model roughness = %3.3e\n",mpre'*K*mpre);
fprintf(fid,"===============================\n");
fprintf(fid,"# of iteration = %i\n",epoch);
fprintf(fid,"# of data = %i\n",idd1+idd4);
fprintf(fid,"# of model parameter = %i\n",M);
fprintf(fid,"Total CPU Time (s) = %f\n",tEnd);
fprintf(fid,"===============================\n");
fprintf(fid,"Evolution of Objective Function\n");
for idi = 1:epoch
    fprintf(fid,"Iteration %i\n",idi);
    fprintf(fid,"Objective function = %3.2e\n",Psi(idi));
    fprintf(fid,"Data misfit = %3.2e\n",Phi(idi));
    fprintf(fid,"Model roughness = %3.2e\n",Stab(idi)/lambda);
    fprintf(fid,"Combined RMS = %3.2e\n",RMS(idi));
end
fclose(fid);

% Export the final model
fid = fopen(strcat(outroot,'_res.txt'),'w');
fprintf(fid,'nel = %10i\n',nel);
fprintf(fid,'nno = %10i\n\n',nno);
fprintf(fid,'EL2NO\n');
for ide = 1:nel
    fprintf(fid,'%12i%12i%12i\n',el2no(ide,1),el2no(ide,2),el2no(ide,3));
end
fprintf(fid,'\n');
fprintf(fid,'NO2YZ\n');
for idn = 1:nno    
    fprintf(fid,'%12.4e%12.4e\n',no2yz(idn,1),no2yz(idn,2));
end
fprintf(fid,'\n');
fprintf(fid,'RESISTIVITY\n');
for ide = 1:nel    
    fprintf(fid,'%12.4e\n',frho(ide));
end

end
