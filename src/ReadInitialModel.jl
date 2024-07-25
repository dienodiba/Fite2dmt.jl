# Dieno Diba 2024
# FITE2DMT
# Read initial resistivity model from input file

function ReadInitialModel(filem0,fix_el,interp_topo)

    drs = readdlm(filem0)
    nel = drs[1,3]
    nno = drs[2,3]
    tmp = only(findall(x->x=="EL2NO",drs[:,1]))+1
    el2no = Int.(drs[tmp:tmp+nel-1,1:3])
    tmp = only(findall(x->x=="NO2YZ",drs[:,1]))+1
    no2yz = 1e+3 * float.(drs[tmp:tmp+nno-1,1:2])
    tmp = only(findall(x->x=="RESISTIVITY",drs[:,1]))+1
    rho0 = float.(drs[tmp:tmp+nel-1,1])

    # Number of parameter M is the number of free resistivity elements
    M = 0
    m0 = zeros(Float64,nel)
    mbin = zeros(Int64,nel)
    m2ee = zeros(Int64,nel)
    Area = zeros(Float64,nel)
    minedge = 0.0
    for ide = 1:nel
        if in(rho0[ide],fix_el) == 0
            M += 1
            mbin[ide] = 1      # 1 for free element, 0 for fixed element
            m0[ide] = log10(rho0[ide])
            m2ee[M] = ide
        end
        edge1 = sqrt((no2yz[el2no[ide,3],1] - no2yz[el2no[ide,2],1])^2 + (no2yz[el2no[ide,3],2] - no2yz[el2no[ide,2],2])^2)
        edge2 = sqrt((no2yz[el2no[ide,1],1] - no2yz[el2no[ide,3],1])^2 + (no2yz[el2no[ide,1],2] - no2yz[el2no[ide,3],2])^2)
        edge3 = sqrt((no2yz[el2no[ide,2],1] - no2yz[el2no[ide,1],1])^2 + (no2yz[el2no[ide,2],2] - no2yz[el2no[ide,1],2])^2)
        tmp = min(edge1,edge2,edge3)
        if ide == 1
            minedge = tmp
        else
            if tmp < minedge
                minedge = tmp
            end
        end
        no_ide = el2no[ide,:]
        yz_ide = no2yz[no_ide,:]
        # Edges
        s1_ide = yz_ide[3,:] - yz_ide[2,:]
        s2_ide = yz_ide[1,:] - yz_ide[3,:]
        s3_ide = yz_ide[2,:] - yz_ide[1,:]
        # Area of the element
        Area[ide]= abs(0.5 * (s2_ide[1]*s3_ide[2] - s2_ide[2]*s3_ide[1]))
    end
    m0 = filter(x->x>0,m0)
    m2ee = filter(x->x>0,m2ee)
    # Assign the first model
    m = copy(m0)
    rho = copy(rho0)
    ymin = minimum(no2yz[:,1])
    ymax = maximum(no2yz[:,1])
    zmin = minimum(no2yz[:,2])
    zmax = maximum(no2yz[:,2])

    # For TM, remove elements and nodes in the air half-space
    rmel = []
    eh2ee = zeros(Int64,nel,1)
    count = 0
    for ide = 1:nel
        tmp = no2yz[el2no[ide,:],:]
        if mean(tmp[:,2]) < interp_topo(mean(tmp[:,1]))
            rmel = vcat(rmel,ide)
        else
            count += 1
            eh2ee[count] = ide
        end
    end
    eh2ee = eh2ee[(eh2ee[:,1] .!= 0),:]
    syel = setdiff(1:nel,rmel)
    el2no_h = el2no[syel,:]
    rho_h = vec(rho[syel,:])
    nel_h = size(el2no_h,1)
    rmno = []
    nh2ne = zeros(Int64,nno)
    tmp = 0
    tmd = unique(el2no_h)
    for idn = 1:nno
        if in(idn,tmd) == 0
            rmno = vcat(rmno,idn)
        else
            tmp += 1
            nh2ne[tmp] = idn
        end
    end
    rmno = filter(x->x!=0,rmno)
    nh2ne = filter(x->x!=0,nh2ne)
    no2yz_h = no2yz[setdiff(1:nno,rmno),:]
    nno_h = size(no2yz_h,1)
    for ide = 1:nel_h
        for idv = 1:3
            el2no_h[ide,idv] = only(findall(nh2ne .== el2no_h[ide,idv]))
        end
    end
    m2eh = zeros(Int64,nel_h)
    tmd = 0
    for ide = 1:nel_h
        if in(rho_h[ide],fix_el) == 0
            tmd += 1
            m2eh[tmd] = ide
        end
    end
    m2eh = filter(x->x!=0,m2eh)

    return nel,nno,el2no,no2yz,M,mbin,m2ee,Area,minedge,m,m0,rho,ymin,ymax,zmin,zmax,eh2ee,syel,el2no_h,no2yz_h,nel_h,nno_h,rho_h,m2eh
end
