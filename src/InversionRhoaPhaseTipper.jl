# Dieno Diba 2024
# FITE2DMT
# Inversion of rhoa phase tipper

function InversionRhoaPhaseTipper(filedat::String, filem0::String, filetopo::String, filestg::String)

    println("=================================================================")
    println("FITE2DMT: Finite Triangular Elements for Two-Dimensional MT")
    println("Dieno Diba, 2023")
    println("=================================================================")

    println("Reading input files and initializing ...")
    
    # Read inversion setting parameters
    lambda,max_iter,outroot,fix_el,damping,hor2ver,static,nodata,wg_rte,wg_pte,wg_rtm,wg_ptm,wg_tyr,wg_tyi,ef_zte,ef_ztm,ef_tzy = ReadInversionSettingParametersRhoaPhaseTipper(filestg)

    # Read observation data of rhoa, phase, tipper
    nst,nper,sta,period,rhe_o,ere_o,phe_o,epe_o,rhh_o,erh_o,phh_o,eph_o,rtz_o,ert_o,itz_o,eit_o,dobs,eobs,N,W,V,T,nnzN,idd1,idd2,idd4 = ReadObservedDataRhoaPhaseTipper(filedat,nodata,wg_rte,wg_pte,wg_rtm,wg_ptm,wg_tyr,wg_tyi,ef_zte,ef_ztm,ef_tzy)

    # Read topography file
    topo = readdlm(filetopo)
    topo = 1e+3 * float.(topo[2:end,:])
    interp_topo = Interpolate(topo[:,1],topo[:,2])

    # Read initial model
    nel,nno,el2no,no2yz,M,mbin,m2ee,Area,minedge,m,m0,rho,ymin,ymax,zmin,zmax,eh2ee,syel,el2no_h,no2yz_h,nel_h,nno_h,rho_h,m2eh = ReadInitialModel(filem0,fix_el,interp_topo)

    # Find boundary nodes for TE mode computation
    zlef_e,zrig_e,ytop_e,ybot_e,ysur_e,nole_e,nori_e,notp_e,nobt_e,nsur_e = FindBoundaryNodesTE(interp_topo,nno,no2yz,ymin,ymax,zmin,zmax,minedge)

    # Find boundary elements for TE mode computation
    elle_e,elri_e,elsur_e = FindBoundaryElementsTE(nel,el2no,no2yz,nole_e,nori_e,nsur_e)

    # Find boundary nodes for TM mode computation
    zlef_h,zrig_h,ysur_h,ybot_h,nole_h,nori_h,nsur_h,nobt_h = FindBoundaryNodesTM(interp_topo,nno_h,no2yz_h,ymin,ymax,zmax,minedge)

    # Find boundary elements for TM mode computation
    elle_h,elri_h,elsur_h = FindBoundaryElementsTM(nel_h,el2no_h,no2yz_h,nole_h,nori_h,nsur_h)

    # Find the elements and nodes beneath stations for TE mode
    ns1e,ns2e,nbse = FindNodesBeneathStationsTE(nst,ysur_e,sta,no2yz,el2no,elsur_e,interp_topo)

    # Find the elements and nodes beneath stations for TM mode
    ns1h,ns2h,nbsh,elsta,tmst = FindNodesbeneathStationsTM(nst,ysur_h,sta,no2yz_h,el2no_h,elsur_h,m2eh,M)

    # Define the boundary nodes (noex), all nodes (noal), interior nodes (noin) for TE and TM modes
    noex_e = vcat(nole_e,nori_e,notp_e,nobt_e)
    noal_e = 1:nno
    noin_e = setdiff(noal_e,noex_e)
    noex_h = vcat(nole_h,nori_h,nsur_h,nobt_h)
    noal_h = 1:nno_h
    noin_h = setdiff(noal_h,noex_h)

    # Construct stiffness matrix TE mode
    De = ConstructStiffnessMatrixTE(nel,nno,el2no,no2yz,Area)

    println("Constructing roughness matrix ...")
    K = ConstructRoughnessMatrix(nel,no2yz,el2no,static,hor2ver,elsta,mbin)
    invK = Threads.@spawn inv(lambda*Matrix(K) + I(M)*damping)

    # Gauss-Newton Iteration
    println("Minimizing objective function")
    Phid1 = zeros(Float64,max_iter)
    Phid2 = zeros(Float64,max_iter)
    Phid4 = zeros(Float64,max_iter)
    RMSd1 = zeros(Float64,max_iter)
    RMSd2 = zeros(Float64,max_iter)
    RMSd4 = zeros(Float64,max_iter)
    Phi = zeros(Float64,max_iter)
    Psi = zeros(Float64,max_iter)
    Stab = zeros(Float64,max_iter)
    RMS = zeros(Float64,max_iter)
    runtime = zeros(Float64,max_iter)
    mpre = zeros(Float64,length(m0))
    dpre = zeros(Float64,length(dobs))
    frho = zeros(Float64,nel)
    epoch = 0

    mu0 = 4 * pi * 1e-7
    t1 = time_ns()
    for itr = 1:max_iter
        epoch += 1
        println("  Iteration ",itr-1," ...")

        dcale1 = zeros(ComplexF64,nper,nst)
        dcale2 = zeros(ComplexF64,nper,nst)
        Aheg1 = zeros(ComplexF64,nper*nst,M)
        Aheg2 = zeros(ComplexF64,nper*nst,M)
        rho_h = rho[syel,:]
        role_h = rho_h[Int.(elle_h[:,1])]
        rori_h = rho_h[Int.(elri_h[:,1])]    
        dcalh = zeros(ComplexF64,nper,nst)
        Ahhg = zeros(ComplexF64,nper*nst,M)

	# Construct the stiffness matrix TM mode
	Dh = ConstructStiffnessMatrixTM(nel_h,nno_h,el2no_h,no2yz_h,rho_h,Area,eh2ee)

        tmpe = ComplexF64.(repeat(real.(De[diagind(De)]),1,nper))
        tmph = ComplexF64.(repeat(real.(Dh[diagind(Dh)]),1,nper))
        for idp = 1:nper
            omega = 2*pi/period[idp]

            # Define ve, vh: vectors of field
            ve = zeros(ComplexF64,nno)
            vh = zeros(ComplexF64,nno_h)

            # Calculate fields at the boundaries
            # - Top boundary
            ve_top = 1
            ve[notp_e] = ve_top * ones(size(notp_e,1),1)
            vh_top = 1
            vh[nsur_h] = vh_top * ones(size(nsur_h,1),1)
            # - Left boundary
            role_e = rho[Int.(elle_e[:,1])]
            ve[nole_e],~ = CalculateFieldAtSideBoundaryTE(role_e,period[idp],zlef_e,ve_top)
            vh[nole_h],~ = CalculateFieldAtSideBoundaryTM(role_h,period[idp],zlef_h,vh_top)
            # - Right boundary
            rori_e = rho[Int.(elri_e[:,1])]
            ve[nori_e],~ = CalculateFieldAtSideBoundaryTE(rori_e,period[idp],zrig_e,ve_top)
            vh[nori_h],~ = CalculateFieldAtSideBoundaryTM(rori_h,period[idp],zrig_h,vh_top)
            # - Bottom boundary
            ve[nobt_e] = transpose(collect(range(ve[nole_e[end]], stop=ve[nori_e[end]], length=size(nobt_e,1))))   
            vh[nobt_h] = transpose(collect(range(vh[nole_h[end]], stop=vh[nori_h[end]], length=size(nobt_h,1))))

            # Assign the diagional elements of De
            for ide = 1:nel
                tmpe[el2no[ide,:],idp] .+= -im*omega*mu0/3/rho[ide]*Area[ide]
            end
            De[diagind(De)] = tmpe[:,idp]

            # Assign the diagional elements of Dh
            for ide = 1:nel_h
                tmph[el2no_h[ide,:],idp] .+= -im*omega*mu0/3*Area[eh2ee[ide]]
            end
            Dh[diagind(Dh)] = tmph[:,idp]

            # Factorization phase TE mode, spawned
            dcDe = Threads.@spawn klu(sparse(@view(De[noin_e,noin_e])))

            # Factorization phase TM mode
            dcDh = klu(sparse(@view(Dh[noin_h,noin_h])))

            # Solving phase TM mode
            vh[noin_h] = dcDh\(-@view(Dh[noin_h,noex_h]) * @view(vh[noex_h]))

            # Calculate TM mode response functions at the stations
            dcalh,ahd,ch = CalculateResponseFunctionsTMRhoaPhase(idp,dcalh,nno_h,nst,no2yz_h,ns1h,ns2h,nbsh,rho_h,elsta,omega,mu0,vh)

            # Construct TM mode jacobian matrix
            uh = zeros(ComplexF64,nno_h,nst)
            uh[noin_h,:] = dcDh\ch[noin_h,:]
            Ahhg = ConstructJacobianTMRhoaPhase(idp,nst,Ahhg,M,el2no_h,el2no,no2yz,m2ee,m2eh,uh,vh,rho_h,tmst,ahd)

            # Solving phase TE mode
            dcDe = fetch(dcDe)
            ve[noin_e] = dcDe\(-@view(De[noin_e,noex_e]) * @view(ve[noex_e]))

            # Calculate TE mode response functions at the stations
            dcale1,dcale2,ce1,ce2,bxv2 = CalculateResponseFunctionsTERhoaPhaseTipper(idp,dcale1,dcale2,nno,nst,no2yz,ns1e,ns2e,nbse,omega,mu0,ve)

            # Construct TE mode jacobian matrix: rhoa, phase, tipper
            ce1[noin_e,:] = dcDe\ce1[noin_e,:]
            ce2[noin_e,:] = dcDe\ce2[noin_e,:]
            Aheg1,Aheg2 = ConstructJacobianTERhoaPhaseTipper(idp,nst,Aheg1,Aheg2,M,el2no,m2ee,omega,mu0,Area,rho,ce1,ce2,ve,bxv2)

        end

        # Merge jacobian matrices
        A = MergeJacobianRhoaPhaseTipper(N,M,nper,nst,Aheg1,Aheg2,Ahhg)

        # Calculate the data misfit, model roughness, and objective function
        dcal,RMSd1,RMSd2,RMSd4,Phid1,Phid2,Phid4,Phi,RMS,Stab,Psi = CalculateObjectiveFunctionRhoaPhaseTipper(RMSd1,RMSd2,RMSd4,Phid1,Phid2,Phid4,Phi,RMS,Stab,Psi,N,itr,nper,nst,dcale1,dcale2,dcalh,dobs,W,T,lambda,m,K,idd1,idd2,idd4)

        # Termination criteria: the objective function increases
        if itr > 1 && Psi[itr] - Psi[itr-1] >= 0
            println("Inversion failed to progress further")
            println("Stopping condition: Psi[",itr-1,"] > Psi[",itr-2,"]")
            println("Try increasing the damping for smaller step size")
            println("The model after ",itr-2," iter is the final model")
            runtime[itr] = time_ns()
            break
        elseif itr > 1 && Psi[itr-1] - Psi[itr] < 0.01*Psi[itr]
            println("Inversion reached a solution (a local minimum of the obj function)")
            println("Stopping condition: Psi[",itr-2,"] - Psi[",itr-1,"] < 0.01*Psi[",itr-1,"]")
            println("The model after ",itr-1," iter is the final model")
            dpre = dcal
            frho = rho
            mpre = m
            runtime[itr] = time_ns()
            break
        end
        # If iteration proceeds (not terminated), save data, model, total misfit
        dpre = dcal
        frho = rho
        mpre = m

        # Update the model when itr < max_iter
        if itr < max_iter
        # Calculate model update with Gauss-Newton data-space method
        invK = fetch(invK)
        m = CalculateModelUpdateDataSpace(invK,A,W,V,dobs,dcal,lambda,K,m,nnzN)
        rho[m2ee] = 10 .^m
        else
            println("Inversion reached the maximum number of iterations")
            println("Try decreasing the damping for larger step size")
            println("or simply increasing the maximum number of iterations")
            println("The model after ",itr-1," iter is the final model")
        end

        runtime[itr] = time_ns()

    end

    t2 = time_ns()

    println("Writing output files")
    # Export Log File
    ExportLogFileRhoaPhaseTipper(outroot,filedat,filem0,filestg,filetopo,epoch,idd1,idd2,idd4,M,t1,t2,Phi,Psi,Stab,lambda,RMS,RMSd1,RMSd2,RMSd4,runtime)
    # Export the final model
    ExportFinalModelResistivity(outroot,nel,nno,el2no,no2yz,frho)
    #Export model responses
    ExportFinalModelRhoaPhaseTipper(outroot,nper,nst,sta,dpre,period,rhe_o,ere_o,phe_o,epe_o,rhh_o,erh_o,phh_o,eph_o,rtz_o,ert_o,itz_o,eit_o)

    return

end
