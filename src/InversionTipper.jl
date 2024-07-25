# Dieno Diba 2024
# FITE2DMT
# Inversion of tipper

function InversionTipper(filedat::String, filem0::String, filetopo::String, filestg::String)

    println("=================================================================")
    println("FITE2DMT: Finite Triangular Elements for Two-Dimensional MT")
    println("Dieno Diba, 2023")
    println("=================================================================")

    println("Reading input files and initializing ...")
    
    # Read inversion setting parameters
    lambda,max_iter,outroot,fix_el,damping,hor2ver,static,nodata,wg_tyr,wg_tyi,ef_tzy = ReadInversionSettingParametersTipper(filestg)

    # Read observation data of tipper
    nst,nper,sta,period,rtz_o,ert_o,itz_o,eit_o,dobs,eobs,N,W,V,T,nnzN,idd2 = ReadObservedDataTipper(filedat,nodata,wg_tyr,wg_tyi,ef_tzy)

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

    # Construct stiffness matrix TE mode
    De = ConstructStiffnessMatrixTE(nel,nno,el2no,no2yz,Area)

    println("Constructing roughness matrix ...")
    K = ConstructRoughnessMatrix(nel,no2yz,el2no,static,hor2ver,elsta,mbin)
    invK = Threads.@spawn inv(lambda*Matrix(K) + I(M)*damping)

    # Gauss-Newton Iteration
    println("Minimizing objective function")
    Phid2 = zeros(Float64,max_iter)
    RMSd2 = zeros(Float64,max_iter)
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

        dcale2 = zeros(ComplexF64,nper,nst)
        Aheg2 = zeros(ComplexF64,nper*nst,M)
        rho_h = rho[syel,:]

        tmpe = ComplexF64.(repeat(real.(De[diagind(De)]),1,nper))
        for idp = 1:nper
            omega = 2*pi/period[idp]

            # Define ve: vector of electric field
            ve = zeros(ComplexF64,nno)

            # Calculate fields at the boundaries
            # - Top boundary
            ve_top = 1
            ve[notp_e] = ve_top * ones(size(notp_e,1),1)
            # - Left boundary
            role_e = rho[Int.(elle_e[:,1])]
            ve[nole_e],~ = CalculateFieldAtSideBoundaryTE(role_e,period[idp],zlef_e,ve_top)
            # - Right boundary
            rori_e = rho[Int.(elri_e[:,1])]
            ve[nori_e],~ = CalculateFieldAtSideBoundaryTE(rori_e,period[idp],zrig_e,ve_top)            
            # - Bottom boundary
            ve[nobt_e] = transpose(collect(range(ve[nole_e[end]], stop=ve[nori_e[end]], length=size(nobt_e,1))))   

            # Assign the diagional elements of De
            for ide = 1:nel
                tmpe[el2no[ide,:],idp] .+= -im*omega*mu0/3/rho[ide]*Area[ide]
            end
            De[diagind(De)] = tmpe[:,idp]

            # Factorization phase TE mode
            dcDe = klu(sparse(@view(De[noin_e,noin_e])))

            # Solving phase TE mode
            ve[noin_e] = dcDe\(-@view(De[noin_e,noex_e]) * @view(ve[noex_e]))

            # Calculate TE mode response functions at the stations
            dcale2,ce2,bxv2 = CalculateResponseFunctionsTETipper(idp,dcale2,nno,nst,no2yz,ns1e,ns2e,nbse,omega,mu0,ve)

            # Construct TE mode jacobian matrix: tipper
            ce2[noin_e,:] = dcDe\ce2[noin_e,:]
            Aheg2 = ConstructJacobianTETipper(idp,nst,Aheg2,M,el2no,m2ee,omega,mu0,Area,rho,ce2,ve,bxv2)

        end

        # Merge jacobian matrices
        A = MergeJacobianTipper(N,M,nper,nst,Aheg2)

        # Calculate the data misfit, model roughness, and objective function
        dcal,RMSd2,Phid2,Phi,RMS,Stab,Psi = CalculateObjectiveFunctionTipper(RMSd2,Phid2,Phi,RMS,Stab,Psi,N,itr,nper,nst,dcale2,dobs,W,T,lambda,m,K,idd2)

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
    ExportLogFileTipper(outroot,filedat,filem0,filestg,filetopo,epoch,idd2,M,t1,t2,Phi,Psi,Stab,lambda,RMS,RMSd2,runtime)
    # Export the final model
    ExportFinalModelResistivity(outroot,nel,nno,el2no,no2yz,frho)
    #Export model responses
    ExportFinalModelTipper(outroot,nper,nst,sta,dpre,period,rtz_o,ert_o,itz_o,eit_o)

    return

end
