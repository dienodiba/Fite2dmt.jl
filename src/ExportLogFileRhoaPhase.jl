# Dieno Diba 2024
# FITE2DMT
# Export log file for rhoa phase inversion

function ExportLogFileRhoaPhase(outroot::Any, filedat::String, filem0::String, filestg::String, filetopo::String, epoch::Integer, idd1::Integer, idd4::Integer, M::Integer, t1::Number, t2::Number, Phi::Vector, Psi::Vector, Stab::Vector, lambda::Number, RMS::Vector, RMSd1::Vector, RMSd4::Vector, runtime::Vector)

    fid = open(string(outroot,"_log.txt"),"w")
    @printf(fid,"FITE2DMT Inversion Report\n")
    @printf(fid,"%s JST\n\n",Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    @printf(fid,"===============================\n")
    @printf(fid,"Input:\n")
    @printf(fid,"%s\n",filedat)
    @printf(fid,"%s\n",filem0)
    @printf(fid,"%s\n",filestg)
    @printf(fid,"%s\n",filetopo)
    @printf(fid,"===============================\n")
    @printf(fid,"# of iteration = %i\n",epoch-1)
    @printf(fid,"# of data = %i\n",idd1+idd4)
    @printf(fid,"# of model parameter = %i\n",M)
    @printf(fid,"Total runtime (s) = %7.1f\n",(t2-t1)/1e+09)
    @printf(fid,"===============================\n")
    @printf(fid,"Convergence history\n")
    for idi = 1:epoch
        @printf(fid,"== Iteration %i\n",idi-1)
        @printf(fid,"Objective function = %3.2e\n",Psi[idi])
        @printf(fid,"Data misfit        = %3.2e\n",Phi[idi])
        @printf(fid,"Model roughness    = %3.2e\n",Stab[idi]/lambda)
        @printf(fid,"Combined RMS       = %3.2e\n",RMS[idi])
        @printf(fid,"RMS_ZTE            = %3.2e\n",RMSd1[idi])
        @printf(fid,"RMS_ZTM            = %3.2e\n",RMSd4[idi])
        @printf(fid,"Runtime (s)        = %3.2e\n",(runtime[idi]-t1)/1e+9)
    end
    close(fid)

    return

end
