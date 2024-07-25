# Dieno Diba 2024
# FITE2DMT
# Export tipper of the final inversion model

function ExportFinalModelTipper(outroot::Any, nper::Integer, nst::Integer, sta::Matrix, dpre::Vector, period::Vector, rtz_o::Matrix, ert_o::Matrix, itz_o::Matrix, eit_o::Matrix)

    rtz_c = zeros(Float64,nper,nst)
    itz_c = zeros(Float64,nper,nst)
    idN = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:2
                idN += 1
                if idj == 1
                    rtz_c[idp,ids] = dpre[idN]
                end
                if idj == 2
                    itz_c[idp,ids] = dpre[idN]
                end
            end
        end
    end
    fid = open(string(outroot,"_cal.txt"),"w")
    @printf(fid,"nstation   = %i\n",nst)
    @printf(fid,"nperiod    = %i\n\n",nper)
    for ids = 1:nst
        @printf(fid,"# %i\n",ids)
        @printf(fid,"y(km) = %f\n",sta[ids,1]*1e-3)
        @printf(fid,"z(km) = %f\n",sta[ids,2]*1e-3)
        print(fid,"      period(s)       tzyr_obs       tzyr_cal       tzyr_err       tzyi_obs       tzyi_cal       tzyi_err\n")
        for idp = 1:nper
            @printf(fid,"%15.5e",period[idp])
            @printf(fid,"%15.5e%15.5e%15.5e",rtz_o[idp,ids],rtz_c[idp,ids],ert_o[idp,ids])
            @printf(fid,"%15.5e%15.5e%15.5e\n",itz_o[idp,ids],itz_c[idp,ids],eit_o[idp,ids])
        end
        @printf(fid,"\n")
    end
    close(fid)

    return

end
