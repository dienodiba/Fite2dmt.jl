# Dieno Diba 2024
# FITE2DMT
# Export rhoa phase of the final inversion model

function ExportFinalModelRhoaPhase(outroot::Any, nper::Integer, nst::Integer, sta::Matrix, dpre::Vector, period::Vector, rhe_o::Matrix, ere_o::Matrix, phe_o::Matrix, epe_o::Matrix, rhh_o::Matrix, erh_o::Matrix, phh_o::Matrix, eph_o::Matrix)

    rhe_c = zeros(Float64,nper,nst)
    phe_c = zeros(Float64,nper,nst)
    rhh_c = zeros(Float64,nper,nst)
    phh_c = zeros(Float64,nper,nst)
    idN = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:4
                idN += 1
                if idj == 1
                    rhe_c[idp,ids] = 10^(2*dpre[idN])
                    rhe_o[idp,ids] = 10^(rhe_o[idp,ids])
                end
                if idj == 2
                    phe_c[idp,ids] = log(10)*dpre[idN]*180/pi
                    phe_o[idp,ids] = phe_o[idp,ids]*180/pi
                end
                if idj == 3
                    rhh_c[idp,ids] = 10^(2*dpre[idN])
                    rhh_o[idp,ids] = 10^(rhh_o[idp,ids])
                end
                if idj == 4
                    phh_c[idp,ids] = log(10)*dpre[idN]*180/pi
                    phh_o[idp,ids] = phh_o[idp,ids]*180/pi
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
        print(fid,"      period(s) roe_obs(Ohm.m) roe_cal(Ohm.m) roe_err(Ohm.m) roh_obs(Ohm.m) roh_cal(Ohm.m) roh_err(Ohm.m)   phe_obs(deg)   phe_cal(deg)   phe_err(deg)   phh_obs(deg)   phh_cal(deg)   phh_err(deg)\n")
        for idp = 1:nper
            @printf(fid,"%15.5e",period[idp])
            @printf(fid,"%15.5e%15.5e%15.5e",rhe_o[idp,ids],rhe_c[idp,ids],ere_o[idp,ids])
            @printf(fid,"%15.5e%15.5e%15.5e",rhh_o[idp,ids],rhh_c[idp,ids],erh_o[idp,ids])
            @printf(fid,"%15.5e%15.5e%15.5e",phe_o[idp,ids],phe_c[idp,ids],epe_o[idp,ids])
            @printf(fid,"%15.5e%15.5e%15.5e\n",phh_o[idp,ids],phh_c[idp,ids],eph_o[idp,ids])
        end
        @printf(fid,"\n")
    end
    close(fid)

    return

end
