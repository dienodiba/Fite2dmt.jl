# Dieno Diba 2024
# FITE2DMT
# Read observation data from input file
# Rhoa, phase, tipper

function ReadObservedDataRhoaPhaseTipper(filedat::String, nodata::Number, wg_rte::Integer, wg_pte::Integer, wg_rtm::Integer, wg_ptm::Integer, wg_tyr::Integer, wg_tyi::Integer, ef_zte::Number, ef_ztm::Number, ef_tzy::Number)

    dat = readdlm(filedat)
    nst = dat[1,3]
    nper = dat[2,3]
    sta = zeros(Float64,nst,2)
    period = zeros(Float64,nper)
    rhe_o = zeros(Float64,nper,nst)
    ere_o = zeros(Float64,nper,nst)
    phe_o = zeros(Float64,nper,nst)
    epe_o = zeros(Float64,nper,nst)
    rhh_o = zeros(Float64,nper,nst)
    erh_o = zeros(Float64,nper,nst)
    phh_o = zeros(Float64,nper,nst)
    eph_o = zeros(Float64,nper,nst)
    rtz_o = zeros(Float64,nper,nst)
    ert_o = zeros(Float64,nper,nst)
    itz_o = zeros(Float64,nper,nst)
    eit_o = zeros(Float64,nper,nst)
    tmp = 2
    for ids = 1:nst
        tmp += 2
        sta[ids,1] = 1e+3 * dat[tmp,3]
        tmp += 1
        sta[ids,2] = 1e+3 * dat[tmp,3]
        tmp += 1
        for idp = 1:nper
            tmp += 1
            period[idp] = dat[tmp,1]
            rhe_o[idp,ids] = log10(dat[tmp,2])
            phe_o[idp,ids] = pi/180*dat[tmp,4]
            if rhe_o[idp,ids] != log10(nodata)
                ere_o[idp,ids] = dat[tmp,3]/dat[tmp,2]/log(10)
                epe_o[idp,ids] = pi/180*dat[tmp,5]
            else
                ere_o[idp,ids] = nodata
                epe_o[idp,ids] = nodata
            end
            rhh_o[idp,ids] = log10(dat[tmp,6])
            phh_o[idp,ids] = pi/180*dat[tmp,8]
            if rhh_o[idp,ids] != log10(nodata)
                erh_o[idp,ids] = dat[tmp,7]/dat[tmp,6]/log(10)
                eph_o[idp,ids] = pi/180*dat[tmp,9]
            else
                erh_o[idp,ids] = nodata
                eph_o[idp,ids] = nodata
            end
            rtz_o[idp,ids] = dat[tmp,10]
            ert_o[idp,ids] = dat[tmp,11]
            itz_o[idp,ids] = dat[tmp,12]
            eit_o[idp,ids] = dat[tmp,13]
        end
    end

    # Data component:
    # - 0.5*log10(rhoa_TE)
    # - 0.5*log10(rhoa_TM)
    # - phase_TE/log(10)
    # - phase_TM/log(10)
    # - real(Tzy)
    # - imag(Tzy)
    # N: Number of data
    N = nst*nper*6
    # Construct vector of observed data considering error floors
    dobs = zeros(Float64,N)
    eobs = zeros(Float64,N)
    idd = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:6
                idd += 1
                if idj == 1
                    dobs[idd] = rhe_o[idp,ids] * 0.5
                    eobs[idd] = ere_o[idp,ids]
                    if eobs[idd] < 2*ef_zte/log(10)
                        eobs[idd] = 2*ef_zte/log(10)
                    end
                end
                if idj == 2
                    dobs[idd] = phe_o[idp,ids] / log(10)
                    eobs[idd] = epe_o[idp,ids]
                    if eobs[idd] < asin(ef_zte)
                        eobs[idd] = asin(ef_zte)
                    end
                end
                if idj == 3
                    dobs[idd] = rtz_o[idp,ids]
                    eobs[idd] = ert_o[idp,ids]
                    if eobs[idd] < ef_tzy
                        eobs[idd] = ef_tzy
                    end
                end
                if idj == 4
                    dobs[idd] = itz_o[idp,ids]
                    eobs[idd] = eit_o[idp,ids]
                    if eobs[idd] < ef_tzy
                        eobs[idd] = ef_tzy
                    end
                end
                if idj == 5
                    dobs[idd] = rhh_o[idp,ids] * 0.5
                    eobs[idd] = erh_o[idp,ids]
                    if eobs[idd] < 2*ef_ztm/log(10)
                        eobs[idd] = 2*ef_ztm/log(10)
                    end
                end
                if idj == 6
                    dobs[idd] = phh_o[idp,ids] / log(10)
                    eobs[idd] = eph_o[idp,ids]
                    if eobs[idd] < asin(ef_ztm)
                        eobs[idd] = asin(ef_ztm)
                    end
                end
            end
        end
    end
    # Diagonal matrix of error inverse for weighting
    W = zeros(Float64,N,N)
    V = zeros(Float64,N,N)
    nnzW = zeros(Int64,N)
    for idN = 1:N
        if eobs[idN] != nodata
            W[idN,idN] = (eobs[idN])^-2
            V[idN,idN] = (eobs[idN])^2
            nnzW[idN] = idN
        end
    end
    nnzW = filter(x->x!=0,nnzW)
    # Count the number of each data component
    T = zeros(Int64,N)
    idd1 = 0
    idd2 = 0
    idd4 = 0
    idN = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:6
                idN += 1
                if idj == 1
                    T[idN] = wg_rte
                    if eobs[idN] != nodata
                        idd1 += wg_rte
                    end
                end
                if idj == 2
                    T[idN] = wg_pte
                    if eobs[idN] != nodata
                        idd1 += wg_pte
                    end
                end
                if idj == 3
                    T[idN] = wg_tyr
                    if eobs[idN] != nodata
                        idd2 += wg_tyr
                    end
                end
                if idj == 4
                    T[idN] = wg_tyi
                    if eobs[idN] != nodata
                        idd2 += wg_tyi
                    end
                end
                if idj == 5
                    T[idN] = wg_rtm
                    if eobs[idN] != nodata
                        idd4 += wg_rtm
                    end
                end
                if idj == 6
                    T[idN] = wg_ptm
                    if eobs[idN] != nodata
                        idd4 += wg_ptm
                    end
                end
            end
        end
    end
    nnzT = findall(T .!= 0)
    nnzN = intersect(nnzT,nnzW)

    return nst,nper,sta,period,rhe_o,ere_o,phe_o,epe_o,rhh_o,erh_o,phh_o,eph_o,rtz_o,ert_o,itz_o,eit_o,dobs,eobs,N,W,V,T,nnzN,idd1,idd2,idd4

end
