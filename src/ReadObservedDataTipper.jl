# Dieno Diba 2024
# FITE2DMT
# Read observation data from input file
# Tipper

function ReadObservedDataTipper(filedat::String, nodata::Number, wg_tyr::Integer, wg_tyi::Integer, ef_tzy::Number)

    dat = readdlm(filedat)
    nst = dat[1,1]
    nper = dat[2,1]
    sta = zeros(Float64,nst,2)
    period = zeros(Float64,nper)
    rtz_o = zeros(Float64,nper,nst)
    ert_o = zeros(Float64,nper,nst)
    itz_o = zeros(Float64,nper,nst)
    eit_o = zeros(Float64,nper,nst)
    tmp = 2
    for ids = 1:nst
        tmp += 1
        sta[ids,1] = 1e+3 * dat[tmp,1]
        tmp += 1
        sta[ids,2] = 1e+3 * dat[tmp,1]
        for idp = 1:nper
            tmp += 1
            period[idp] = dat[tmp,1]
            rtz_o[idp,ids] = dat[tmp,2]
            ert_o[idp,ids] = dat[tmp,3]
            itz_o[idp,ids] = dat[tmp,4]
            eit_o[idp,ids] = dat[tmp,5]
        end
    end

    # Data component:
    # - real(Tzy)
    # - imag(Tzy)
    # N: Number of data
    N = nst*nper*2
    # Construct vector of observed data considering error floors
    dobs = zeros(Float64,N)
    eobs = zeros(Float64,N)
    idd = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:2
                idd += 1
                if idj == 1
                    dobs[idd] = rtz_o[idp,ids]
                    eobs[idd] = ert_o[idp,ids]
                    if eobs[idd] < ef_tzy
                        eobs[idd] = ef_tzy
                    end
                end
                if idj == 2
                    dobs[idd] = itz_o[idp,ids]
                    eobs[idd] = eit_o[idp,ids]
                    if eobs[idd] < ef_tzy
                        eobs[idd] = ef_tzy
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
    idd2 = 0
    idN = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:2
                idN += 1
                if idj == 1
                    T[idN] = wg_tyr
                    if eobs[idN] != nodata
                        idd2 += wg_tyr
                    end
                end
                if idj == 2
                    T[idN] = wg_tyi
                    if eobs[idN] != nodata
                        idd2 += wg_tyi
                    end
                end
            end
        end
    end
    nnzT = findall(T .!= 0)
    nnzN = intersect(nnzT,nnzW)

    return nst,nper,sta,period,rtz_o,ert_o,itz_o,eit_o,dobs,eobs,N,W,V,T,nnzN,idd2

end
