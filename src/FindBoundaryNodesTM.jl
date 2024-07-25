# Dieno Diba 2024
# FITE2DMT
# Find nodes in the boundaries and air-earth interface for TM mode

function FindBoundaryNodesTM(interp_topo::Any, nno_h::Integer, no2yz_h::Matrix, ymin::Number, ymax::Number, zmax::Number, minedge::Number)

    zlef_h = []
    zrig_h = []
    ysur_h = []
    ybot_h = []
    nole_h = []
    nori_h = []
    nsur_h = []
    nobt_h = []
    for idn = 1:nno_h
        if no2yz_h[idn,1] == ymin
            zlef_h = vcat(zlef_h,no2yz_h[idn,2])
            nole_h = vcat(nole_h,idn)
        elseif no2yz_h[idn,1] == ymax
            zrig_h = vcat(zrig_h,no2yz_h[idn,2])
            nori_h = vcat(nori_h,idn)
        elseif no2yz_h[idn,2] == zmax
            ybot_h = vcat(ybot_h,no2yz_h[idn,1])
            nobt_h = vcat(nobt_h,idn)
        elseif abs(no2yz_h[idn,2] - interp_topo(no2yz_h[idn,1])) < minedge/2
            ysur_h = vcat(ysur_h,no2yz_h[idn,1])
            nsur_h = vcat(nsur_h,idn)
        end
    end
    nole_h = Int.(nole_h[sortperm(zlef_h)])
    zlef_h = float.(zlef_h[sortperm(zlef_h)])
    nori_h = Int.(nori_h[sortperm(zrig_h)])
    zrig_h = float.(zrig_h[sortperm(zrig_h)])
    nobt_h = Int.(nobt_h[sortperm(ybot_h)])
    ybot_h = float.(ybot_h[sortperm(ybot_h)])
    nsur_h = Int.(nsur_h[sortperm(ysur_h)])
    ysur_h = float.(ysur_h[sortperm(ysur_h)])

    return zlef_h,zrig_h,ysur_h,ybot_h,nole_h,nori_h,nsur_h,nobt_h 

end
