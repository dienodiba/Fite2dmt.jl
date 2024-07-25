# Dieno Diba 2024
# FITE2DMT
# Find nodes in the boundaries and air-earth interface for TE mode

function FindBoundaryNodesTE(interp_topo::Any, nno::Integer, no2yz::Matrix, ymin::Number, ymax::Number, zmin::Number, zmax::Number, minedge::Number)

    zlef_e = []
    zrig_e = []
    ytop_e = []
    ybot_e = []
    ysur_e = []
    nole_e = []
    nori_e = []
    notp_e = []
    nobt_e = []
    nsur_e = []
    for idn = 1:nno
        if no2yz[idn,1] == ymin
            zlef_e = vcat(zlef_e,no2yz[idn,2])
            nole_e = vcat(nole_e,idn)
        elseif no2yz[idn,1] == ymax
            zrig_e = vcat(zrig_e,no2yz[idn,2])
            nori_e = vcat(nori_e,idn)
        elseif no2yz[idn,2] == zmin
            ytop_e = vcat(ytop_e,no2yz[idn,1])
            notp_e = vcat(notp_e,idn)
        elseif no2yz[idn,2] == zmax
            ybot_e = vcat(ybot_e,no2yz[idn,1])
            nobt_e = vcat(nobt_e,idn)
        end
        if abs(no2yz[idn,2] - interp_topo(no2yz[idn,1])) < minedge/2
            ysur_e = vcat(ysur_e,no2yz[idn,1])
            nsur_e = vcat(nsur_e,idn)
        end
    end
    nole_e = Int.(nole_e[sortperm(zlef_e)])
    zlef_e = float.(zlef_e[sortperm(zlef_e)])
    nori_e = Int.(nori_e[sortperm(zrig_e)])
    zrig_e = float.(zrig_e[sortperm(zrig_e)])
    notp_e = Int.(notp_e[sortperm(ytop_e)])
    ytop_e = float.(ytop_e[sortperm(ytop_e)])
    nobt_e = Int.(nobt_e[sortperm(ybot_e)])
    ybot_e = float.(ybot_e[sortperm(ybot_e)])
    nsur_e = Int.(nsur_e[sortperm(ysur_e)])
    ysur_e = float.(ysur_e[sortperm(ysur_e)])

    return zlef_e,zrig_e,ytop_e,ybot_e,ysur_e,nole_e,nori_e,notp_e,nobt_e,nsur_e
end
