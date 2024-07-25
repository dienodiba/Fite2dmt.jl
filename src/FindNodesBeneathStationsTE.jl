# Dieno Diba 2024
# FITE2DMT
# Find the nodes beneath each station for TE mode

function FindNodesBeneathStationsTE(nst::Integer, ysur_e::Vector, sta::Matrix, no2yz::Matrix, el2no::Matrix, elsur_e::Vector, interp_topo::Any)

    ns1e = zeros(Int64,nst)
    ns2e = zeros(Int64,nst)
    nbse = zeros(Int64,nst)
    for ids = 1:nst
        (~,itp1) = findmin(abs.(ysur_e .- sta[ids,1]))
        itp1 = itp1[1]
        if ysur_e[itp1] > sta[ids,1]
            itp2 = itp1 - 1
        else
            itp2 = itp1 + 1
        end
        tmp = findall(no2yz[:,1] .== ysur_e[itp1])
        (~,tmd) = findmin(abs.(no2yz[tmp,2] .- sta[ids,2]))
        ns1e[ids] = tmp[tmd]
        tmp = findall(no2yz[:,1] .== ysur_e[itp2])
        (~,tmd) = findmin(abs.(no2yz[tmp,2] .- sta[ids,2]))
        ns2e[ids] = tmp[tmd]
        @threads for ide in elsur_e
            if length(intersect([ns1e[ids] ns2e[ids]],el2no[ide,:])) == 2 &&
                only(no2yz[setdiff(el2no[ide,:],[ns1e[ids] ns2e[ids]]),2]) > interp_topo(no2yz[setdiff(el2no[ide,:],[ns1e[ids] ns2e[ids]]),1])
                nbse[ids] = only(setdiff(el2no[ide,:],[ns1e[ids] ns2e[ids]]))
            end
        end
    end

  return ns1e,ns2e,nbse

end
