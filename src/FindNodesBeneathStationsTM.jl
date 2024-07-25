# Dieno Diba 2024
# FITE2DMT
# Find the nodes beneath each station for TM mode

function FindNodesbeneathStationsTM(nst::Integer, ysur_h::Vector, sta::Matrix, no2yz_h::Matrix, el2no_h::Matrix, elsur_h::Vector, m2eh::Vector, M::Integer)

    ns1h = zeros(Int64,nst)
    ns2h = zeros(Int64,nst)
    nbsh = zeros(Int64,nst)
    elsta = zeros(Int64,nst)
    for ids = 1:nst
        (~,itp1) = findmin(abs.(ysur_h .- sta[ids,1]))
        itp1 = itp1[1]
        if ysur_h[itp1] > sta[ids,1]
            itp2 = itp1 - 1
        else
            itp2 = itp1 + 1
        end
        tmp = findall(no2yz_h[:,1] .== ysur_h[itp1])
        (~,tmd) = findmin(abs.(no2yz_h[tmp,2] .- sta[ids,2]))
        ns1h[ids] = tmp[tmd]
        tmp = findall(no2yz_h[:,1] .== ysur_h[itp2])
        (~,tmd) = findmin(abs.(no2yz_h[tmp,2] .- sta[ids,2]))
        ns2h[ids] = tmp[tmd]
        @threads for ideh in elsur_h
            if length(intersect([ns1h[ids] ns2h[ids]],el2no_h[ideh,:])) == 2
                nbsh[ids] = only(setdiff(el2no_h[ideh,:],[ns1h[ids] ns2h[ids]]))
                elsta[ids] = ideh
                break
            end
        end
    end
    tmst = zeros(Int64,M)
    @threads for idM = 1:M
        if in(m2eh[idM],elsta) == 1
            tmst[idM] = only(findall(elsta .== m2eh[idM]))
        end
    end

    return ns1h,ns2h,nbsh,elsta,tmst

end
