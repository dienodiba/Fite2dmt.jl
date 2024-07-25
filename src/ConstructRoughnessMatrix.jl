# Dieno Diba 2024
# FITE2DMT
# Construction of roughness matrix

# Type of roughness operator: penalizing elements sharing an edge, index only, no weighting between adjacent elements

function ConstructRoughnessMatrix(nel::Integer, no2yz::Matrix, el2no::Matrix, static::Integer, hor2ver::Number, elsta::Vector, mbin::Vector)

    L = zeros(Float64,nel,nel)
    tmpa = el2no[:,1]
    tmpb = el2no[:,2]
    tmpc = el2no[:,3]
    for ide = 1:nel
        if static == 1 && in(ide,elsta) == 1
            continue
        end
        #edge1 = sqrt((no2yz[el2no[ide,3],1] - no2yz[el2no[ide,2],1])^2 + (no2yz[el2no[ide,3],2] - no2yz[el2no[ide,2],2])^2)
        #edge2 = sqrt((no2yz[el2no[ide,1],1] - no2yz[el2no[ide,3],1])^2 + (no2yz[el2no[ide,1],2] - no2yz[el2no[ide,3],2])^2)
        #edge3 = sqrt((no2yz[el2no[ide,2],1] - no2yz[el2no[ide,1],1])^2 + (no2yz[el2no[ide,2],2] - no2yz[el2no[ide,1],2])^2)
        #dis = mean([edge1,edge2,edge3])
        tmy = (no2yz[el2no[ide,1],1] + no2yz[el2no[ide,2],1] + no2yz[el2no[ide,3],1]) / 3
        tmz = (no2yz[el2no[ide,1],2] + no2yz[el2no[ide,2],2] + no2yz[el2no[ide,3],2]) / 3
        tmp = el2no[ide,:]
        va1 = findall(tmpa .== tmp[1])
        vb1 = findall(tmpb .== tmp[1])
        vc1 = findall(tmpc .== tmp[1])
        vu1 = unique(vcat(va1,vb1,vc1))
        va2 = findall(tmpa .== tmp[2])
        vb2 = findall(tmpb .== tmp[2])
        vc2 = findall(tmpc .== tmp[2])
        vu2 = unique(vcat(va2,vb2,vc2))
        va3 = findall(tmpa .== tmp[3])
        vb3 = findall(tmpb .== tmp[3])
        vc3 = findall(tmpc .== tmp[3])
        vu3 = unique(vcat(va3,vb3,vc3))
        sheg = zeros(Int64,3)
        count = 0
        for ideg = 1:3
            if ideg == 1; a = vu1; b = vu2; end
            if ideg == 2; a = vu2; b = vu3; end
            if ideg == 3; a = vu3; b = vu1; end
            if isempty(setdiff(intersect(a,b),ide)) == 0
                sheg[ideg] = only(setdiff(intersect(a,b),ide))
                if static == 1 && in(sheg[ideg],elsta) == 1
                    sheg[ideg] = 0
                    continue
                end
                if mbin[sheg[ideg]] == 1
                    count += 1
                    tmyn = (no2yz[el2no[sheg[ideg],1],1] + no2yz[el2no[sheg[ideg],2],1] + no2yz[el2no[sheg[ideg],3],1]) / 3
                    tmzn = (no2yz[el2no[sheg[ideg],1],2] + no2yz[el2no[sheg[ideg],2],2] + no2yz[el2no[sheg[ideg],3],2]) / 3
                    L[ide,ide] += -1 / sqrt(hor2ver*(tmy-tmyn)^2 + (tmz-tmzn)^2) * 1e+3
                    L[ide,sheg[ideg]] = 1 / sqrt(hor2ver*(tmy-tmyn)^2 + (tmz-tmzn)^2) * 1e+3
                end
            end
        end
        #L[ide,ide] = -1.0
        #L[ide,filter(x->x>0,sheg)] .= 1/count
        #L[ide,ide] = -1.0 / (dis * 1e-3)
        #L[ide,filter(x->x>0,sheg)] .= 1/count / (dis * 1e-3)
    end
    tmp = findall(mbin .== 1)
    L = L[tmp,tmp]
    K = sparse(L') * sparse(L)

    return K
end
