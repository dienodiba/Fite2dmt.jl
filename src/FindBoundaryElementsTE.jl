# Dieno Diba 2024
# FITE2DMT
# Find elements in the boundaries and air-earth interface for TE mode

function FindBoundaryElementsTE(nel::Integer, el2no::Matrix, no2yz::Matrix, nole_e::Vector, nori_e::Vector, nsur_e::Vector)

    elle_e = zeros(Float64,nel,2)
    elri_e = zeros(Float64,nel,2)
    elsur_e = zeros(Int64,nel)
    for ide = 1:nel
        if length(intersect(Set(nole_e),Set(el2no[ide,:]))) == 2
            pt = findall(in(nole_e),el2no[ide,:])
            tmp = mean([no2yz[el2no[ide,pt[1]],2] no2yz[el2no[ide,pt[2]],2]])
            elle_e[ide,:] = [float(ide) tmp]
        end
        if length(intersect(Set(nori_e),Set(el2no[ide,:]))) == 2
            pt = findall(in(nori_e),el2no[ide,:])
            tmp = mean([no2yz[el2no[ide,pt[1]],2] no2yz[el2no[ide,pt[2]],2]])
            elri_e[ide,:] = [float(ide) tmp]
        end
        if length(intersect(Set(nsur_e),Set(el2no[ide,:]))) == 2
            elsur_e[ide] = ide
        end
    end
    elle_e = float.(elle_e[sortperm(elle_e[:,2]),:])
    elle_e = elle_e[(elle_e[:,1] .!= 0),:]
    elri_e = float.(elri_e[sortperm(elri_e[:,2]),:])
    elri_e = elri_e[(elri_e[:,1] .!= 0),:]
    elsur_e = filter(x->x!=0,elsur_e)

    return elle_e, elri_e, elsur_e

end
