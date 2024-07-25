# Dieno Diba 2024
# FITE2DMT
# Find elements in the boundaries and air-earth interface for TM mode

function FindBoundaryElementsTM(nel_h::Integer, el2no_h::Matrix, no2yz_h::Matrix, nole_h::Vector, nori_h::Vector, nsur_h::Vector)

    elle_h = zeros(Float64,nel_h,2)
    elri_h = zeros(Float64,nel_h,2)
    elsur_h = zeros(Int64,nel_h)
    tmd = 0
    for ide = 1:nel_h
        if length(intersect(Set(nole_h),Set(el2no_h[ide,:]))) == 2
            pt = findall(in(nole_h),el2no_h[ide,:])
            tmp = mean([no2yz_h[el2no_h[ide,pt[1]],2] no2yz_h[el2no_h[ide,pt[2]],2]])
            elle_h[ide,:] = [ide tmp]
        end
        if length(intersect(Set(nori_h),Set(el2no_h[ide,:]))) == 2
            pt = findall(in(nori_h),el2no_h[ide,:])
            tmp = mean([no2yz_h[el2no_h[ide,pt[1]],2] no2yz_h[el2no_h[ide,pt[2]],2]])
            elri_h[ide,:] = [ide tmp]
        end
        if length(intersect(Set(nsur_h),Set(el2no_h[ide,:]))) == 2
            elsur_h[ide] = ide
        end
    end
    elle_h = float.(elle_h[sortperm(elle_h[:,2]),:])
    elle_h = elle_h[(elle_h[:,1] .!= 0),:]
    elri_h = float.(elri_h[sortperm(elri_h[:,2]),:])
    elri_h = elri_h[(elri_h[:,1] .!= 0),:]
    elsur_h = filter(x->x!=0,elsur_h)

    return elle_h,elri_h,elsur_h

end
