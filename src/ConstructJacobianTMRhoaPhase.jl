# Dieno Diba 2024
# FITE2DMT
# Construct jacobian matrix for TM mode apparent resistivity and phase

function ConstructJacobianTMRhoaPhase(idp::Integer, nst::Integer, Ahhg::Matrix, M::Integer, el2no_h::Matrix, el2no::Matrix, no2yz::Matrix, m2ee::Vector, m2eh::Vector, uh::Matrix, vh::Vector, rho_h::Matrix, tmst::Vector, ahd::Vector)

    for idM = 1:M
        jiw = el2no_h[m2eh[idM],:]
        yz = no2yz[el2no[m2ee[idM],:],:]
        # Edges
        s1 = yz[3,:] - yz[2,:]
        s2 = yz[1,:] - yz[3,:]
        s3 = yz[2,:] - yz[1,:]
        # Area of the element ide
        Ar = abs(0.5 * (s2[1]*s3[2] - s2[2]*s3[1]))
        # Gradient of test functions
        grad_phi1 = [-s1[2],s1[1]]/(2 * Ar)
        grad_phi2 = [-s2[2],s2[1]]/(2 * Ar)
        grad_phi3 = [-s3[2],s3[1]]/(2 * Ar)
        grad_phi = hcat(grad_phi1,grad_phi2,grad_phi3)
        Ahhg[(idp-1)*nst+1:idp*nst,idM] = -transpose(@view(uh[jiw,:])) * (-grad_phi'*(grad_phi*Ar*@view(vh[jiw]))) * rho_h[m2eh[idM]]
        if tmst[idM] != 0
            Ahhg[(idp-1)*nst+tmst[idM],idM] += ahd[tmst[idM]] * rho_h[m2eh[idM]]
        end
    end

    return Ahhg

end
