# Dieno Diba 2024
# FITE2DMT
# Construct jacobian matrix for TE mode tipper

function ConstructJacobianTETipper(idp::Integer, nst::Integer, Aheg2::Matrix, M::Integer, el2no::Matrix, m2ee::Vector, omega::Number, mu0::Number, Area::Vector, rho::Vector, ce2::Matrix, ve::Vector, bxv2::Any)

    for idM = 1:M
        kiw = @views el2no[m2ee[idM],:]
        pDe = im*omega*mu0*Area[m2ee[idM]]/3/rho[m2ee[idM]]^2
        Aheg2[(idp-1)*nst+1:idp*nst,idM] = (bxv2.^-2) .* (-transpose(@view(ce2[kiw,:]))) * (pDe*ve[kiw]) * rho[m2ee[idM]] * log(10)
    end

    return Aheg2

end
