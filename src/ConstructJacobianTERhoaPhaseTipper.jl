# Dieno Diba 2024
# FITE2DMT
# Construct jacobian matrix for TE mode apparent resistivity, phase, and tipper

function ConstructJacobianTERhoaPhaseTipper(idp::Integer, nst::Integer, Aheg1::Matrix, Aheg2::Matrix, M::Integer, el2no::Matrix, m2ee::Vector, omega::Number, mu0::Number, Area::Vector, rho::Vector, ce1::Matrix, ce2::Matrix, ve::Vector, bxv2::Any)

    for idM = 1:M
        kiw = @views el2no[m2ee[idM],:]
        pDe = im*omega*mu0*Area[m2ee[idM]]/3/rho[m2ee[idM]]^2
        Aheg1[(idp-1)*nst+1:idp*nst,idM] = -transpose(@view(ce1[kiw,:])) * (pDe*@view(ve[kiw])) * rho[m2ee[idM]]
        Aheg2[(idp-1)*nst+1:idp*nst,idM] = (bxv2.^-2) .* (-transpose(@view(ce2[kiw,:]))) * (pDe*ve[kiw]) * rho[m2ee[idM]] * log(10)
    end

    return Aheg1,Aheg2

end
