# Dieno Diba 2024
# FITE2DMT
# Calculate data misfit, model roughness, and objective function
# Tipper inversion

function CalculateObjectiveFunctionTipper(RMSd2::Vector, Phid2::Vector, Phi::Vector, RMS::Vector, Stab::Vector, Psi::Vector, N::Integer, itr::Integer, nper::Integer, nst::Integer, dcale2::Matrix, dobs::Vector, W::Matrix, T::Vector, lambda::Number, m::Vector, K::Any, idd2::Integer)

    dcal = zeros(Float64,N)
    idN = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:2
                idN += 1
                if idj == 1
                    dcal[idN] = real(dcale2[idp,ids])
                    Phid2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                end
                if idj == 2
                    dcal[idN] = imag(dcale2[idp,ids])
                    Phid2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                end
            end
        end
    end
    # Calculate total data misfit
    Phi[itr] = Phid2[itr]
    # Calculate model rouhgness
    Stab[itr] = lambda * m' * K * m
    # Calculate total objective function
    Psi[itr] = Phi[itr] + Stab[itr]
    # Calculate combined RMS
    RMS[itr] = sqrt((RMSd2[itr])/(idd2))
    # Calculate individual RMS
    RMSd2[itr] = sqrt(RMSd2[itr]/idd2)

    return dcal,RMSd2,Phid2,Phi,RMS,Stab,Psi
end
