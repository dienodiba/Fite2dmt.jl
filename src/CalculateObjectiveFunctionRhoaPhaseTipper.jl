# Dieno Diba 2024
# FITE2DMT
# Calculate data misfit, model roughness, and objective function
# Rhoa phase tipper inversion

function CalculateObjectiveFunctionRhoaPhaseTipper(RMSd1::Vector, RMSd2::Vector, RMSd4::Vector, Phid1::Vector, Phid2::Vector, Phid4::Vector, Phi::Vector, RMS::Vector, Stab::Vector, Psi::Vector, N::Integer, itr::Integer, nper::Integer, nst::Integer, dcale1::Matrix, dcale2::Matrix, dcalh::Matrix, dobs::Vector, W::Matrix, T::Vector, lambda::Number, m::Vector, K::Any, idd1::Integer, idd2::Integer, idd4::Integer)

    dcal = zeros(Float64,N)
    idN = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:6
                idN += 1
                if idj == 1
                    dcal[idN] = real(dcale1[idp,ids])
                    Phid1[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd1[itr] += (2*dcal[idN] - 2*dobs[idN])^2*W[idN,idN]*T[idN]
                end
                if idj == 2
                    dcal[idN] = imag(dcale1[idp,ids])
                    Phid1[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd1[itr] += (log(10)*dcal[idN] - log(10)*dobs[idN])^2*W[idN,idN]*T[idN]
                end
                if idj == 3
                    dcal[idN] = real(dcale2[idp,ids])
                    Phid2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                end
                if idj == 4
                    dcal[idN] = imag(dcale2[idp,ids])
                    Phid2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd2[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                end
                if idj == 5
                    dcal[idN] = real(dcalh[idp,ids])
                    Phid4[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd4[itr] += (2*dcal[idN] - 2*dobs[idN])^2*W[idN,idN]*T[idN]
                end
                if idj == 6
                    dcal[idN] = imag(dcalh[idp,ids])
                    Phid4[itr] += (dcal[idN] - dobs[idN])^2*W[idN,idN]*T[idN]
                    RMSd4[itr] += (log(10)*dcal[idN] - log(10)*dobs[idN])^2*W[idN,idN]*T[idN]
                end
            end
        end
    end
    # Calculate total data misfit
    Phi[itr] = Phid1[itr] + Phid2[itr] + Phid4[itr]
    # Calculate model rouhgness
    Stab[itr] = lambda * m' * K * m
    # Calculate total objective function
    Psi[itr] = Phi[itr] + Stab[itr]
    # Calculate combined RMS
    RMS[itr] = sqrt((RMSd1[itr]+RMSd2[itr]+RMSd4[itr])/(idd1+idd2+idd4))
    # Calculate individual RMS
    RMSd1[itr] = sqrt(RMSd1[itr]/idd1)
    RMSd2[itr] = sqrt(RMSd2[itr]/idd2)
    RMSd4[itr] = sqrt(RMSd4[itr]/idd4)

    return dcal,RMSd1,RMSd2,RMSd4,Phid1,Phid2,Phid4,Phi,RMS,Stab,Psi
end
