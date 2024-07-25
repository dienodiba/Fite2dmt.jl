# Dieno Diba 2024
# FITE2DMT
# Merge jacobian matrix for rhoa phase inversion

function MergeJacobianRhoaPhase(N::Integer, M::Integer, nper::Integer, nst::Integer, Aheg1::Matrix, Ahhg::Matrix)

    A = zeros(Float64,N,M)
    idN = 0
    id1 = 0
    id4 = 0
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:4
                idN += 1
                if idj == 1
                    id1 += 1
                    A[idN,:] = real(Aheg1[id1,:])
                end
                if idj == 2
                    A[idN,:] = imag(Aheg1[id1,:])
                end
                if idj == 3
                    id4 += 1
                    A[idN,:] = real(Ahhg[id4,:])
                end
                if idj == 4
                    A[idN,:] = imag(Ahhg[id4,:])
                end
            end
        end
    end

    return A

end
