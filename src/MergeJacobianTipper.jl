# Dieno Diba 2024
# FITE2DMT
# Merge jacobian matrix for tipper inversion

function MergeJacobianTipper(N::Integer, M::Integer, nper::Integer, nst::Integer, Aheg2::Matrix)

    A = zeros(Float64,N,M)
    idN = 0
    id2 = 0    
    for idp = 1:nper
        for ids = 1:nst
            for idj = 1:2
                idN += 1
                if idj == 1
                    id2 += 1
                    A[idN,:] = real(Aheg2[id2,:])
                end
                if idj == 2
                    A[idN,:] = imag(Aheg2[id2,:])
                end
            end
        end
    end

    return A

end
