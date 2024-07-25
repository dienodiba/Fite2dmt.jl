# Dieno Diba 2024
# FITE2DMT
# Calculate model update with Gauss-newton data-space method

function CalculateModelUpdateDataSpace(invK::Matrix, A::Matrix, W::Matrix, V::Matrix, dobs::Vector, dcal::Vector ,lambda::Number, K::Any, m::Vector, nnzN::Vector)

    Grd = A[nnzN,:]'*(W[nnzN,nnzN]*(dobs[nnzN]-dcal[nnzN])) - lambda*K*m
    DS_A = V[nnzN,nnzN] + A[nnzN,:]*invK*A[nnzN,:]'
    DS_b = A[nnzN,:]*invK*Grd
    m += invK*Grd - invK*(A[nnzN,:]'*ldiv!(lu(DS_A),DS_b))

    return m

end
