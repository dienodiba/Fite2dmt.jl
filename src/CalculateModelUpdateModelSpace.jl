# Dieno Diba 2024
# FITE2DMT
# Calculate model update with Gauss-newton model-space method

function (m::Vector, A::Matrix, W::Matrix, lambda::Number, K::Any, damping::Number, Psi::Number, M::Integer, dobs::Vector, dcal::Vector, nnzN::Vector)

    MS_A = A[nnzN,:]'*W[nnzN,nnzN]*A[nnzN,:] + lambda*K + damping*Psi*I(M)
    MS_b = -A[nnzN,:]'*W[nnzN,nnzN]*(dobs[nnzN]-dcal[nnzN]) + lambda*K*m
    m -= MS_A\MS_b

    return m

end
