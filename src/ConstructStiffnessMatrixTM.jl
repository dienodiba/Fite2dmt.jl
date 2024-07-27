# Dieno Diba 2024
# FITE2DMT
# Construct stiffness matrix (FEM coefficient matrix) TM mode

function ConstructStiffnessMatrixTM(nel_h::Integer, nno_h::Integer,  el2no_h::Matrix, no2yz_h::Matrix, rho_h::Matrix, Area::Vector, eh2ee::Vector)

    Dh = zeros(ComplexF64,nno_h,nno_h)
    for ide = 1:nel_h
        # Nodes and their coordinates for the element ide
        no_ide = el2no_h[ide,:]
        yz_ide = no2yz_h[no_ide,:]
        # Edges
        s1_ide = yz_ide[3,:]-yz_ide[2,:]
        s2_ide = yz_ide[1,:]-yz_ide[3,:]
        s3_ide = yz_ide[2,:]-yz_ide[1,:]
        # Area of the element ide
        Ar_ide = Area[eh2ee[ide]]
        # Gradient of test functions
        grad_phi1 = [-s1_ide[2],s1_ide[1]]/(2*Ar_ide)
        grad_phi2 = [-s2_ide[2],s2_ide[1]]/(2*Ar_ide)
        grad_phi3 = [-s3_ide[2],s3_ide[1]]/(2*Ar_ide)
        grad_phi = hcat(grad_phi1,grad_phi2,grad_phi3)
        # Compute all the integrals for this particular element
        Dh[no_ide,no_ide] += -rho_h[ide] * grad_phi'*grad_phi * Ar_ide
    end

    return Dh

end
