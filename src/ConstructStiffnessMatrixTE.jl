# Dieno Diba 2024
# FITE2DMT
# Construct stiffness matrix (FEM coefficient matrix) TE mode

function ConstructStiffnessMatrixTE(nel::Integer, nno::Integer,  el2no::Matrix, no2yz::Matrix, Area::Vector)

    De = zeros(ComplexF64,nno,nno)
    for ide = 1:nel
        # Nodes and their coordinates for the element ide
        no_ide = el2no[ide,:]
        yz_ide = no2yz[no_ide,:]
        # Edges
        s1_ide = yz_ide[3,:] - yz_ide[2,:]
        s2_ide = yz_ide[1,:] - yz_ide[3,:]
        s3_ide = yz_ide[2,:] - yz_ide[1,:]
        # Gradient of test functions
        grad_phi1 = [-s1_ide[2],s1_ide[1]]/(2 * Area[ide])
        grad_phi2 = [-s2_ide[2],s2_ide[1]]/(2 * Area[ide])
        grad_phi3 = [-s3_ide[2],s3_ide[1]]/(2 * Area[ide])
        grad_phi = hcat(grad_phi1,grad_phi2,grad_phi3)
        # Compute all the integrals for this particular element
        Ar_ide = Area[ide]
        De[no_ide,no_ide] += -(grad_phi'*grad_phi) * Ar_ide
    end

    return De

end
