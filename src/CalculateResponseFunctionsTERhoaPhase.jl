# Dieno Diba 2024
# FITE2DMT
# Calculate TE mode apparent resistivity, phase at the stations

function CalculateResponseFunctionsTERhoaPhase(idp::Integer, dcale1::Matrix, nno::Integer, nst::Integer, no2yz::Matrix, ns1e::Vector, ns2e::Vector, nbse::Vector, omega::Number, mu0::Number, ve::Vector)

    # ae, be: coefficients to interpolate and/or differentiate fields
    # ce: required for pseudo-forward computation
    ce1 = zeros(ComplexF64,nno,nst)
    for ids = 1:nst
        a = sqrt((no2yz[ns1e[ids],1]-no2yz[nbse[ids],1])^2 + (no2yz[ns1e[ids],2]-no2yz[nbse[ids],2])^2)
        b = sqrt((no2yz[ns2e[ids],1]-no2yz[nbse[ids],1])^2 + (no2yz[ns2e[ids],2]-no2yz[nbse[ids],2])^2)
        c = sqrt((no2yz[ns1e[ids],1]-no2yz[ns2e[ids],1])^2 + (no2yz[ns1e[ids],2]-no2yz[ns2e[ids],2])^2)
        C = acos((a^2 + b^2 - c^2)/(2 * a * b))
        h = a*b/c * sin(C)
        d = a * cos(asin(h/a))
        bY = zeros(ComplexF64,3)
        slope = abs(atan((no2yz[ns2e[ids],2] - no2yz[ns1e[ids],2])/(no2yz[ns2e[ids],1] - no2yz[ns1e[ids],1])))
        if no2yz[ns1e[ids],2] < no2yz[ns2e[ids],2]
            egv = h * sin(slope) / sin(pi - slope - acos(h/b))
            egw = h * sin(acos(h/b)) / sin(pi - slope - acos(h/b))
            bY[3] = -1/im/omega/mu0/egw * (1 - egv/b)
            bY[2] = -1/im/omega/mu0/egw * egv/b
        else
            egv = h * sin(slope) / sin(pi - slope - acos(h/a))
            egw = h * sin(acos(h/a)) / sin(pi - slope - acos(h/a))
            bY[3] = -1/im/omega/mu0/egw * (1 - egv/a)
            bY[1] = -1/im/omega/mu0/egw * egv/a
        end
        if d > c    # obtuse triangle beneath the iyy-ref
            bY[1] = bY[1] + -1 * -1/im/omega/mu0/egw * 0.5
            bY[2] = bY[2] + -1 * -1/im/omega/mu0/egw * 0.5
        else        # acute triangle beneath the iyy-ref
            bY[1] = bY[1] + -1 * -1/im/omega/mu0/egw * (1 - d/c)
            bY[2] = bY[2] + -1 * -1/im/omega/mu0/egw * d/c
        end
        # rho-app and phase
        ae1 = zeros(3)
        ae1[1] = 0.5
        ae1[2] = 0.5
        ce1[[ns1e[ids],ns2e[ids],nbse[ids]],ids] = ae1/(transpose(ae1)*ve[[ns1e[ids],ns2e[ids],nbse[ids]]]) -
            bY/(transpose(bY)*ve[[ns1e[ids],ns2e[ids],nbse[ids]]])
        dcale1[idp,ids] = log10(sqrt(1/omega/mu0)*(transpose(ae1)*ve[[ns1e[ids],ns2e[ids],nbse[ids]]])/
            (transpose(bY)*ve[[ns1e[ids],ns2e[ids],nbse[ids]]]))
    end

    return dcale1,ce1

end
