# Dieno Diba 2024
# FITE2DMT
# Calculate TM mode apparent resistivity and phase at the stations

function CalculateResponseFunctionsTMRhoaPhase(idp::Integer, dcalh::Matrix, nno_h::Integer, nst::Integer, no2yz_h::Matrix, ns1h::Vector, ns2h::Vector, nbsh::Vector, rho_h::Matrix, elsta::Vector, omega::Number, mu0::Number, vh::Vector)

    # ah, bh: coefficients to interpolate and/or differentiate fields
    # ch: required for pseudo-forward computation
    ch = zeros(ComplexF64,nno_h,nst)
    ahd = zeros(ComplexF64,nst)
    for ids = 1:nst
        a = sqrt((no2yz_h[ns1h[ids],1]-no2yz_h[nbsh[ids],1])^2 + (no2yz_h[ns1h[ids],2]-no2yz_h[nbsh[ids],2])^2)
        b = sqrt((no2yz_h[ns2h[ids],1]-no2yz_h[nbsh[ids],1])^2 + (no2yz_h[ns2h[ids],2]-no2yz_h[nbsh[ids],2])^2)
        c = sqrt((no2yz_h[ns1h[ids],1]-no2yz_h[ns2h[ids],1])^2 + (no2yz_h[ns1h[ids],2]-no2yz_h[ns2h[ids],2])^2)
        C = acos((a^2 + b^2 - c^2)/(2 * a * b))
        h = a*b/c * sin(C)
        ah = zeros(3)
        bh = zeros(3)
        bh[1] = 0.5
        bh[2] = 0.5
        slope = abs(atan((no2yz_h[ns2h[ids],2] - no2yz_h[ns1h[ids],2])/(no2yz_h[ns2h[ids],1] - no2yz_h[ns1h[ids],1])))
        if no2yz_h[ns1h[ids],2] < no2yz_h[ns2h[ids],2]
            egv = h * sin(slope) / sin(pi - slope - acos(h/b))
            egw = h * sin(acos(h/b)) / sin(pi - slope - acos(h/b))
            ah[3] = rho_h[elsta[ids]]/egw * (1 - egv/b)
            ah[1] = rho_h[elsta[ids]]/egw * (egv/b - 1)
        else
            egv = h * sin(slope) / sin(pi - slope - acos(h/a))
            egw = h * sin(acos(h/a)) / sin(pi - slope - acos(h/a))
            ah[3] = rho_h[elsta[ids]]/egw * (1 - egv/a)
            ah[1] = rho_h[elsta[ids]]/egw * (egv/a - 1)
        end
        dcalh[idp,ids] = log10(sqrt(1/omega/mu0) * (transpose(ah)*vh[[ns1h[ids],ns2h[ids],nbsh[ids]]]) /
         (transpose(bh)*vh[[ns1h[ids],ns2h[ids],nbsh[ids]]]))
        ch[[ns1h[ids],ns2h[ids],nbsh[ids]],ids] = ah/(transpose(ah)*vh[[ns1h[ids],ns2h[ids],nbsh[ids]]]) -
         bh/(transpose(bh)*vh[[ns1h[ids],ns2h[ids],nbsh[ids]]])
        ahd[ids] = transpose(ah/rho_h[elsta[ids]]/(transpose(ah)*vh[[ns1h[ids],ns2h[ids],nbsh[ids]]]))*vh[[ns1h[ids],ns2h[ids],nbsh[ids]]]
    end

    return dcalh,ahd,ch

end
