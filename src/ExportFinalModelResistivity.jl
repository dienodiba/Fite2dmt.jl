# Dieno Diba 2024
# FITE2DMT
# Export final inversion model into a txt file

function ExportFinalModelResistivity(outroot::Any, nel::Integer, nno::Integer, el2no::Matrix, no2yz::Matrix, frho::Vector)

    fid = open(string(outroot,"_res.txt"),"w")
    @printf(fid,"nel = %10i\n",nel)
    @printf(fid,"nno = %10i\n\n",nno)
    @printf(fid,"EL2NO\n")
    for ide = 1:nel
        @printf(fid,"%12i%12i%12i\n",el2no[ide,1],el2no[ide,2],el2no[ide,3])
    end
    @printf(fid,"\n")
    @printf(fid,"NO2YZ\n")
    for idn = 1:nno
        @printf(fid,"%12.4e%12.4e\n",1e-3*no2yz[idn,1],1e-3*no2yz[idn,2])
    end
    @printf(fid,"\n")
    @printf(fid,"RESISTIVITY\n")
    for ide = 1:nel
        @printf(fid,"%12.4e\n",frho[ide])
    end
    close(fid)

    return

end 
