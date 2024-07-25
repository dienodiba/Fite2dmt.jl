# Dieno Diba 2024
# FITE2DMT
# Read inversion setting parameters from input file
# Rhoa, phase, tipper inversion

function ReadInversionSettingParametersTipper(filestg::String)

    # required keywords
    stg = readdlm(filestg)
    lambda = stg[only(findall(x->x=="TRADEOFF",stg[:,1]))+1,1]
    max_iter = stg[only(findall(x->x=="MAXITER",stg[:,1]))+1,1]
    outroot = stg[only(findall(x->x=="OUTFILE",stg[:,1]))+1,1]
    fix_n = stg[only(findall(x->x=="FIXBLOCK",stg[:,1]))+1,1]
    fix_el = stg[only(findall(x->x=="FIXBLOCK",stg[:,1]))+2,1:fix_n]
    # optional keywords with default values
    wg_tyr = 1
    if isempty(findall(x->x=="USE_TYR",stg[:,1])) == 0
        damping = stg[only(findall(x->x=="USE_TYR",stg[:,1]))+1,1]
    end
    wg_tyi = 1
    if isempty(findall(x->x=="USE_TYI",stg[:,1])) == 0
        damping = stg[only(findall(x->x=="USE_TYI",stg[:,1]))+1,1]
    end
    damping = 0.1
    if isempty(findall(x->x=="DAMPING",stg[:,1])) == 0
        damping = stg[only(findall(x->x=="DAMPING",stg[:,1]))+1,1]
    end
    hor2ver = 1
    if isempty(findall(x->x=="HOR2VER",stg[:,1])) == 0
        hor2ver = stg[only(findall(x->x=="HOR2VER",stg[:,1]))+1,1]
    end
    static = 0
    if isempty(findall(x->x=="STATIC",stg[:,1])) == 0
        static = stg[only(findall(x->x=="STATIC",stg[:,1]))+1,1]
    end
    nodata = 99999
    if isempty(findall(x->x=="NODATA",stg[:,1])) == 0
        static = stg[only(findall(x->x=="NODATA",stg[:,1]))+1,1]
    end
    ef_tzy = 1e-10
    if isempty(findall(x->x=="ERF_TZY",stg[:,1])) == 0
        ef_tzy = stg[only(findall(x->x=="ERF_TZY",stg[:,1]))+1,1]
    end

    return lambda,max_iter,outroot,fix_el,damping,hor2ver,static,nodata,wg_tyr,wg_tyi,ef_tzy

end
