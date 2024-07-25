# Dieno Diba 2024
# FITE2DMT
# Read inversion setting parameters from input file
# Rhoa, phase inversion

function ReadInversionSettingParametersRhoaPhase(filestg::String)

    # required keywords
    stg = readdlm(filestg)
    lambda = stg[only(findall(x->x=="TRADEOFF",stg[:,1]))+1,1]
    max_iter = stg[only(findall(x->x=="MAXITER",stg[:,1]))+1,1]
    outroot = stg[only(findall(x->x=="OUTFILE",stg[:,1]))+1,1]
    fix_n = stg[only(findall(x->x=="FIXBLOCK",stg[:,1]))+1,1]
    fix_el = stg[only(findall(x->x=="FIXBLOCK",stg[:,1]))+2,1:fix_n]
    # optional keywords with default values
    wg_rte = 1
    if isempty(findall(x->x=="USE_RTE",stg[:,1])) == 0
        damping = stg[only(findall(x->x=="USE_RTE",stg[:,1]))+1,1]
    end
    wg_pte = 1
    if isempty(findall(x->x=="USE_PTE",stg[:,1])) == 0
        damping = stg[only(findall(x->x=="USE_PTE",stg[:,1]))+1,1]
    end
    wg_rtm = 1
    if isempty(findall(x->x=="USE_RTM",stg[:,1])) == 0
        damping = stg[only(findall(x->x=="USE_RTM",stg[:,1]))+1,1]
    end
    wg_ptm = 1
    if isempty(findall(x->x=="USE_PTM",stg[:,1])) == 0
        damping = stg[only(findall(x->x=="USE_PTM",stg[:,1]))+1,1]
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
    ef_zte = 1e-10
    if isempty(findall(x->x=="ERF_ZTE",stg[:,1])) == 0
        ef_zte = stg[only(findall(x->x=="ERF_ZTE",stg[:,1]))+1,1]
    end
    ef_ztm = 1e-10
    if isempty(findall(x->x=="ERF_ZTM",stg[:,1])) == 0
        ef_ztm = stg[only(findall(x->x=="ERF_ZTM",stg[:,1]))+1,1]
    end

    return lambda,max_iter,outroot,fix_el,damping,hor2ver,static,nodata,wg_rte,wg_pte,wg_rtm,wg_ptm,ef_zte,ef_ztm

end
