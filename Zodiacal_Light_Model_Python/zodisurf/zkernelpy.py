# If any old scripts try to call zkernelpy.zkernel(), they hit this wrapper, 
# which instantly redirects the math over to your ultra-fast zkernel_fast.py engine.
def zkernel(data, a, phase_type='kelsall', L2=False, indxpar=None, losinfo=False,
            no_colcorr=True, solar_irr=None, new_iso_comp=False, iso_comp_only=False):
    """
    Compatibility wrapper: forwards all requests to the fast vectorized kernel.
    """
    if losinfo or indxpar is not None:
        raise NotImplementedError("losinfo/partials: use zodisurf.zkernelpy_original.")
    
    from zodisurf.zkernel_fast import zkernel_fast   # imported here to avoid a circular import
    
    return zkernel_fast(data, a, phase_type=phase_type, L2=L2, no_colcorr=no_colcorr,
                        solar_irr=solar_irr, new_iso_comp=new_iso_comp,
                        iso_comp_only=iso_comp_only)