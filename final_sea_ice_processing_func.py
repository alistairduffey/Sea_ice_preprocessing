
def get_SIa_SI_e(ds, areacell_path=None, 
                 model=None, grid_type='ocean', 
                 var=None, siconc_in_percent=True,
                 hemisphere='NH', UKESM_MASS_data=False):
    """ function takes in a dataset containing siconc/siconca,
    and returns a timeseries dataset of sea ice area and extent 
    by hemisphere.  
    
    INPUTS
    ds: a dataset containing sea ice concentration
    areacell_path (optional): file path to cell areas .nc file
    model (optional): if no areacell_path, must specify a model to look up grid areas file on JASMSIN
    grid_type: 'ocean' (recommended) or 'atmosphere'. The only reason to use 'atmosphere' is to..
    .. compare the two options as the pre-computed cmip6 regridding does not conserve sea ice area or extent 
    var (optional): 'siconc' or 'siconca', but this is worked out from grid_type if not given
    siconc_in_percent: assumed true. This is not the case for direct MASS imported UKESM data
    hemisphere: 'NH' or 'SH', selects hemipshere for outputs.
    UKESM_MASS_data: True for data imported from MASS without CMIP6 CEDA processsing. Selecting..
    .. True means minor differences which break coordinate matching are accounted for. 
    """

    # make sure our naming conventions are correct, 
    # and that we have a nominal latitude with which to 
    # select hemisphere (nominal lat is not used otherwise) 
    if grid_type=='ocean':
        ds = replace_x_y_nominal_lat_lon(correct_lon(rename_cmip6(ds)))
    else:
        ds = rename_cmip6(ds)
                     
    # set up some dicts
    area_strings = {'ocean':'Ofx/areacello', 
                'atmosphere':'fx/areacella'}
                     
    var_strings = {'ocean':'siconc',
                   'atmosphere':'siconca'}
    cell_area_var_strings = {'ocean':'areacello',
                             'atmosphere':'areacella'}
    hemisphere_slices = {'NH':[0, 90],
                         'SH':[-90, 0]}
                     
    # select variable names:
    if not var:        
        var = var_strings[grid_type]
    cell_area_var = cell_area_var_strings[grid_type]
        
    # find grid cell area filepath (on jasmin):
    if not areacell_path:
        folder = glob.glob('/badc/cmip6/data/CMIP6/CMIP/*/{m}/*/*/{a}/*/latest/'.format(
                    m=model, a=area_strings[grid_type]))[0]
        file = os.listdir(folder)[0]
        areacell_path = folder+file
       
    # read in cell_areas
    # if data comes with tarea (MASS does) then use that:
    if 'tarea' in list(ds.data_vars):
        cell_areas = ds.tarea.to_dataset().rename({'tarea':'areacello'})
    else:
        cell_areas = rename_cmip6(xr.open_dataset(areacell_path))
        if grid_type=='ocean':
            cell_areas = replace_x_y_nominal_lat_lon(correct_lon(cell_areas))
    
    # convert si concentration into fraction if in percent
    if siconc_in_percent:
        ds[var] = ds[var]*0.01

    # define unit conversion m2 to M_km2:
    unit_conv = 10**12

    # select hemipshere and calc area and extent
    lats = hemisphere_slices[hemisphere]                 
    ds = ds.sel(y=slice(lats[0], lats[1]))
    cell_areas = cell_areas.sel(y=slice(lats[0], lats[1])) 

    # need to cheat slightly for the UKESM runs that are direct exports from MASS
    # this code accounts for v slight differences in nominal lat/lon as calculated by xmip
    # use with caution!! 

    if UKESM_MASS_data:
        # check validity of replacing coords:
        if np.max((ds.y - cell_areas.y) > 0.1):
            print('WARNING: replacing coords in ds with significantly different ones from cell areas.')
        if np.max((ds.x - cell_areas.x) > 0.1):
            print('WARNING: replacing coords in ds with significantly different ones from cell areas.')
        ds['x'] = cell_areas.x
        ds['y'] = cell_areas.y

    SI_area = ds[var]*cell_areas[cell_area_var]/unit_conv
    SI_area = SI_area.sum(dim=['x', 'y']).to_dataset(name='si_area_{}_Mkm2'.format(hemisphere))

    SI_extent_mask = xr.where(ds[var]>0.15, 1, 0)
    SI_extent = SI_extent_mask*cell_areas[cell_area_var]/unit_conv
    SI_extent = SI_extent.sum(dim=['x', 'y']).to_dataset(name='si_extent_{}_Mkm2'.format(hemisphere))                 

    out_ds = xr.merge([SI_area, SI_extent])
    out_ds['Grid_type']=grid_type
    return out_ds
