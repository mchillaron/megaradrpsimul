import teareduce as tea

# How to obtain a value from any keyword:
# print(CCDregions.regions_kernel["bias_topCCD"]["slice2d"])

# -------------------------------------------------------------------------------------------------

overscan3 = tea.SliceRegion2D("[51:4146, 2057:2156]", mode='fits')  # middle region
overscan3_bottom = tea.SliceRegion2D("[51:4146, 2057:2106]", mode='fits')  # middle region
overscan3_top = tea.SliceRegion2D("[51:4146, 2107:2156]", mode='fits')  # middle region

# -------------------------------------------------------------------------------------------------

topCCD_full = tea.SliceRegion2D("[1:4196, 2107:4212]", mode='fits') # with 8 upper pixels and 50 on each side
bottomCCD_full = tea.SliceRegion2D("[1:4196, 1:2106]", mode='fits') # with 8 lower pixels and 50 on each side

# -------------------------------------------------------------------------------------------------

# cosmic rays cleaning:
regions_cosmicrays = {
    "bias_topCCD": {
        'slice2d': tea.SliceRegion2D("[51:4146, 2107:4204]", mode='fits'), # without 8 upper pixels
        'median_size': (9,9), 
        'tsigma_peak': 5, 
        'tsigma_tail': 3
    },
    "bias_bottomCCD": {
        'slice2d': tea.SliceRegion2D("[51:4146, 9:2106]", mode='fits'), # without 8 lower pixels
        'median_size': (9,9), 
        'tsigma_peak': 5, 
        'tsigma_tail': 3
    }
}

# -------------------------------------------------------------------------------------------------
# Smoothing filters:
regions_kernel = {
    "bias_topCCD": {
        'slice2d': tea.SliceRegion2D("[51:4146, 2107:4212]", mode='fits'), # including 8 upper pixels
        'num_filters_SG': 2,
        'axis1_SG': 1,
        'size1_SG': 11,
        'pol_order1_SG': 1, 
        'axis2_SG': 0,
        'size2_SG': 11, 
        'pol_order2_SG': 2
    },
    "bias_bottomCCD": {
        'slice2d': tea.SliceRegion2D("[51:4146, 1:2106]", mode='fits'), # including 8 lower pixels 
        'num_filters_SG': 2,
        'axis1_SG': 1,
        'size1_SG': 11,
        'pol_order1_SG': 1, 
        'axis2_SG': 0,
        'size2_SG': 11, 
        'pol_order2_SG': 2
    },
    "overscan1_top": {
        'slice2d': tea.SliceRegion2D("[1:50, 2107:4212]", mode='fits'),
        'num_filters_SG': 2,
        'axis1_SG': 1,
        'size1_SG': 11,
        'pol_order1_SG': 1, 
        'axis2_SG': 0,
        'size2_SG': 11, 
        'pol_order2_SG': 1
    },
    "overscan1_bottom": {
        'slice2d': tea.SliceRegion2D("[1:50, 1:2106]", mode='fits'), 
        'num_filters_SG': 0,
        'median_size': (101, 1)
    },
    "overscan5_top": {
        'slice2d': tea.SliceRegion2D("[4147:4196, 2107:4212]", mode='fits'), 
        'num_filters_SG': 0,
        'median_size': (101, 1)
    },
    "overscan5_bottom": {
        'slice2d': tea.SliceRegion2D("[4147:4196, 1:2106]", mode='fits'),
        'num_filters_SG': 2,
        'axis1_SG': 1,
        'size1_SG': 11,
        'pol_order1_SG': 1, 
        'axis2_SG': 0,
        'size2_SG': 11, 
        'pol_order2_SG': 1
    }       
}
# -------------------------------------------------------------------------------------------------

