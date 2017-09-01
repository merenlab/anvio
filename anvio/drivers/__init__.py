import anvio.drivers as drivers

driver_modules = {}
driver_modules['phylogeny'] = {"default":  drivers.FastTree,
                               "fasttree": drivers.FastTree}

driver_modules['psa'] = {"default": drivers.muscle,
                         "muscle": drivers.muscle}
