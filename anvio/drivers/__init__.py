from anvio.drivers.fasttree import FastTree

driver_modules = {}
driver_modules['phylogeny'] = {"default":  FastTree,
                               "fasttree": FastTree}
