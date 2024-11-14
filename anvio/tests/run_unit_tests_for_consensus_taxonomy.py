#!/usr/bin/env python
# some tests for SCG taxonomy string processing

import argparse

import anvio.terminal as terminal
import anvio.taxonomyops as taxonomyops
import anvio.taxonomyops.scg as scgtaxonomyops

levels_of_taxonomy = [
    "t_domain",
    "t_phylum",
    "t_class",
    "t_order",
    "t_family",
    "t_genus",
    "t_species",
]

c = scgtaxonomyops.PopulateContigsDatabaseWithSCGTaxonomy(
    argparse.Namespace(skip_sanity_check=True), run=terminal.Run(verbose=False)
)

p = scgtaxonomyops.SCGTaxonomyEstimatorSingle(
    argparse.Namespace(skip_sanity_check=True, skip_init=True),
    run=terminal.Run(verbose=False),
)

cX = lambda: c.get_consensus_hit(scg_raw_hits)
cT = lambda level: cX()[level]


def pX(scg_dict):
    for i in scg_dict:
        scg_dict[i]["tax_hash"] = taxonomyops.HASH(scg_dict[i])
    return p.get_consensus_taxonomy(scg_dict)


pT = lambda level: pX(scg_dict)[level]


#########################################

scg_raw_hits = [
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
]

assert cT("t_species") == "G x"
assert cT("t_genus") == "F"

#########################################

scg_raw_hits = [
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "T x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "T x",
    },
]

assert cT("t_species") == None
assert cT("t_genus") == "F"

#########################################

scg_raw_hits = [
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "W",
        "t_species": "T x",
    },
]

assert cT("t_species") == None
assert cT("t_genus") == None
assert cT("t_family") == "E"

#########################################

scg_raw_hits = [
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 99.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 99.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "W",
        "t_species": "T x",
    },
]

assert cT("t_species") == "G x"

#########################################

scg_raw_hits = [
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": None,
    },
    {
        "percent_identity": 100.0,
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": None,
    },
]

assert cT("t_species") == "G x"


#########################################

scg_dict = {
    1: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    3: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "T x",
    },
}

assert pT("t_species") == None
assert pT("t_genus") == "F"

#########################################

scg_dict = {
    1: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    2: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    3: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "W",
        "t_species": "T x",
    },
}

assert pT("t_species") == "G x"

#########################################

scg_dict = {
    1: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    2: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "F",
        "t_species": "G x",
    },
    3: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "J",
        "t_genus": "X",
        "t_species": "J x",
    },
    4: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "J",
        "t_genus": "X",
        "t_species": "J x",
    },
    5: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "W",
        "t_species": "T x",
    },
    6: {
        "t_domain": "A",
        "t_phylum": "B",
        "t_class": "C",
        "t_order": "D",
        "t_family": "E",
        "t_genus": "W",
        "t_species": "T x",
    },
}

assert pT("t_species") == None
assert pT("t_genus") == None
assert pT("t_family") == None
assert pT("t_order") == "D"
