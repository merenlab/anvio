var named_layers = {
    '__parent__': {
        'pretty_name': 'Parent',
    },
    'taxonomy': {
        'pretty_name': 'Taxonomy',
    },
    'num_genes': {
        'height': 0,
        'color': '#414141',
        'norm': 'none',
        'pretty_name': 'Number of genes',
    },
    'avg_gene_length': {
        'height': 0,
        'color': '#414141',
        'norm': 'none',
        'pretty_name': 'Avg. gene length',
    },
    'ratio_coding': {
        'height': 0,
        'color': '#414141',
        'norm': 'none',
        'pretty_name': 'Ratio coding region',
    },
    'ratio_hypothetical': {
        'height': 0,
        'color': '#072c59',
        'norm': 'none',
        'pretty_name': 'Ratio hypothetical',
    },
    'ratio_with_tax': {
        'height': 80,
        'color': '#555907',
        'norm': 'none',
        'pretty_name': 'Ratio w/ taxonomy',
    },
    'tax_accuracy': {
        'height': 0,
        'color': '#40012e',
        'norm': 'none',
        'pretty_name': 'Accuracy of taxonomy',
    },
    'length': {
        'height': 80,
        'color': '#414141',
        'norm': 'none',
        'pretty_name': 'Length',
        'type': 'bar',
    },
    'gc_content': {
        'height': 180,
        'color': '#004a0f',
        'norm': 'none',
        'pretty_name': 'GC-content',
        'type': 'bar',
    },
    'num_genes_per_kb': {
        'height': 180,
        'color': '#414141',
        'norm': 'none',
        'pretty_name': 'Num genes per kbp',
    },
    'num_genomes_pc_has_hits': {
        'height': 90,
        'color': '#002e4f',
        'norm': 'none',
        'pretty_name': 'Num genomes',
    },
    'num_genes_in_pc': {
        'height': 90,
        'color': '#002e4f',
        'norm': 'sqrt',
        'pretty_name': 'Num genes in PC',
    },
    'total_length': {
        'height': 180,
        'color': '#121261',
        'norm': 'none',
        'pretty_name': 'Total length',
    },
    'normalized_coverage': {
        'height': 0,
        'color': '#414141',
        'norm': 'log',
        'pretty_name': 'Normalized coverage',
    },
    'std_coverage': {
        'height': 80,
        'color': '#414141',
        'norm': 'log',
        'pretty_name': 'Coverage STD',
    },
    'portion_covered': {
        'height': 0,
        'color': '#616161',
        'norm': 'none',
        'pretty_name': 'Percent covered',
    },
    'abundance': {
        'height': 0,
        'color': '#818181',
        'norm': 'none',
        'pretty_name': 'Abundance',
    },
    'mean_coverage': {
        'height': 300,
        'color': '#141414',
        'norm': 'log',
        'pretty_name': 'Mean coverage',
    },
    'max_normalized_ratio': {
        'height': 0,
        'color': '#141414',
        'norm': 'log',
        'pretty_name': 'Max-normalized ratio',
    },
    'relative_abundance': {
        'height': 0,
        'color': '#141414',
        'norm': 'log',
        'pretty_name': 'Relative abundance',
    },
    'variability': {
        'height': 180,
        'color': '#4a000f',
        'norm': 'none',
        'pretty_name': 'Variablity',
    },
    'percent_completion': {
        'height': 200,
        'color': '#004400',
        'norm': 'none',
        'min': 0,
        'max': 100,
        'max_disabled': false,
        'min_disabled': false,
        'pretty_name': 'Completion',
        'type': 'bar',
    },
    'percent_redundancy': {
        'height': 200,
        'color': '#440000',
        'norm': 'none',
        'min': 0,
        'max': 100,
        'max_disabled': false,
        'min_disabled': false,
        'pretty_name': 'Redundancy',
        'type': 'bar',
    },
    'bin_name': {
        'type': 'text',
        'pretty_name': 'Bin name',
    },
    'blank_view': {
        'height': 0,
        'color': '#FFFFFF',
        'norm': 'none',
        'pretty_name': '_',
    },
};

pretty_names = {
    'tnf-cov': 'Seq. Composition + Diff. Coverage',
    'cov': 'Differential coverage',
    'tnf': 'Sequence composition',
    'tnf-splits': 'Sequence composition (w/independent splits)'
};

function getPrettyName(name)
{

    if (['none', 'sqrt', 'log', 'bar', 'intensity'].indexOf(name) >= 0){
        return name;
    }

    if (name in named_layers){
        if ('pretty_name' in named_layers[name]){
            return named_layers[name]['pretty_name']
        }
    }

    if (name in pretty_names){
        return pretty_names[name]
    }

    name = name.replace(/_/g, " ").replace(/-/g, " ");
    name = name.charAt(0).toUpperCase() + name.slice(1);

    return name;
}

function getNamedLayerDefaults(layer, attribute, default_value)
{
    if (typeof default_value == "string" && default_value.charAt(0) != '#'){
        default_value = getPrettyName(default_value)
    }

    /* Some ad-hoc manipulation of special hmmx_ split hmm layers */
    if (layer.substring(0, 5) == "hmmx_") {
        if (attribute == 'height') return '30';
        if (attribute == 'norm')   return 'none';
        if (attribute == 'color')  return '#882222'
    }

    /* Some ad-hoc manipulation of special hmms_ single hmm layers */ 
    if (layer.substring(0, 5) == "hmms_"){
        if (attribute == 'height') return '150';
        if (attribute == 'norm')   return 'sqrt';
        if (attribute == 'color')  return '#882222'
    }

    if (layer in named_layers)
    {
        if (attribute in named_layers[layer])
        {
            return named_layers[layer][attribute];
        }
        else if (attribute == 'color')
        {
            // layer exists but no color
            return randomColor();
        }
        return default_value;
    }
    return default_value;
}
