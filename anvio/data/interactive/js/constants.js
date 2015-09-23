var named_layers = {
	'__parent__': {
		'pretty_name': 'Parent',
	},
	'taxonomy': {
		'pretty_name': 'Taxonomy',
	},
	'num_genes': {
		'height': 80,
		'color': '#414141',
		'norm': 'sqrt',
        'pretty_name': 'Number of genes',
	},
	'avg_gene_length': {
		'height': 0,
		'color': '#414141',
		'norm': 'sqrt',
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
	},
	'gc_content': {
		'height': 180,
		'color': '#004a0f',
		'norm': 'none',
		'pretty_name': 'GC-content',
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
};

function getNamedLayerDefaults(layer, attribute, default_value)
{
    if (typeof default_value == "string" && default_value.charAt(0) != '#'){
        // make the default value a bit prettier, if possible.
        default_value = default_value.replace(/_/g, " ");
        default_value = default_value.charAt(0).toUpperCase() + default_value.slice(1);
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
