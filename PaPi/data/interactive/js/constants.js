var named_layers = {
	'num_genes': {
		'height': 80,
		'color': '#414141',
		'norm': 'sqrt',
	},
	'avg_gene_length': {
		'height': 80,
		'color': '#414141',
		'norm': 'sqrt',
	},
	'ratio_coding': {
		'height': 80,
		'color': '#414141',
		'norm': 'none',
	},
	'ratio_hypothetical': {
		'height': 80,
		'color': '#072c59',
		'norm': 'none',
	},
	'ratio_with_tax': {
		'height': 80,
		'color': '#555907',
		'norm': 'none',
	},
	'tax_accuracy': {
		'height': 80,
		'color': '#40012e',
		'norm': 'none',
	},
	'length': {
		'height': 80,
		'color': '#414141',
		'norm': 'none',
	},
	'GC_content': {
		'height': 180,
		'color': '#004a0f',
		'norm': 'none',
	},

	'normalized_coverage': {
		'height': 120,
		'color': '#414141',
		'norm': 'log',
	},
	'std_coverage': {
		'height': 60,
		'color': '#414141',
		'norm': 'log',
	},
	'portion_covered': {
		'height': 60,
		'color': '#616161',
		'norm': 'none',
	},
	'abundance': {
		'height': 60,
		'color': '#818181',
		'norm': 'none',
	},

	'variability': {
		'height': 90,
		'color': '#4a000f',
		'norm': 'none',
	},
};

function getNamedLayerDefaults(layer, attribute, default_value)
{
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
