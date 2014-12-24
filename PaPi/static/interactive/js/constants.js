var named_layers = {
    'length': {
        'height': 100,
        'color': '#808080',
        'norm': 'none',
    },
    'GC_content': {
        'height': 100,
        'color': '#115511',
        'norm': 'none',
    },
    'num_genes': {
        'height': 80,
        'color': '#551111',
        'norm': 'none',
    },
    'num_tax_calls': {
        'height': 80,
        'color': '#552222',
        'norm': 'none',
    },
    'num_function_calls': {
        'height': 80,
        'color': '#553333',
        'norm': 'none',
    },
    'tax_accuracy': {
        'height': 80,
        'color': '#554444',
        'norm': 'none',
    }
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