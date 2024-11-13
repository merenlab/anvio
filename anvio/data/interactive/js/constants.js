/**
 * Constants used from various interface functions
 *
 *  Authors: A. Murat Eren <a.murat.eren@gmail.com>
 *           Isaac Fink <iafink@uchicago.edu>
 *           Ozcan Esen
 *
 * Copyright 2015-2021, The anvi'o project (http://anvio.org)
 *
 * Anvi'o is a free software. You can redistribute this program
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with anvi'o. If not, see <http://opensource.org/licenses/GPL-3.0>.
 *
 * @license GPL-3.0+ <http://opensource.org/licenses/GPL-3.0>
 */

var COG_categories = {
    'A': '[A] RNA processing and modification',
    'B': '[B] Chromatin Structure and dynamics',
    'C': '[C] Energy production and conversion',
    'D': '[D] Cell cycle control and mitosis',
    'E': '[E] Amino Acid metabolism and transport',
    'F': '[F] Nucleotide metabolism and transport',
    'G': '[G] Carbohydrate metabolism and transport',
    'H': '[H] Coenzyme metabolis',
    'I': '[I] Lipid metabolism',
    'J': '[J] Translation',
    'K': '[K] Transcription',
    'L': '[L] Replication and repair',
    'M': '[M] Cell wall/membrane/envelop biogenesis',
    'N': '[N] Cell motility',
    'O': '[O] Post-translational modification, protein turnover, chaperone functions',
    'P': '[P] Inorganic ion transport and metabolism',
    'Q': '[Q] Secondary Structure',
    'T': '[T] Signal Transduction',
    'U': '[U] Intracellular trafficing and secretion',
    'V': '[V] Defense mechanisms',
    'W': '[W] Extracellular structures',
    'Y': '[Y] Nuclear structure',
    'Z': '[Z] Cytoskeleton',
    'R': '[R] General Functional Prediction only',
    'S': '[S] Function Unknown'
}

var default_COG_colors = {
  'A': '#98d67c',
  'B': '#446af3',
  'C': '#7e6500',
  'D': '#ff80ec',
  'E': '#68a300',
  'F': '#bc0066',
  'G': '#acd441',
  'H': '#c9a1ff',
  'I': '#eb9800',
  'J': '#0297f0',
  'K': '#e74d22',
  'L': '#01c0ab',
  'M': '#d1002e',
  'N': '#00a45c',
  'O': '#7b3979',
  'P': '#d6c777',
  'Q': '#ff9bc6',
  'T': '#455d00',
  'U': '#ff8241',
  'V': '#125e40',
  'W': '#c35b00',
  'Y': '#9ecea1',
  'Z': '#933215',
  'R': '#ff9a92',
  'S': '#81402e'
}

var KEGG_categories = {
  'C': 'Carbohydrate metabolism',
  'E': 'Energy metabolism',
  'L': 'Lipid metabolism',
  'N': 'Nucleotide metabolism',
  'A': 'Amino acid metabolism',
  'G': 'Glycan biosynthesis and metabolism',
  'V': 'Metabolism of cofactors and vitamins',
  'T': 'Metabolism of terpenoids and polyketides',
  'S': 'Biosynthesis of other secondary metabolites',
  'X': 'Xenobiotics biodegradation and metabolism'
}

var default_KEGG_colors = {
  'C': '#0000ee',
  'E': '#9933cc',
  'L': '#009999',
  'N': '#ff0000',
  'A': '#ff9933',
  'G': '#3399ff',
  'V': '#ff6699',
  'T': '#00cc33',
  'S': '#cc3366',
  'X': '#ccaa99'
}

var default_source_colors = {
  'None': '#808080',
  'Function': '#008000',
  'tRNA': '#226ab2',
  'rRNA': '#b22222'
}

// default colors for user-supplied functional annotations
var custom_cag_colors = Object.values(default_COG_colors).concat(Object.values(default_KEGG_colors));

var named_functional_sources = {
    'EGGNOG_BACT': {
        'accession_decorator': (function (d) {
                                    return '<a href="http://www.uniprot.org/uniprot/?query=' + d + '&sort=score" target="_blank">' + d + '</a>';
                                }),

    },

    'COG14_FUNCTION': {
        'accession_decorator': (function (d) {
                                    var cogs = d.split(', ').map((function (c){return '<a href="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=' + c +'" target=_"blank">' + c + '</a>';}));
                                    return cogs.join(', ');
                                }),
    },

    'COG14_CATEGORY': {
        'annotation_decorator': (function (d) {
                                    var cogs = d.split(', ').map((function (c){if (c in COG_categories) {return COG_categories[c];} else {return c;}}));
                                    return cogs.join('; ');
                                }),
    },

    'COG20_FUNCTION': {
        'accession_decorator': (function (d) {
                                    var cogs = d.split(', ').map((function (c){return '<a href="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=' + c +'" target=_"blank">' + c + '</a>';}));
                                    return cogs.join(', ');
                                }),
    },

    'COG20_CATEGORY': {
        'annotation_decorator': (function (d) {
                                    var cogs = d.split(', ').map((function (c){if (c in COG_categories) {return COG_categories[c];} else {return c;}}));
                                    return cogs.join('; ');
                                }),
    },


    'KOfam': {
        'accession_decorator': (function (d) {
                                    var kos = d.split(', ').map((function (c){return '<a href="https://www.genome.jp/dbget-bin/www_bget?ko:' + c +'" target=_"blank">' + c + '</a>';}));
                                    return kos.join(', ');
                                }),
    },

    'KEGG_Module': {
        'accession_decorator': (function (d) {
                                    var modules = d.split(', ').map((function (c){return '<a href="https://www.genome.jp/kegg-bin/show_module?' + c +'" target=_"blank">' + c + '</a>';}));
                                    return modules.join(', ');
                                }),
    },

    'Pfam': {
        'accession_decorator': (function (d) {
                                    var pfams = d.split(', ').map((function (c){return '<a href="https://pfam.xfam.org/family/' + c +'" target=_"blank">' + c + '</a>';}));
                                    return pfams.join(', ');
                                }),
    },

    'KEGG_PATHWAYS': {
        'annotation_decorator': (function (d) {
                                    var maps = d.split(', ').map((function (m){return '<a href="http://www.genome.jp/dbget-bin/www_bget?' + m +'" target=_"blank">' + m + '</a>';}));
                                    return maps.join(', ');
                                }),
    },

    'KEGG_Module': {
        'accession_decorator': (function (d) {
                                    var maps = d.split(', ').map((function (m){return '<a href="http://www.genome.jp/dbget-bin/www_bget?' + m +'" target=_"blank">' + m + '</a>';}));
                                    return maps.join(', ');
                                }),
    },

    'KEGG_BRITE': {
        'accession_decorator': (function (d) {
                                    var maps = d.split(', ').map((function (m){return '<a href="https://www.genome.jp/brite/' + m +'" target=_"blank">' + m + '</a>';}));
                                    return maps.join(', ');
                                }),
    },

    'GO_TERMS': {
        'annotation_decorator': (function (d) {
                                    var gos = d.split(', ').map((function (g){return '<a href="http://amigo.geneontology.org/amigo/term/' + g +'" target=_"blank">' + g + '</a>';}));
                                    return gos.join(', ');
                                }),
    },
}


function getPrettyFunctionsString(fstring, source) {
    if (!fstring)
        return ["-"];
    else
        farray = fstring.split('!!!');

        if ((source == null) || source == '' || source == 'None')
            return farray.join(' / ');
        else {
            var dec = new Array();

            for (f in farray)
                dec.push(decorateAccession(source, farray[f]))

            return dec.join(' / ')
        }
}


function decorateAccession(source, accession_id){
    if (!accession_id)
        return accession_id;

    if (source in named_functional_sources){
        if ('accession_decorator' in named_functional_sources[source]){
            return named_functional_sources[source]['accession_decorator'](accession_id);
        }
    }

    return accession_id;
}


function decorateAnnotation(source, annotation){
    if (source in named_functional_sources){
        if ('annotation_decorator' in named_functional_sources[source]){
            return named_functional_sources[source]['annotation_decorator'](annotation);
        }
    }

    return annotation;
}


var named_layers = {
    '__parent__': {
        'pretty_name': 'Parent',
    },
    'taxonomy': {
        'pretty_name': 'Taxonomy',
        'type': 'color',
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
    'num_gene_clusters': {
        'height': 220,
        'color': '#661111',
        'norm': 'none',
        'pretty_name': 'Num gene clusters',
        'type': 'bar',
    },
    'num_gene_clusters_raw': {
        'height': 220,
        'color': '#661111',
        'norm': 'none',
        'pretty_name': 'Num gene clusters (Raw)',
        'type': 'bar',
    },
    'singleton_gene_clusters': {
        'height': 220,
        'color': '#661111',
        'norm': 'none',
        'pretty_name': 'Singleton gene clusters',
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
    'num_genomes_gene_cluster_has_hits': {
        'height': 180,
        'color': '#002e4f',
        'norm': 'none',
        'pretty_name': 'Num contributing genomes',
    },
    'num_genes_in_gene_cluster': {
        'height': 180,
        'color': '#002e4f',
        'norm': 'sqrt',
        'pretty_name': 'Num genes in GC',
    },
    'num_genomes_pc_has_hits': {
        'height': 180,
        'color': '#002e4f',
        'norm': 'none',
        'pretty_name': 'Num contributing genomes',
    },
    'num_genes_in_pc': {
        'height': 180,
        'color': '#002e4f',
        'norm': 'sqrt',
        'pretty_name': 'Num genes in GC',
    },
    'functional_homogeneity_index': {
        'height': 180,
        'color': '#3D9970',
        'norm': 'none',
        'min': 0,
        'max': 1,
        'type': 'bar',
        'pretty_name': 'Func. Homogeneity Ind.',
    },
    'geometric_homogeneity_index': {
        'height': 180,
        'color': '#3D8870',
        'norm': 'none',
        'min': 0,
        'max': 1,
        'type': 'bar',
        'pretty_name': 'Geo. Homogeneity Ind.',
    },
    'combined_homogeneity_index': {
        'height': 180,
        'color': '#326b59',
        'norm': 'none',
        'min': 0,
        'max': 1,
        'type': 'bar',
        'margin': 90,
        'pretty_name': 'Comb. Homogeneity Ind.',
    },
    'AAI_min': {
        'height': 180,
        'color': '#9e009e',
        'norm': 'none',
        'min': 0,
        'max': 1,
        'type': 'bar',
        'pretty_name': 'Min AAI',
    },
    'AAI_max': {
        'height': 180,
        'color': '#9e009e',
        'norm': 'none',
        'min': 0,
        'max': 1,
        'type': 'bar',
        'pretty_name': 'Max AAI',
    },
    'AAI_avg': {
        'height': 180,
        'color': '#590059',
        'norm': 'none',
        'min': 0,
        'max': 1,
        'type': 'bar',
        'pretty_name': 'Avg AAI',
    },
    'Gene_cluster_type_LLR': {
        'height': 180,
        'color': '#00AA00',
        'color-start': '#AA0000',
        'norm': 'none',
        'type': 'intensity',
        'pretty_name': 'Gene Cluster Type (LLR)',
    },
    'Gene_cluster_type': {
        'height': 180,
        'pretty_name': 'Gene Cluster Type (Class)',
        'margin': 30,
    },
    'max_num_paralogs': {
        'height': 180,
        'color': '#002e4f',
        'norm': 'none',
        'min': 1,
        'margin': 90,
        'pretty_name': 'Max num paralogs',
    },
    'SCG': {
        'height': 180,
        'color': '#4f1111',
        'norm': 'none',
        'margin': 90,
        'pretty_name': 'SCG Clusters',
    },
    'total_length': {
        'height': 180,
        'color': '#414141',
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
    'detection': {
        'height': 0,
        'color': '#616161',
        'norm': 'none',
        'pretty_name': 'Detection',
    },
    'ECGs_and_EAGs!EAG;ECG;NA': {
        'height': 400,
        'pretty_name': 'ECGs and EAGs',
    },
    'EAG_ECG_ratio': {
        'height': 400,
        'pretty_name': 'EAG Ratio',
        'color': '#940000',
        'norm': 'none',
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
        'color': '#000077',
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
    'matching_domain': {
        'type': 'text',
        'pretty_name': 'Matching domain',
    },
    'blank_view': {
        'height': 0,
        'color': '#FFFFFF',
        'norm': 'none',
        'pretty_name': '_',
    },
    'num_gene_clusters_in_psgc': {
        'height': 180,
        'color': '#002e4f',
        'norm': 'none',
        'pretty_name': 'Num gene clusters in PSGC',
    },
    'num_genes_in_psgc': {
        'height': 180,
        'color': '#002e4f',
        'norm': 'none',
        'pretty_name': 'Num genes in PSGC',
    },
};

named_category_colors = {
    'KNOWN'           : '#00AA00',
    'UNKNOWN'         : '#F0F0F0',
    'ECG'             : '#00AA00',
    'EAG'             : '#AA0000',
    'NA'              : '#F0F0F0',
    'TSC'             : '#e38181',
    'TSA'             : '#0000AA',
    'TNC'             : '#00AA00',
    'TNA'             : '#00d1ca',
    'NaN'             : '#73727a',
    'K'               : '#233B43',
    'KWP'             : '#556C74',
    'GU'              : '#65ADC2',
    'EU'              : '#E84646',
    'SINGL'           : '#BCC8CC',
    'DISC'            : '#BCC8CC',
    'CORE'            : '#00AA00',
    'ACCESSORY'       : '#AA0000',
};

function getNamedCategoryColor(name)
{
    if ((name == null) || name == '' || name == 'None')
        return '#FFFFFF';
    else {
        if (name in named_category_colors)
            return named_category_colors[name];
        else
            return randomColor();
    }
}

pretty_names = {
    'tnf-cov': 'Seq. Composition + Diff. Coverage',
    'cov': 'Differential coverage',
    'tnf': 'Sequence composition',
    'tnf-splits': 'Sequence composition (w/independent splits)'
};

function getPrettyLayerTitle(layer_title) {
    if (layer_title in named_layers && 'pretty_name' in named_layers[layer_title]) {
        layer_title = named_layers[layer_title]['pretty_name'];
    } else if(layer_title.substring(0, 5) == "hmmx_") {
        layer_title = layer_title.replace(/hmmx_/g, "").replace(/_/g, " ");
    } else if(layer_title.substring(0, 5) == "hmms_") {
        layer_title = layer_title.replace(/hmms_/g, "").replace(/_/g, " ");
    } else {
        layer_title = layer_title.replace(/_/g, " ");
    }

    if (layer_title.indexOf('!') > -1 )
    {
        layer_title = layer_title.split('!')[0];
    }

    return layer_title;
}

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

function getClusteringPrettyName(name)
{
    var name_parts = name.split(':').map(getPrettyName);

    // build -> $clustering (D: $distance; L: $linkage)
    return name_parts[0] + ' (D: ' + name_parts[1] + '; L: ' + name_parts[2] + ')';
}

function getNamedLayerDefaults(layer, attribute, default_value, group)
{
    if (typeof default_value == "string" && default_value.charAt(0) != '#'){
        default_value = getPrettyName(default_value)
    }

    if (typeof group !== 'undefined' && group.startsWith('ANI_')) {
        if (group === 'ANI_percentage_identity' || group === 'ANI_full_percentage_identity' || group === 'ANI_ani') {
            if (attribute == 'min') return 0.7;
            if (attribute == 'max') return 1;
        }

        if (attribute == 'height') return '180';
        if (attribute == 'color')  return '#FF0000';
        if (attribute == 'color-start')  return '#F2F2F2';
        if (attribute == 'type')   return 'intensity';
    }

    if (typeof group !== 'undefined' && group.startsWith('SourMash_')) {
        if (attribute == 'height') return '180';
        if (attribute == 'color')  return '#FF006F';
        if (attribute == 'color-start')  return '#F2F2F2';
        if (attribute == 'type')   return 'intensity';
    }

    /* Some ad-hoc manipulation of special hmmx_ split hmm layers */
    if (layer.substring(0, 5) == "hmmx_") {
        if (attribute == 'height') return '30';
        if (attribute == 'norm')   return 'none';
        if (attribute == 'color')  return '#882222';
    }

    /* Some ad-hoc manipulation of special hmms_ single hmm layers */
    if (layer.substring(0, 5) == "hmms_"){
        if (attribute == 'type') return 'intensity';
        if (attribute == 'height') return '150';
        if (attribute == 'norm')   return 'sqrt';

        if (layer.substring(0, 13) == "hmms_Transfer"){
            console.log(layer, attribute);
            if (attribute == 'color-start')  return '#bfd9f3';
            if (attribute == 'color')  return '#226ab2';
        }
        else if (layer.substring(0, 14) == "hmms_Ribosomal"){
            if (attribute == 'color-start')  return '#FFDDDD';
            if (attribute == 'color')  return '#882222';
        }
        else {
            if (attribute == 'color')  return '#444444';
            if (attribute == 'color-start')  return '#DDDDDD';
        }
    }

    if (layer.substring(0, 6) == "motif_") {
        if (attribute == 'norm')   return 'none';
        if (attribute == 'color')  return '#222288';
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
