import os

from anvio.errors import ConfigError

dir_path = os.path.dirname(os.path.abspath(__file__))
J = lambda x: os.path.join(dir_path, x)

miscellaneous_hmm_sources = {
        'POLYMERASES': {
            'DNA_Polymerase_Type_B': {'path': J('POLYMERASES/DNA_Polymerase_Type_B'),
                                      'info': 'A universal gene that is shared across all domains of life and viruses.'},
            'RNA_Polymerase_Type_A': {'path': J('POLYMERASES/RNA_Polymerase_Type_A'),
                                      'info': "Just like the DNA_Polymerase_Type_B, but easier to assemble from metagenomes compared to it (Meren's personal observation -- a grain of salt may be appropriate here)."},
            'RNA_Polymerase_Type_B': {'path': J('POLYMERASES/RNA_Polymerase_Type_B'),
                                      'info': 'Just like the RNA_Polymerase_Type_B, but different (often assemble poorly compared to RNA_Polymerase_Type_A -- very volatile insights here, so do your own research (but not like that)).'},
        },
        'rRNAs': {
            'Ribosomal_RNA_5S':  {'path': J('rRNAs/Ribosomal_RNA_5S'),
                                  'info': 'The rRNA 5S. Unfortunately it does not work. Feel free to fix it please :)'},
            'Ribosomal_RNA_12S': {'path': J('rRNAs/Ribosomal_RNA_12S'),
                                  'info': 'Just like the rRNA 5S. Bad model. Needs fixing.'},
            'Ribosomal_RNA_23S': {'path': J('rRNAs/Ribosomal_RNA_23S'),
                                  'info': 'rRNA 23S. The subunit that is commonly found in Bacteria and Archaea. It works, but does not compete well with its cousin, 16S rRNA model.'},
            'Ribosomal_RNA_28S': {'path': J('rRNAs/Ribosomal_RNA_28S'),
                                  'info': 'rRNA 28S. The subunit that is commonly found in Eukaryotic genomes. This model works, but does not compete well with its cousin, 18S rRNA model.'},
        },
}

def display_miscellaneous_models(run=None):
    if not run:
        from anvio.terminal import Run
        run = Run()

    misc_dir = os.path.abspath(os.path.dirname(__file__))

    run.warning(None, header='Miscellaneous HMM models', lc='green')
    run.info_single("These models are not run by default, but shipped with anvi'o due to their "
                    "relevance. If you wish to run any of these models on your contigs-db file, "
                    "You can simply copy-paste the full path shown below after the `-H` paramter.",
                    nl_after=1, level=0)

    for category in sorted(miscellaneous_hmm_sources.keys()):
        models = miscellaneous_hmm_sources[category]
        run.info_single(f"{category} ({len(models)} available)", mc='blue', nl_after=1)
        for model_name in sorted(models.keys()):
            model_data = models[model_name]
            model_path = os.path.abspath(model_data.get('path', ''))
            model_info = model_data.get('info', '')

            run.info_single(f"'{model_name}'. {model_info if model_info else ''} Please use either of these options "
                            f"to run this model on your contigs-db:", level=2)
            run.info_single(f"          anvi-run-hmms -c contigs-db -H {model_path}", mc='red', level=0, cut_after=0, nl_before=1, pretty_indentation=False)
            run.info_single(f"          anvi-run-hmms -c contigs-db -M {model_name}", mc='red', level=0, cut_after=0, nl_after=1, pretty_indentation=False)

        run.info_single('', nl_after=0, level=0)


def resolve_miscellaneous_model(key, run=None):
    """Resolve a user-provided miscellaneous model key to its directory path.

    Accepted formats: 'CATEGORY:MODEL_NAME' or 'MODEL_NAME' (if unique across categories)
    """

    if not run:
        from anvio.terminal import Run
        run = Run()

    if ':' in key:
        category, model = [p.strip() for p in key.split(':', 1)]
        if not category or not model:
            raise ConfigError(f"Malformed miscellaneous model key '{key}'. Use 'CATEGORY:MODEL_NAME'.")
        if category not in miscellaneous_hmm_sources:
            raise ConfigError(f"No miscellaneous category named '{category}'. Available categories: "
                              f"{', '.join(sorted(miscellaneous_hmm_sources.keys()))}.")
        if model not in miscellaneous_hmm_sources[category]:
            raise ConfigError(f"Category '{category}' has no model named '{model}'. Available models: "
                              f"{', '.join(sorted(miscellaneous_hmm_sources[category].keys()))}.")
        model_dict = miscellaneous_hmm_sources[category][model]
    else:
        matches = []
        for cat, models in miscellaneous_hmm_sources.items():
            if key in models:
                matches.append((cat, models[key]))
        if not matches:
            raise ConfigError(f"No miscellaneous model named '{key}' :( Use '--list-miscellaneous-models' to see options.")
        if len(matches) > 1:
            raise ConfigError(f"Model name '{key}' is ambiguous across categories ({', '.join([m[0] for m in matches])}) "
                              f"(congratulations). Please use 'CATEGORY:{key}' to be specific :)")
        category, model_dict = matches[0]
        key = f"{category}:{key}"

    model_path = os.path.abspath(model_dict.get('path', ''))
    if not os.path.isdir(model_path):
        raise ConfigError(f"The resolved miscellaneous model path does not exist: '{model_path}' :/")

    run.info('Found a miscellaneous model', key)
    run.info('In the HMM profile directory', model_path)

    return model_path
