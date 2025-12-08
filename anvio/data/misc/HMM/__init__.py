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

    misc_sources = miscellaneous_hmm_sources
    misc_dir = os.path.abspath(os.path.dirname(__file__))

    run.warning(None, header='Miscellaneous HMM models', lc='green')
    run.info_single("These models are not run by default, but shipped with anvi'o due to their "
                    "relevance. If you wish to run any of these models on your contigs-db file, "
                    "You can simply copy-paste the full path shown below after the `-H` paramter.",
                    nl_after=1, level=0)

    for category in sorted(misc_sources.keys()):
        models = misc_sources[category]
        run.info_single(f"{category} ({len(models)} available)", mc='blue', nl_after=1)
        for model_name in sorted(models.keys()):
            model_data = models[model_name]
            model_path = os.path.abspath(model_data.get('path', ''))
            model_info = model_data.get('info', '')

            run.info_single(f"{model_name}. {model_info if model_info else ''} Please use the following "
                            f"full path for this model:", level=2)
            run.info_single(f"          {model_path}", mc='red', level=0, cut_after=0, nl_before=1, nl_after=1, pretty_indentation=False)

        run.info_single('', nl_after=0, level=0)
