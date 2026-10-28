"""Workflow discovery helpers."""


def get_workflow_module_dict():
    """Return the mapping of workflow names to workflow classes."""
    from anvio.workflows.contigs import ContigsDBWorkflow
    from anvio.workflows.metagenomics import MetagenomicsWorkflow
    from anvio.workflows.pangenomics import PangenomicsWorkflow
    from anvio.workflows.phylogenomics import PhylogenomicsWorkflow
    from anvio.workflows.trnaseq import TRNASeqWorkflow
    from anvio.workflows.ecophylo import EcoPhyloWorkflow
    from anvio.workflows.sra_download import SRADownloadWorkflow

    workflows_dict = {'contigs': ContigsDBWorkflow,
                      'metagenomics': MetagenomicsWorkflow,
                      'pangenomics': PangenomicsWorkflow,
                      'phylogenomics': PhylogenomicsWorkflow,
                      'trnaseq': TRNASeqWorkflow,
                      'ecophylo': EcoPhyloWorkflow,
                      'sra_download': SRADownloadWorkflow}

    return workflows_dict
