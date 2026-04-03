This artifact stores the data downloaded by %(anvi-setup-prostt5)s and is essential for running %(anvi-pan-genome)s. It includes the ProstT5 model weights, which are required by Foldseek to efficiently perform searches for protein structural similarities.

As detailed in the Foldseek documentation, this data consists of the pre-trained ProstT5 model that accelerates protein structure search tasks. The ProstT5 model is crucial for ensuring accurate results when using Foldseek in your anvi'o workflows.

By default, the ProstT5 model weights are stored in anvio/data/misc/PROSTT5/weights, but users can specify a custom path during setup by using the --prostt5-weight-dir parameter in %(anvi-setup-foldseek)s.