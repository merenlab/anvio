
This program, like other anvi-setup commands, prepares your environment by downloading and configuring the ProstT5 model required for %(anvi-pan-genome --pan-mode structure-informed)s. ProstT5 is essential for Foldseek to perform efficient and accurate searches for protein structural similarities. You only need to run this setup once.

By executing this command, the necessary ProstT5 model weights will be downloaded and stored in the %(foldseek-model-data)s artifact, ensuring that Foldseek can function optimally in your anvi'o workflows.

Setting up the ProstT5 model is simple:

{{ codestart }}
anvi-setup-prostt5
{{ codestop }}

When running this program, you can provide a path to store your ProstT5 model in. The default path is `anvio/data/misc/PROSTT5/weights`; if you use a custom path, you will have to provide it to %(anvi-pan-genome)s with the same parameter. Here is an example run:


{{ codestart }}
anvi-setup-prostt5 --prostt5-weight-dir path/to/directory 
{{ codestop }}

If you want to overwrite any data that you have already downloaded (for example if you suspect something went wrong in the download), add the `--reset` flag: 

{{ codestart }}
anvi-setup-prostt5  --prostt5-weight-dir path/to/directory \ 
                        --reset
{{ codestop }}

