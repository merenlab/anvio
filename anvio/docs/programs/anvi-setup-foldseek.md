
This program, like other anvi-setup commands, prepares your environment by downloading and configuring the Foldseek ProstT5 model required for %(anvi-display-pan --structure)s. ProstT5 is essential for Foldseek to perform efficient and accurate searches for protein structural similarities. You only need to run this setup once.

By executing this command, the necessary ProstT5 model weights will be downloaded and stored in the %(foldseek-model-data)s artifact, ensuring that Foldseek can function optimally in your anvi'o workflows.

Setting up the Foldseek ProstT5 model is simple:

It's easy as 1-2-3:

{{ codestart }}
anvi-setup-foldseek
{{ codestop }}

When running this program, you can provide a path to store your Foldseek ProstT5 model in. The default path is `anvio/data/misc/PROSTT5/weights`; if you use a custom path, you will have to provide it to %(anvi-run-interacdome)s with the same parameter. Here is an example run:


{{ codestart }}
anvi-setup-foldseek --foldseek-weight-dir path/to/directory 
{{ codestop }}

If you want to overwrite any data that you have already downloaded (for example if you suspect something went wrong in the download), add the `--reset` flag: 

{{ codestart }}
anvi-setup-foldseek  --foldseek-weight-dir path/to/directory \ 
                        --reset
{{ codestop }}

