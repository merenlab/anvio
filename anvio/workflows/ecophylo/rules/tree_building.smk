if M.run_iqtree == True:

    rule iqtree:
        """Calculate a phylogenetic tree using iqtree"""
        input:
            fasta=rules.remove_sequences_with_X_percent_gaps.output.fasta,
        output:
            tree=os.path.join(dirs_dict["TREES"], "{group}", "{group}.iqtree"),
            done=touch(os.path.join(dirs_dict["TREES"], "{group}", "{group}-tree.done")),
        log:
            rule_log("iqtree", "iqtree_{group}"),
        threads: M.T("iqtree")
        params:
            outfile=os.path.join(dirs_dict["TREES"], "{group}"),
            model=M.get_param_value_from_config(["iqtree", "-m"]),
            additional_params=M.get_param_value_from_config(
                ["iqtree", "additional_params"]
            ),
        shell:
            "iqtree -s {input} -nt AUTO -m {params.model} -pre {params.outfile}  -T AUTO {params.additional_params} >> {log} 2>&1"

elif M.run_fasttree == True:

    rule fasttree:
        """Want to go faster?? Then Fasttree"""
        input:
            fasta=rules.remove_sequences_with_X_percent_gaps.output.fasta,
        output:
            tree=os.path.join(dirs_dict["TREES"], "{group}", "{group}.nwk"),
            done=touch(os.path.join(dirs_dict["TREES"], "{group}", "{group}-tree.done")),
        log:
            rule_log("fasttree", "fasttree_{group}"),
        threads: M.T("fasttree")
        params:
            additional_params=M.get_param_value_from_config(
                ["fasttree", "additional_params"]
            ),
            tree=os.path.join(dirs_dict["TREES"], "{group}", "{group}.nwk"),
        shell:
            "FastTree -fastest {input} 2> {log} 1> {params.tree}"


rule rename_tree_tips:
    """Add "_split_00001" suffix to tree tips names to bind with profileDB"""
    input:
        tree=os.path.join(dirs_dict["TREES"], "{group}", "{group}-tree.done"),
    output:
        done=os.path.join(dirs_dict["TREES"], "{group}", "{group}_combined.done"),
        tree=os.path.join(dirs_dict["TREES"], "{group}", "{group}_renamed.nwk"),
        fasta=os.path.join(dirs_dict["TREES"], "{group}", "{group}_renamed.faa"),
        fasta_all=os.path.join(dirs_dict["TREES"], "{group}", "{group}_renamed_all.faa"),
    log:
        rule_log("rename_tree_tips", "rename_tree_tips_{group}"),
    threads: M.T("rename_tree_tips")
    params:
        fasttree=os.path.join(dirs_dict["TREES"], "{group}", "{group}.nwk"),
        iqtree=os.path.join(dirs_dict["TREES"], "{group}", "{group}.iqtree"),
        fasta=os.path.join(
            dirs_dict["MSA"], "{group}", "{group}_aligned_trimmed_filtered.fa"
        ),
        fasta_all=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"], "{group}", "{group}-all.faa"
        ),
    run:
        # rename tree tips
        from ete3 import Tree

        def add_split_name_to_tree_tips(tree, outfile):
            """This function adds a string to the end of tree tip names in a newick file

            Parameters
            ==========
            tree: newick file

            outfile: str
                path to returned newick tree
            """

            t = Tree(tree)

        for leaf in t:
            leaf.name = leaf.name + "_split_00001"
        t.write(format=1, outfile=outfile)
        if M.run_iqtree == True:
            add_split_name_to_tree_tips(tree=params.iqtree, outfile=output.tree)
        elif M.run_fasttree == True:
            add_split_name_to_tree_tips(tree=params.fasttree, outfile=output.tree)

        # rename fasta headers to match tree
        def add_split_name_to_fasta_headers(fasta, outfile):
            """This function adds a string "_split_00001\n" to the end of fasta file headers

            Parameters
            ==========
            fasta: str
                path to fasta file

            outfile: str
                path returned fasta
            """

            fasta = open(str(fasta))

        newfasta = open(outfile, "a")
        for line in fasta:
            if line.startswith(">"):
                newname = line.rstrip("\n") + "_split_00001\n"
                newfasta.write(newname)
            else:
                newfasta.write(line)
        fasta.close()
        newfasta.close()
        add_split_name_to_fasta_headers(params.fasta, output.fasta)
        add_split_name_to_fasta_headers(params.fasta_all, output.fasta_all)
        shell("touch {output.done}")
