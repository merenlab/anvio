output for parser 'myrast_cmdline_dont_use' was generated using:

    svr_assign_to_dna_using_figfams < ../contigs.fa > svr_assign_to_dna_using_figfams.txt

output for parser 'myrast_cmdline' was generated using:

    svr_call_pegs < ../contigs.fa > svr_call_pegs.txt
    svr_assign_using_figfams -otu -all -blastOutput 1 -reliability 20 < svr_call_pegs.txt > svr_assign_using_figfams.txt
