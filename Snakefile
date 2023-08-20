#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene
# --reference-name {params.ref}
#       """
configfile: "config/configfile.yaml"
A_OR_B = ["b"]

rule all:
    input:
        expand("results/{a_or_b}/graphs/cumulative_sum_syn.png", a_or_b=A_OR_B)

rule branch_from_root:
    message:
        """Adding mutations to root sequence to reconstruct all branches"""
    input:
        root = "data/rsv_{a_or_b}_root_sequence.json",
        sequences = "data/{a_or_b}_sequences.fasta",
        tree = "data/{a_or_b}_tree.nwk",
        treejson = "data/rsv_{a_or_b}_genome.json"
    output:
        reconstructed_seq = "results/{a_or_b}/reconstructed_sequences.fasta"
    shell:
        """
        python3 scripts/reconstruct_from_root.py \
        --input-root {input.root} \
        --input-tree-json {input.treejson} \
        --sequences {input.sequences} \
        --input-tree {input.tree} \
        --output {output.reconstructed_seq} 
        """

rule align_to_ref:
    message:
        """Aligning reconstructed branches and tips pairwise"""
    input:
        reconstructed_seq = rules.branch_from_root.output.reconstructed_seq,
        reference = "config/{a_or_b}reference.fasta"
    output:
        aligned = "results/{a_or_b}/pairwise_G.fasta"
    shell:
        """
        nextalign run -j 4 \
        --reference {input.reference} \
        --output-fasta {output.aligned} \
          {input.reconstructed_seq}
        """

rule just_duplication:
    message:
        """Extracting duplication from alignment based on given parameters"""
    input:
        reconstructed_seq = rules.align_to_ref.output.aligned
    params:
        start = lambda w: config["just_dupl"]["start"].get(w.a_or_b),
        end = lambda w: config["just_dupl"]["end"].get(w.a_or_b)
    output:
        only_dupl = "results/{a_or_b}/only_duplication.fasta"
    shell:
        """
        python3 scripts/duplication.py \
        --input {input.reconstructed_seq} \
        --start {params.start} \
        --end {params.end} \
        --output {output.only_dupl}
        """

rule align_G:
    message:
        """Aligning duplications - MSA"""
    input:
        only_dupl = rules.just_duplication.output.only_dupl
    output:
        aligned_dupl = "results/{a_or_b}/aligned_duplication.fasta"
    params:
        ref = lambda w: config["ref"].get(w.a_or_b)
    shell:
        """
        augur align --nthreads 4 \
        --sequences {input.only_dupl} \
        --output {output.aligned_dupl} \
        --reference-name {params.ref}
        """

rule graphs:
    input:
        data = rules.align_G.output.aligned_dupl,
        tree_ = rules.branch_from_root.input.tree
    output:
        graph_cumsum_syn = "results/{a_or_b}/graphs/cumulative_sum_syn.png"
    params:
        length = lambda w: config["graphs"].get(w.a_or_b),
        tsv = "results/{a_or_b}/graphs/cumulative_sum_"
    shell:
        """
        python3 scripts/reconstruct_and_graph.py \
        --input {input.data} \
        --tree {input.tree_} \
        --length {params.length} \
        --output {output.graph_cumsum_syn} \
        --tsv {params.tsv} \
        """