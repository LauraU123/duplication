#this workflow reconstructs the mutations in the RSV A and RSV B duplicated regions of the G gene
# --reference-name {params.ref}
#       """
configfile: "config/configfile.yaml"
A_OR_B = ["a"]


build_dir = 'results'

rule all:
    input:
        expand("results/{a_or_b}/{build}/graphs/cumulative_sum_syn.png", a_or_b=A_OR_B, build = config.get("buildstorun", ['genome']))


'''
This part of the workflow expects input files
            sequences = "data/sequences.fasta"
            metadata = "data/metadata.tsv"
'''

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta"
    output:
        sequence_index = build_dir + "/{a_or_b}/{build_name}/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule newreference:
    message:
        """
        Making new reference
        """
    input:
        oldreference = "config/{a_or_b}reference.gbk"
    output:
        newreferencegbk = build_dir + "/{a_or_b}/{build_name}/newreference.gbk",
        newreferencefasta = build_dir + "/{a_or_b}/{build_name}/newreference.fasta",
        greference = build_dir + "/{a_or_b}/{build_name}/greference.fasta"
    params:
        gene = lambda w: w.build_name,
        newreference = build_dir + "/{a_or_b}/{build_name}/newreference",
        oldreference = 'config/{a_or_b}reference',
        greference = build_dir + "/{a_or_b}/{build_name}/greference"
    shell:
        """
        python scripts/newreference.py \
            --greference {params.greference} \
            --reference {params.oldreference} \
            --output {params.newreference} \
            --gene {params.gene}
        """


rule filter:
    message:
        """
        Aligning sequences to {input.reference}
            - gaps relative to reference are considered real
        """
    input:
        sequences = "data/{a_or_b}/sequences.fasta",
        reference = "config/{a_or_b}reference.gbk",
        metadata = "data/{a_or_b}/metadata.tsv",
        sequence_index = rules.index_sequences.output
    output:
    	sequences = build_dir + "/{a_or_b}/{build_name}/filtered.fasta"
    params:
    	group_by = config["filter"]["group_by"],
    	min_length = lambda w: config["filter"]["min_length"].get(w.build_name, 10000),
    	subsample_max_sequences = config["filter"]["subsample_max_sequences"],
    	strains = "config/dropped_strains.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --exclude {params.strains} \
            --subsample-max-sequences {params.subsample_max_sequences} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = "config/{a_or_b}reference.fasta"
    output:
        alignment = build_dir + "/{a_or_b}/{build_name}/sequences.aligned.fasta",
        insertionsfile = build_dir + "/{a_or_b}/{build_name}/insertions.csv"
    threads: 4
    shell:
        """
        nextalign run -j {threads}\
            --reference {input.reference} \
            --output-fasta {output.alignment} \
            --output-insertions {output.insertionsfile} \
            {input.sequences}
        """

rule cut_G_out:
    input:
        oldalignment = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        newalignment = build_dir + "/{a_or_b}/{build_name}/aligned_without_G.fasta"
    shell:
        """
        python scripts/cut_G_out.py \
            --oldalignment {input.oldalignment} \
            --newalignment {output.newalignment} \
            --reference {input.reference} \
        """
        
rule cut:
    input:
        oldalignment = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        newalignment = build_dir + "/{a_or_b}/{build_name}/newalignment.fasta"
    shell:
        """
        python scripts/cut.py \
            --oldalignment {input.oldalignment} \
            --newalignment {output.newalignment} \
            --reference {input.reference} \
        """

rule realign:
    input:
        newalignment = rules.cut.output.newalignment,
        reference = rules.newreference.output.greference
    output:
        realigned = build_dir + "/{a_or_b}/{build_name}/realigned.fasta"
    threads: 4
    shell:
        """
        augur align --nthreads {threads} \
            --sequences {input.newalignment} \
            --reference-sequence {input.reference} \
            --output {output.realigned}
        """

rule alignment_for_tree:
    input:
        realigned = rules.realign.output.realigned,
        original = rules.align.output.alignment,
        reference = "config/{a_or_b}reference.gbk"
    output:
        aligned_for_tree = build_dir + "/{a_or_b}/{build_name}/alignment_for_tree.fasta"
    params:
        gene = lambda w: w.build_name
    shell:
        """
        python scripts/align_for_tree.py \
            --realign {input.realigned} \
            --original {input.original} \
            --reference {input.reference} \
            --output {output.aligned_for_tree} \
            --gene {params.gene}
        """

rule tree:
    message: "Building tree with max likelihood method using GTR model"
    input:
        alignment = rules.cut_G_out.output.newalignment
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree_raw.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.tree.input.alignment,
        metadata = rules.filter.input.metadata
    output:
        tree = build_dir + "/{a_or_b}/{build_name}/tree.nwk",
        node_data = build_dir + "/{a_or_b}/{build_name}/branch_lengths.json"
    params:
    	coalescent = config["refine"]["coalescent"],
    	clock_filter_iqd = config["refine"]["clock_filter_iqd"],
    	date_inference = config["refine"]["date_inference"]
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --date-confidence \
            --timetree \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message:
        """
        Reconstructing ancestral sequences and mutations
          - inferring ambiguous mutations
        """
    input:
        tree = rules.refine.output.tree,
        alignment = rules.alignment_for_tree.output.aligned_for_tree
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/nt_muts.json"
    params:
    	inference = config["ancestral"]["inference"]
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = rules.newreference.output.newreferencegbk
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/aa_muts.json",
        aa_data = build_dir + "/{a_or_b}/{build_name}/alignedG.fasta"
    params:
    	alignment_file_mask = build_dir + "/{a_or_b}/{build_name}/aligned%GENE.fasta"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
            --alignment-output {params.alignment_file_mask}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.filter.input.metadata
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/traits.json"
    log:
        "logs/{a_or_b}/traits_{build_name}_rsv.txt"
    params:
    	columns = config["traits"]["columns"]
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

def get_node_data(w):
    node_data = [rules.refine.output.node_data,
                    rules.traits.output.node_data,
                    rules.ancestral.output.node_data,
                    rules.translate.output.node_data]
    return node_data

rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.traits.input.metadata,
        node_data = get_node_data,
        auspice_config = config["files"]["auspice_config"],
        description = config["description"]
    output:
        auspice_json =  build_dir + "/{a_or_b}/{build_name}/tree.json",
        root_sequence = build_dir + "/{a_or_b}/{build_name}/tree_root-sequence.json"
    params:
    	title = lambda w: f"RSV-{w.a_or_b.upper()} phylogeny"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --title {params.title:q} \
            --description {input.description} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --output {output.auspice_json}
        """

rule branch_from_root:
    message:
        """Adding mutations to root sequence to reconstruct all branches"""
    input:
        root = rules.export.output.root_sequence,
        sequences = rules.index_sequences.input,
        tree = rules.refine.output.tree,
        treejson = rules.export.output.auspice_json
    output:
        reconstructed_seq = "results/{a_or_b}/{build_name}/reconstructed_sequences.fasta"
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
        aligned = "results/{a_or_b}/{build_name}/pairwise_G.fasta"
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
        only_dupl = "results/{a_or_b}/{build_name}/only_duplication.fasta"
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
        aligned_dupl = "results/{a_or_b}/{build_name}/aligned_duplication.fasta"
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
        graph_cumsum_syn = "results/{a_or_b}/{build_name}/graphs/cumulative_sum_syn.png"
    params:
        length = lambda w: config["graphs"].get(w.a_or_b),
        tsv = "results/{a_or_b}/{build_name}/graphs/cumulative_sum_"
    shell:
        """
        python3 scripts/reconstruct_and_graph.py \
        --input {input.data} \
        --tree {input.tree_} \
        --length {params.length} \
        --output {output.graph_cumsum_syn} \
        --tsv {params.tsv} \
        """