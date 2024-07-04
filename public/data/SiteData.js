const {Command, Parameter} = require('./Structure');

// ---- DNA ---- //

const DNACategories = [
    { title: "File Conversion", entries: ['fasta-to-fastq', 'fast5-to-fastq', 'bcl2fastq', 'sam-to-bam']},
    { title: "File Processing", entries: ['sort-bam', 'index-bam', 'fastqc', 'trimming', 'alignment', 'mark-or-remove-duplicates']},
    { title: "Statistics", entries: ['add-or-replace-read-groups', 'bam-index-stats', 'flag-stats']},
    { title: "Summary and Graphs", entries: ['alignment-summary', 'gc-bias-summary', 'insert-size-summary']},
    // { title: "Metrics", entries: ['', '', '']}, // CBTT
    { title: "Sequencing", entries: ['create-sequence-dictionary', 'sequence-depth', 'sequence-coverage']},
];

let DNACommands = new Map([
    // File Conversion
    ["fasta-to-fastq", new Command("Convert FastA to FastQ", false, 'fasta-to-fastq')],
    ["fast5-to-fastq", new Command("Convert Fast5 to FastQ", false, 'fast5-to-fastq')],
    ["bcl2fastq", new Command("Convert BCL to FastQ", false, 'bcl2fastq')],
    ["sam-to-bam", new Command("Convert SAM to BAM", false, 'sam-to-bam')],

    // File Processing
    ["sort-bam", new Command("Sort BAM", false, 'sort-bam')],
    ["index-bam", new Command("Index BAM", false, 'index-bam')],
    ["fastqc", new Command("Run FastQC", false, 'fastqc')],
    ["trimming", new Command("Trimming", false, 'trimming')],
    ["alignment", new Command("Alignment", false, 'alignment')],
    ["mark-or-remove-duplicates", new Command("Mark Or Remove Duplicates", false, 'mark-or-remove-duplicates')],

    // Statistics
    ["add-or-replace-read-groups", new Command("Add or Replace Read Groups", false, 'add-or-replace-read-groups')],
    ["bam-index-stats", new Command("Extract BAM Index Stats", false, 'bam-index-stats')],
    ["flag-stats", new Command("Examine Flag Stats", false, 'flag-stats')],

    // Summary and Graphs
    ["alignment-summary", new Command("Generate Alignment Summary", false, 'alignment-summary')],
    ["gc-bias-summary", new Command("Generate GC Bias Summary", false, 'gc-bias-summary')],
    ["insert-size-summary", new Command("Generate Insert Size Summary", false, 'insert-size-summary')],

    // Sequencing
    ["create-sequence-dictionary", new Command("Create Sequence Dictionary", false, 'create-sequence-dictionary')],
    ["sequence-depth", new Command("Examine Sequencing Depth", false, 'sequence-depth')],
    ["sequence-coverage", new Command("Examine Sequencing Coverage", false, 'sequence-coverage')],
]);

let DNAParameters = new Map([
    // File Conversion, no parameters
    ["fasta-to-fastq", []],
    ["fast5-to-fastq", []],
    ["bcl2fastq", []],
    ["sam-to-bam", []],

    // File Processing
    ["sort-bam", []], // none
    ["index-bam", []], // none
    ["fastqc", []], // none
    ["trimming", [
        new Parameter("Perform Adapter Trim", false, 'checkbox', [], 'adapter_trim', false),
            new Parameter("Adapter Trim: Adapter File", true, 'upload', [], 'adapter_trim_file', ''),
        new Parameter("Perform Read Length Trim", false, 'checkbox', [], 'read_length_trim', false),
            new Parameter("Read Length Trim: Minimum Length", true, 'number', [], 'minlen', ''),
        new Parameter("Perform Sliding Window Trim", false, 'checkbox', [], 'sliding_window_trim', false),
            new Parameter("Sliding Window Trim: Start", true,  'number', [], 'window_start', ''),
            new Parameter("Sliding Window Trim: End", true,  'number', [], 'window_end', ''),
        new Parameter("Perform Leading Trim", false, 'checkbox', [], 'leading_trim', false),
            new Parameter("Leading Trim: Input", true, 'number', [], 'leading', ''),
        new Parameter("Perform Trailing Trim", false, 'checkbox', [], 'trailing_trim', false),
            new Parameter("Trailing Trim: Input", true, 'checkbox', [], 'trailing', false),
    ]],
    ["alignment", [
        new Parameter("Alignment Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
        new Parameter("Alignment Mode", false, 'select', ['BWA', 'Bowtie', 'Bowtie2'], 'mode', 'BWA'),
    ]],
    ["mark-or-remove-duplicates", [
        new Parameter("Remove Duplicates", false, 'select', ['Do Not Remove', 'Remove All', 'Remove Sequencing Duplicates Only'], 'remove_dupes', 'Do Not Remove'),
    ]],

    // Statistics
    ["add-or-replace-read-groups", [
        new Parameter("Overwrite Existing Read Groups", false, 'checkbox', [], 'editMode', false),
        new Parameter("New Read Groups Line", false, 'text', [], 'newReadGroupLine', ''),
    ]],
    ["bam-index-stats", []], // none
    ["flag-stats", []], // none

    // Summary and Graphs
    ["alignment-summary", [
        new Parameter("Alignment Summary Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
        new Parameter("Create Alignment Graph", false, 'checkbox', [], 'isVisual', false),
    ]], 
    ["gc-bias-summary", [
        new Parameter("GC Bias Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
        new Parameter("Create Bias Graph", false, 'checkbox', [], 'isVisual', false),
    ]],
    ["insert-size-summary", [
        new Parameter("Create Size Graph", false, 'checkbox', [], 'isVisual', false),
    ]],

    // Sequencing
    ["create-sequence-dictionary", [
        new Parameter("Sequence Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
    ]],
    ["sequence-depth", []], // none
    ["sequence-coverage", [
        new Parameter("Create Sequence Graph", false, 'checkbox', [], 'isVisual', false),
    ]]
]);

// ---- RNA ---- //

const RNACategories = [
    { title: "File Conversion", entries: ['fasta-to-fastq', 'fast5-to-fastq', 'bcl2fastq', 'sam-to-bam']},
    { title: "Pre Processing", entries: ['sort-bam', 'index-bam', 'trimming', 'alignment-mapping']},
    { title: "Quality Control", entries: ['fastqc', 'fastq-screen', 'alignmentqc', 'rseqc']},
    { title: "Analysis", entries: ['quantification', 'normal-diff', 'quality-yield-metrics', 'rseq-metrics', 'flagstats-rna', 'coverage-rna']},
];

let RNACommands = new Map([
    // File Conversion
    ["fasta-to-fastq", new Command("Convert FastA to FastQ", false, 'fasta-to-fastq')],
    ["fast5-to-fastq", new Command("Convert Fast5 to FastQ", false, 'fast5-to-fastq')],
    ["bcl2fastq", new Command("Convert BCL to FastQ", false, 'bcl2fastq')],
    ["sam-to-bam", new Command("Convert SAM to BAM", false, 'sam-to-bam')],

    // PreProcessing
    ["sort-bam", new Command("Sort BAM", false, 'sort-bam')],
    ["index-bam", new Command("Index BAM", false, 'index-bam')],
    ["trimming", new Command("Trimming", false, 'trimming')],
    ["alignment-mapping", new Command("Aligment / Mapping", false, 'alignment-mapping')],

    // Quality Control
    ["fastqc", new Command("FastQC", false, 'fastqc')],
    ["fastq-screen", new Command("FastQScreen", false, 'fastq-screen')],
    ["alignmentqc", new Command("Alignment QC", false, 'alignmentqc')],
    ["rseqc", new Command("RSeQC", false, 'rseqc')],

    // Analysis
    ["quantification", new Command("Quantification", false, 'quantification')],
    ["normal-diff", new Command("Normalization and Differential Expression", false, 'normal-diff')],
    ["quality-yield-metrics", new Command("Collect Quality Yield Metrics", false, 'quality-yield-metrics')],
    ["rseq-metrics", new Command("Collect RNA RSeq Metrics", false, 'rseq-metrics')],
    ["flagstats-rna", new Command("Flagstat", false, 'flagstats-rna')],
    ["coverage-rna", new Command("Coverage", false, 'coverage-rna')],
]);

let RNAParameters = new Map([
    // File Conversion, no parameters
    ["fasta-to-fastq", []],
    ["fast5-to-fastq", []],
    ["bcl2fastq", []],
    ["sam-to-bam", []],

    // Pre Processing
    ["sort-bam", []], // none
    ["index-bam", []], // none
    ["trimming", [
        new Parameter("Perform Adapter Trim", false, 'checkbox', [], 'adapter_trim', false),
            new Parameter("Adapter Trim: Adapter File", true, 'upload', [], 'adapter_trim_file', ''),
        new Parameter("Perform Read Length Trim", false, 'checkbox', [], 'read_length_trim', false),
            new Parameter("Read Length Trim: Minimum Length", true, 'number', [], 'minlen', ''),
        new Parameter("Perform Sliding Window Trim", false, 'checkbox', [], 'sliding_window_trim', false),
            new Parameter("Sliding Window Trim: Start", true,  'number', [], 'window_start', ''),
            new Parameter("Sliding Window Trim: End", true,  'number', [], 'window_end', ''),
        new Parameter("Perform Leading Trim", false, 'checkbox', [], 'leading_trim', false),
            new Parameter("Leading Trim: Input", true, 'number', [], 'leading', ''),
        new Parameter("Perform Trailing Trim", false, 'checkbox', [], 'trailing_trim', false),
            new Parameter("Trailing Trim: Input", true, 'checkbox', [], 'trailing', false),
    ]],
    ["alignment-mapping", [
        new Parameter("Version", false, 'select', ['Bowtie2', 'TopHat', 'Star', 'Salmon', 'HISAT2'], 'version', 'Bowtie2'),
    ]],

    // Quality Control
    ["fastqc", []], // none
    ["fastq-screen", []],
    ["alignmentqc", []],
    ["rseqc", []],

    // Analysis
    ["quantification", [
        new Parameter("Version", false, 'select', ['htseq-count', 'featureCounts'], 'version', 'htseq-count'),
    ]],
    ["normal-diff", []],
    ["quality-yield-metrics", []],
    ["rseq-metrics", []],
    ["flagstats-rna", []],
    ["coverage-rna", []],

]);
    
module.exports = {DNACategories, DNACommands, DNAParameters, RNACategories, RNACommands, RNAParameters};