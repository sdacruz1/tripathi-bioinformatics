const {Command, Parameter, Downloadable, Executable} = require('./Structure');

// ---- DNA ---- //

const DNACategories = [
    { title: "File Conversion", entries: ['fasta-to-fastq', 'fast5-to-fastq', 'bcl2fastq', 'sam-to-bam']},
    { title: "File Processing", entries: ['sort-bam', 'index-bam', 'fastqc', 'trimming', 'alignment', 'mark-or-remove-duplicates']},
    { title: "Statistics", entries: ['add-or-replace-read-groups', 'bam-index-stats', 'flag-stats']},
    { title: "Summary and Graphs", entries: ['alignment-summary', 'gc-bias-summary', 'insert-size-summary']},
    // { title: "Metrics", entries: ['', '', '']}, // CBTT
    { title: "Sequencing", entries: ['create-sequence-dictionary', 'sequence-depth', 'sequence-coverage']},
];

let DNAExecutables = new Map([
    ["fastqc", new Executable(
        '',
        'fastqc',
        [['uploaded_file', '']],
        "fastqc <uploaded_file>",
        [new Downloadable(false, false, "FastQC Results", "CBTT")]
    )],
    ["fasta-to-fastq", new Executable(
        'convertedFastQ.fq',
        'seqtk',
        [['char_to_replace_quality_scores', '#'], ['in_fa', ''], ['out_fq', '']],
        "seqtk seq -F <char_to_replace_quality_scores> <in_fa> > <out_fq>",
        [new Downloadable(false, false, "FastQ from FastA", "FastAConvertedToFastQ.fq")]
    )],
    ["sam-to-bam", new Executable(
        'AlignedSAM.sam',
        'samtools',
        [['output_file', 'ConvertedToBAM.bam']],
        "samtools view -b <main_file> > <output_file>",
        [new Downloadable(false, false, "Converted BAM File", "AlignedSAM.sam")]
    )],
    ["sort-bam-file", new Executable(
        'ConvertedToBAM.bam',
        'samtools',
        [['output_file', 'ConvertedToBAM.bam']],
        "samtools sort <main_file> -o <output_file>",
        [new Downloadable(false, false, "Sorted BAM", "SortedBAM.bam")]
    )],
    ["index-bam-file", new Executable(
        'SortedBAM.bam',
        'samtools',
        [],
        "samtools index <main_file>",
        [new Downloadable(false, false, "BAM Index", "SortedBAM.bam.bai"),
         new Downloadable(false, false, "BAM Index", "SortedBAM.bam.csi")
        ]
    )],
    ["add-or-replace-read-groups", new Executable(
        'SortedBAM.bam',
        'samtools',
        [['new_readgroup_line', '@RG\tID:sample1\tSM:sample1\tLB:library1\tPL:illumina'], ['edit_mode', ''],['output_file', 'RG_bam.bam']],
        "samtools addreplacerg -r <new_readgroup_line> -m <edit_mode> -u -o <output_file> <main_file>",
        [new Downloadable(false, false, "Read Groups Added/Replaced BAM", "RG_bam.bam")]
    )],
    ["bam-index-stats", new Executable(
        'SortedBAM.bam',
        'samtools',
        [],
        "samtools idxstats <main_file>",
        [new Downloadable(false, false, "BAM Index Stats", "BAM_Index_Stats.txt")]
    )],
    ["alignment-summary", new Executable(
        'output/SortedBAM.bam',
        'picard',
        [['ref_genome', 'hgch38_index'], ['output_file', 'AlignmentSummary.txt']],
        "java -jar ../picard/picard.jar CollectAlignmentSummaryMetrics -I usr/<main_file> -O usr/output/<output_file> -R usr/refs/<ref_genome>/<ref_genome>.fna",
        [new Downloadable(false, false, "Alignment Summary", "AlignmentSummary.txt")]
    )],
    ["gc-bias-summary", new Executable(
        'output/SortedBAM.bam',
        'picard',
        [['ref_genome', 'hgch38_index'], ['output_metrics_txt', 'GC_BIAS_Metrics.txt'], ['output_chart', 'GC_BIAS_OutputChart.txt'], ['output_summary', 'GC_BIAS_SummaryOutput.txt']],
        "java -jar ../picard/picard.jar CollectGcBiasMetrics -I usr/<main_file> -O usr/output/<output_metrics_txt> -CHART usr/output/<output_chart> -S usr/output/<output_summary> -R usr/refs/<ref_genome>/<ref_genome>.fna",
        [new Downloadable(false, false, "GC Bias Metrics", "GC_BIAS_Metrics.txt"),
            new Downloadable(false, true, "GC Bias Output Chart", "GC_BIAS_OutputChart.pdf"),
            new Downloadable(false, false, "GC Bias Summary Output", "GC_BIAS_SummaryOutput.txt")
        ]
    )],
    ["insert-size-summary", new Executable(
        'output/SortedBAM.bam',
        'picard',
        [['output_raw_data', 'Insert_Size_RawData.txt'], ['output_histogram_name', 'Insert_Size_Histogram.pdf']],
        "java -jar ../picard/picard.jar CollectInsertSizeMetrics -I usr/<main_file> -O <output_raw_data> -H <output_histogram_name> -M 0.5",
        [new Downloadable(false, false, "Insert Size Raw Data", "Insert_Size_RawData.txt"),
            new Downloadable(false, false, "Insert Size Histogram", "Insert_Size_Histogram.pdf")
        ]
    )],
    ["create-sequence-dictionary", new Executable(
        '',
        'samtools',
        [['ref_genome', 'hgch38_index'], ['output_file', 'SeqDict.txt']],
        "samtools dict usr/refs/<ref_genome>/<ref_genome>.fna -o usr/output/<output_file>",
        [new Downloadable(false, false, "Sequence Dictionary", "SeqDict.txt")]
    )],
    ["flag-stats", new Executable(
        'SortedBAM.bam',
        'samtools',
        [],
        "samtools flagstat <main_file>",
        [new Downloadable(false, false, "Flag Stats", "Flag_Stats.txt")]
    )],
    ["sequence-depth", new Executable(
        'SortedBAM.bam',
        'samtools',
        ['output_file', 'Seq_Depth.txt'],
        "samtools depth -o <output_file> <main_file>",
        [new Downloadable(false, false, "Sequence Depth Data", "Seq_Depth.txt")]
    )],
    ["sequence-coverage", new Executable(
        'SortedBAM.bam',
        'samtools',
        ['output_file', 'Seq_Coverage.txt'],
        "samtools coverage -o <output_file> -m <main_file>",
        [new Downloadable(false, false, "Sequence Coverage Histogram", "Seq_Coverage_Histogram.pdf"),
            new Downloadable(false, false, "Sequence Coverage Data", "Seq_Coverage_Data.txt")
        ]
    )],
]);

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
    ["fasta-to-fastq", [], [
        new Downloadable(false, false, "FastQ From FastA", "../output/convertedFastQ.fq")
    ]],
    ["fast5-to-fastq", [], [
        new Downloadable(false, false, "FastQ From Fast5", "../output/FastQFromFast5.fastq")
    ]],
    ["bcl2fastq", [], [
        new Downloadable(false, false, "BCL2FastQ Output", "bcl2fastq_Output")
    ]],
    ["sam-to-bam", [], [
        new Downloadable(false, false, "Converted BAM File", "AlignedSAM.sam")
    ]],

    // File Processing
    ["sort-bam", [] [
        new Downloadable(false, false, "Sorted BAM", "SortedBAM.bam")
    ]],
    ["index-bam", [], [
        new Downloadable(false, false, "Indexed BAM", "SortedBAM.bam.bai"),
        new Downloadable(false, false, "Indexed BAM", "SortedBAM.bam.csi")
    ]], // none
    ["fastqc", [], [
        // new Downloadable() // CBTT
    ]],
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
    ], [
        new Downloadable(false, false, "Trimming Paired Read 1", "output_paired_Read1.fastq.gz"),
        new Downloadable(false, false, "Trimming Paired Read 2", "output_paired_Read2.fastq.gz"),
        new Downloadable(false, false, "Trimming Unpaired Read 1", "output_unpaired_Read1.fastq.gz"),
        new Downloadable(false, false, "Trimming Unpaired Read 2", "output_unpaired_Read2.fastq.gz"),
    ]],
    ["alignment", [
        new Parameter("Alignment Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
        new Parameter("Alignment Mode", false, 'select', ['BWA', 'Bowtie', 'Bowtie2'], 'mode', 'BWA'),
    ], [
        new Downloadable(false, false, "Aligned SAM File", "AlignedSAM.sam")
    ]],
    ["mark-or-remove-duplicates", [
        new Parameter("Remove Duplicates", false, 'select', ['Do Not Remove', 'Remove All', 'Remove Sequencing Duplicates Only'], 'remove_dupes', 'Do Not Remove'),
    ], [
        new Downloadable(false, false, "Marked Duplicates BAM", "MarkedDuplicatesBAM.bam"),
        new Downloadable(false, false, "Duplicate Metrics", "DuplicateMetrics.txt")
    ]],

    // Statistics
    ["add-or-replace-read-groups", [
        new Parameter("Overwrite Existing Read Groups", false, 'checkbox', [], 'editMode', false),
        new Parameter("New Read Groups Line", false, 'text', [], 'newReadGroupLine', ''),
    ], [
        new Downloadable(false, false, "Read Groups Added/Replaced BAM", "RG_bam.bam")
    ]],
    ["bam-index-stats", [], [
        new Downloadable(false, false, "BAM Index Stats", "BAM_Index_Stats.txt")
    ]],
    ["flag-stats", [], [
        new Downloadable(false, false, "Flag Stats", "Flag_Stats.txt")
    ]],

    // Summary and Graphs
    ["alignment-summary", [
        new Parameter("Alignment Summary Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
        new Parameter("Create Alignment Graph", false, 'checkbox', [], 'isVisual', false),
    ], [
        new Downloadable(false, false, "Alignment Summary", "AlignmentSummary.txt")
    ]], 
    ["gc-bias-summary", [
        new Parameter("GC Bias Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
        new Parameter("Create Bias Graph", false, 'checkbox', [], 'isVisual', false),
    ], [
        new Downloadable(false, false, "GC Bias Metrics", "GC_BIAS_Metrics.txt"),
        new Downloadable(false, false, "GC Bias Output Chart", "GC_BIAS_OutputChart.pdf"),
        new Downloadable(false, false, "GC Bias Summary Output", "GC_BIAS_SummaryOutput.txt"),
    ]],
    ["insert-size-summary", [
        new Parameter("Create Size Graph", false, 'checkbox', [], 'isVisual', false),
    ], [
        new Downloadable(false, false, "Insert Size Raw Data", "Insert_Size_RawData.txt"),
        new Downloadable(false, false, "Insert Size Histogram", "Insert_Size_Histogram.pdf"),
    ]],

    // Sequencing
    ["create-sequence-dictionary", [
        new Parameter("Sequence Reference Genome", false, 'select', ['Human', 'Mouse', 'Ecoli', 'HIV', 'Pig', 'Staphylococcus_aureus'], 'ref_genome', 'Human'),
    ]],
    ["sequence-depth", [], [
        new Downloadable(false, false, "Sequence Depth Data", "Seq_Depth.txt")
    ]], // none
    ["sequence-coverage", [
        new Parameter("Create Sequence Graph", false, 'checkbox', [], 'isVisual', false),
    ], [
        new Downloadable(false, false, "Sequence Coverage Histogram", "Seq_Coverage.txt"),
        new Downloadable(false, true, "Sequence Coverage Data", "Seq_Coverage.txt"),
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