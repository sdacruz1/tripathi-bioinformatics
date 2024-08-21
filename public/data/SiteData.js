const { Command, Parameter, Downloadable, Executable } = require('./Structure');

// ---- DNA ---- //

const DNACategories = [
    { title: "File Conversion", entries: ['fasta-to-fastq', 'fast5-to-fastq', 'bcl2fastq', 'sam-to-bam'] },
    { title: "File Processing", entries: ['sort-bam', 'index-bam', 'fastqc', 'trimming', 'alignment', 'mark-or-remove-duplicates'] },
    { title: "Statistics", entries: ['add-or-replace-read-groups', 'bam-index-stats', 'flag-stats'] },
    { title: "Summary and Graphs", entries: ['alignment-summary', 'gc-bias-summary', 'insert-size-summary'] },
    { title: "Sequencing", entries: ['create-sequence-dictionary', 'sequence-depth', 'sequence-coverage'] },
];

const refGenomeOptions = [['Human', 'hgch38_index'], ['Mouse', 'House_Mouse_GRCm39'], ['Ecoli', 'Ecoli'], ['HIV', 'HIV'], ['Pig', 'Pig'], ['Staphylococcus_aureus', 'Staphylococcus_aureus']];

let DNAExecutables = new Map([
    ["fasta-to-fastq", new Executable(
        false,
        "Convert FastA to FastQ",
        'dummy', // use uploaded
        'seqtk',
        ['', ''],
        [],
        "seqtk seq -F # <main_file> > usr/output/convertedFastQ.fq",
        [new Downloadable(false, false, "FastQ from FastA", "convertedFastQ.fq")]
    )],
    ["fast5-to-fastq", new Executable(
        false,
        "Convert Fast5 to FastQ",
        'dummy', // use uploaded
        'fast5-to-fastq',
        ['output_dir', 'Fast5ConvertedToFastQ'],
        [],
        "/usr/uploads/<main_file> > usr/output/Fast5ConvertedToFastQ",
        [new Downloadable(false, false, "FastQ from Fast5", "Fast5ConvertedToFastQ")]
    )],
    ["bcl2fastq", new Executable(
        false,
        "Convert BCL to FastQ",
        'dummy', // use uploaded
        'illumina-bcl2fastq',
        ['output_dir', 'BCLConvertedToFastQ'],
        [],
        "/bin/bash -c bcl2fastq --runfolder-dir /usr/uploads/<main_file> --output-dir /usr/output/BCLConvertedToFastQ --processing-threads 1 --no-lane-splitting",
        [new Downloadable(false, false, "BCL2FastQ Conversion Output", "BCLConvertedToFastQ")]
    )],
    ["fastqc", new Executable(
        false,
        "Run FastQC",
        'convertedFastQ.fq',
        'fastqc',
        ['output_dir', 'FastQC_Output'],
        [],
        "fastqc <main_file> -o usr/output/FastQC_Output",
        [new Downloadable(false, false, "FastQC Results", "FastQC_Output")]
    )],
    ["trimming", new Executable(
        false,
        "Trimming",
        '',
        'picard',
        ['paired', ['output_files', ['usr/output/output_paired_Read1.fastq.gz usr/output/output_unpaired_Read1.fastq.gz usr/output/output_paired_Read2.fastq.gz usr/output/output_unpaired_Read2.fastq.gz', 'usr/output/output_single_trim.fastq.gz']]],
        [new Parameter("Perform Adapter Trim", 'select', [['Yes', 'ILLUMINACLIP:<adapter_filepath>2:30:2'], ['No', '']], 'adapter_trim', ''),
        new Parameter("Perform Read Length Trim", 'select', [['Yes', 'MINLEN:'], ['No', '']], 'read_length_trim', ''),
        new Parameter("Read Length Trim: Minimum Length", 'number', [], 'minlen', ''),
        new Parameter("Perform Sliding Window Trim", 'select', [['Yes', 'SLIDINGWINDOW:<window_start>:<window_end>'], ['No', '']], 'sliding_window_trim', ''),
        new Parameter("Window Start", 'number', [], 'window_start', ''),
        new Parameter("Window End", 'number', [], 'window_end', ''),
        new Parameter("Perform Leading Trim", 'select', [['Yes', 'LEADING:'], ['No', '']], 'leading_trim', ''),
        new Parameter("Leading Trim: Input", 'number', [], 'leading_int', ''),
        new Parameter("Perform Trailing Trim", 'select', [['Yes', 'TRAILING:'], ['No', '']], 'trailing_trim', ''),
        new Parameter("Trailing Trim: Input", 'number', [], 'trailing_int', ''),
        ],
        "java -jar picard.jar <paired_or_single> <main_file> <output_files> <adapter_trim> <read_length_trim><minlen> <sliding_window_trim> <leading_trim><leading_int> <trailing_trim><trailing_int>",
        [new Downloadable(false, false, "Trimming Paired Read 1", "output_paired_Read1.fastq.gz"),
        new Downloadable(false, false, "Trimming Unpaired Read 1", "output_unpaired_Read1.fastq.gz"),
        new Downloadable(false, false, "Trimming Paired Read 2", "output_paired_Read2.fastq.gz"),
        new Downloadable(false, false, "Trimming Unpaired Read 2", "output_unpaired_Read2.fastq.gz"),
        new Downloadable(false, false, "Single Trimming Output", "output_single_trim.fastq.gz"),
        ]
    )],
    ["alignment-bwamem", new Executable(
        false,
        "Alignment BWA",
        '',
        'bwamem',
        ['', ''],
        [new Parameter("Alignment BWA Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "mem usr/refs/<ref_genome> <main_file_read1> <main_file_read2> > usr/output/AlignedSAM.sam",
        [new Downloadable(false, false, "Aligned SAM File: BWA", "AlignedSAM.sam")]
    )],
    ["alignment-bowtie", new Executable(
        false,
        "Alignment Bowtie",
        '',
        'bowtie',
        ['', ''],
        [new Parameter("Alignment Bowtie Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "bowtie -x usr/refs/<ref_genome>/<ref_genome>.fna --12 r <main_file_read1>, <main_file_read2> -S usr/output/AlignedSAM.sam",
        [new Downloadable(false, false, "Aligned SAM File: Bowtie", "AlignedSAM.sam")]
    )],
    ["alignment-bowtie2", new Executable(
        false,
        "Alignment Bowtie2",
        '',
        'bowtie2',
        ['', ''],
        [new Parameter("Alignment Bowtie2 Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "bowtie2 -x usr/refs/<ref_genome>/<ref_genome>.fna -U r <main_file_read1>, <main_file_read2> -S usr/output/AlignedSAM.sam",
        [new Downloadable(false, false, "Aligned SAM File: Bowtie2", "AlignedSAM.sam")]
    )],
    ["sam-to-bam", new Executable(
        false,
        "Convert SAM to BAM",
        'AlignedSAM.sam',
        'samtools',
        ['', ''],
        [],
        "samtools view -b <main_file> > ConvertedToBAM.bam",
        [new Downloadable(false, false, "Converted BAM File", "ConvertedToBAM.bam")]
    )],
    ["sort-bam-file", new Executable(
        false,
        "Sort BAM File",
        'ConvertedToBAM.bam',
        'samtools',
        ['', ''],
        [],
        "samtools sort <main_file> -o SortedBAM.bam",
        [new Downloadable(false, false, "Sorted BAM", "SortedBAM.bam")]
    )],
    ["index-bam-file", new Executable(
        false,
        "Index BAM File",
        'SortedBAM.bam',
        'samtools',
        ['readout', "BAM_Index_Stats.txt"],
        [],
        "samtools index <main_file>",
        [new Downloadable(false, false, "BAM Index", "SortedBAM.bam.bai"),
        new Downloadable(false, false, "BAM Index", "SortedBAM.bam.csi")
        ]
    )],
    ["mark-or-remove-duplicates", new Executable(
        false,
        "Mark Or Remove Duplicates",
        'SortedBAM.bam',
        'picard',
        ['', ''],
        [new Parameter("Remove Duplicates", 'select', [['Do Not Remove', ''], ['Remove All', '--REMOVE_DUPLICATES true'], ['Remove Sequencing Duplicates Only', '--REMOVE_SEQUENCING_DUPLICATES true']], 'remove_dupes', '')],
        "java -jar picard.jar MarkDuplicates -I <main_file> -O usr/output/MarkedDuplicatesBAM.bam -M usr/output/DuplicateMetrics.txt <remove_dupes>",
        [new Downloadable(false, false, "Marked Duplicates BAM", "MarkedDuplicatesBAM.bam"),
        new Downloadable(false, false, "Duplicate Metrics", "DuplicateMetrics.txt")
        ]
    )],
    ["add-or-replace-read-groups", new Executable(
        false,
        "Add or Replace Read Groups",
        'SortedBAM.bam',
        'samtools',
        ['', ''],
        [new Parameter("Overwrite Existing Read Groups", 'select', [['Overwrite All', 'overwrite_all'], ['Orphan Only', 'orphan_only']], 'edit_mode', 'overwrite_all'),
        new Parameter("New Read Groups Line", 'text', [], 'new_readgroup_line', '@RG\tID:sample1\tSM:sample1\tLB:library1\tPL:illumina'),
        ],
        "samtools addreplacerg -r <new_readgroup_line> -m <edit_mode> -u -o usr/output/RG_bam.bam <main_file>",
        [new Downloadable(false, false, "Read Groups Added/Replaced BAM", "RG_bam.bam")]
    )],
    ["bam-index-stats", new Executable(
        false,
        "Extract BAM Index Stats",
        'SortedBAM.bam',
        'samtools',
        ['', ''],
        [],
        "samtools idxstats <main_file>",
        [new Downloadable(false, false, "BAM Index Stats", "BAM_Index_Stats.txt")]
    )],
    ["alignment-summary", new Executable(
        false,
        "Generate Alignment Summary",
        'output/SortedBAM.bam',
        'picard',
        ['', ''],
        [new Parameter("Alignment Summary Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "java -jar ../picard/picard.jar CollectAlignmentSummaryMetrics -I usr/<main_file> -O usr/output/AlignmentSummary.txt -R usr/refs/<ref_genome>/<ref_genome>.fna",
        [new Downloadable(false, false, "Alignment Summary", "AlignmentSummary.txt")]
    )],
    ["gc-bias-summary", new Executable(
        false,
        "Generate GC Bias Summary",
        'output/SortedBAM.bam',
        'picard',
        ['', ''],
        [new Parameter("GCBias Summary Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "java -jar ../picard/picard.jar CollectGcBiasMetrics -I usr/<main_file> -O usr/output/GC_BIAS_Metrics.txt -CHART usr/output/GC_BIAS_OutputChart.txt -S usr/output/GC_BIAS_SummaryOutput.txt -R usr/refs/<ref_genome>/<ref_genome>.fna",
        [new Downloadable(false, false, "GC Bias Metrics", "GC_BIAS_Metrics.txt"),
        new Downloadable(false, true, "GC Bias Output Chart", "GC_BIAS_OutputChart.pdf"),
        new Downloadable(false, false, "GC Bias Summary Output", "GC_BIAS_SummaryOutput.txt")
        ]
    )],
    ["insert-size-summary", new Executable(
        false,
        "Generate Insert Size Summary",
        'output/SortedBAM.bam',
        'picard',
        ['', ''],
        [],
        "java -jar ../picard/picard.jar CollectInsertSizeMetrics -I usr/<main_file> -O usr/output/Insert_Size_RawData.txt -H usr/output/Insert_Size_Histogram.pdf -M 0.5",
        [new Downloadable(false, false, "Insert Size Raw Data", "Insert_Size_RawData.txt"),
        new Downloadable(false, false, "Insert Size Histogram", "Insert_Size_Histogram.pdf")
        ]
    )],
    ["create-sequence-dictionary", new Executable(
        false,
        "Create Sequence Dictionary",
        '',
        'samtools',
        ['', ''],
        [new Parameter("Sequence Dictionary Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "samtools dict usr/refs/<ref_genome>/<ref_genome>.fna -o usr/output/SeqDict.txt",
        [new Downloadable(false, false, "Sequence Dictionary", "SeqDict.txt")]
    )],
    ["flag-stats", new Executable(
        false,
        "Examine Flag Stats",
        'SortedBAM.bam',
        'samtools',
        ['readout', "Flag_Stats.txt"],
        [],
        "samtools flagstat <main_file>",
        [new Downloadable(false, false, "Flag Stats", "Flag_Stats.txt")]
    )],
    ["sequence-depth", new Executable(
        false,
        "Examine Sequencing Depth",
        'SortedBAM.bam',
        'samtools',
        ['', ''],
        [],
        "samtools depth -o usr/output/Seq_Depth.txt <main_file>",
        [new Downloadable(false, false, "Sequence Depth Data", "Seq_Depth.txt")]
    )],
    ["sequence-coverage", new Executable(
        false,
        "Examine Sequencing Coverage",
        'SortedBAM.bam',
        'samtools',
        ['', ''],
        [],
        "samtools coverage -o usr/output/Seq_Coverage.txt -m <main_file>",
        [new Downloadable(false, false, "Sequence Coverage Histogram", "Seq_Coverage_Histogram.pdf"),
        new Downloadable(false, false, "Sequence Coverage Data", "Seq_Coverage_Data.txt")
        ]
    )],
]);

// ---- RNA ---- //

// const RNACategories = [
//     { title: "File Conversion", entries: ['fasta-to-fastq', 'fast5-to-fastq', 'bcl2fastq', 'sam-to-bam'] },
//     { title: "Pre Processing", entries: ['sort-bam', 'index-bam', 'trimming', 'alignment-mapping'] },
//     { title: "Quality Control", entries: ['fastqc', 'fastq-screen', 'alignmentqc', 'rseqc'] },
//     { title: "Analysis", entries: ['quantification', 'normal-diff', 'quality-yield-metrics', 'rseq-metrics', 'flagstats-rna', 'coverage-rna'] },
// ];

// let RNAExecutables = new Map([
//     ["fastqscreen", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // will be an input fastq
//         'fastqscreen',
//         ['', ''],
//         [['fastq_screen_config', ''], ['output_dir', '']],
//         "fastq_screen --conf <fastq_screen_config> <main_file> --outdir <output_dir>",
//         [new Downloadable(false, false, "FastQScreen Results", "CBTT")]
//     )],
//     ["bowtie2-alignment", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // will be an input reads.fq
//         'bowtie2',
//         ['', ''],
//         [['index', ''], ['output_sam', 'Bowtie_Alignment_RNA.sam']],
//         "bowtie2 -x <index> -U <main_file> -S <output_sam>",
//         [new Downloadable(false, false, "Bowtie2 Alignment Results", "Bowtie_Alignment_RNA.sam")]
//     )],
//     ["star-alignment", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // will be an input reads.fq
//         'star',
//         ['', ''],
//         [['ref_genome', ''], ['read1_fastq', ''], ['read2_fastq', ''], ['output_prefix', '']],
//         "STAR --genomeDir usr/refs/<ref_genome>/<ref_genome>.fna --readFilesIn <read1_fastq> <read2_fastq> --outFileNamePrefix <output_prefix>",
//         [new Downloadable(false, false, "Star Alignment Results", "CBTT")]
//     )],
//     ["salmon-alignment", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // will be an input reads.fq
//         'salmon',
//         ['', ''],
//         [['quant', ''], ['index', ''], ['read1_fq', ''], ['read2_fq', ''], ['output_dir', '']],
//         "salmon <quant> -i <index> -l A -1 <read1_fq> -2 <read2_fq> -o <output_dir>",
//         [new Downloadable(false, false, "Salmon Alignment Results", "CBTT")]
//     )],
//     ["hisat2-alignment", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // will be an input reads.fq
//         'hisat2',
//         ['', ''],
//         [['index', ''], ['read1_fq', ''], ['read2_fq', ''], ['output_sam', 'hisat_output.sam']],
//         "hisat2 -x <index> -1 <read1_fq> -2 <read2_fq> -S usr/output/<output_sam>",
//         [new Downloadable(false, false, "HISAT2 Alignment Results", "hisat_output.sam")]
//     )],
//     ["htseq-alignment", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // bam file? how many
//         'htseq',
//         ['', ''],
//         [['bam', ''], ['gene_id', ''], ['aligned_reads_bam', ''], ['genes_gtf', ''], ['counts_txt', 'htseq_counts.txt']],
//         "htseq-count -f <bam> -s no -t exon -i <gene_id> <aligned_reads_bam> <genes_gtf> > <counts_txt>",
//         [new Downloadable(false, false, "Htseq-counts Alignment Results", "htseq_counts.txt")]
//     )],
//     ["featureCounts-alignment", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // will be an input reads.fq
//         'featureCount',
//         ['', ''],
//         [['genes_gtf', ''], ['counts_txt', 'feature_counts.txt'], ['aligned_reads_bam', '']],
//         "featureCounts -a <genes_gtf> -o <counts_txt> <aligned_reads_bam>",
//         [new Downloadable(false, false, "Feature-counts Alignment Results", "feature_counts.txt")]
//     )],
//     ["rna-coverage-r", new Executable(
//         "fasta-to-fastq",
//         false,
//         '', // will be an input CBTT
//         'r',
//         ['', ''],
//         [['reads', ''], ['transcripts', '']],
//         "library(GenomicAlignments) coverage <- coverageByTranscript(<reads>, <transcripts>)",
//         [new Downloadable(false, false, "RNA R Coverage Results", "CBTT")]
//     )],
//     ["rna-coverage-samtools", new Executable(
//         "fasta-to-fastq",
//         false,
//         'SortedBAM.bam', // will be an input bam
//         'samtools',
//         ['', ''],
//         [['exons_bed', ''], ['coverage_txt', 'rna_coverage.txt']],
//         "samtools coverage -b <exons_bed> <main_file> > coverage_txt",
//         [new Downloadable(false, false, "RNA Samtools Coverage Results", "rna_coverage.txt")]
//     )],
// ]);

// let RNACommands = new Map([
//     // File Conversion
//     ["fasta-to-fastq", new Command("Convert FastA to FastQ", false, 'fasta-to-fastq')],
//     ["fast5-to-fastq", new Command("Convert Fast5 to FastQ", false, 'fast5-to-fastq')],
//     ["bcl2fastq", new Command("Convert BCL to FastQ", false, 'bcl2fastq')],
//     ["sam-to-bam", new Command("Convert SAM to BAM", false, 'sam-to-bam')],

//     // PreProcessing
//     ["sort-bam", new Command("Sort BAM", false, 'sort-bam')],
//     ["index-bam", new Command("Index BAM", false, 'index-bam')],
//     ["trimming", new Command("Trimming", false, 'trimming')],
//     ["alignment-mapping", new Command("Aligment / Mapping", false, 'alignment-mapping')],

//     // Quality Control
//     ["fastqc", new Command("FastQC", false, 'fastqc')],
//     ["fastq-screen", new Command("FastQScreen", false, 'fastq-screen')],
//     ["alignmentqc", new Command("Alignment QC", false, 'alignmentqc')],
//     ["rseqc", new Command("RSeQC", false, 'rseqc')],

//     // Analysis
//     ["quantification", new Command("Quantification", false, 'quantification')],
//     ["normal-diff", new Command("Normalization and Differential Expression", false, 'normal-diff')],
//     ["quality-yield-metrics", new Command("Collect Quality Yield Metrics", false, 'quality-yield-metrics')],
//     ["rseq-metrics", new Command("Collect RNA RSeq Metrics", false, 'rseq-metrics')],
//     ["flagstats-rna", new Command("Flagstat", false, 'flagstats-rna')],
//     ["coverage-rna", new Command("Coverage", false, 'coverage-rna')],
// ]);

// let RNAParameters = new Map([
//     // File Conversion, no parameters
//     ["fasta-to-fastq", []],
//     ["fast5-to-fastq", []],
//     ["bcl2fastq", []],
//     ["sam-to-bam", []],

//     // Pre Processing
//     ["sort-bam", []], // none
//     ["index-bam", []], // none
//     ["trimming", [
//         new Parameter("Perform Adapter Trim", false, 'checkbox', [], 'adapter_trim', false),
//         new Parameter("Adapter Trim: Adapter File", true, 'upload', [], 'adapter_trim_file', ''),
//         new Parameter("Perform Read Length Trim", false, 'checkbox', [], 'read_length_trim', false),
//         new Parameter("Read Length Trim: Minimum Length", true, 'number', [], 'minlen', ''),
//         new Parameter("Perform Sliding Window Trim", false, 'checkbox', [], 'sliding_window_trim', false),
//         new Parameter("Sliding Window Trim: Start", true, 'number', [], 'window_start', ''),
//         new Parameter("Sliding Window Trim: End", true, 'number', [], 'window_end', ''),
//         new Parameter("Perform Leading Trim", false, 'checkbox', [], 'leading_trim', false),
//         new Parameter("Leading Trim: Input", true, 'number', [], 'leading', ''),
//         new Parameter("Perform Trailing Trim", false, 'checkbox', [], 'trailing_trim', false),
//         new Parameter("Trailing Trim: Input", true, 'checkbox', [], 'trailing', false),
//     ]],
//     ["alignment-mapping", [
//         new Parameter("Version", false, 'select', ['Bowtie2', 'TopHat', 'Star', 'Salmon', 'HISAT2'], 'version', 'Bowtie2'),
//     ]],

//     // Quality Control
//     ["fastqc", []], // none
//     ["fastq-screen", []],
//     ["alignmentqc", []],
//     ["rseqc", []],

//     // Analysis
//     ["quantification", [
//         new Parameter("Version", false, 'select', ['htseq-count', 'featureCounts'], 'version', 'htseq-count'),
//     ]],
//     ["normal-diff", []],
//     ["quality-yield-metrics", []],
//     ["rseq-metrics", []],
//     ["flagstats-rna", []],
//     ["coverage-rna", []],

// ]);

// module.exports = { DNACategories, DNACommands, DNAParameters, RNACategories, RNACommands, RNAParameters };
module.exports = { DNAExecutables };