const { Command, Parameter, Downloadable, Executable } = require('./Structure');

let FastAToFastQ = new Executable(
    false,
    "Convert FastA to FastQ",
    'dummy', // use uploaded
    'seqtk',
    ['', ''],
    [],
    "seqtk seq -F # <main_file> > usr/output/convertedFastQ.fq",
    [new Downloadable(false, false, "FastQ from FastA", "convertedFastQ.fq")]
);
let Fast5ToFastQ = new Executable(
    false,
    "Convert Fast5 to FastQ",
    'dummy', // use uploaded
    'fast5-to-fastq',
    ['output_dir', 'Fast5ConvertedToFastQ'],
    [],
    "/usr/uploads/<main_file> > usr/output/Fast5ConvertedToFastQ",
    [new Downloadable(false, false, "FastQ from Fast5", "Fast5ConvertedToFastQ")]
);
let BCL2FastQ = new Executable(
    false,
    "Convert BCL to FastQ",
    'dummy', // use uploaded
    'illumina-bcl2fastq',
    ['output_dir', 'BCLConvertedToFastQ'],
    [],
    "/bin/bash -c bcl2fastq --runfolder-dir /usr/uploads/<main_file> --output-dir /usr/output/BCLConvertedToFastQ --processing-threads 1 --no-lane-splitting",
    [new Downloadable(false, false, "BCL2FastQ Conversion Output", "BCLConvertedToFastQ")]
);
let FastQC = new Executable(
    false,
    "Run FastQC",
    'convertedFastQ.fq',
    'fastqc',
    ['output_dir', 'FastQC_Output'],
    [],
    "fastqc <main_file> -o usr/output/FastQC_Output",
    [new Downloadable(false, false, "FastQC Results", "FastQC_Output.zip")]
);
let Trimming = new Executable(
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
);
let SamToBam = new Executable(
    false,
    "Convert SAM to BAM",
    'AlignedSAM.sam',
    'samtools',
    ['', ''],
    [],
    "samtools view -b <main_file> > usr/output/ConvertedToBAM.bam",
    [new Downloadable(false, false, "Converted BAM File", "ConvertedToBAM.bam")]
);
let SortBamFile = new Executable(
    false,
    "Sort BAM File",
    'ConvertedToBAM.bam',
    'samtools',
    ['', ''],
    [],
    "sort <main_file> -o usr/output/SortedBAM.bam",
    [new Downloadable(false, false, "Sorted BAM", "SortedBAM.bam")]
);

let IndexBamFile = new Executable(
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
);

// ---- DNA ---- //

const DNACategories = [
    { title: "File Conversion", entries: ['fasta-to-fastq', 'fast5-to-fastq', 'bcl2fastq', 'sam-to-bam'] },
    { title: "File Processing", entries: ['sort-bam-file', 'index-bam-file', 'fastqc', 'trimming', 'alignment-bwamem', 'alignment-bowtie', 'alignment-bowtie2', 'mark-or-remove-duplicates'] },
    { title: "Statistics", entries: ['add-or-replace-read-groups', 'bam-index-stats', 'flag-stats'] },
    { title: "Summary and Graphs", entries: ['alignment-summary', 'gc-bias-summary', 'insert-size-summary'] },
    { title: "Sequencing", entries: ['create-sequence-dictionary', 'sequence-depth', 'sequence-coverage'] },
];

const refGenomeOptions = [['Human', 'hgch38_index'], ['Mouse', 'House_Mouse_GRCm39'], ['Ecoli', 'Ecoli'], ['HIV', 'HIV'], ['Pig', 'Pig'], ['Staphylococcus_aureus', 'Staphylococcus_aureus']];

let DNAExecutables = new Map([
    ["fasta-to-fastq", FastAToFastQ],
    ["fast5-to-fastq", Fast5ToFastQ],
    ["bcl2fastq", BCL2FastQ],
    ["fastqc", FastQC],
    ["trimming", Trimming],
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
    ["sam-to-bam", SamToBam],
    ["sort-bam-file", SortBamFile],
    ["index-bam-file", IndexBamFile],
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
        // new Parameter("New Read Groups Line", 'text', [], 'new_readgroup_line', '@RG\\tID:sample1\\tSM:sample1\\tLB:library1\\tPL:illumina'), NOTE: putting special characters like \ breaks the JSON. If i stringify it first it works, but then whenever i parse it bck the issue will persist again
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

const RNACategories = [
    { title: "File Conversion", entries: ['fasta-to-fastq', 'fast5-to-fastq', 'bcl2fastq', 'sam-to-bam'] },
    { title: "Pre Processing", entries: ['sort-bam', 'index-bam', 'trimming'] },
    { title: "Quality Control", entries: ['fastqc', 'fastq-screen'] },
    { title: "Alignment", entries: ['alignment-bowtie2', 'star-alignment', 'salmon-alignment', 'hisat2-alignment', 'htseq-alignment', 'featureCounts-alignment'] },
    { title: "Coverage", entries: ['rna-coverage-samtools', 'build-index'] },
    // { title: "Quality Control", entries: ['fastqc', 'fastq-screen', 'alignmentqc', 'rseqc'] },
    // { title: "Analysis", entries: ['quantification', 'normal-diff', 'quality-yield-metrics', 'rseq-metrics', 'flagstats-rna', 'coverage-rna'] },
];
let RNAExecutables = new Map([
    ["fasta-to-fastq", FastAToFastQ],
    ["fast5-to-fastq", Fast5ToFastQ],
    ["bcl2fastq", BCL2FastQ],
    ["sam-to-bam", SamToBam],
    ["sort-bam", SortBamFile],
    ["index-bam", IndexBamFile],
    ["trimming", Trimming],
    ["fastqc", FastQC],
    ["fastq-screen", new Executable(
        false,
        "FastQ Screen",
        'convertedFastQ.fq',
        'fastqscreen',
        ['output_dir', 'FastQScreen_Output'],
        [new Parameter("FastQ Screen Configurations", 'select', [['', '']], 'fastq_screen_config', '')], //CBTT
        "fastq_screen --conf <fastq_screen_config> <main_file> --outdir usr/output/FastQScreen_Output",
        [new Downloadable(false, false, "FastQScreen Results", "FastQScreen_Output.zip")]
    )],
    ["alignment-bowtie2", new Executable(
        false,
        "Bowtie2 Alignment",
        '',
        'bowtie2',
        ['', ''],
        [new Parameter("Bowtie2 Alignment Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "bowtie2 -x usr/refs/<ref_genome>/<ref_genome>.fna -U r <main_file_read1>, <main_file_read2> -S usr/output/AlignedSAM.sam",
        [new Downloadable(false, false, "Bowtie 2 Alignment Results", "AlignedSAM.sam")]
    )],
    ["star-alignment", new Executable(
        false,
        "Star Alignment",
        '', // will be an input reads.fq
        'star',
        ['', ''],
        [new Parameter("Star Alignment Reference Genome", 'select', refGenomeOptions, 'ref_genome', 'hgch38_index')],
        "STAR --genomeDir usr/refs/<ref_genome>/<ref_genome>.fna --readFilesIn <main_file_read1> <main_file_read2> --outFileNamePrefix <output_prefix>", // CBTT, choose a static output prefix
        [new Downloadable(false, false, "Star Alignment Results", "CBTT")]
    )],
    ["salmon-alignment", new Executable(
        false,
        "Salmon Alignment",
        '', // will be an input reads.fq
        'salmon',
        ['output_dir', 'Salmon_Alignment_Output'],
        [new Parameter("Salmon Quantity", 'text', [], 'quant', ''),
            new Parameter("Salmon Index", 'text', [], 'index', ''),
        ],
        "salmon <quant> -i <index> -l A -1 <main_file_read1> -2 <main_file_read2> -o usr/output/Salmon_Alignment_Output",
        [new Downloadable(false, false, "Salmon Alignment Results", "Salmon_Alignment_Output.zip")]
    )],
    ["hisat2-alignment", new Executable(
        false,
        "HISAT2 Alignment",
        '', // will be an input reads.fq
        'hisat2',
        ['', ''],
        [new Parameter("HISAT2 Index", 'text', [], 'index', ''),],
        "hisat2 -x <index> -1 <main_file_read1> -2 <main_file_read2> -S usr/output/HITSAT_Alignment_Output.sam",
        [new Downloadable(false, false, "HISAT2 Alignment Results", "hisat_output.sam")]
    )],
    ["htseq-alignment", new Executable(
        false,
        "HTSeq Alignment",
        'Aligned_Reads_BAM.bam', // aligned_reads_bam file
        'htseq',
        ['', ''],
        [new Parameter("Reverse?", 'select', [['yes', 'yes'], ['no', 'no']], 'reverse_or_not', 'no'),
            new Parameter("__On Type", 'select', [['exon', 'exon'], ['intron', 'intron']], 'on_type', 'exon'),
            new Parameter("HTSeq Gene ID", 'text', [], 'gene_id', '664792'),
            new Parameter("HTSeq Genes GTF", 'select', [['Human', 'GRCh38.p14_genomic'], ['Mouse', 'House_Mouse_GRCm39'], ['Ecoli', 'Ecoli'], ['HIV', 'HIV'], ['Pig', 'Pig'], ['Staphylococcus_aureus', 'Staphylococcus_aureus']], 'genes_gtf', 'House_Mouse_GRCm39'),
        ],
        "htseq-count -f bam -s <reverse_or_not> -t <on_type> -i <gene_id> /<main_file> /usr/refs/RNA/<genes_gtf>/<genes_gtf>.gtf > /usr/output/htseq_counts.txt",
        [new Downloadable(false, false, "Htseq-counts Alignment Results", "htseq_counts.txt")]
    )],
    ["featureCounts-alignment", new Executable(
        false,
        "FeatureCounts Alignment",
        'Aligned_Reads_BAM.bam', // aligned_reads_bam file
        'featurecounts:2.0.6',
        ['', ''],
        [new Parameter("FeatureCounts Genes GTF", 'select', [['Human', 'GRCh38.p14_genomic'], ['Mouse', 'House_Mouse_GRCm39'], ['Ecoli', 'Ecoli'], ['HIV', 'HIV'], ['Pig', 'Pig'], ['Staphylococcus_aureus', 'Staphylococcus_aureus']], 'genes_gtf', 'House_Mouse_GRCm39')],
        "featureCounts -a /usr/refs/RNA/<genes_gtf>/<genes_gtf>.gtf -o /usr/output/feature_counts.txt /<main_file>",
        [new Downloadable(false, false, "Feature-counts Alignment Results", "feature_counts.txt")]
    )],
    // ["rna-coverage-r", new Executable(
    //     false,
    //     "RNA Coverage (R)",
    //     '', // will be an input CBTT
    //     'r',
    //     ['', ''],
    //     [new Parameter("CBTT", 'text', [], 'reads', ''), // CBTT?
    //         new Parameter("RNA Coverage Transcripts", 'text', [], 'transcripts', '')], //CBTT
    //     "library(GenomicAlignments) coverage <- coverageByTranscript(<reads>, <transcripts>)",
    //     [new Downloadable(false, false, "RNA R Coverage Results", "CBTT")]
    // )],
    ["rna-coverage-samtools", new Executable(
        false,
        "RNA Coverage Samtools",
        'SortedBAM.bam', // will be an input bam
        'samtools',
        ['', ''],
        [new Parameter("CBTT", 'text', [], 'exons_bed', '')], // CBTT
        "samtools coverage -b <exons_bed> <main_file> > usr/output/rna_coverage_txt",
        [new Downloadable(false, false, "RNA Samtools Coverage Results", "rna_coverage.txt")]
    )],
    ["build-index", new Executable(
        true,
        "Make Index Bowtie2",
        '',
        'bowtie2',
        ['', ''],
        [],
        "bowtie2-build /usr/refs/House_Mouse_GRCm39/House_Mouse_GRCm39.fna /usr/output/bowtie2_index/house_mouse",
        [new Downloadable(false, false, "RNA Samtools Coverage Results", "rna_coverage.txt")]
    )],
]);

module.exports = { DNACategories, DNAExecutables, RNACategories, RNAExecutables };