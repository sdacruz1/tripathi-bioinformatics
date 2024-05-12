var express = require('express');
const multer = require('multer');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const bodyParser = require('body-parser');
const Docker = require('dockerode');

const docker = new Docker();
var router = express.Router();


router.use(bodyParser.urlencoded({ extended: true }));
router.use(bodyParser.json());

// Use express.urlencoded middleware
router.use(express.urlencoded({ extended: true }));

// Serve uploaded images statically
router.use('/images', express.static('images'));
router.use('/output', express.static('output'));
router.use('/ref_genomes', express.static('ref_genomes'));
router.use('/uploads', express.static('uploads'));

// Set up Multer for file uploads
const storage = multer.diskStorage({
  destination: function (req, file, cb) {
    cb(null, 'uploads/'); // Specify the folder where uploaded files will be stored
  },
  filename: function (req, file, cb) {
    cb(null, Date.now() + '-' + file.originalname); // Set the file name to be unique
  },
});
const upload = multer({ storage: storage });

//#region  ... Constants ...

// Categories
let files = ["File Processing",
  [
    // name, toggle, options, file
    ["Convert to FastQ", false, [], []],

    ["Trimming", false, [
      // name, input type, options, selected
      ["Trim Type", "select", ["Adapter Trim", "Read Length Trim", "Quality Score Trim", "Duplicate Trim"], []]
    ], []],

    ["Alignment", false, [
      ["Type of Alignment", "select", ["BWA mem", "Bowtie", "Bowtie 2"], []],
      ["Number of Mismatches", "number", [], [0]],
      ["Discard Unpaired", "checkbox", [], []],
    ], []],

    ["Demultiplex FastQ", false, [], []],

    ["Convert to BAM File", false, [], []],

    ["Cleanup BAM File", false, [
      ["Sort", "checkbox", [], []],
      ["Index", "checkbox", [], []],
    ], []],
  ]];
let stats = ["Statistics",
  [
    ["Add or Replace Read Groups", false, [], []],
    ["Bam Index Stats", false, [], []],
  ]];
let summary = ["Summary",
  [
    ["Alignment Summary", false, [], []],
    ["GC Bias Summary", false, [], []],
    ["Insert Size Summary", false, [], []],
  ]];
let graphs = ["Graphs",
  [
    ["Alignment Graph", false, [], []],
    ["GC Bias Graphs", false, [], []],
    ["Insert Size Graphs", false, [], []],
  ]];
let metrics = ["Metrics",
  [
    ["Quality Yield Metrics", false, [], []],
    ["Whole Genome Sequencing Metrics", false, [], []],
    ["Targeted PCR Metrics", false, [], []],
  ]];
let cleanup = ["Cleanup and Sequencing",
  [
    ["Create Sequence Dictionary", false, [], []],
    ["Mark Duplicates", false, [], []],
    ["Sort Bam File", false, [], []],
    ["Flag Stats", false, [], []],
    ["Sequencing Depth", false, [], []],
    ["Sequencing Coverage", false, [], []],
  ]];

let categories = [files, stats, summary, graphs, metrics, cleanup];

// Program Mode and Input Type
let mode = "";
let inputType = "";

// Stored Outputs
let infoSteps;          // An array that determines which file processing steps will be available in the timeline
let uploadedFilePath = "";   // The path to the original file that was uploaded or the directory, if it was a BCL folder
let uploadedFileType;   // A string representing the original uploaded file type
let fastQConversion;    // The path to the fastQ file result of any conversions (BCL, FastA, Fast5)
let fastQCResults;      // The path to the directory containing the results of running FastQC on the uploaded file
let BAMFile;            // The path to the BAM file result of any conversions
let downloadable_content = [  // An object array holding the result of the steps in the actual timeline
  // {
  //     enabled: true,
  //     has_visual_component: false,
  //     label: "Test Object BAM",
  //     content: "../output/SortedBAM.bam"
  // },
  // {
  //     enabled: false,
  //     has_visual_component: true,
  //     label: "Test Object 2",
  //     content: "../images/Checked_Circle.png"
  // },
];

//#endregion

//#region  ... GET and POST Views ... 

/* Home Page */
router.get('/', function (req, res, next) {
  // Render the 'home' view
  res.render('home', { toolbar_index: 1 });
});

/* File Info */
router.post('/file-information', function (req, res, next) {
  // Set the mode
  mode = req.body.mode;
  res.render('file-info', { toolbar_index: 2});
});

/* DNA Goalposts */
router.post('/dna-goalposts', function (req, res, next) {
  // Store the uploaded file and any conversions
  inputType = req.body.inputType;
  infoSteps = JSON.parse(req.body.infoSteps);
  // uploadedFilePath = req.body.uploadedFilePath;
  uploadedFileType = req.body.uploadedFileType;
  fastQConversion = req.body.fastQConversion;
  fastQCResults = req.body.fastQCResults;

  res.render('dna-goalposts', { toolbar_index: 3, categories, infoSteps });
});

/* DNA Pipeline */
router.post('/dna-pipeline', function (req, res, next) {
  categories = JSON.parse(decodeURIComponent(req.body.categories || '[]'));

  res.render('dna-pipeline', { toolbar_index: 4, categories, infoSteps });
});

/* Running Page */
router.post('/running', function (req, res, next) {
  const categories = JSON.parse(decodeURIComponent(req.body.categories || '[]'));

  res.render('running', { toolbar_index: 5, categories, uploadedFilePath });
});

/* Running Page */
router.get('/running', function (req, res, next) {
  // Render the 'output' view
  res.render('running', { categories, uploadedFilePath, toolbar_index: 5 });
});

/* Output Page */
router.post('/output', function (req, res, next) {
  // Render the 'output' view
  res.render('output', { downloadable_content, toolbar_index: 5 });
});

/* Output Page */
router.get('/output', function (req, res, next) {
  // Render the 'output' view
  res.render('output', { downloadable_content, toolbar_index: 5 });
});

/* Test Page */
router.get('/test-commands', function (req, res, next) {
  res.render('test-commands', { toolbar_index: 1 });
});

//#endregion

//#region  ... Utility ...

// Store Files
// Uploads an array of files. If there is only one, just returns the path
// If there are multiple, puts them in a directory and sends that
router.post('/store-files', upload.array('files'), (req, res) => {
  const uploadedFiles = req.files;

  // Check if the files were uploaded
  if (!uploadedFiles || uploadedFiles.length === 0) {
    return res.status(400).send('No files were uploaded.');
  }

  // If it's just one file, send the new path to that
  if (uploadedFiles.length === 1) {
    uploadedFilePath = path.join(__dirname, '..', 'uploads', uploadedFiles[0].filename);
    return res.status(200).send({ storedPath: uploadedFilePath });
  }

  // Create an empty directory to hold them in
  const BCLdirectory = '/uploads/BCLFiles-' + Date.now();
  if (!fs.existsSync(BCLdirectory)) {
    fs.mkdirSync(BCLdirectory);
  }

  // Move each uploaded file to the empty directory
  uploadedFiles.forEach(file => {
    const destinationPath = path.join(BCLdirectory, file.originalname);
    fs.renameSync(file.path, destinationPath);
  });

  uploadedFilePath = BCLdirectory;
  // Send the path to the directory
  return res.status(200).send({ storedPath: uploadedFilePath });
});

// Route to delete all files in the 'uploads' folder
router.get('/cleanup', (req, res) => {
  const uploadsPath = 'uploads';

  // Delete all files in the 'uploads' folder
  fs.readdir(uploadsPath, (err, files) => {
    if (err) {
      console.error(`Error reading 'uploads' folder: ${err.message}`);
      res.status(500).send('Error cleaning up files');
      return;
    }

    files.forEach((file) => {
      const filePath = `${uploadsPath}/${file}`;

      fs.unlink(filePath, (err) => {
        if (err) {
          console.error(`Error deleting file ${filePath}: ${err.message}`);
        } else {
          console.log(`File deleted: ${filePath}`);
        }
      });
    });

    res.send('All files in the "uploads" folder deleted');
  });
});

// Run BCL2FastQ
// Expects the contents of a BCL File
router.post('/bcl2fastq' , upload.array('files'), (req, res) => {

  // =============== Handle the BCL File Folder ===============
  const uploadedFiles = req.files;

  // Check if the files were uploaded
  if (!uploadedFiles || uploadedFiles.length === 0) {
    return res.status(400).send('No files were uploaded.');
  }

  // Create an empty directory to hold them in
  const BCLdirectory = './uploads/BCLFiles-' + Date.now();
  if (!fs.existsSync(BCLdirectory)) {
    fs.mkdirSync(BCLdirectory);
  }

  // Move each uploaded file to the empty directory
  uploadedFiles.forEach(file => {
    const destinationPath = path.join(BCLdirectory, file.originalname);
    fs.renameSync(file.path, destinationPath);
  });

  // =============== Run BCL To FastQ ===============

  // Create an empty directory to hold the output
  const OutputDirectory = './uploads/bclFile';
  if (!fs.existsSync(OutputDirectory)) {
    fs.mkdirSync(OutputDirectory);
  }

  // Spawn the process
  const BCL2FastQCommand = path.join(__dirname, '..', 'bio_modules', 'FastQC.app', 'Contents', 'MacOS', 'bcl2fastq');
  const BCL2FastQArgs = [' --runfolder-dir ' + BCLdirectory + ' --output-dir ' + OutputDirectory
                      + ' --sample-sheet ' + path.join(BCLdirectory, 'SampleSheet.csv') + ' --barcode-mismatches 1'];
  const runBCL2FastQ = spawn(BCL2FastQCommand, BCL2FastQArgs);

  let outputData = '';

  // Handle the terminal response
  runBCL2FastQ.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runBCL2FastQ.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runBCL2FastQ.on('exit', (code) => {
    console.log(`BCL2FastQ process exited with code ${code}`);
  });

  // Respond to the client
  res.json({ fastQfile:  "nope"});
  // NOTE: recongifgure run-fastqc to take the file name

});

// Run FastQC
// Expects a single FastQ File
router.post('/fastqc', upload.none(), (req, res) => {
  const uploadedFile = req.file;

  // Check if a file was uploaded
  if (!uploadedFile) {
    res.status(400).send('No file uploaded');
    return;
  }

  // Run FastQC
  const FastQCCommand = path.join(__dirname, '..', 'bio_modules', 'FastQC.app', 'Contents', 'MacOS', 'fastqc');
  const FastQArgs = [uploadedFile.path];

  const runFastQC = spawn(FastQCCommand, FastQArgs);

  let outputData = '';

  runFastQC.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runFastQC.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runFastQC.on('exit', (code) => {
    console.log(`FastQC process exited with code ${code}`);
  });

  // Return the results as a JSON
  const trimmed = false; // set this based on the result of the QC
  const demultiplexed = false;

  res.json({ trimmed, demultiplexed });

});

// ======= COMMANDS ======= //

router.post('/trimming', upload.none(), (req, res) => {

  // Variables
  let path_to_Fastq_Read1 = '';
  let path_to_Fastq_Read2 = '';
  let output_paired_Read1 = path.join(__dirname, '..', 'output', 'output_paired_Read1.fastq.gz');
  let output_unpaired_Read1 = path.join(__dirname, '..', 'output', 'output_unpaired_Read1.fastq.gz');
  let output_paired_Read2 = path.join(__dirname, '..', 'output', 'output_paired_Read2.fastq.gz');
  let output_unpaired_Read2 = path.join(__dirname, '..', 'output', 'output_unpaired_Read2.fastq.gz');

  // Choices
  let adapt = [];
  let read_length = [];
  let quality_score = [];

  if (req.adapter_trim) {
    adapt = [('ILLUMINACLIP:' + AdapterFile.path + ':2:30:2')];
  }

  if (req.read_length_trim) {
    read_length = ['MINLEN: ' + req.minlen];
  }

  if (!(req.quality_score_trim.length === 0)) {
    if (req.quality_score_trim === 'window') {
      quality_score = ['SLIDINGWINDOW:' + req.windowSize + ':' + req.quality];
    } else if (req.quality_score_trim === 'leading') {
      quality_score = ['LEADING:' + req.quality];
    } else if (req.quality_score_trim === 'trailing') {
      quality_score = ['TRAILING:' + req.quality];
    }
  }

  if (!req.adapter_trim && !req.read_length_trim && (req.quality_score_trim.length === 0)) {
    console.log("Error! At least one trim mode must be specified");
    return;
  }


  // Run Command
  // java -jar trimmomatic-0.39.jar PE <Read1 Fastq(.gz)> <Read2 Fastq(.gz)> <desired name of paired output file R1>
  // <desired name of unpaired output file R1> <desired name of paired output file R2> <desired name of unpaired output file R2>

  // -- Adapter Trim
  // ILLUMINACLIP:<adapter fasta file>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>
  //      Seed mismatches: max mismatch count that will still count as  matchduring initial seed alignment step (ie if set to 1 and there is a 1 bp difference between the adapter and the sequence it will count it as a match and identify it as a potential adapter; default to 2)
  //      Palindrome Clip threshhold: how accurate the match between two reads for palindrome read alignment (default to 30)
  //      Simple Clip threshhold: max mismatch count that will still count as  match (ie if set to 1 and there is a 1 bp difference between the adapter and the sequence it will count it as a match and cut default to 2)

  // -- Read Length Trim
  // MINLEN: <minimum length a read must be to not be cut>

  // -- Quality Score Trim
  // SLIDINGWINDOW:<windowSize>:<RequiredQuality>   ---- looks at average quality score of a window of specified length. If the average is below the specific quality score it will cut from 3' end until the average is above the threshold. Once it is above, the window will move one position to create a new window, starting the process again
  // LEADING:<quality>   ---- cuts reads until it reaches a read above the quality score specified starting from the 5' end
  // TRAILING:<quality>  ---- cuts reads until it reaches a read above the quality score specified starting from the 3' end

  // .push(...array) will append the contents of the array onto the trimArgs array, and will append nothing if empty, so that we can choose 1 or more trim modes
  const TrimCommand = 'java';
  const TrimArgs = ['-jar', path.join(__dirname, '..', 'bio_modules', 'trimmomatic-0.39.jar'), 'PE', path_to_Fastq_Read1, path_to_Fastq_Read2,
                  output_paired_Read1, output_unpaired_Read1, output_paired_Read2, output_unpaired_Read2].push(...adapt).push(...read_length).push(...quality_score);

  const runTrim = spawn(TrimCommand, TrimArgs);

  // Handle Response
  let outputData = '';

  runTrim.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runTrim.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runTrim.on('exit', (code) => {
    console.log(`runTrim process exited with code ${code}`);

    // Prepare the data you want to send back as JSON
    const responseData = {
      status: 'success',
      message: 'Trimming completed successfully',
      output: outputData, // Include any relevant output data
    };

    // Make a downloadable_content entry for each file that exists
    if (fs.existsSync(output_paired_Read1)) {
      downloadable_content.push({
        enabled: false,
        has_visual_component: false,
        label: 'Trimming Paired Read 1',
        content: 'output_paired_Read1.fastq.gz',
      });
    }

    if (fs.existsSync(output_unpaired_Read1)) {
      downloadable_content.push({
        enabled: false,
        has_visual_component: false,
        label: 'Trimming Unpaired Read 1',
        content: 'output_unpaired_Read1.fastq.gz',
      });
    }

    if (fs.existsSync(output_paired_Read2)) {
      downloadable_content.push({
        enabled: false,
        has_visual_component: false,
        label: 'Trimming Paired Read 2',
        content: 'output_paired_Read2.fastq.gz',
      });
    }

    if (fs.existsSync(output_unpaired_Read2)) {
      downloadable_content.push({
        enabled: false,
        has_visual_component: false,
        label: 'Trimming Unpaired Read 2',
        content: 'output_unpaired_Read2.fastq.gz',
      });
    }

    // Send the JSON response
    res.json(responseData);
  });

  

});

router.post('/alignment', upload.none(), (req, res) => {

  // Variables
  let reference_file_name = '';
  let Read1_FastQ_file = '';
  let Read2_FastQ_file = '';
  let name_of_output_file = path.join(__dirname, '..', 'output', 'AlignedSAM.sam');

  // Run Command
  let AlignmentCommand =' ';
  let AlignmentArgs = [];

  // --- BWA MEM
  // bwa mem <reference file name> <Read1 FastQ file (paired file from trimming)> <Read2 FastQ file (paired file from trimming)> > <name of output file>.sam
  if (req.type === 'bwa') {
    AlignmentCommand = path.join(__dirname, '..', 'bio_modules', 'bwa');
    AlignmentArgs = ['mem', reference_file_name, Read1_FastQ_file, Read2_FastQ_file, '>', name_of_output_file];
  }
  // --- BOWTIE
  // bowtie [options]* -x <ebwt> {-1 <m1> -2 <m2> | --12 <r> | --interleaved <i> | <s>} -S <output sam file>
  else if (req.type === 'bowtie') {
    AlignmentCommand = path.join(__dirname, '..', 'bio_modules', 'bowtie');
    AlignmentArgs = ['-x', req.ref, '--12', 'r', uploadedFiles[0].path + ', ' + uploadedFiles[1].path , '-S', name_of_output_file];
  } 
  // --- BOWTIE2
  // bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | b <bam>} -S [<sam>]
  else {
    AlignmentCommand = path.join(__dirname, '..', 'bio_modules', 'bowtie2');
    // NOTE: Can come back later and add the option to run this on a fastQ file if I do the m1 m2 option instead
    AlignmentArgs = ['-x', req.ref, 'b', uploadedFiles[0].path, '-S', name_of_output_file];
  } 
  
  const runAlignment = spawn(AlignmentCommand, AlignmentArgs);

  // Handle Response
  let outputData = '';

  runAlignment.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runAlignment.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runAlignment.on('exit', (code) => {
    console.log(`runAlignment process exited with code ${code}`);

    // Make a downloadable_content entry for it
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Aligned SAM File: ' + req.type,
      content: 'AlignedSAM.sam',
    });
  
    // Return the results as a JSON
    // Need to find and return the path of the output: Aligned SAM file in same directory
    res.json();
  });


});

router.post('/sam-to-bam', upload.none(), (req, res) => {
  // Variables
  let name_of_output_file_bam = path.join(__dirname, '..', 'output', 'ConvertedToBAM.bam');

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'AlignedSAM.sam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // samtools view -bS <path to sam file> > <name of output file>.bam
  const ConvertToBamCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const ConvertToBamArgs = ['view', '-bS', mainFilePath, '>', name_of_output_file_bam];

  const runConvertToBam = spawn(ConvertToBamCommand, ConvertToBamArgs);

  // Handle Response
  let outputData = '';

  runConvertToBam.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runConvertToBam.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runConvertToBam.on('exit', (code) => {
    console.log(`runConvertToBam process exited with code ${code}`);

    // Make a downloadable_content entry for it
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Converted BAM File',
      content: 'AlignedSAM.sam',
    });
  
    // Return the results as a JSON
    // Need to find and return the path of the output: Aligned SAM file in same directory
    res.json();
  });

});

router.post('/sort-bam-file', upload.none(), (req, res) => {
  // Variables
  let outputFile = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'ConvertedToBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    console.log(uploadedFilePath);
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // samtools sort <path to bam file> -o <name of output file>
  const SortBamCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const SortBamArgs = ['sort', mainFilePath, '-o', outputFile];

  const runSortBam = spawn(SortBamCommand, SortBamArgs);

  // Handle Response
  let outputData = '';

  runSortBam.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runSortBam.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runSortBam.on('exit', (code) => {
    console.log(`runSortBam process exited with code ${code}`);

    // Make a downloadable_content entry for it
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Sorted BAM',
      content: 'SortedBAM.bam'
    });

    // Return the results as a JSON
    // Need to find and return the path of the output: Aligned SAM file in same directory
    res.json({ });
  });

});

router.post('/index-bam-file', upload.none(), (req, res) => {

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // samtools index <path to sorted bam file>
  const IndexBamCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const IndexBamArgs = ['index', mainFilePath];

  const runIndexBam = spawn(IndexBamCommand, IndexBamArgs);

  // Handle Response
  let outputData = '';

  runIndexBam.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runIndexBam.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runIndexBam.on('exit', (code) => {
    console.log(`runIndexBam process exited with code ${code}`);

    let indexPath = '';
    // Check if it made a .bai or .csi index file
    if (fs.existsSync(path.join(__dirname, '..', 'output', 'SortedBAM.bam.bai'))) {
      indexPath = 'SortedBAM.bam.bai';
    } else {
      indexPath = 'SortedBAM.bam.csi';
    }

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'BAM Index',
      content: indexPath,
    });

    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/add-or-replace-read-groups', upload.none(), (req, res) => {
  // Specific Variables
  // let newReadGroupLine = '@RG\tID:sample1\tSM:sample1\tLB:library1\tPL:illumina';
  // let editMode = 'overwrite_all';

  let editMode = req.body.editMode; // overwrite_all or orphan_only depending on user selected mode
  let newReadGroupLine = req.body.newReadGroupLine; // This is the string that will be replaced/added to the file
  let outputFile = path.join(__dirname, '..', 'output', 'RG_bam.bam');

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Add or Replace Read Groups
  const RGCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const RGArgs = ['addreplacerg', '-r', newReadGroupLine, '-m', editMode, '-u', '-o', outputFile, mainFilePath];
  console.log(RGArgs);

  const runRG = spawn(RGCommand, RGArgs);

  let outputData = '';

  runRG.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runRG.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runRG.on('exit', (code) => {
    console.log(`Add or Replace Read Groups process exited with code ${code}`);

    // Make a downloadable_content entry for it
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Read Groups Added/Replaced BAM',
      content: 'RG_bam.bam',
    });
  
    res.json({ });
  });

});

router.post('/bam-index-stats', upload.none(), (req, res) => {
  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // samtools idxstats in.bam
  const BamIndexStatsCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const BamIndexStatsArgs = ['idxstats', mainFilePath];

  const runBamIndexStats = spawn(BamIndexStatsCommand, BamIndexStatsArgs);

  // Handle Response
  let outputData = '';

  runBamIndexStats.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runBamIndexStats.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runBamIndexStats.on('exit', (code) => {
    console.log(`runIndexBam process exited with code ${code}`);

    // Write outputData to a text file
    fs.writeFile(path.join(__dirname, '..', 'output', 'BAM_Index_Stats.txt'), outputData, (err) => {
      if (err) {
          console.error('Error writing file:', err);
          return;
      }
      console.log('stdout data written to BAM_Index_Stats.txt');
    });

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'BAM Index Stats',
      content: 'BAM_Index_Stats.txt',
    });

    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/alignment-summary', upload.none(), (req, res) => {
  // Variables
  let chosenRef = '';
  chosenRef = req.body.ref;
  let path_to_reference_fastA = path.join(__dirname, '..', 'ref_genomes', chosenRef, chosenRef + '.fna');

  let output_file_name_txt = path.join(__dirname, '..', 'output', 'AlignmentSummary.txt');

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // java -jar picard.jar CollectAlignmentSummaryMetrics R=<path to reference fasta> I=<path the sorted BAM file> O=<output file name.txt>
  const AlignmentDataArgs = ['-jar', path.join(__dirname, '..', 'bio_modules', 'picard.jar'), 'CollectAlignmentSummaryMetrics', 'R=', path_to_reference_fastA, 'I=', mainFilePath, 'O=', output_file_name_txt];

  const runAlignmentData = spawn('java', AlignmentDataArgs);

  // Handle Response
  let outputData = '';

  runAlignmentData.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runAlignmentData.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runAlignmentData.on('exit', (code) => {
    console.log(`runAlignmentData process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Alignment Summary',
      content: 'AlignmentSummary.txt',
    });

    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/gc-bias-summary', upload.none(), (req, res) => {
  // Variables
  let output_GC_bias_metrics_txt = path.join(__dirname, '..', 'output', 'GC_BIAS_Metrics.txt');
  let GC_bias_outputchart_pdf = path.join(__dirname, '..', 'output', 'GC_BIAS_OutputChart.pdf');
  let GC_Bias_summary_output_txt = path.join(__dirname, '..', 'output','GC_BIAS_SummaryOutput.txt');

  let chosenRef = '';
  chosenRef = req.body.ref;
  let path_to_reference_fastA = path.join(__dirname, '..', 'uploads', 'ref_genomes', chosenRef);

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // java -jar picard.jar CollectGcBiasMetrics -I <path to sorted BAM> -O <output GC bias metrics.txt> -CHART <GC bias ouputchart.pdf> -S <GC Bias summary output.txt> -R <reference fasta>
  const GCBiasDataCommand = 'java';
  const GCBiasDataArgs = ['-jar', path.join(__dirname, '..', 'bio_modules', 'picard.jar'), 'CollectGcBiasMetrics', '-I', mainFilePath, '-O', output_GC_bias_metrics_txt, '-CHART' + GC_bias_outputchart_pdf, '-S', GC_Bias_summary_output_txt + '-R', + path_to_reference_fastA];

  const runGCBiasData = spawn(GCBiasDataCommand, GCBiasDataArgs);

  // Handle Response
  let outputData = '';

  runGCBiasData.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runGCBiasData.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runGCBiasData.on('exit', (code) => {
    console.log(`runGCBiasData process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'GC Bias Metrics',
      content: 'GC_BIAS_Metrics.txt'
    });

    downloadable_content.push({
      enabled: false,
      has_visual_component: true,
      label: 'GC Bias Output Chart',
      content: 'GC_BIAS_OutputChart.pdf',
    });

    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'GC Bias Summary Output',
      content: 'GC_BIAS_SummaryOutput.txt',
    });


    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/insert-size-data', upload.none(), (req, res) => {
  // Variables
  let output_raw_data_txt = path.join(__dirname, '..', 'output', 'Insert_Size_RawData.txt');
  let output_histogram_name_pdf = path.join(__dirname, '..', 'output', 'Insert_Size_Histogram.pdf');

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // java -jar picard.jar CollectInsertSizeMetrics -I <path to sorted bam> -O <output raw data.txt> -H <output histogram name.pdf> M=.5
  const InsertSizeDataCommand = 'java';
  const InsertSizeDataArgs = ['-jar', path.join(__dirname, '..', 'bio_modules', 'picard.jar'), 'CollectInsertSizeMetrics', '-I', mainFilePath,
                          '-O', output_raw_data_txt, '-H' + output_histogram_name_pdf, 'M=.5'];

  const InsertSizeData = spawn(InsertSizeDataCommand, InsertSizeDataArgs);

  // Handle Response
  let outputData = '';

  InsertSizeData.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  InsertSizeData.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  InsertSizeData.on('exit', (code) => {
    console.log(`InsertSizeData process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Insert Size Raw Data',
      content: 'Insert_Size_RawData.txt',
    });
  
    downloadable_content.push({
      enabled: false,
      has_visual_component: true,
      label: 'Insert Size Histogram',
      content: 'Insert_Size_Histogram.pdf',
    });
  
  
    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/create-seq-dict', upload.none(), (req, res) => {
  // Variables
  let chosenRef = '';
  chosenRef = req.body.ref;
  let path_to_reference_fastA = path.join(__dirname, '..', 'uploads', 'ref_genomes', chosenRef);

  let outputFile = path.join(__dirname, '..', 'output', 'SeqDict.txt');

  // Run Command
  // samtools dict ref.fasta|ref.fasta.gz -o <output file name>
  const SeqDictCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const SeqDictArgs = ['idxstats', path_to_reference_fastA.path, '-o', outputFile];

  const runSeqDict = spawn(SeqDictCommand, SeqDictArgs);

  // Handle Response
  let outputData = '';

  runSeqDict.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runSeqDict.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runSeqDict.on('exit', (code) => {
    console.log(`SeqDict process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Sequence Dictionary',
      content: 'SeqDict.txt',
    });
  
    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/flag-stats', upload.none(), (req, res) => {
  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // samtools flagstat in.sam|in.bam|in.cram
  const FlagStatsCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const FlagStatsArgs = ['flagstat', mainFilePath];

  const runFlagStats = spawn(FlagStatsCommand, FlagStatsArgs);

  // Handle Response
  let outputData = '';

  runFlagStats.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runFlagStats.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runFlagStats.on('exit', (code) => {
    console.log(`FlagStats process exited with code ${code}`);

    // Write outputData to a text file
    fs.writeFile(path.join(__dirname, '..', 'output', 'Flag_Stats.txt'), outputData, (err) => {
      if (err) {
          console.error('Error writing file:', err);
          return;
      }
      console.log('stdout data written to Flag_Stats.txt');
    });

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Flag Stats',
      content: 'Flag_Stats.txt',
    });
  
    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/seq-depth', upload.none(), (req, res) => {
  let output_file_name = path.join(__dirname, '..', 'output', 'Seq_Depth.txt');

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // samtools depth -o FILE in.sam|in.bam
  const SeqDepthCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const SeqDepthArgs = ['depth', '-o', output_file_name, mainFilePath];

  const runSeqDepth = spawn(SeqDepthCommand, SeqDepthArgs);

  // Handle Response
  let outputData = '';

  runSeqDepth.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runSeqDepth.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runSeqDepth.on('exit', (code) => {
    console.log(`SeqDepth process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Sequence Depth Data',
      content: 'Seq_Depth.txt',
    });
  
    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/seq-coverage', upload.none(), (req, res) => {
  let output_file_name = path.join(__dirname, '..', 'output', 'Seq_Coverage.txt');
  let isVisual = req.body.visual ? '-m' : '';

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // samtools coverage -o <output file name> [-m] in.sam|in.bam
  const SeqCovCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const SeqCovArgs = ['coverage', '-o', output_file_name, isVisual, mainFilePath];

  const runSeqCov = spawn(SeqCovCommand, SeqCovArgs);

  // Handle Response
  let outputData = '';

  runSeqCov.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runSeqCov.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runSeqCov.on('exit', (code) => {
    console.log(`SeqCoverage process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: req.visual,
      label: req.visual ? 'Sequence Coverage Histogram' : 'Sequence Coverage Data',
      content: 'Seq_Coverage.txt',
    });

    // Return the results as a JSON
    res.json({ });
  });

});

router.post('/mark-remove-duplicates', upload.none(), (req, res) => {
  let output_bam_file_name = path.join(__dirname, '..', 'output', 'MarkedDuplicatesBAM.bam');
  let output_metrics_file_name = path.join(__dirname, '..', 'output', 'DuplicateMetrics.txt');
  let removeDupes = '';
  if (req.body.remove === 'yes') {
    removeDupes = '--REMOVE_DUPLICATES';
  } else if (req.body.remove === 'select') {
    removeDupes = '--REMOVE_SEQUENCING_DUPLICATES';
  }
  let removeDupHelper = (removeDupes === '') ? '' : 'true';

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    mainFilePath = uploadedFilePath;
    // We need the file, so check if a file was uploaded
    if (!fs.existsSync(uploadedFilePath)) {
      res.status(400).send('No file uploaded');
      return;
    }
  }

  // Run Command
  // NOTE: Bam file must be sorted
  // java -jar picard.jar MarkDuplicates -I <input bam file> -O <output bam with marked duplicates> -M <output metrics for marked duplicates> [--REMOVE_SEQUENCING_DUPLICATES | REMOVE_SEQUENCING_DUPLICATES  <true]
  const MarkRemDupCommand = 'java';
  const MarkRemDupArgs = ['-jar', path.join(__dirname, '..', 'bio_modules', 'picard.jar'), 'MarkDuplicates', '-I', mainFilePath, '-O', output_bam_file_name, '-M', output_metrics_file_name, removeDupes, removeDupHelper];

  const runMarkRemDup = spawn(MarkRemDupCommand, MarkRemDupArgs);

  // Handle Response
  let outputData = '';

  runMarkRemDup.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runMarkRemDup.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json();
  });

  runMarkRemDup.on('exit', (code) => {
    console.log(`Mark/Remove Duplicates process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Marked Duplicates BAM',
      content: 'MarkedDuplicatesBAM.bam',
    });
  
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Duplicate Metrics',
      content: 'DuplicateMetrics.txt',
    });
  
    // Return the results as a JSON
    res.json({ });
  });

});


/// SCREM

router.get('/testing', function (req, res) {
  res.render('testing');
});

router.post('/test-docker', upload.none(), async (req, res) => {
  try {
  // 'CollectAlignmentSummaryMetrics R=' + path.join(__dirname, '..', 'ref_genomes', 'Ecoli',  + 'Ecoli.fna') + ' I=' + path.join(__dirname, '..', 'output', 'SortedBAM.bam') + ' O=' + path.join(__dirname, '..', 'output', 'AlignmentSummary.txt');

  // 'java -jar ../picard/picard.jar MarkIlluminaAdapters -h'
    // let r_path = '../refs/Ecoli/Ecoli.fna';
    let i_path = '../output/SortedBAM.bam';
    let o_path = '../output/AlignmentSummary.txt';
    // const output = await runPicardCommand('java -jar ../picard/picard.jar CollectAlignmentSummaryMetrics -R ' + r_path + ' -I ' + i_path + ' -O ' + o_path);
    const output = await runPicardCommand('java -jar ../picard/picard.jar CollectAlignmentSummaryMetrics -I ' + i_path + ' -O ' + o_path);

    console.log('Output:', output);
    res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }
});


async function runPicardCommand(command) {
  const container = await docker.createContainer({
    Image: 'picard-tools-image', // Replace with your Docker image name
    Cmd: ['/bin/bash', '-c', command],
    AttachStdout: true,
    AttachStderr: true,
    Tty: false, // Important: Set Tty to false to prevent the container from running indefinitely
    HostConfig: {
      Mounts: [
        {
          Type: 'bind',
          Source: path.join(__dirname, '..', 'output'), // Replace with the absolute path to the output folder on your host machine
          Target: '/usr/output',
        },
        {
          Type: 'bind',
          Source: path.join(__dirname, '..', 'ref_genomes'), // Replace with the absolute path to the output folder on your host machine
          Target: '/usr/refs',
        },
      ],
    },
  });

  await container.start();

  const outputPromise = new Promise((resolve, reject) => {
    container.wait((err, data) => {
      if (err) {
        reject(err);
      } else {
        console.log('Raw Data:', data); // Log the raw data object
        resolve(data.toString());
      }
    });
  });

  let output;
  try {
    output = await outputPromise;
  } catch (err) {
    console.error('Error while waiting for container output:', err);
    throw new Error('Error waiting for container output');
  } finally {
    try {
      // Get container logs after command execution
      const logs = await container.logs({ stdout: true, stderr: true });
      console.log('Container Logs:', logs.toString());
    } catch (logErr) {
      console.error('Error getting container logs:', logErr);
    }

    container.remove(); // Always remove the container after use
  }

  return output;
}



//#endregion

module.exports = router;