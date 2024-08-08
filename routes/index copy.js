var express = require('express');
const multer = require('multer');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const bodyParser = require('body-parser');
const Docker = require('dockerode');

const { DNACategories, DNACommands, DNAParameters, RNACategories, RNACommands, RNAParameters } = require('../public/data/SiteData');
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

// Program Mode and Input Type
let mode = "";
let inputType = "";

// Stored Outputs
let infoSteps;          // An array that determines which file processing steps will be available in the timeline
let uploadedFilePath = "";   // The path to the original file that was uploaded or the directory, if it was a BCL folder
let uploadedPairedPath = "";
let uploadedAdapterPath = "";
let uploadedGenomeFile = "";
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
  {
      enabled: false,
      has_visual_component: true,
      label: "GC Bias Summary Graph",
      content: "../output/GC_Bias_Summary.png"
  },
  {
    enabled: false,
    has_visual_component: true,
    label: "Insert Size Summary Graph",
    content: "../output/Insert_Size_Graph.png"
  },
];

//#endregion

//#region  ... GET and POST Views ... 

router.get('/test', function(req, res) {
  res.render('test');
});

/* Home Page */
router.get('/', function (req, res, next) {
  // Render the 'home' view
  res.render('home', { toolbar_index: 1 });
});

/* File Info */
router.post('/file-information', function (req, res, next) {
  // Set the mode
  mode = req.body.mode;
  res.render('file-info', { toolbar_index: 2 });
});

/* DNA Goalposts */
router.post('/dna-goalposts', function (req, res, next) {
  // Store the uploaded file and any conversions
  inputType = req.body.inputType;
  infoSteps = JSON.parse(req.body.infoSteps);
  uploadedFilePath = req.body.uploadedFilePath;
  uploadedPairedPath = req.body.uploadedPairedPath;
  uploadedAdapterPath = req.body.uploadedAdapterPath;
  uploadedFileType = req.body.uploadedFileType;
  fastQConversion = req.body.fastQConversion;
  fastQCResults = req.body.fastQCResults;
  uploadedGenomeFile = req.body.uploadedGenomeFile;

  let Categories;
  let Commands;

  if (mode == 'DNA') {
    Categories = DNACategories;
    Commands = DNACommands;
  } else {
    Categories = RNACategories;
    Commands = RNACommands;
  }

  res.render('dna-goalposts', { toolbar_index: 3, Categories, Commands, infoSteps });
});

/* DNA Pipeline */
router.post('/dna-pipeline', function (req, res, next) {
  let temp = JSON.parse(decodeURIComponent(req.body.Commands || '[]'));

  if (mode == 'DNA') {
    // Clear the existing commands
    DNACommands.clear();

    // Convert the plain object to a Map and then iterate over it
    const dnaCommandsMap = new Map(temp);
    for (const [key, value] of dnaCommandsMap) {
      DNACommands.set(key, value);
    }

    Categories = DNACategories;
    Commands = DNACommands;
    Parameters = DNAParameters;
  } else {
    // Clear the existing commands
    RNACommands.clear();

    // Convert the plain object to a Map and then iterate over it
    const rnaCommandsMap = new Map(temp);
    for (const [key, value] of rnaCommandsMap) {
      RNACommands.set(key, value);
    }

    Categories = RNACategories;
    Commands = RNACommands;
    Parameters = RNAParameters;
  }

  res.render('dna-pipeline', { toolbar_index: 4, Categories, Commands, Parameters, infoSteps });
});

/* Running Page */
router.post('/running', function (req, res, next) {
  let tempCOM = JSON.parse(decodeURIComponent(req.body.Commands || '[]'));
  let tempPAR = JSON.parse(decodeURIComponent(req.body.Parameters || '[]'));

  if (mode == 'DNA') {
    // Clear the existing commands
    DNACommands.clear();
    DNAParameters.clear();

    // Convert the plain object to a Map and then iterate over it
    const dnaCommandsMap = new Map(tempCOM);
    for (const [key, value] of dnaCommandsMap) {
      DNACommands.set(key, value);
    }
    const dnaParametersMap = new Map(tempPAR);
    for (const [key, value] of dnaParametersMap) {
      DNAParameters.set(key, value);
    }

    Commands = DNACommands;
    Parameters = DNAParameters;
  } else {
    // Clear the existing commands
    RNACommands.clear();
    RNAParameters.clear();

    // Convert the plain object to a Map and then iterate over it
    const rnaCommandsMap = new Map(tempCOM);
    for (const [key, value] of rnaCommandsMap) {
      RNACommands.set(key, value);
    }
    const rnaParametersMap = new Map(tempPAR);
    for (const [key, value] of rnaParametersMap) {
      RNAParameters.set(key, value);
    }
    
    Commands = RNACommands;
    Parameters = RNAParameters;
  }

  res.render('running', { toolbar_index: 5, Commands, Parameters, uploadedFilePath });
});

/* Output Page */
router.post('/output', function (req, res, next) {
  // Render the 'output' view
  res.render('output', { downloadable_content, toolbar_index: 5 });
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
    return res.status(200).send({ storedPath: [uploadedFilePath] });
  }

  // If it's exactly two files, send the new paths to those
  if (uploadedFiles.length === 1) {
    uploadedFilePath = path.join(__dirname, '..', 'uploads', uploadedFiles[0].filename);
    uploadedFilePath2 = path.join(__dirname, '..', 'uploads', uploadedFiles[1].filename);
    return res.status(200).send({ storedPath: [uploadedFilePath, uploadedFilePath2] });
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
  return res.status(200).send({ storedPath: [uploadedFilePath] });
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
router.post('/bcl2fastq', upload.array('files'), async (req, res) => {

  // =============== Handle the BCL File Folder ===============
  const uploadedFiles = req.files;

  // Check if the files were uploaded
  if (!uploadedFiles || uploadedFiles.length === 0) {
    return res.status(400).send('No files were uploaded.');
  }

  // Create an empty directory to hold them in
  const BCLdirectory = path.join(__dirname, 'uploads', 'BCLFiles-' + Date.now());
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
  const OutputDirectory = path.join(__dirname, 'output', 'bcl2fastq_Output');
  if (!fs.existsSync(OutputDirectory)) {
    fs.mkdirSync(OutputDirectory);
  }

  // Run Command
  try {
    const output = await runBCLCommand(BCLdirectory);

    console.log('Output:', output);
    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'BCL2FastQ Output',
      content: 'bcl2fastq_Output',
    });

    // Return the results as a JSON
    res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }

});

// ======= COMMANDS ======= //

router.post('/fast5-to-fastq', upload.array('files'), async (req, res) => {

  // =============== Handle the fast5 File Folder ===============
  const uploadedFiles = req.files;

  // Check if the files were uploaded
  if (!uploadedFiles || uploadedFiles.length === 0) {
    return res.status(400).send('No files were uploaded.');
  }

  // Create an empty directory to hold them in
  const folderName = 'Fast5Files-' + Date.now();
  const Fast5directory = path.join(__dirname, 'uploads', folderName);
  if (!fs.existsSync(Fast5directory)) {
    fs.mkdirSync(Fast5directory);
  }

  // Move each uploaded file to the empty directory
  uploadedFiles.forEach(file => {
    const destinationPath = path.join(Fast5directory, file.originalname);
    fs.renameSync(file.path, destinationPath);
  });

  // =============== Run BCL To FastQ ===============

  // Create an empty directory to hold the output
  const OutputFastQ = path.join(__dirname, 'output', 'FastQFromFast5.fastq');

  // Run Command
  try {
    const output = await runFast5ToFastQ(folderName);

    console.log('Output:', output);
    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'FastQ From Fast5',
      content: OutputFastQ,
    });

    // Return the results as a JSON
    res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }

});

// Run FastQC
// Expects a single FastQ File
router.post('/fastqc', upload.none(), (req, res) => {
  // Check if a file was uploaded
  if (!fs.existsSync(uploadedFilePath)) {
    res.status(400).send('No file uploaded');
    return;
  }

  // Run FastQC
  const FastQCCommand = path.join(__dirname, '..', 'bio_modules', 'FastQC.app', 'Contents', 'MacOS', 'fastqc');
  const FastQArgs = [uploadedFilePath];

  const runFastQC = spawn(FastQCCommand, FastQArgs);

  let outputData = '';

  runFastQC.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runFastQC.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json({});
    return;
  });

  runFastQC.on('exit', (code) => {
    console.log(`FastQC process exited with code ${code}`);
  });

  // Return the results as a JSON
  const trimmed = false; // set this based on the result of the QC
  const demultiplexed = false;

  res.json({ trimmed, demultiplexed });

});


router.post('/fasta-to-fastq', upload.none(), (req, res) => {

  let convertedOutputFQ = path.join(__dirname, '..', 'output', 'convertedFastQ.fq')

  // We need the file, so check if a file was uploaded
  if (!fs.existsSync(uploadedFilePath)) {
    res.status(400).send('No file uploaded');
    return;
  }

  // Run Command
  // ./seqtk seq -F '<char to replace quality scores>' in.fa > out.fq
  const AtoQCommand = path.join(__dirname, '..', 'bio_modules', 'seqtk');
  const AtoQArgs = ['seq', '-F', '#', uploadedFilePath, '>', convertedOutputFQ];

  const runAtoQ = spawn(AtoQCommand, AtoQArgs);

  // Handle Response
  let outputData = '';

  runAtoQ.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runAtoQ.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json({});
    return;
  });

  runAtoQ.on('exit', (code) => {
    console.log(`runAtoQ process exited with code ${code}`);

    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'FastQ from FastA',
      content: 'convertedFastQ.fq',
    });

    // Return the results as a JSON
    res.json({});
    return;
  });

});

router.post('/trimming', upload.none(), (req, res) => {
  console.log('ummmm');


  // Variables
  let path_to_Single_FastQ = uploadedFilePath;
  let single_trim_output = path.join(__dirname, '..', 'output', 'single_trim_output.fastq.gz');
  let path_to_Fastq_Read1 = uploadedFilePath;
  let path_to_Fastq_Read2 = uploadedPairedPath;
  let output_paired_Read1 = path.join(__dirname, '..', 'output', 'output_paired_Read1.fastq.gz');
  let output_unpaired_Read1 = path.join(__dirname, '..', 'output', 'output_unpaired_Read1.fastq.gz');
  let output_paired_Read2 = path.join(__dirname, '..', 'output', 'output_paired_Read2.fastq.gz');
  let output_unpaired_Read2 = path.join(__dirname, '..', 'output', 'output_unpaired_Read2.fastq.gz');
  let adapter_filepath = '';

  // Choices
  let adapt = [];
  let read_length = [];
  let window = [];
  let leading = [];
  let trailing = [];

  console.log('before first check');


  if (req.body.adapter_trim) {
    adapt.push(`ILLUMINACLIP:${adapter_filepath}:2:30:2`); // note: needs adapter file to exist
  }
  console.log('after first check');


  if (req.body.read_length_trim) {
    read_length.push(`MINLEN:${req.body.minlen}`);
  }

  console.log('after second check');

  if (req.body.window) {
    window.push(`SLIDINGWINDOW:${req.body.window[0]}:${req.body.window[1]}`);
  }
  if (req.body.leading) {
    leading.push(`LEADING:${req.body.leading}`);
  }
  if (req.body.trailing) {
    trailing.push(`TRAILING:${req.body.trailing}`);
  }

  console.log('after third check');


  if (!req.body.adapter_trim && !req.body.read_length_trim && !req.body.quality_score_trim) {
    console.log("Error! At least one trim mode must be specified");
    return res.status(400).json({
      status: 'error',
      message: 'At least one trim mode must be specified'
    });
  }

  console.log('okayyy');

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

  // Paired End or Single End depending what was given
  let TrimArgs = [];
  if (uploadedPairedPath) {
    TrimArgs = [
      '-jar', path.join(__dirname, '..', 'bio_modules', 'Trimmomatic-0.39', 'trimmomatic-0.39.jar'),
      'PE', path_to_Fastq_Read1, path_to_Fastq_Read2,
      output_paired_Read1, output_unpaired_Read1, output_paired_Read2, output_unpaired_Read2,
      ...adapt, ...read_length, ...window, ...leading, ...trailing
    ];
  } else {
    TrimArgs = [
      '-jar', path.join(__dirname, '..', 'bio_modules', 'Trimmomatic-0.39', 'trimmomatic-0.39.jar'),
      'SE', path_to_Single_FastQ, single_trim_output,
      ...adapt, ...read_length, ...window, ...leading, ...trailing
    ];
  }


  const runTrim = spawn(TrimCommand, TrimArgs);

  // Handle Response
  let outputData = '';
  let errorData = '';

  runTrim.stdout.on('data', (data) => {
    outputData += data.toString();
    console.log(`stdout: ${data}`);
  });

  runTrim.stderr.on('data', (data) => {
    errorData += data.toString();
    console.error(`stderr: ${data}`);
  });

  runTrim.on('exit', (code) => {
    console.log(`runTrim process exited with code ${code}`);

    if (code !== 0) {
      res.status(500).json({
        status: 'error',
        message: 'Trimming process failed',
        error: errorData,
      });
      return;
    }

    // Prepare the data you want to send back as JSON
    const responseData = {
      status: 'success',
      message: 'Trimming completed successfully',
      output: outputData,
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

  // Request should have: reference file option, choice of mode

  // ref finder
  let chosenRef = '';
  chosenRef = req.body.ref_genome;
  if (chosenRef = 'Human') { chosenRef = 'hgch38_index'; }
  if (chosenRef = 'House_Mouse') { chosenRef = 'House_Mouse_GRCm39'; }
  let path_to_reference_fastA = path.join(__dirname, '..', 'ref_genomes', chosenRef, chosenRef + '.fna');

  // Variables
  let name_of_output_file = path.join(__dirname, '..', 'output', 'AlignedSAM.sam');

  // Run Command
  let AlignmentCommand = '';
  let AlignmentArgs = [];

  // --- BWA MEM
  // bwa mem <reference file name> <Read1 FastQ file (paired file from trimming)> <Read2 FastQ file (paired file from trimming)> > <name of output file>.sam
  if (req.body.type === 'bwa') {
    AlignmentCommand = path.join(__dirname, '..', 'bio_modules', 'bwa');
    AlignmentArgs = ['mem', path.join(__dirname, '..', 'ref_genomes', chosenRef), uploadedFilePath, uploadedPairedPath, '>', name_of_output_file];
  }
  // --- BOWTIE
  // bowtie [options]* -x <ebwt> {-1 <m1> -2 <m2> | --12 <r> | --interleaved <i> | <s>} -S <output sam file>
  else if (req.body.type === 'bowtie') {
    AlignmentCommand = path.join(__dirname, '..', 'bio_modules', 'bowtie');
    AlignmentArgs = ['-x', path_to_reference_fastA, '--12', 'r', uploadedFilePath + ', ' + uploadedPairedPath, '-S', name_of_output_file];
  }
  // --- BOWTIE2
  // bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | b <bam>} -S [<sam>]
  else {
    AlignmentCommand = path.join(__dirname, '..', 'bio_modules', 'bowtie2');
    AlignmentArgs = ['-x', path.join(__dirname, '..', 'ref_genomes', chosenRef), '-U', 'r', uploadedFilePath, + ', ' + uploadedPairedPath, '-S', name_of_output_file];
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
    res.json({});
    return;
  });

  runAlignment.on('exit', (code) => {
    console.log(`runAlignment process exited with code ${code}`);

    // Make a downloadable_content entry for it
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Aligned SAM File: ' + req.body.type,
      content: 'AlignedSAM.sam',
    });

    // Return the results as a JSON
    // Need to find and return the path of the output: Aligned SAM file in same directory
    res.json({});
    return;
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
  // samtools view -b <path to sam file> > <name of output file>.bam
  const ConvertToBamCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const ConvertToBamArgs = ['view', '-b', mainFilePath, '>', name_of_output_file_bam];

  const runConvertToBam = spawn(ConvertToBamCommand, ConvertToBamArgs);

  // Handle Response
  let outputData = '';

  runConvertToBam.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runConvertToBam.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json({});
    return;
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
    res.json({});
    return;
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
    res.json({});
    return;
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
    res.json({});
    return;
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
    res.json({});
    return;
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
    res.json({});
    return;
  });

});

router.post('/add-or-replace-read-groups', upload.none(), (req, res) => {
  // Specific Variables
  // let newReadGroupLine = '@RG\tID:sample1\tSM:sample1\tLB:library1\tPL:illumina';
  // let editMode = 'overwrite_all';

  let editMode = req.body.editMode ? 'overwrite_all' : 'orphan_only'; // overwrite_all or orphan_only depending on user selected mode
  let newReadGroupLine = req.body.newReadGroupLine; // This is the string that will be replaced/added to the file
  let outputFile = path.join(__dirname, '..', 'output', 'RG_bam.bam');

  // Main File Check
  let mainFilePath = path.join(__dirname, '..', 'output', 'SortedBAM.bam');

  if (!fs.existsSync(mainFilePath)) {
    console.log("where is file");
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
    res.json({});
    return;
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

    res.json({});
    return;
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
    res.json({});
    return;
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
    res.json({});
    return;
  });

});

router.post('/alignment-summary', upload.none(), async (req, res) => {
  // Variables
  let chosenRef = '';
  chosenRef = req.body.ref_genome;
  if (chosenRef = 'Human') { chosenRef = 'hgch38_index'; }
  if (chosenRef = 'House_Mouse') { chosenRef = 'House_Mouse_GRCm39'; }

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
  try {
    // let r_path = '../refs/Human/hgch38_index.fna';
    let r_path = '../refs/' + chosenRef + '/' + chosenRef + '.fna';
    let i_path = '../output/SortedBAM.bam'; // CBTT this negates my main file check
    let o_path = '../output/AlignmentSummary.txt';
    const output = await runPicardCommand('java -jar ../picard/picard.jar CollectAlignmentSummaryMetrics -I ' + i_path + ' -O ' + o_path + ' -R ' + r_path);

    console.log('Output:', output);
    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Alignment Summary',
      content: 'AlignmentSummary.txt',
    });

    // Return the results as a JSON
    res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }

  // runAlignmentData.on('exit', (code) => {
  //   console.log(`runAlignmentData process exited with code ${code}`);
  // });

});

router.post('/gc-bias-summary', upload.none(), async (req, res) => {
  // Variables
  let output_GC_bias_metrics_txt = path.join('..', 'output', 'GC_BIAS_Metrics.txt');
  let GC_bias_outputchart_pdf = path.join('..', 'output', 'GC_BIAS_OutputChart.pdf');
  let GC_Bias_summary_output_txt = path.join('..', 'output', 'GC_BIAS_SummaryOutput.txt');

  let chosenRef = '';
  chosenRef = req.body.ref_genome;
  if (chosenRef = 'Human') { chosenRef = 'hgch38_index'; }
  if (chosenRef = 'House_Mouse') { chosenRef = 'House_Mouse_GRCm39'; }

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
  try {
    let r_path = '../refs/' + chosenRef + '/' + chosenRef + '.fna';
    let i_path = '../output/SortedBAM.bam'; // CBTT this negates my main file check
    let command = 'java -jar ../picard/picard.jar CollectGcBiasMetrics -I ' + i_path + ' -O ' + output_GC_bias_metrics_txt + ' -CHART ' + GC_bias_outputchart_pdf + ' -S ' + GC_Bias_summary_output_txt + ' -R ' + r_path;
    if (!req.body.isVisual) {
      command = 'java -jar ../picard/picard.jar CollectGcBiasMetrics -I ' + i_path + ' -O ' + output_GC_bias_metrics_txt + ' -S ' + GC_Bias_summary_output_txt + ' -R ' + r_path;
    }  
    const output = await runPicardCommand(command);

    console.log('Output:', output);
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
    res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }

});

router.post('/insert-size-summary', upload.none(), async (req, res) => {
  // Variables
  let output_raw_data_txt = path.join('..', 'output', 'Insert_Size_RawData.txt');
  let output_histogram_name_pdf = path.join('..', 'output', 'Insert_Size_Histogram.pdf');

  // Main File Check
  let mainFilePath = path.join('..', 'output', 'SortedBAM.bam');

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
  try {
    // let r_path = '../refs/Ecoli/Ecoli.fna';
    let i_path = '../output/SortedBAM.bam'; // CBTT this negates my main file check
    const output = await runPicardCommand('java -jar ../picard/picard.jar CollectInsertSizeMetrics -I ' + i_path + ' -O ' + output_raw_data_txt + ' -H ' + output_histogram_name_pdf + ' -M 0.5');

    console.log('Output:', output);
    // Make a downloadable_content entry for the results
    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Insert Size Raw Data',
      content: 'Insert_Size_RawData.txt',
    });

    downloadable_content.push({
      enabled: false,
      has_visual_component: false,
      label: 'Insert Size Histogram',
      content: 'Insert_Size_Histogram.pdf',
    });

    // Return the results as a JSON
    res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }

});

router.post('/create-sequence-dictionary', upload.none(), (req, res) => {
  // Variables
  let chosenRef = '';
  chosenRef = req.body.ref_genome;
  if (chosenRef = 'Human') { chosenRef = 'hgch38_index'; }
  if (chosenRef = 'House_Mouse') { chosenRef = 'House_Mouse_GRCm39'; }
  let path_to_reference_fastA = path.join(__dirname, '..', 'ref_genomes', chosenRef, chosenRef + '.fna');

  let outputFile = path.join(__dirname, '..', 'output', 'SeqDict.txt');

  // Run Command
  // samtools dict ref.fasta|ref.fasta.gz -o <output file name>
  const SeqDictCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const SeqDictArgs = ['dict', path_to_reference_fastA, '-o', outputFile];

  const runSeqDict = spawn(SeqDictCommand, SeqDictArgs);

  // Handle Response
  let outputData = '';

  runSeqDict.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runSeqDict.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json({});
    return;
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
    res.json({});
    return;
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
    res.json({});
    return;
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
    res.json({});
    return;
  });

});

router.post('/sequence-depth', upload.none(), (req, res) => {
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
    res.json({});
    return;
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
    res.json({});
    return;
  });

});

router.post('/sequence-coverage', upload.none(), (req, res) => {
  let output_file_name = path.join(__dirname, '..', 'output', 'Seq_Coverage.txt');
  // let isVisual = req.body.visual ? '-m' : '';

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
  const SeqCovArgs = req.body.visual ? ['coverage', '-o', output_file_name, '-m', mainFilePath] : ['coverage', '-o', output_file_name, mainFilePath];

  const runSeqCov = spawn(SeqCovCommand, SeqCovArgs);

  // Handle Response
  let outputData = '';

  runSeqCov.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  runSeqCov.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
    res.json({});
    return;
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
    res.json({});
    return;
  });

});

router.post('/mark-or-remove-duplicates', upload.none(), async (req, res) => {
  let output_bam_file_name = path.join('..', 'output', 'MarkedDuplicatesBAM.bam');
  let output_metrics_file_name = path.join('..', 'output', 'DuplicateMetrics.txt');
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
  // java -jar picard.jar MarkDuplicates -I <input bam file> -O <output bam with marked duplicates> -M <output metrics for marked duplicates> [--REMOVE_SEQUENCING_DUPLICATES | REMOVE_SEQUENCING_DUPLICATES  <true]
  try {
    // let r_path = '../refs/Ecoli/Ecoli.fna';
    let i_path = '../output/SortedBAM.bam'; // CBTT this negates my main file check
    let command = 'java -jar ../picard/picard.jar MarkDuplicates -I ' + i_path + ' -O ' + output_bam_file_name + ' -M ' + output_metrics_file_name;
    command += removeDupes + removeDupHelper;
    const output = await runPicardCommand(command);

    console.log('Output:', output);
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
    res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }

});


// async function runCommand(checkFile, environment, command, downloadables) {
//   // Main File Check (We need a file. Check if we have it. If not, use the user uploaded file.)
//   let mainFilePath = path.join(__dirname, '..', 'output', checkFile);
//   if (!fs.existsSync(mainFilePath)) {
//     mainFilePath = uploadedFilePath;
//     // We need the file, so check if a file was uploaded
//     if (!fs.existsSync(uploadedFilePath)) {
//       res.status(400).send('No file uploaded');
//       return;
//     }
//   }

//   // Run the command in a Docker environment
//   try {
//     const output = await runDockerCommand(environment, command);

//     // output will hold the output from stdout and stderr
//     // NOTE: how to deal with things that read stdout into a file
//     console.log('Output:', output);
    
//     // Make downloadable content entries for the results
//     downloadables.forEach(downloadable => {
//       // check if the result was successfully created
//       if (fs.existsSync(path.join(__dirname, '..', 'output', downloadable.content))) {
//         // if so, push it to the downloadable list
//         downloadable_content.push(downloadable);
//       }
//     })

//     // Return the results as a JSON
//     res.status(200).send(JSON.stringify(output)); // Convert output to JSON string
//   } catch (error) {
//     console.error('Error:', error.message);
//     res.status(500).send('Error executing Docker command');
//   }
// }

// async function runDockerCommand(image, command) { 
//   const container = await docker.createContainer({
//     Image: image,
//     Cmd: command,
//     AttachStdout: true,
//     AttachStderr: true,
//     Tty: false, // Important: Set Tty to false to prevent the container from running indefinitely
//     HostConfig: {
//       Mounts: [
//         {
//           Type: 'bind',
//           Source: path.join(__dirname, '..', 'output'),
//           Target: '/usr/output',
//         },
//         {
//           Type: 'bind',
//           Source: path.join(__dirname, '..', 'ref_genomes'),
//           Target: '/usr/refs',
//         },
//         {
//           Type: 'bind',
//           Source: path.join(__dirname, '..', 'uploads'),
//           Target: '/usr/uploads',
//         },
//       ],
//     },
//   });

//   await container.start();

//   const outputPromise = new Promise((resolve, reject) => {
//     container.wait((err, data) => {
//       if (err) {
//         reject(err);
//       } else {
//         console.log('Raw Data:', data); // Log the raw data object
//         resolve(data.toString());
//       }
//     });
//   });

//   let output;
//   try {
//     output = await outputPromise;
//   } catch (err) {
//     console.error('Error while waiting for container output:', err);
//     throw new Error('Error waiting for container output');
//   } finally {
//     try {
//       // Get container logs after command execution
//       const logs = await container.logs({ stdout: true, stderr: true });
//       console.log('Container Logs:', logs.toString());
//     } catch (logErr) {
//       console.error('Error getting container logs:', logErr);
//     }

//     container.remove(); // Always remove the container after use
//   }

//   return output;
// }


/// SCREM

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

async function runBCLCommand(inputDirectory) {
  const container = await docker.createContainer({
    Image: 'illumina-bcl2fastq', // Replace with your Docker image name
    Cmd: ['/bin/bash', '-c', 'bcl2fastq --runfolder-dir /usr/uploads --output-dir /usr/output --processing-threads 1 --no-lane-splitting'], // the options are for RAM
    AttachStdout: true,
    AttachStderr: true,
    Tty: false, // Important: Set Tty to false to prevent the container from running indefinitely
    HostConfig: {
      Mounts: [
        {
          Type: 'bind',
          Source: path.join(__dirname, 'output', 'bcl2fastq_Output'), // Replace with the absolute path to the output folder on your host machine
          Target: '/usr/output',
        },
        {
          Type: 'bind',
          Source: inputDirectory, // Replace with the absolute path to the output folder on your host machine
          Target: '/usr/uploads',
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

async function runFast5ToFastQ(inputDirectory) {
  const container = await docker.createContainer({
    Image: 'fast5-to-fastq', // Replace with your Docker image name
    Cmd: ['/usr/uploads/' + folderName, '>', 'usr/output/FastQFromFast5.fastq'],
    AttachStdout: true,
    AttachStderr: true,
    Tty: false, // Important: Set Tty to false to prevent the container from running indefinitely
    HostConfig: {
      Mounts: [
        {
          Type: 'bind',
          Source: path.join(__dirname, 'output'), // Replace with the absolute path to the output folder on your host machine
          Target: '/usr/output',
        },
        {
          Type: 'bind',
          Source: inputDirectory, // Replace with the absolute path to the output folder on your host machine
          Target: '/usr/uploads',
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

//#region ... Execute Commands ...

/** Example usage:
* const pairs = [["name", "Jack"], ["age", "twenty"]];
* const inputString = "fillForm <name> -addFile -useWarnings <age>";
* const result = substitutePlaceholders(pairs, inputString);
* console.log(result); // Output: "fillForm Jack -addFile -useWarnings twenty"
**/
function substitutePlaceholders(pairs, inputString) {
  // Loop through each pair in the array
  pairs.forEach(pair => {
    const placeholder = `<${pair[0]}>`; // Create the placeholder format
    const value = pair[1]; // Get the value to substitute

    // Replace all occurrences of the placeholder in the input string
    inputString = inputString.split(placeholder).join(value);
  });

  return inputString;
}

router.post('/run-command', upload.none(), async (req, res) => {
  // Get the info needed
  // let Executable = DNAExecutables[req.body.command];
  let Executable = JSON.parse(req.body.executable);

  // -- Main File Check --
  // Find the output file we need to run on from the database
  if (Executable.checkFile != '') {
    let mainFilePath = path.join(__dirname, '..', 'output', Executable.checkFile);
    // If the file we think we need does not exist, use the user uploaded file
    if (!fs.existsSync(mainFilePath)) {
      mainFilePath = uploadedFilePath;
      // If there is no user uploaded file, we have an error
      if (!fs.existsSync(uploadedFilePath)) {
        res.status(400).send('No file uploaded');
        return;
      }
    }
  }

  // -- Build the Command --
  let allVars = (Executable.checkFile == '') ? Executable.variables : [['mainFile', mainFilePath]].push(Executable.variables);
  let commandString = substitutePlaceholders(allVars, Executable.command);

  // -- Run the Command in a Docker Environment --
  try {
    const output = await runDockerCommand(Executable.environment, commandString);
    // NOTE: how to deal with things that read stdout into a file?
    console.log('Output:', output);

    // Make downloadable content entries for the results
    Executable.downloadables.forEach(downloadable => {
      // Check if the result was successfully created by the command
      if (fs.existsSync(path.join(__dirname, '..', 'output', downloadable.path))) {
        // if so, push it to the downloadable list
        downloadable_content.push(downloadable);
      }
    })

    // Return the results as a JSON
    res.status(200).send(JSON.stringify(output));
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send('Error executing Docker command');
  }
});

async function runDockerCommand(image, command) { 
  const container = await docker.createContainer({
    Image: image,
    Cmd: command.split(' '),
    AttachStdout: true,
    AttachStderr: true,
    Tty: false, // Important: Set Tty to false to prevent the container from running indefinitely
    HostConfig: {
      Mounts: [
        {
          Type: 'bind',
          Source: path.join(__dirname, '..', 'output'),
          Target: '/usr/output',
        },
        {
          Type: 'bind',
          Source: path.join(__dirname, '..', 'ref_genomes'),
          Target: '/usr/refs',
        },
        {
          Type: 'bind',
          Source: path.join(__dirname, '..', 'uploads'),
          Target: '/usr/uploads',
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



//#endregion

module.exports = router;