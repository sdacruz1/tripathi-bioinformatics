var express = require('express');
const multer = require('multer');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
var router = express.Router();

// Use express.urlencoded middleware
router.use(express.urlencoded({ extended: true }));

// Serve uploaded images statically
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
let uploadedFilePath;   // The path to the original file that was uploaded or the directory, if it was a BCL folder
let uploadedFileType;   // A string representing the original uploaded file type
let fastQConversion;    // The path to the fastQ file result of any conversions (BCL, FastA, Fast5)
let fastQCResults;      // The path to the directory containing the results of running FastQC on the uploaded file
let BAMFile;            // The path to the BAM file result of any conversions
let Outputs;            // An object array holding the result of the steps in the actual timeline

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
  infoSteps = req.body.infoSteps;
  uploadedFilePath = req.body.uploadedFilePath;
  uploadedFileType = req.body.uploadedFileType;
  // fastQConversion = req.body.fastQConversion;
  // fastQCResults = req.body.fastQCResults;

  res.render('dna-goalposts', { toolbar_index: 3, categories });
});

/* DNA Pipeline */
router.post('/dna-pipeline', function (req, res, next) {
  categories = JSON.parse(decodeURIComponent(req.body.categories || '[]'));

  res.render('dna-pipeline', { toolbar_index: 4, categories });
});

/* Running Page */
router.post('/running', function (req, res, next) {
  const categories = JSON.parse(decodeURIComponent(req.body.categories || '[]'));

  res.render('running', { toolbar_index: 5, categories });
});

//#endregion

//#region  ... Utility ...

// Run BCL2FastQ
// Expects the contents of a BCL File
router.post('/run-bcl2fastq' , upload.array('files'), (req, res) => {

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
  const OutputDirectory = './uploads/output';
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
router.post('/run-fastqc', upload.single('file'), (req, res) => {
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
  });

  runFastQC.on('exit', (code) => {
    console.log(`FastQC process exited with code ${code}`);
  });

  // Return the results as a JSON
  const trimmed = false; // set this based on the result of the QC
  const demultiplexed = false;

  res.json({ trimmed, demultiplexed });

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

//#endregion

module.exports = router;
