var express = require('express');
const multer = require('multer');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
var router = express.Router();

// Use express.urlencoded middleware
router.use(express.urlencoded({ extended: true }));

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

const toolbar = ['Select a Workflow', 'Upload Files', 'Select Goalposts', 'Edit Timeline', 'Run Timeline'];

/* GET home page. */
router.get('/', function (req, res, next) {
  // Render the 'home' view with the myArray
  res.render('home', { toolbar, index: 1 });
});

router.get('/file-info', function (req, res, next) {
  res.render('file-info', {toolbar, index: 2});
});

router.get('/dna-goalposts', function (req, res, next) {
  // Categories
  const files = ["File Processing",
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

  const stats = ["Statistics",
    [
      ["Add or Replace Read Groups", false, [], []],
      ["Bam Index Stats", false, [], []],
    ]];

  const summary = ["Summary",
    [
      ["Alignment Summary", false, [], []],
      ["GC Bias Summary", false, [], []],
      ["Insert Size Summary", false, [], []],
    ]];

  const graphs = ["Graphs",
    [
      ["Alignment Graph", false, [], []],
      ["GC Bias Graphs", false, [], []],
      ["Insert Size Graphs", false, [], []],
    ]];

  const metrics = ["Metrics",
    [
      ["Quality Yield Metrics", false, [], []],
      ["Whole Genome Sequencing Metrics", false, [], []],
      ["Targeted PCR Metrics", false, [], []],
    ]];

  const cleanup = ["Cleanup and Sequencing",
    [
      ["Create Sequence Dictionary", false, [], []],
      ["Mark Duplicates", false, [], []],
      ["Sort Bam File", false, [], []],
      ["Flag Stats", false, [], []],
      ["Sequencing Depth", false, [], []],
      ["Sequencing Coverage", false, [], []],
    ]];

  const categories = [files, stats, summary, graphs, metrics, cleanup];

  res.render('dna-goalposts', { toolbar, index: 3, categories });
});

router.post('/dna-pipeline', function (req, res, next) {
  const categoriesString = req.body.categories || '[]';
  const categories = JSON.parse(decodeURIComponent(categoriesString));

  res.render('dna-pipeline', {toolbar, index:4, categories });
});

router.post('/running', function (req, res, next) {
  const categoriesString = req.body.categories || '[]';
  const categories = JSON.parse(decodeURIComponent(categoriesString));

  res.render('running', { toolbar, index: 5, categories });
});

router.post('/qc-check-2', (req, res) => {
  // stuff here
  res.redirect('dna-goalposts');
});

router.post('/file-information', function (req, res, next) {
  const mode = req.body.mode;
  res.render('file-info', {toolbar, index: 2});
});

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

router.post('/upload-test', upload.array('files'), (req, res) => {
  const uploadedFile = req.files[0];
  const secondFile = req.files[1];

  // Check if a file was uploaded
  if (!uploadedFile) {
    res.status(400).send('No file uploaded');
    return;
  }

  if (req.files.length > 2) {
    res.status(400).send('Too many files uploaded');
    return;
  }

  // Run FastQC

  // return wether the file has adapters or not

});

// Handle file upload
router.post('/upload', upload.single('file'), (req, res) => {
  // Access the uploaded file information
  const uploadedFile = req.file;

  // Check if a file was uploaded
  if (!uploadedFile) {
    res.status(400).send('No file uploaded');
    return;
  }

  // Access file details
  const fileName = uploadedFile.originalname;
  const filePath = uploadedFile.path;

  // // Delete the uploaded file after processing (can't display it if i do this)
  // fs.unlink(filePath, (err) => {
  //   if (err) {
  //     console.error(`Error deleting file: ${err.message}`);
  //   } else {
  //     console.log(`File deleted: ${filePath}`);
  //   }
  // });

  res.send(`<img src="/${filePath}" alt="Processed Image">`);
  // res.send(`File uploaded: ${filePath}`);
});

// Serve uploaded images statically
router.use('/uploads', express.static('uploads'));


router.post('/', (req, res) => {
  const username = req.body.username;
  const password = req.body.password;

  res.render('index', { title: password });
})

router.post('/run-picard-tools', (req, res) => {
  //// I CAN RUN PICARD TOOLS
  // const picardCommand = 'java';
  // const picardArgs = ['-jar', path.join(__dirname, '..', 'bio_modules', 'picard.jar'), 'CommandLineTool' ];

  //// I CAN RUN FASTQC
  // const picardCommand = path.join(__dirname, '..', 'bio_modules', 'FastQC.app', 'Contents', 'MacOS', 'fastqc');
  // const picardArgs = ['-version' ];

  //// I CAN RUN BWA
  // const picardCommand = path.join(__dirname, '..', 'bio_modules', 'bwa');
  // const picardArgs = ['index' ];

  //// I CAN RUN SAMTOOLS BABEY
  const picardCommand = path.join(__dirname, '..', 'bio_modules', 'samtools');
  const picardArgs = ['help'];

  const picardProcess = spawn(picardCommand, picardArgs);

  let outputData = '';

  picardProcess.stdout.on('data', (data) => {
    outputData += data;
    console.log(`stdout: ${data}`);
  });

  picardProcess.stderr.on('data', (data) => {
    console.error(`stderr: ${data}`);
  });

  picardProcess.on('exit', (code) => {
    console.log(`Picard Tools process exited with code ${code}`);
    res.send(`Picard Tools process exited with code ${code}`);
  });
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



module.exports = router;
