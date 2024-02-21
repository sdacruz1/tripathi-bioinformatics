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

/* GET home page. */
router.get('/', function (req, res, next) {
  res.render('home');
});

router.get('/file-info', function (req, res, next) {
  res.render('file-info');
});

router.get('/test', function (req, res, next) {
  res.render('test');
});

router.get('/dna-goalposts', (req, res) => {
  // Category titles
  const categories = ["File Processing", "Statistics", "Summary", "Graphs", "Metrics", "Cleanup and Sequencing"];
  // Options by category
  const fileCategory = ["Trimming", "Alignment", "BAM File"];
  const statsCategory = ["Add or Replace Read Groups", "Bam Index Stats"];
  const summaryCategory = ["Alignment Summary", "GC Bias Summary", "Insert Size Summary"];
  const graphsCategory = ["Alignment Graph", "GC Bias Graphs", "Insert Size Graphs"];
  const metricsCategory = ["Quality Yield Metrics", "Whole Genome Sequencing Metrics", "Targeted PCR Metrics"];
  const cleanupCategory = ["Create Sequence Dictionary", "Mark Duplicates", "Sort Bam File", "Flag Stats", "Sequencing Depth", "Sequencing Coverage"];
  // Assembled options by category
  const categoryOptions = [fileCategory, statsCategory, summaryCategory, graphsCategory, metricsCategory, cleanupCategory];
  // Checkbox status
  const checkboxStatus = [false, false, false, false, false, false, false, false, false, false, false, false];

  res.render('dna-goalposts', { categories, categoryOptions, checkboxStatus });
});

router.post('/qc-check-2', (req, res) => {
  // stuff here
  res.redirect('dna-goalposts');
})

router.post('/file-information', function (req, res, next) {
  const mode = req.body.mode;

  // Perform different actions based on the buttonType
  if (mode === 'custom') {
    // res.send('<h1>Custom</h1>');
    res.render('file-info');
    // res.render('index', {title: mode});
  } else if (mode === 'recommended') {
    // res.send('<h1>Recommended</h1>');
    res.render('file-info');
    // res.render('index', {title: mode});
  } else {
    res.send('<h1>Unknown button type</h1>');
  }
});

router.post('/main-file-upload', upload.array('files'), (req, res) => {
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

  const hasAdapters = false; // set this based on the result of the QC

  // return whether the file has adapters or not
  res.json({ hasAdapters });

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
