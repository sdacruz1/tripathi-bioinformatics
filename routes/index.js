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
router.get('/', function(req, res, next) {
  res.render('home');
});

router.get('/file-info', function(req, res, next) {
  res.render('file-info');
});

router.get('/test', function(req, res, next) {
  res.render('test');
});

router.post('/file-information', function(req, res, next) {
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

router.post('/upload-test', upload.single('file'), (req, res) => {
  const uploadedFile = req.file;

  // Check if a file was uploaded
  if (!uploadedFile) {
    res.status(400).send('No file uploaded');
    return;
  }

  // Access file details
  const fileName = uploadedFile.originalname;
  const fileSize = uploadedFile.size;

  // Log the file name on the server
  console.log(`Uploaded file name: ${fileName}`);
  console.log(`File Size: ${fileSize} bytes`);

  // You can perform additional processing or respond to the client here

  // Respond with information about the uploaded file
  res.send(`File Size: ${fileSize} bytes`);
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
  const picardArgs = ['help' ];

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
