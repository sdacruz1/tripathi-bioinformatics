//#region ... Imports ...
var express = require('express');
const multer = require('multer');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const bodyParser = require('body-parser');
const Docker = require('dockerode');
const docker = new Docker();
var router = express.Router();

const { DNACategories, DNACommands, DNAParameters, RNACategories, RNACommands, RNAParameters } = require('../public/data/SiteData');

//#endregion

//#region ... Middleware ...

router.use(bodyParser.urlencoded({ extended: true }));
router.use(bodyParser.json());
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

//#endregion

//#region  ... Constants ...

// Program Mode and Input Type
let mode = "";
let inputType = "";

// Stored Outputs
let infoSteps;          // An array that determines which file processing steps will be available in the timeline
let uploadedFilePath = "";   // The path to the original file that was uploaded or the directory, if it was a BCL folder
let uploadType = '';
let uploadedPairedPath = "";
let uploadedAdapterPath = "";
let uploadedGenomeFile = "";
let uploadedFileType;   // A string representing the original uploaded file type
let fastQConversion;    // The path to the fastQ file result of any conversions (BCL, FastA, Fast5)
let fastQCResults;      // The path to the directory containing the results of running FastQC on the uploaded file
let BAMFile;            // The path to the BAM file result of any conversions
let downloadable_content = [  // An object array holding the result of the steps in the actual timeline
  // {  -- EXAMPLE --
  //     enabled: true,
  //     has_visual_component: false,
  //     label: "Test Object BAM",
  //     content: "../output/SortedBAM.bam"
  // },
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
  uploadType ='single';

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
  if (uploadedFiles.length === 2) {
    uploadedFilePath = path.join(__dirname, '..', 'uploads', uploadedFiles[0].filename);
    uploadedFilePath2 = path.join(__dirname, '..', 'uploads', uploadedFiles[1].filename);
    return res.status(200).send({ storedPath: [uploadedFilePath, uploadedFilePath2] });
  }

  // Create an empty directory to hold them in
  const UploadedDirectory = 'UploadedFiles-' + Date.now();
  if (!fs.existsSync(path.join(__dirname, '..', 'uploads', UploadedDirectory))) {
    fs.mkdirSync(path.join(__dirname, '..', 'uploads', UploadedDirectory));
  }

  // Move each uploaded file to the empty directory
  uploadedFiles.forEach(file => {
    const destinationPath = path.join(path.join(__dirname, '..', 'uploads', UploadedDirectory), file.originalname);
    fs.renameSync(file.path, destinationPath);
  });

  uploadedFilePath = UploadedDirectory;
  uploadType ='multiple';

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

//#endregion

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

  // -- Special Case: Requires an output directory --
  if(Executable.special[0] == 'output_dir') {
    // Create an empty directory to hold the output
    const OutputDirectory = path.join(__dirname, 'output', Executable.special[1]);
    if (!fs.existsSync(OutputDirectory)) {
      fs.mkdirSync(OutputDirectory);
    }
  }

  // -- Build the Command --
  let allVars = (Executable.checkFile == '') ? Executable.variables : [['mainFile', mainFilePath]].push(Executable.variables);
  let commandString = substitutePlaceholders(allVars, Executable.command);

  // -- Run the Command in a Docker Environment --
  try {
    const output = await runDockerCommand(Executable.environment, commandString);
    console.log('Output:', output);

    // -- Special Case: Writes to stdout, have to record the output manually --
    if(Executable.special[0] == 'readout') {
      // Write output data to a text file
      fs.writeFile(path.join(__dirname, '..', 'output', Executable.special[1]), output, (err) => {
        if (err) {
          console.error('Error writing file:', err);
          return;
        }
        console.log('stdout data written to ' + Executable.special[1]);
      });
    }

    // -- Make downloadable content entries for the results --
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

module.exports = router;