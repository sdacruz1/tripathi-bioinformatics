//#region ... Imports ...

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

//#endregion

//#region  ... GET and POST Views ... 

/* Home Page */
router.get('/', function (res) {
  res.render('home', { toolbar_index: 1 });
});

/* File Info */
router.post('/file-information', function (req, res) {
  mode = req.body.mode;
  res.render('file-info', { toolbar_index: 2 });
});

/* DNA Goalposts */
router.post('/dna-goalposts', function (req, res) {
  infoSteps = JSON.parse(req.body.infoSteps);
  uploadedFilePath = req.body.uploadedFilePath;
  uploadedPairedPath = req.body.uploadedPairedPath;
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
router.post('/dna-pipeline', function (req, res) {
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
router.post('/running', function (req, res) {
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
router.post('/output', function (res) {
  res.render('output', { downloadable_content, toolbar_index: 5 });
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

router.post('/run-command', async (req, res) => {
  // Get the info needed
  let Executable = DNAExecutables[req.body.command];

  // -- Main File Check --
  // Find the output file we need to run on from the database
  if (mainFilePath != '') {
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
  let allVars = (mainFilePath == '') ? Executable.variables : [['mainFile', mainFilePath]].push(Executable.variables);
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
    Cmd: command,
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

