var express = require('express');
const { spawn } = require('child_process');
const path = require('path');
var router = express.Router();

/* GET home page. */
router.get('/', function(req, res, next) {
  res.render('index', { title: 'Express' });
});

// Use express.urlencoded middleware
router.use(express.urlencoded({ extended: true }));

router.post('/', (req, res) => {
  const username = req.body.username;
  const password = req.body.password;

  res.render('index', { title: password });
})

router.post('/run-picard-tools', (req, res) => {
  const picardCommand = 'java';
  const picardArgs = ['-jar', path.join(__dirname, '..', 'bio_modules', 'picard.jar'), '-h'];

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



module.exports = router;
