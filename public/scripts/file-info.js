// The variables that will be set in this view and passed to the server
let infoSteps = [0, 0, 0, 0, 0, 0];
let uploadedFilePath = "";
let fastQConversion = "";
let fastQCResults = "";
let uploadedFileType = "";
let trimmed = false;
let demultiplexed = false;
let uploadedAdapterFile = "";

let input_type = '';
let uploadedFile;
let uploadedFile2;

//#region ... REFERENCES ... 

// Get references to the buttons and the RNA Genome Uploader
const RNA_Genome_Uploader = document.getElementById('RNA_Genome_Uploader');

// Get referneces to the sections
const fileSection = document.getElementById('file_uploader_section');
const adapterSection = document.getElementById('adapter_section');
const rnaSection = document.getElementById('rna_section');

// Get references to the upload button and the file information text
const Main_File_Uploader = document.getElementById('Main_File_Uploader');
const Main_File_Uploader_Second = document.getElementById('Main_File_Uploader-Second');
const File1Label = document.getElementById('file1label');
const File2Label = document.getElementById('file2label');
const typeDropdown = document.getElementById('directoryOption');

// Get references to the Adapter objects
hasAdapters = false;
const ACHeading = document.getElementById('adapterCheckHeading');
const adapterFileInput = document.getElementById('Adapter_File_Uploader');
const adapterUploadButton = document.getElementById('adapterUploadButton');

// Get a reference to the Genome objects
const genomeTypeInfo = document.getElementById('genomeTypeInfo');
const genomeSection = document.getElementById('genome_section');
const genomeFileInput = document.getElementById('RNA_Genome_Uploader');

// Get a reference to the continue button
const continue_button = document.getElementById('file_info_continue_button');

//#endregion

/** @ select_input_type
 * 
 * Sets the input type to DNA or RNA
 *
 * @param {string} type - The type to set it to. 
 */
function select_input_type(type) {
    input_type = type;
    Main_File_Uploader.disabled = false;
    fileSection.hidden = false;
    if (type == 'RNA') {
        rnaSection.hidden = false;
        RNA_Genome_Uploader.disabled = false;
    } else {
        rnaSection.hidden = true;
    }
};

// Function to change main file input type based on dropdown selection
function changeInputType() {
    if (typeDropdown.value === 'BCL') {
        Main_File_Uploader.setAttribute('webkitdirectory', 'true');
        Main_File_Uploader_Second.disabled = true;
        File1Label.innerText = "Upload A BCL Folder";
        File2Label.innerText = "";
    } else if (typeDropdown.value === 'Paired') {
        Main_File_Uploader.removeAttribute('webkitdirectory');
        Main_File_Uploader_Second.disabled = false;
        File1Label.innerText = "Upload FastQ File Read 1";
        File2Label.innerText = "Upload FastQ File Read 2";
    } else {
        Main_File_Uploader.removeAttribute('webkitdirectory');
        Main_File_Uploader_Second.disabled = true;
        File1Label.innerText = "Upload A " + typeDropdown.value + " File";
        File2Label.innerText = "";
    }
}

/** @ processMainFile
 * 
 * This function is called when the main data file is uploaded.
 * It stores the file, determines its type, runs FastQC if applicable, then unlocks the next step.
 *
 */
function processMainFile() {

    // If it's a directory, this should be a BCL folder.
    if (typeDropdown.value === 'BCL') {
        /** This function will run BCL2FastQ on the file.
         * Then it will run FastQC on the resulting file.
         * It will return:
         *      - The path to the folder holding the original files
         *      - The path to the folder holding the converted FastQ files
         *      - The path to the folder holding the results of FastQC
         *      - A boolean (trimmed) determining wether the adapter file is trimmed or not
         *      - A boolean (demultiplexed) determining wether the fastQ file is demultiplexed or not
        */
        BCLresults = handleBCL(Main_File_Uploader.files);
        
        // Set the variables needed later
        infoSteps           = [1, 1, 1, 1, 1, 1];
        uploadedFileType    = 'BCL';
        uploadedFilePath    = BCLresults.originalFile;
        fastQConversion     = BCLresults.bcl2fastq_results;
        fastQCResults       = BCLresults.fastqc_results;
        trimmed             = BCLresults.trimmed;
        demultiplexed       = BCLresults.demultiplexed;

        // NOTE: Do I need to worry about adapter trimming steps here?
        // Nothing more to do for DNA
        if (input_type === 'DNA') {
            continue_button.disabled = false;
        } else {
            // For RNA, unhide the RNA section and enable the RNA uploader
            // NOTE: Return here for RNA
        }
        return;
    }
    // If the file somehow doesn't exist, abort process
    uploadedFile = Main_File_Uploader.files[0];
    if (!uploadedFile) {
        return;
    }

    // Store it
    storeFiles(Main_File_Uploader.files)
        .then(fileStorage => {
            uploadedFilePath = fileStorage.storedPath[0];
            if (fileStorage.storedPath.length > 1) {
                uploadedPairedPath = fileStorage.storedPath[1];
            }
        })
        .catch(error => {
            // Handle any errors
            console.error('Error analyzing file:', error);
        });

    // If it's not a directory, it's a file
    // Determine the type
    switch (typeDropdown.value) {
        case 'BAM':
            // No processing, Cleaning BAM is required
            infoSteps = [0, 0, 0, 0, 0, 2];
            uploadedFileType = 'BAM';
            fastQConversion = "";
            fastQCResults = "";
            trimmed = false;
            demultiplexed = false;

            // if input type is DNA enable the continue button
            if (input_type === 'DNA') {
                continue_button.disabled = false;
            } else {
                // For RNA, unhide the RNA section and enable the RNA uploader
                // NOTE: Return here for RNA
            }
            return;
        case 'SAM':
            // No processing, Convert to BAM is required, Clean BAM is optional
            infoSteps = [0, 0, 0, 0, 2, 1];
            uploadedFileType = 'SAM';
            fastQConversion = "";
            fastQCResults = "";
            trimmed = false;
            demultiplexed = false;

            // if input type is DNA enable the continue button
            if (input_type === 'DNA') {
                continue_button.disabled = false;
            } else {
                // For RNA, unhide the RNA section and enable the RNA uploader
                // NOTE: Return here for RNA
            }
            return;
        case 'FastQ':
        case 'Paired':
            // Leave all of the steps on
            infoSteps = [1, 1, 1, 1, 1, 1];
            uploadedFileType = 'FastQ';
            fastQConversion = "";
            uploadedFilePath = uploadedFile.path;

            // run fastQC
            runFastQC(uploadedFile)
                .then(analysis => {
                    fastQCResults = analysis.results;
                    trimmed = analysis.trimmed;
                    demultiplexed = analysis.demultiplexed;

                    // check for demultiplexing
                    if (!analysis.demultiplexed) {
                        // still multiplexed? Set Demultiplexing to required
                        infoSteps[3] = 2;
                    }
                })
                .catch(error => {
                    // Handle any errors
                    console.error('Error analyzing file:', error);
                });
            break;
        case 'FastA':
        case 'Fast5':
            // NOTE: Nothing yet, come back to this after meeting
            break;
        default:
            // provide error
            console.log('Provided file is not one of the supported file types.');
    }

    // Set the adapter text/variable and uploader enabled settings accordingly
    adapterSection.hidden = false;
    ACHeading.textContent = 'If applicable, upload an adapter file.'
    adapterFileInput.disabled = false;
    // if input type is DNA and there are no adapters, enable the continue button
    if (input_type === 'DNA') {
        continue_button.disabled = false;
    } else {
        genomeSection.hidden = false;
        genomeFileInput.disabled = false;
    }

}

/** @ handleBCL
 * This function takes in a BCL File folder, converts it to fastQ, runs fastQC, and returns trimming and multiplexing information
 * Returns:
 *      - The path to the folder holding the original files
 *      - The path to the folder holding the converted FastQ files
 *      - The path to the folder holding the results of FastQC
 *      - A boolean (trimmed) determining wether the adapter file is trimmed or not
 *      - A boolean (demultiplexed) determining wether the fastQ file is demultiplexed or not
 * */
function handleBCL (BCLFiles) {
    const formData = new FormData();
    for (const file of BCLFiles) {
        formData.append('files', file);
    }

    // Run Bcl2FastQ
    return new Promise((resolve, reject) => {
    fetch('/run-bcl2fastq', {
        method: 'POST',
        body: formData,
    })
        .then(response => response.json())
        .then(data => {
            resolve();
            // // Run fastQC on the returned file
            // runFastQC(data.fastQFile)
            //     .then(analysis => {
                
                // resolve({ originalFile: data.BCLFolder,
                //         bcl2fastq_results: data.results,
                //         fastqc_results: analysis.results,
                //         trimmed: analysis.trimmed,
                //         demultiplexed: analysis.demultiplexed });
                // })

                // .catch(error => {
                //     // Handle any errors
                //     console.error('Error analyzing file:', error);
                // });
        })
        .catch(error => {
            console.error('Error:', error);
            reject(error);
        });
    });
}

// @ runFastQC
// This function takes in a fastQ file and runs fastQC on that file
// Returns:
//      - a string that is the path to the directory holding the results of the FastQC process
//      - a boolean that indicates whether the adapter is trimmed
//      - a boolean that indicates wether the fastQ file is demultiplexed
function runFastQC (fastQFile) {
    const formData = new FormData();
    formData.append('file', fastQFile);

    return new Promise((resolve, reject) => {
        // Make an asynchronous request to the server using the Fetch API
        fetch('/run-fastqc', {
            method: 'POST',
            body: formData,
        })
            .then(response => response.json()) // Assume the server sends JSON data
            .then(data => {
                // return the results, the Adapter Trim bool and the multiplexed bool based on the response
                resolve({ results: data.results, trimmed: data.trimmed, demultiplexed: data.demultiplexed });
            })
            .catch(error => {
                console.error('Error:', error);
                reject(error);
            });
    });
}

// If an adapter file is uploaded...
// nothing to do!
// Send it to the backend
function uploadAdapterFiles() {
    let adapterFile = adapterFileInput.files[0];
    if (!adapterFile) {
        return;
    }
    storeFiles(adapterFileInput.files)
        .then(fileStorage => {
            uploadedAdapterFile = fileStorage.storedPath[0];
        })
        .catch(error => {
            // Handle any errors
            console.error('Error analyzing file:', error);
        });
}
// and set up a listener to activate the continue button if input type is DNA
adapterUploadButton.addEventListener('click', function () {
    if (input_type === 'DNA') {
        continue_button.disabled = false;
    }
});

/** @ storeFiles
 * 
 * Takes in the files to be stored
 * 
 * Returns:
 *      - a path to the files in storage
*/
function storeFiles(files) {
    const fileFormData = new FormData();

    for (const file of files) {
        fileFormData.append('files', file);
    }

    return new Promise((resolve, reject) => {
        // Make an asynchronous request to the server using the Fetch API
        fetch('/store-files', {
            method: 'POST',
            body: fileFormData,
        })
            .then(response => response.json()) // Assume the server sends JSON data
            .then(data => {
                // Store the path to it
                resolve({storedPath: data.storedPath});
            })
            .catch(error => {
                console.error('Error:', error);
                reject(error);
            });
    });
}

// NOTE: Come back to this after RNA
// If a genome file is uploaded...
function processGenomeFile() {
    if (!genomeFileInput.files[0]) {
        return;
    }
    storeFiles(genomeFileInput.files)
        .then(fileStorage => {
            uploadedGenomeFile = fileStorage.storedPath[0];
            continue_button.disabled = false;
        })
        .catch(error => {
            // Handle any errors
            console.error('Error analyzing file:', error);
        });
}

// Redirects the submit call and adds the variables to pass to the server
function submitContinueForm() {

    document.getElementById('infoStepsInput').value = JSON.stringify(infoSteps);
    document.getElementById('uploadedFileInput').value = uploadedFilePath;
    document.getElementById('uploadedFileInput2').value = uploadedPairedPath;
    document.getElementById('uploadedAdapter').value = uploadedAdapterFile;
    document.getElementById('uploadedFileTypeInput').value = uploadedFileType;
    document.getElementById('inputTypeInput').value = input_type;
    document.getElementById('fastQConversionInput').value = fastQConversion;
    document.getElementById('fastQCResultsInput').value = fastQCResults;
    document.getElementById('trimmedInput').value = trimmed;
    document.getElementById('demultiplexedInput').value = demultiplexed;
    document.getElementById('uploadedAdapterFileInput').value = uploadedAdapterFile;
    document.getElementById('uploadedGenomeFileInput').value = uploadedGenomeFile;

    // Submit the form
    document.getElementById('continueForm').submit();
}