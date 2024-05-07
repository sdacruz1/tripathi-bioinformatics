// ==== Helper Functions ==== //

// Takes a list of indexes and checks if any of them correspond to our selected commands to be run
function containedInCommandList(commandIndexes) {
    commandIndexes.foreach(function(index) {
        if (commandList.contains(index)) { // CBTT to check if this is a real JS function lol
            return true;
        }
    });
    return false;
}

// runs the Command given if it was selected by the user.
// returns true if this is the last command or beyond, and false otherwise
function runCommand(commandIndex, command) {
    if (commandList[commandIndex]) {
        MakeRequest(command);
    }

    if (lastCommandIndex >= commandIndex) {
        return true;
    }
    return false;
}

// Make Request
const MakeRequest = (command, file = null) => {
    const formData = new FormData();
    command.formDataArray.forEach(element => {
        formData.append(element[0], element[1]);
    });

    setStatus(command.commandIndex[0], command.commandIndex[1], 'running');

    fetch(command.commandToRun, {
        method: 'POST',
        body: formData,
    })
        .then(response => response.json()) // Assume the server sends JSON data
        .then(data => {
            // Handle Response
            console.log(title + ' response:', data);
            setStatus(command.commandIndex[0], command.commandIndex[1], 'done');
        })
        .catch(error => {
            console.error('Error with ' + command.title + ' :', error);
        });
};

// ==== Variables and Data Structures ==== //

const Commands = {
    Store_BCL: {
        title: 'Store_BCL',
        commandToRun: '/store-files',
        formDataArray: [],
        commandIndex: [0, 1]
    },
    BCL2FastQ: {
        title: 'BCL2FastQ',
        commandToRun: '/bcl2fastq',
        formDataArray: [],
        commandIndex: [0, 1]
   },
    FastA_To_Fast5: {
        title: 'FastA_To_Fast5',
        commandToRun: '/fasta-to-fast5',
        formDataArray: [],
        commandIndex: [0, 1]
   },
    Fast5_To_FastQ: {
        title: 'Fast5_To_FastQ',
        commandToRun: '/fast5-to-fastq',
        formDataArray: [],
        commandIndex: [0, 1]
   },
    FastQC: {
        title: 'FastQC',
        commandToRun: '/fastqc',
        formDataArray: [],
        commandIndex: [0, 1]
   },
    Trimming: {
        title: 'Trimming',
        commandToRun: '/trimming',
        formDataArray: [],
        commandIndex: [0, 1]
    },
    Alignment: {
        title: 'Alignment',
        commandToRun: '/alignment',
        formDataArray: [],
        commandIndex: [0, 2]
    },
    FastQ_To_Sam: {
        title: 'FastQ_To_Sam',
        commandToRun: '/fastq-to-sam',
        formDataArray: [],
        commandIndex: [0, 1]
    },
    Sam_To_Bam: {
        title: 'Convert to BAM File',
        commandToRun: '/convert-to-bam',
        formDataArray: [],
        commandIndex: [0, 4]
    },
    Sort_BAM_File: {
        title: 'Sort',
        commandToRun: '/sort-bam',
        formDataArray: [],
        commandIndex: [0, 5]
    },
    Index_BAM_File: {
        title: 'Index',
        commandToRun: '/index-bam',
        formDataArray: [],
        commandIndex: [0, 5]
    },
    Add_Or_Replace_Read_Groups: {
        title: 'Add_Or_Replace_Read_Groups',
        commandToRun: '/add-or-replace-read-groups',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Bam_Index_Stats: {
        title: 'Index Size',
        commandToRun: '/insert-size-data',
        formDataArray: [],
        commandIndex: [2, 2]
    },
    Alignment_Summary: {
        title: 'Alignment',
        commandToRun: '/alignment-data',
        formDataArray: [],
        commandIndex: [2, 0]
    },
    GC_Bias_Summary: {
        title: 'GC Bias',
        commandToRun: '/gc-bias-data',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Insert_Size_Summary: {
        title: 'Insert_Size_Summary',
        commandToRun: '/insert-size-summary',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Create_Sequence_Dictionary: {
        title: 'Create_Sequence_Dictionary',
        commandToRun: '/create-sequence-dictionary',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Flag_Stats: {
        title: 'Flag_Stats',
        commandToRun: '/flag-stats',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Sequence_Depth: {
        title: 'Sequence_Depth',
        commandToRun: '/sequence-depth',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Sequence_Coverage: {
        title: 'Sequence_Coverage',
        commandToRun: '/sequence-coverage',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Add_Or_Replace_Read_Groups: {
        title: 'Add_Or_Replace_Read_Groups',
        commandToRun: '/add-or-replace-read-groups',
        formDataArray: [],
        commandIndex: [2, 1]
    },
    Mark_Or_Remove_Duplicates: {
        title: 'Mark_Or_Remove_Duplicates',
        commandToRun: '/mark-or-remove-duplicates',
        formDataArray: [],
        commandIndex: [2, 1]
    },
};


// optional variables that will hold filepaths.
let FastQFile = '';
let Fast5File = '';
let SAMFile = '';
let BAMFile = '';

// variables to track mandatory conversions
let requiresTrimming = containedInCommandList([1, 2, 3, 4]); // CBTT to put in the actual command indexes
let requiresAlignment = false; // CBTT to put in the actual command indexes
let reqiresSAM = false; // CBTT to put in the actual command indexes
let reqiresBAM = false; // CBTT to put in the actual command indexes
let reqiresSortedBAM = false; // CBTT to put in the actual command indexes
let reqiresIndexedBAM = false; // CBTT to put in the actual command indexes

// otherwise necessary variables
let originalInput = ''; // This will be the file(s) I was given by the user on the upload page
let firstCommandIndex = 0; // This will indicate the first 'true', AKA activated step the user selected.
switch(firstCommandIndex) {
    case 0: // BCL
        let BCL_folder = store_BCL(originalInput);
        FastQFile = BCL2FastQ(BCL_folder);
        // note the lack of a 'break;' here. This is intentional. I want it to keep going through commands sequentially.
    case 1: // FastA
        if (firstCommandIndex != 0) { // don't want this step for BCL case
            Fast5File = MakeRequest(Commands.FastA_To_Fast5, originalInput);
        }
    case 2: // Fast5
        if (firstCommandIndex != 0) {// don't want this step for BCL case
            Fast5File = (Fast5File == '') ? originalInput : Fast5File;
            FastQFile = MakeRequest(Commands.Fast5_To_FastQ, originalInput);
        }
    case 3: // FastQ
        FastQFile = (FastQFile == '') ? originalInput : FastQFile;
        fastQC(FastQFile); // COMMAND
    case 4: // Trimming
        // This is the first optional command. 
        if (commandList[4] || requiresTrimming) { // run only if command is selected OR if it will be needed later
            FastQFile = MakeRequest(Commands.Trimming, FastQFile); // COMMAND, CBTT
        }
    case 5: // Alignment
        if (commandList[5] || requiresAlignment) { // run only if command is selected
            FastQFile = MakeRequest(Commands.Alignment, FastQFile); // CBTT, needs alignment type and info and stuff
        }
    case 6: // SAM
        // either they are converting to SAM, or they gave me a SAM file
        // CBTT: if they don't want this step at all, what then?
        SAMFile = (commandList[6] || reqiresSAM) ? MakeRequest(Commands.FastQ_To_Sam, FastQFile) : originalInput; // COMMAND, CBTT
    case 7: // BAM
        // either they are converting to BAM, or they gave me a BAM file
        // CBTT: if they don't want this step at all, what then?
        BAMFile = (commandList[7] || reqiresBAM) ? MakeRequest(Commands.Sam_To_Bam, SAMFile) : originalInput; // COMMAND, CBTT

    // That covers all of the steps to do with conversion.
    
    default:
        // Sort BAM
        if (commandList[8] || reqiresSortedBAM) {
            BAMFile = MakeRequest(Commands.Sort_BAM_File, BAMFile);
        }
        // Index BAM
        if (commandList[9] || reqiresIndexedBAM) {
            MakeRequest(Commands.Index_BAM_File, BAMFile);
        }
        // Add or Replace Read Groups
        if (runCommand(10, Commands.Add_Or_Replace_Read_Groups)) { break; }
        // Bam Index Stats
        if (runCommand(11, Commands.Bam_Index_Stats)) { break; }
        // Alignment Summary
        if (runCommand(12, Commands.Alignment_Summary)) { break; }
        // GC Bias Summary
        if (runCommand(13, Commands.GC_Bias_Summary)) { break; }
        // Insert Size Summary
        if (runCommand(14, Commands.Insert_Size_Summary)) { break; }
        // Create Sequence Dictionary
        if (runCommand(15, Commands.Create_Sequence_Dictionary)) { break; }
        // Flag Stats
        if (runCommand(16, Commands.Flag_Stats)) { break; }
        // Sequence Depth
        if (runCommand(17, Commands.Sequence_Depth)) { break; }
        // Sequence Coverage
        if (runCommand(18, Commands.Sequence_Coverage)) { break; }
        // Mark or Remove Duplicates
        if (runCommand(19, Commands.Mark_Or_Remove_Duplicates)) { break; }
}
