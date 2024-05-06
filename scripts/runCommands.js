// Helper Functions

// Takes a list of indexes and checks if any of them correspond to our selected commands to be run
function containedInCommandList(commandIndexes) {
    commandIndexes.foreach(function(index) {
        if (commandList.contains(index)) { // CBTT to check if this is a real JS function lol
            return true;
        }
    });
    return false;
}

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
        if (firstCommandIndex == 1) { // don't want this step for BCL case
            Fast5File = ConvertToFast5(originalInput);
        }
    case 2: // Fast5
        if (firstCommandIndex == 1 || firstCommandIndex == 2) {// don't want this step for BCL case
            if (Fast5File == '') {
                Fast5File = originalInput;
            }
            FastQFile = Fast5ToFastQ(Fast5File);
        }
    case 3: // FastQ
        if (FastQFile == '') {
            FastQFile = originalInput;
        }
        fastQC(FastQFile); // COMMAND
    case 4: // Trimming
        // This is the first optional command. 
        if (commandList[4] || requiresTrimming) { // run only if command is selected OR if it will be needed later
            FastQFile = Trimming(trim_type, FastQFile); // COMMAND, CBTT
        }
    case 5: // Alignment
        if (commandList[5] || requiresAlignment) { // run only if command is selected
            FastQFile = Alignment(align_type, FastQFile); // COMMAND, CBTT
        }
    case 6: // SAM
        // either they are converting to SAM, or they gave me a SAM file
        // CBTT: if they don't want this step at all, what then?
        SAMFile = (commandList[6] || reqiresSAM) ? ConvertToSam(FastQFile) : originalInput; // COMMAND, CBTT
    case 7: // BAM
        // either they are converting to BAM, or they gave me a BAM file
        // CBTT: if they don't want this step at all, what then?
        BAMFile = (commandList[7] || reqiresBAM) ? ConvertToBam(SAMFile) : originalInput; // COMMAND, CBTT

    // That covers all of the steps to do with conversion.
    
    case 8: // Sort BAM
        if (commandList[8] || reqiresSortedBAM) {
            BAMFile = SortBAM(BAMFile); // COMMAND, CBTT
        }
    case 9: // Index BAM
        if (commandList[9] || reqiresIndexedBAM) {
            IndexBAM(BAMFile); // COMMAND, CBTT, don't need to store BAM File this time I think, it will be used later
        }
    case 10: // Add or Replace Read Groups
        if (commandList[10]) {
            MakeRequest(Commands.Add_Or_Replace_Read_Groups);
        }
    case 11: // Bam Index Stats
        if (commandList[11]) {
            MakeRequest(Commands.Bam_Index_Stats);
        }
    case 12: // Alignment Summary
        if (commandList[12]) {
            MakeRequest(Commands.Alignment_Dummary);
        }
    case 13: // GC Bias Summary
        if (commandList[13]) {
            MakeRequest(Commands.GC_Bias_Summary);
        }
    case 14: // Insert Size Summary
        if (commandList[14]) {
            MakeRequest(Commands.Insert_Size_Summary);
        }
    case 15: // Create Sequence Dictionary
        if (commandList[15]) {
            MakeRequest(Commands.Create_Sequence_Dictionary);
        }
    case 16: // Flag Stats
        if (commandList[16]) {
            MakeRequest(Commands.Add_Or_Replace_Read_Groups);
        }
    case 17: // Sequence Depth
        if (commandList[17]) {
            MakeRequest(Commands.Add_Or_Replace_Read_Groups);
        }
    case 18: // Sequence Coverage
        if (commandList[18]) {
            MakeRequest(Commands.Add_Or_Replace_Read_Groups);
        }
    case 19: // Mark or Remove Duplicates
        if (commandList[19]) {
            MakeRequest(Commands.Add_Or_Replace_Read_Groups);
        }

}
