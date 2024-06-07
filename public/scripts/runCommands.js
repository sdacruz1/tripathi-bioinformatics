// ==== Helper Functions ==== //

// Set Status
function setStatus(tag, status) {
    const imageElement = document.getElementById(`status_img_${tag}`);
    const textElement = document.getElementById(`status_text_${tag}`);

    if (imageElement === null || textElement === null) {
        return;
    }

    switch (status) {
        case 'running':
            imageElement.src = "../images/Hourglass.png";
            textElement.textContent = 'Running';
            break;
        case 'done':
            imageElement.src = '../images/Checked_Circle.png';
            textElement.textContent = 'Done';
            break;
        default:
            // Handle unknown status
            break;
    }
}

// Make Request
const MakeRequest = (command) => {
    return new Promise((resolve) => {
        let formData = new FormData();
        DNAParameters.get(command.serverCall).forEach(parameter => {
            formData.append(parameter.serverTag, parameter.value);
        });

        setStatus(command.serverCall, 'running');

        fetch('/' + command.commandToRun, {
            method: 'POST',
            body: formData,
        })
            .then(response => {
                if (!response.ok) {
                // If response is not OK, throw an error with the status text
                throw new Error(`HTTP error! Status: ${response.status} ${response.statusText}`);
                }
                // Assume the server sends JSON data
                return response.json();
            })
            .then(data => {
                // Handle Response
                console.log(command.title + ' response:', data);
                setStatus(command.serverCall, 'done');
                resolve(10);
            })
            .catch(error => {
                // Handle the error more gracefully
                console.error('Error with ' + command.title + ' :', error);
                setStatus(command.serverCall, 'error');
                resolve(0);
            });
    });
};

// Call RunCommands the first time
document.addEventListener('DOMContentLoaded', function() {
    // The DOM content has been fully loaded, including stylesheets, images, etc.
    console.log('View content has finished rendering');
    // Call your function here
    runCommands();
});

// ==== Variables and Data Structures ==== //

const DNACommands = JSON.parse('<%- JSON.stringify(DNACommands) %>');
const DNAParameters = JSON.parse('<%- JSON.stringify(DNAParameters) %>');

async function runCommands() {
    try {
        DNACommands.forEach(async command => {
            if (command.isEnabled) {
                await MakeRequest(command);
            }
        })
    } catch (error) {
        console.error('Error fetching data:', error);
    }
}