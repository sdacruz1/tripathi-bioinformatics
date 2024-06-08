import { Command, Parameter } from '../data/Structure.js';  // Adjust the path as needed

const orm = document.getElementById('goalpostorm');
const dnaCommandsData = JSON.parse(orm.dataset.dnaCommands);
const DNACommands = new Map(dnaCommandsData);

// Now you can use DNACommands in your script
console.log(DNACommands);
                    const form = document.getElementById('goalpostForm');
                    // let DNACommands = new Map(JSON.parse('<%- JSON.stringify(Array.from(DNACommands.entries())) %>'));
                    // let DNACategories = JSON.parse('<%- JSON.stringify(DNACategories) %>');
                    // const infoSteps = JSON.parse('<%- JSON.stringify(infoSteps) %>');

                    form.addEventListener('change', function (event) {
                        if (event.target.type === 'checkbox') {
                            // Set the correct array element accordingly
                            console.log ("BEFORE");
                            console.log(DNACommands);
                            console.log(DNACommands.get(event.target.name));
                            
                            changingEntry = DNACommands.get(event.target.name);
                            DNACommands.set(event.target.name, changingEntry.setEnabled(event.target.checked));

                            console.log ("AFTER");
                            console.log(DNACommands);
                            console.log(DNACommands.get(event.target.name));

                        }
                    });

                    form.addEventListener('submit', function(event) {
                        // prevent submission for now
                        event.preventDefault();

                        console.log(DNACommands);

                        // send the DNACommands array as a hidden input, to make sure it updates BTS
                        const sendCommands = document.querySelector('[name="DNACommands"]');
                        sendCommands.value = JSON.stringify(Array.from(DNACommands.entries()));

                        // Submit the form
                        form.submit();
                    });

                    // JavaScript to handle the toggling
                    function toggleCollapsible(sectionId) {
                        const collapsibleContent = document.getElementById(sectionId);
                        collapsibleContent.classList.toggle('active');
                    }