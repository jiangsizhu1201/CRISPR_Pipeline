 function openTab(event, tabName) {
    var i, tabcontent, tabbuttons;
    tabcontent = document.getElementsByClassName("tab-content");
    tabbuttons = document.getElementsByClassName("tab-button");

    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].classList.remove("active");
    }

    for (i = 0; i < tabbuttons.length; i++) {
        tabbuttons[i].classList.remove("active");
    }

    document.getElementById(tabName).classList.add("active");
    event.currentTarget.classList.add("active");
}

// Open the first tab by default
document.addEventListener('DOMContentLoaded', (event) => {
    document.querySelector('.tab-button').click();
});

function addBatch(tabName) {
    var batchName = prompt("Enter Batch Name:");
    if (batchName) {
        var batchesContainer = document.getElementById(tabName + '-batches');

        var batchContainer = document.createElement('div');
        batchContainer.className = 'batch-container';

        var batchHeader = document.createElement('div');
        batchHeader.className = 'batch-header';

        var batchTitle = document.createElement('h4');
        batchTitle.textContent = 'Batch: ' + batchName;

        var removeBatchBtn = document.createElement('button');
        removeBatchBtn.className = 'remove-batch-button';
        removeBatchBtn.textContent = 'Remove Batch';
        removeBatchBtn.onclick = function() {
            batchesContainer.removeChild(batchContainer);
        };

        batchHeader.appendChild(batchTitle);
        batchHeader.appendChild(removeBatchBtn);

        batchContainer.appendChild(batchHeader);

        var sequencesContainer = document.createElement('div');
        sequencesContainer.className = 'sequences-container';

        batchContainer.appendChild(sequencesContainer);

        var addSequenceBtn = document.createElement('button');
        addSequenceBtn.className = 'add-sequence-button';
        addSequenceBtn.textContent = 'Add Sequence';
        addSequenceBtn.onclick = function() {
            addSequence(sequencesContainer, tabName, batchName);
        };

        batchContainer.appendChild(addSequenceBtn);

        // Initially add one sequence
        addSequence(sequencesContainer, tabName, batchName);

        batchesContainer.appendChild(batchContainer);

        // Add covariate section
        var covariatesContainer = document.createElement('div');
        covariatesContainer.className = 'covariates-container';

        var addCovariateBtn = document.createElement('button');
        addCovariateBtn.className = 'add-covariate-button';
        addCovariateBtn.textContent = 'Add Covariate';
        addCovariateBtn.onclick = function() {
            addCovariate(covariatesContainer);
        };

        batchContainer.appendChild(covariatesContainer);
        batchContainer.appendChild(addCovariateBtn);

        batchesContainer.appendChild(batchContainer);
    }
}

function addSequence(container, tabName, batchName) {
    var inputGroup = document.createElement('div');
    inputGroup.className = 'input-group';

    var r1Input = document.createElement('input');
    r1Input.type = 'text';
    r1Input.placeholder = 'Enter R1 file path...';
    r1Input.name = tabName.toLowerCase() + '_r1[' + batchName + '][]';

    var r2Input = document.createElement('input');
    r2Input.type = 'text';
    r2Input.placeholder = 'Enter R2 file path...';
    r2Input.name = tabName.toLowerCase() + '_r2[' + batchName + '][]';

    var removeBtn = document.createElement('button');
    removeBtn.textContent = 'Remove';
    removeBtn.onclick = function() {
        container.removeChild(inputGroup);
    };

    inputGroup.appendChild(r1Input);
    inputGroup.appendChild(r2Input);
    inputGroup.appendChild(removeBtn);

    container.appendChild(inputGroup);
}

function addCovariate(covariatesContainer) {
    var covariateName = prompt("Enter Covariate Name:");
    if (covariateName) {
        var covariateBtn = document.createElement('button');
        covariateBtn.className = 'covariate-button';
        covariateBtn.textContent = covariateName;

        covariateBtn.onclick = function() {
            // Remove the covariate when clicked
            covariatesContainer.removeChild(covariateBtn);
        };

        // Add the covariate button to the container
        covariatesContainer.appendChild(covariateBtn);
    }
}

// Function to count the number of rows in the table
function updateRowCount(tabName) {
    var tableBody = document.getElementById(tabName.toLowerCase() + '-table').getElementsByTagName('tbody')[0];
    var rowCount = tableBody.rows.length;
    document.getElementById(tabName.toLowerCase() + '-row-count').textContent = rowCount; // Update the span with the row count
}

// Load File Functionality
function loadFile(tabName) {
    var fileInput = document.getElementById(tabName.toLowerCase() + '-file-input');
    var file = fileInput.files[0];

    if (file) {
        Papa.parse(file, {
            header: true,
            delimiter: detectDelimiter(file.name),
            complete: function(results) {
                populateTable(results.data, tabName);
                updateRowCount(tabName);
            }
        });
    } else {
        alert("Please select a file.");
    }
}

// Detect Delimiter Based on File Extension
function detectDelimiter(fileName) {
    if (fileName.endsWith('.tsv')) {
        return '\t';
    } else {
        return ','; // Default to CSV
    }
}

// Toggle Paste Area Visibility
function togglePasteArea(tabName) {
    var pasteArea = document.getElementById(tabName + '-paste-area');
    var parseButton = document.getElementById(tabName + '-parse-paste-button');

    if (pasteArea.style.display === 'none') {
        pasteArea.style.display = 'block';
        parseButton.style.display = 'inline-block';
    } else {
        pasteArea.style.display = 'none';
        parseButton.style.display = 'none';
    }
}

// Parse Pasted Data
function parsePastedData(tabName) {
    var data = document.getElementById(tabName + '-paste-area').value;

    Papa.parse(data, {
        header: true,
        delimiter: '\t',
        complete: function(results) {
            populateTable(results.data, tabName);
            updateRowCount(tabName);
        }
    });
}

// Populate Table with Data
function populateTable(data, tabName) {
    var tableBody;
    var requiredFields;

    if (tabName === 'Hash') {
        tableBody = document.querySelector('#hash-table tbody');
        requiredFields = ['HTO_ID', 'HTO_sequences'];
    } else if (tabName === 'Guides') {
        tableBody = document.querySelector('#guides-table tbody');
        requiredFields = ['sgRNA_ID', 'sgRNA_sequences', 'Target_name', 'chr', 'start', 'end'];
    } else if (tabName === 'Pairs') {
        tableBody = document.querySelector('#pairs-table tbody');
        requiredFields = ['guide_id', 'gene_name', 'intended_target_name', 'pair_type'];
    }  
    else {
        return;
    }
    // Check for required fields
    var firstRow = data[0];
    var missingFields = requiredFields.filter(field => !(field in firstRow));

    if (missingFields.length > 0) {
        alert("Error: Missing required fields: " + missingFields.join(', '));
        return;
    }

    tableBody.innerHTML = ''; // Clear existing data

    data.forEach(function(row) {
        var tr = document.createElement('tr');

        requiredFields.forEach(function(field) {
            var td = document.createElement('td');
            td.textContent = row[field] || '';
            tr.appendChild(td);
        });

        tableBody.appendChild(tr);
    });
    updateRowCount(tabName);
}

// Export Table Data
function exportTable(tableId) {
    var table = document.getElementById(tableId);
    var rows = Array.from(table.querySelectorAll('tr'));
    var csvContent = "";

    rows.forEach(function(row) {
        var cols = Array.from(row.querySelectorAll('th, td'));
        var rowData = cols.map(function(col) {
            return '"' + col.textContent.replace(/"/g, '""') + '"';
        }).join(',');
        csvContent += rowData + "\r\n";
    });

    var blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    var link = document.createElement("a");

    if (navigator.msSaveBlob) { // IE 10+
        navigator.msSaveBlob(blob, tableId + ".csv");
    } else {
        link.href = URL.createObjectURL(blob);
        link.setAttribute('download', tableId + '.csv');
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }
}

// Show/Hide Inference Options
function showInferenceOption(option) {
    var distanceInput = document.getElementById('distance-input');
    var predefinedPairsSection = document.getElementById('predefined-pairs-section');

    if (option === 'by_distance') {
        distanceInput.style.display = 'block';
        predefinedPairsSection.style.display = 'none';
    } else if (option === 'in_trans') {
        distanceInput.style.display = 'none';
        predefinedPairsSection.style.display = 'none';
    } else if (option === 'predefined_pairs') {
        distanceInput.style.display = 'none';
        predefinedPairsSection.style.display = 'block';
    }
}

function showEvaluationOption(option) {
    document.getElementById('user_central_nodes_section').style.display = option === 'user_central_nodes' ? 'block' : 'none';
    document.getElementById('central_nodes_num_section').style.display = option === 'central_nodes_num' ? 'block' : 'none';
}

// Disables all other fields
function toggleFields() {
    // For Perturbation Inference Tab
    const inferenceMethod = document.getElementById('inference_method').value;
    const isSceptrePerturbation = inferenceMethod === 'perturbo';

    const fieldsToDisablePerturbation = ['moi', 'side', 'grna_integration_strategy', 'resampling_approximation', 'control_group', 'resampling_mechanism', 'formula_object'];

    fieldsToDisablePerturbation.forEach(id => {
        const field = document.getElementById(id);
        field.disabled = isSceptrePerturbation;
        field.style.backgroundColor = isSceptrePerturbation ? '#e0e0e0' : '';
    });

    // For Cell Guide Assignment Settings
    const assignmentMethod = document.getElementById('assignment_method').value;
    const isSceptreAssignment = assignmentMethod === 'sceptre';

    const thresholdField = document.getElementById('THRESHOLD');
    thresholdField.disabled = isSceptreAssignment;
    thresholdField.style.backgroundColor = isSceptreAssignment ? '#e0e0e0' : '';
}

// Initialize the fields when the page loads
toggleFields();

// Attach event listeners to both dropdowns
document.getElementById('inference_method').addEventListener('change', toggleFields);
document.getElementById('assignment_method').addEventListener('change', toggleFields);


function collectBatchSequences(tabName, configData) {
    var batchesContainer = document.getElementById(tabName + '-batches');
    var batchContainers = batchesContainer.getElementsByClassName('batch-container');

    for (var i = 0; i < batchContainers.length; i++) {
        var batchContainer = batchContainers[i];
        var batchName = batchContainer.querySelector('.batch-header h4').textContent.replace('Batch: ', '');

        // Collect covariates for this batch
        var covariatesContainer = batchContainer.querySelector('.covariates-container');
        var covariateButtons = covariatesContainer.getElementsByClassName('covariate-button');
        var covariates = [];

        for (var k = 0; k < covariateButtons.length; k++) {
            var covariateName = covariateButtons[k].textContent;
            covariates.push(covariateName);
        }

        // If covariates exist, append them to the batchName
        if (covariates.length > 0) {
            batchName += ', ' + covariates.join(', ');
        }

        console.log("Final batchName with covariates: ", batchName); // Debugging

        // Collect sequences for this batch
        var sequences = batchContainer.getElementsByClassName('input-group');

        for (var j = 0; j < sequences.length; j++) {
            var sequenceGroup = sequences[j];
            var read1 = sequenceGroup.children[0].value;
            var read2 = sequenceGroup.children[1].value;

            configData.push({
                tab_name: tabName,
                variable: 'sequence',
                variable_value: '',
                batch_name: batchName,
                read1: read1,
                read2: read2
            });
        }
    }
}

// Export Configuration Functionality
function exportConfig() {
    var configData = [];
    var tabs = ['GeneralConfig', 'scRNA', 'Guides', 'Hash', 'CellGuideAssignment', 'PerturbationInference', 'InferencePairs', 'QualityFilters'];

    // Collect data from General Config Tab
    var runName = document.getElementById('run_name').value;
    var runDescription = document.getElementById('run_description').value;
    var transcriptome = document.getElementById('transcriptome').value;
    var seqspecsDirectory = document.getElementById('seqspecs_directory').value;
    var genomeDownloadPath = document.getElementById('genome_download_path').value;
    var genomeLocalPath = document.getElementById('genome_local_path').value;

    configData.push({ tab_name: 'GeneralConfig', variable: 'run_name', variable_value: runName });
    configData.push({ tab_name: 'GeneralConfig', variable: 'run_description', variable_value: runDescription });
    configData.push({ tab_name: 'GeneralConfig', variable: 'transcriptome', variable_value: transcriptome });
    configData.push({ tab_name: 'GeneralConfig', variable: 'seqspecs_directory', variable_value: seqspecsDirectory });
    configData.push({ tab_name: 'GeneralConfig', variable: 'genome_download_path', variable_value: genomeDownloadPath });
    configData.push({ tab_name: 'GeneralConfig', variable: 'genome_local_path', variable_value: genomeLocalPath });

    // Collect data from scRNA Tab
    var scrnaSeqspecYaml = document.getElementById('scRNA-seqspec-yaml').value;
    var gtfUrl = document.getElementById('gtf_url').value;

    configData.push({ tab_name: 'scRNA', variable: 'scRNA_seqspec_yaml', variable_value: scrnaSeqspecYaml });
    configData.push({ tab_name: 'scRNA', variable: 'gtf_url', variable_value: gtfUrl });

    // Collect sequences from scRNA batches
    collectBatchSequences('scRNA', configData);

    // Collect data from Guides Tab
    var guidesSeqspecYaml = document.getElementById('Guides-seqspec-yaml').value;
    configData.push({ tab_name: 'Guides', variable: 'Guides_seqspec_yaml', variable_value: guidesSeqspecYaml });

    // Collect sequences from Guides batches
    collectBatchSequences('Guides', configData);

    // Collect data from Hash Tab
    var hashSeqspecYaml = document.getElementById('Hash-seqspec-yaml').value;
    configData.push({ tab_name: 'Hash', variable: 'Hash_seqspec_yaml', variable_value: hashSeqspecYaml });

    // Collect sequences from Hash batches
    collectBatchSequences('Hash', configData);

    // Collect data from Cell Guide Assignment Tab
    var assignmentMethod = document.getElementById('assignment_method').value;
    var threshold = document.getElementById('THRESHOLD').value;

    configData.push({ tab_name: 'CellGuideAssignment', variable: 'assignment_method', variable_value: assignmentMethod });
    configData.push({ tab_name: 'CellGuideAssignment', variable: 'THRESHOLD', variable_value: threshold });

    // Collect data from Perturbation Inference Tab
    var inferenceMethod = document.getElementById('inference_method').value;
    var moi = document.getElementById('moi').value
    var side = document.getElementById('side').value;
    var grnaIntegrationStrategy = document.getElementById('grna_integration_strategy').value;
    var resamplingApproximation = document.getElementById('resampling_approximation').value;
    var controlGroup = document.getElementById('control_group').value;
    var resamplingMechanism = document.getElementById('resampling_mechanism').value;
    var formulaObject = document.getElementById('formula_object').value;

    configData.push({ tab_name: 'PerturbationInference', variable: 'inference_method', variable_value: inferenceMethod });
    configData.push({ tab_name: 'PerturbationInference', variable: 'moi', variable_value: moi });
    configData.push({ tab_name: 'PerturbationInference', variable: 'side', variable_value: side });
    configData.push({ tab_name: 'PerturbationInference', variable: 'grna_integration_strategy', variable_value: grnaIntegrationStrategy });
    configData.push({ tab_name: 'PerturbationInference', variable: 'resampling_approximation', variable_value: resamplingApproximation });
    configData.push({ tab_name: 'PerturbationInference', variable: 'control_group', variable_value: controlGroup });
    configData.push({ tab_name: 'PerturbationInference', variable: 'resampling_mechanism', variable_value: resamplingMechanism });
    configData.push({ tab_name: 'PerturbationInference', variable: 'formula_object', variable_value: formulaObject });

    // Collect data from Inference Pairs Tab
    var inferenceOption = document.querySelector('input[name="inference_option"]:checked').value;
    configData.push({ tab_name: 'InferencePairs', variable: 'inference_option', variable_value: inferenceOption });

    var distanceFromCenter = document.getElementById('distance_from_center').value;
        configData.push({ tab_name: 'InferencePairs', variable: 'distance_from_center', variable_value: distanceFromCenter });

    // Collect data from Quality Filters Tab
    var minGenes = document.getElementById('min_genes').value;
    var minCells = document.getElementById('min_cells').value;
    var pctMito = document.getElementById('pct_mito').value;

    configData.push({ tab_name: 'QualityFilters', variable: 'min_genes', variable_value: minGenes });
    configData.push({ tab_name: 'QualityFilters', variable: 'min_cells', variable_value: minCells });
    configData.push({ tab_name: 'QualityFilters', variable: 'pct_mito', variable_value: pctMito });

    // Collect data from Evaluation Tab
    var userNodes = document.getElementById('user_central_nodes').value;
    var numNodes = document.getElementById('central_nodes_num').value;

    configData.push({ tab_name: 'Evaluation', variable: 'user_central_nodes', variable_value: userNodes });
    configData.push({ tab_name: 'Evaluation', variable: 'central_nodes_num', variable_value: numNodes });

    // Generate CSV Content
    var csvContent = "tab_name,variable,variable_value,batch_name,read1,read2\n";
    configData.forEach(function(item) {
        csvContent += '"' + item.tab_name + '","' + item.variable + '","' + (item.variable_value || '') + '","' + (item.batch_name || '') + '","' + (item.read1 || '') + '","' + (item.read2 || '') + '"\n';
    });

    // Download CSV
    var blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    var link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.setAttribute('download', 'configuration.csv');
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}