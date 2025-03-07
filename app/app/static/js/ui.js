// ui.js
import { loadImageDataApi, startAnalysis, checkTaskStatus } from "./api.js";
import {
    getAnalysisTypeWithoutSuffix,
    getBrowserAnalysisTypeWithoutSuffix,
    isNumeric,
    constructImageUrl,
} from "./utils.js";

let previousPlantSelection = [];

/**
 * Populate the file type dropdown with available file extensions.
 * @param {string} baseUrl The base URL.
 */
export function populateFileTypeDropdown(baseUrl) {
    const sessionId = sessionStorage.getItem("session_id");
    fetch(`${baseUrl}/api/list_file_extensions`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: sessionId }),
    })
        .then((response) => response.json())
        .then((response) => {
            if (response.status === "success") {
                const selectElement = document.getElementById("fileTypeSelect");
                if (!selectElement) return;
                selectElement.innerHTML = "";
                // Create the "All Results" option
                const allOption = document.createElement("option");
                allOption.value = "all";
                allOption.text = "All Results";
                selectElement.appendChild(allOption);
                // Add unique extensions
                const uniqueExtensions = new Set(response.extensions);
                uniqueExtensions.forEach((extension) => {
                    const option = document.createElement("option");
                    option.value = extension;
                    option.text = extension;
                    selectElement.appendChild(option);
                });
                // If ".html" is present, select it by default
                if (uniqueExtensions.has(".html")) {
                    selectElement.value = ".html";
                    updateFileListBasedOnExtension(baseUrl);
                } else {
                    selectElement.selectedIndex = 0;
                }
            }
        })
        .catch((error) => {
            console.error("Error fetching file extensions:", error);
        });
}

/**
 * Update the file list UI based on the selected extension.
 * @param {string} baseUrl The base URL.
 */
export function updateFileListBasedOnExtension(baseUrl) {
    const selectElement = document.getElementById("fileTypeSelect");
    if (!selectElement) return;
    const selectedExtension = selectElement.value;
    const sessionId = sessionStorage.getItem("session_id");
    fetch(`${baseUrl}/api/list_filtered_files`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
            session_id: sessionId,
            extension: selectedExtension,
        }),
    })
        .then((response) => response.json())
        .then((response) => {
            const analysisListEl = document.getElementById("analysisList");
            if (!analysisListEl) return;
            if (response.status === "success") {
                let fileListHtml = '<ul class="list-group">';
                response.files.forEach((file) => {
                    const filePath = `${BASE_URL}/${OUTPUT_DIR}/${sessionId}/${file}`;
                    fileListHtml += `<li class="list-group-item"><a href="${filePath}" target="_blank">${file}</a></li>`;
                });
                fileListHtml += "</ul>";
                analysisListEl.innerHTML = fileListHtml;
                $("#downloadResultsButton").prop("disabled", false);
                $("#clearFilesButton").prop("disabled", false);
            } else {
                analysisListEl.innerHTML = "<p>No files found.</p>";
            }
        })
        .catch((error) => {
            console.error("Error fetching filtered files:", error);
            const analysisListEl = document.getElementById("analysisList");
            if (analysisListEl) {
                analysisListEl.innerHTML = "<p>Error retrieving files.</p>";
            }
        });
}  

/**
 * Load image data and display images.
 * @param {string} baseUrl The base URL.
 * @param {string} plant The plant identifier.
 * @param {string} category The image category.
 */
export function loadImageDataUI(baseUrl, plant, category) {
    loadImageDataApi(baseUrl, plant, category)
        .then((response) => {
            if (response.status === "success") {
                displayImages(response.figures);
            } else {
                alert("Failed to load images: " + response.message);
            }
        })
        .catch((error) => {
            alert("Error loading images: " + error);
        });
}

/**
 * Display images in the gallery.
 * @param {Array<Object>} figures Array of image data.
 */
export function displayImages(figures) {
    const gallery = $("#imageDiv");
    gallery.empty();
    figures.sort((a, b) => {
        if (!a.original || !b.original) {
            console.error("Undefined original property in:", a, b);
            return 0;
        }
        return a.original.localeCompare(b.original);
    });
    figures.forEach((imageData) => {
        const thumb = $("<a>")
            .attr({
                href: imageData.original,
                "data-fancybox": "gallery",
                "data-caption":
                    "Click the right or left side of the image to move to the next or previous image.",
            })
            .append(
                $("<img>")
                    .addClass("img-fluid m-2 lightbox-image")
                    .attr("src", imageData.thumbnail)
                    .css({ width: "100%", height: "auto", display: "block" })
            );
        gallery.append(thumb);
    });
    adjustThumbnailSizeUI();
    gallery.show();
}

/**
 * Adjust thumbnail sizes based on window width.
 */
export function adjustThumbnailSizeUI() {
    const gallery = $("#imageDiv");
    const thumbnails = $("#imageDiv a");
    let thumbnailWidth;
    if (window.innerWidth < 576) {
        thumbnailWidth = 55;
    } else if (window.innerWidth < 768) {
        thumbnailWidth = 50;
    } else if (window.innerWidth < 992) {
        thumbnailWidth = 33.33;
    } else {
        thumbnailWidth = 25;
    }
    thumbnails.css("width", thumbnailWidth + "%");
}

/**
 * Initialize the browser iframe.
 * @param {string} baseUrl The base URL.
 */
export function startBrowserUI(baseUrl) {
    const initialSpeciesFilter = $("#plant option:selected")
        .map(function () {
            return $(this).text();
        })
        .get();
    $("#loadingBarBrowser").show();
    $("#browserResults").hide();
    $("#startBrowser").prop("disabled", true);
    const speciesFilterParam = encodeURIComponent(
        JSON.stringify(initialSpeciesFilter)
    );
    $("#browserIframe").attr(
        "src",
        `${baseUrl}/api/browser?species_filter=${speciesFilterParam}`
    );
    $("#browserIframe").on("load", function () {
        $("#loadingBarBrowser").hide();
        $("#browserResults").fadeIn();
        $("#browserAnalysisCard").fadeIn();
        $("#startBrowser").prop("disabled", false);
    });
}

/**
 * Initialize the info iframe.
 * @param {string} baseUrl The base URL.
 */
export function startInfoUI(baseUrl) {
    const selectedPlant = $("#plant").val();
    const initialSpeciesFilter = $("#plant option:selected")
        .map(function () {
            return $(this).text();
        })
        .get();
    $("#loadingBarInfo").show();
    $("#infoDiv").hide();
    const speciesFilterParam = encodeURIComponent(
        JSON.stringify(initialSpeciesFilter)
    );
    $("#infoIframe").attr(
        "src",
        `${baseUrl}/api/general_info?species_filter=${speciesFilterParam}`
    );
    if (selectedPlant.length === 1) {
        $("#imageCategoryDiv").fadeIn();
    }
}

/**
 * Handle plant selection changes and update UI only if the selection is new.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function handleSelectionChange(baseUrl, sessionId) {
    // Get the current selection as an array of selected values
    const currentSelection = $("#plant option:selected")
        .map(function () {
            return $(this).val();
        })
        .get();

    // Compare with the previous selection
    if (
        previousPlantSelection.length === currentSelection.length &&
        previousPlantSelection.join(",") === currentSelection.join(",")
    ) {
        // The same selection was made – no action needed.
        $("#plant").prop("disabled", false);
        return;
    }
    // Save the current selection as the new reference
    previousPlantSelection = currentSelection;

    // If at least one entry is selected:
    if (currentSelection.length > 0) {
        selectPlants();
        $("#browserAnalysisCard, #browserAnalysisResults, #browserResults, #results").hide();
        startInfoUI(baseUrl);
        startBrowserUI(baseUrl);
    } else {
        $("#infoDiv, #settingsCard, #downstreamCard, #browserCard, #results, #browserResults, #browserAnalysisCard, #browserAnalysisResults, #imageCategoryDiv").hide();
    }

    if (currentSelection.length > 1) {
        $("#imageCategoryDiv").hide();
        $("#plotTypePlotly").prop("disabled", true);
        $("#plotTypeCytoscape").prop("checked", true);
    } else {
        $("#plotTypePlotly").prop("disabled", false);
    }

    $("#plant").prop("disabled", false);
    updatePlantSelection();
}

/**
 * Update plant selection and load images if only one plant is selected.
 */
export function updatePlantSelection() {
    let selectedPlant = $("#plant").val();
    const gallery = $("#results");
    if (typeof selectedPlant === "string") {
        selectedPlant = [selectedPlant];
    }
    console.log("Plant change event", selectedPlant);
    gallery.empty();
    if (selectedPlant.length === 1) {
        const selectedCategory = $("#imageCategory").val();
        loadImageDataUI(BASE_URL, selectedPlant[0], selectedCategory);
    }
}

/**
 * Update analysis input fields based on the selected analysis type.
 */
export function updateInputFields() {
    const shapeLabel = document.getElementById("useShapesLabel");
    const selectedAnalysis = getAnalysisTypeWithoutSuffix(
        document.getElementById("analysisType").value
    );
    let hasParameters = false;
    $(".general-analysis-description").hide();
    $("#generalHelpIcon")
        .show()
        .off("click")
        .on("click", function () {
            $("#description_" + selectedAnalysis).toggle();
        });
    $("#nTopGroup, #nTopPercentGroup, #thresholdGroup, #positiveGroup, #maxPvalGroup, #minFeGroup, #minDepthGroup, #useShapesGroup, #useColorsGroup, #interSpeciesOnlyGroup").hide();
    shapeLabel.innerHTML =
        selectedAnalysis === "plot_correlation_network"
            ? "Use Shapes for Modules"
            : "Use Shapes for Species";
    switch (selectedAnalysis) {
        case "plot_co_expression_network":
            $("#thresholdGroup").show();
            hasParameters = true;
            break;
        case "highest_expr_genes":
            $("#nTopGroup").show();
            hasParameters = true;
            break;
        case "plot_correlation_network":
            $("#thresholdGroup, #useShapesGroup, #useColorsGroup").show();
            hasParameters = true;
            break;
        case "plot_species_correlation_network":
            $("#thresholdGroup, #maxPvalGroup, #positiveGroup, #useShapesGroup, #useColorsGroup, #interSpeciesOnlyGroup").show();
            hasParameters = true;
            break;
        case "plot_module_goea":
            $("#maxPvalGroup, #minFeGroup, #minDepthGroup, #nTopPercentGroup").show();
            hasParameters = true;
            break;
        case "plot_jaccard_tissues":
            $("#nTopGroup, #thresholdGroup, #useShapesGroup, #interSpeciesOnlyGroup").show();
            hasParameters = true;
            break;
        case "plot_tissues_corr":
            $("#nTopGroup").show();
            hasParameters = true;
            break;
        case "plot_jaccard_modules":
            $("#thresholdGroup, #useShapesGroup, #useColorsGroup, #interSpeciesOnlyGroup").show();
            hasParameters = true;
            break;
        default:
            hasParameters = false;
    }
    if (hasParameters) {
        $("#wholeParametersToggle").show();
    } else {
        $("#wholeParametersToggle").hide();
    }
}

/**
 * Update browser analysis input fields based on the selected browser analysis type.
 */
export function updateBrowserInputFields() {
    const selectedAnalysis = getBrowserAnalysisTypeWithoutSuffix(
        document.getElementById("browserAnalysisType").value
    );
    const browserButton = $("#runBrowserAnalysis");
    let hasParameters = false;
    $(".analysis-description").hide();
    $("#browserHelpIcon")
        .show()
        .off("click")
        .on("click", function () {
            $("#description_" + selectedAnalysis).toggle();
        });
    browserButton.prop("disabled", false);
    $("#nTopPercentGroupBrowser, #nTopGroupBrowser, #thresholdGroupBrowser, #maxPvalGroupBrowser, #minFeGroupBrowser, #minDepthGroupBrowser, #useShapesGroupBrowser, #useShapesSpeciesGroupBrowser, #useColorsGroupBrowser, #interSpeciesOnlyGroupBrowser, #minOrthosGroupBrowser, #highlightListGroupBrowser, #maxNeighborsGroupBrowser, #forceDetailedViewGroupBrowser, #minClusterSizeGroupBrowser").hide();
    switch (selectedAnalysis) {
        case "plot_co_expression_network":
            $("#thresholdGroupBrowser, #highlightListGroupBrowser, #useShapesSpeciesGroupBrowser, #useColorsGroupBrowser, #maxNeighborsGroupBrowser, #forceDetailedViewGroupBrowser, #minClusterSizeGroupBrowser").show();
            hasParameters = true;
            break;
        case "plot_go_terms":
            $("#maxPvalGroupBrowser, #minFeGroupBrowser, #minDepthGroupBrowser, #nTopPercentGroupBrowser").show();
            hasParameters = true;
            break;
        case "plot_filtered_jaccard_modules":
            $("#thresholdGroupBrowser, #useShapesSpeciesGroupBrowser, #useColorsGroupBrowser, #interSpeciesOnlyGroupBrowser, #minOrthosGroupBrowser").show();
            hasParameters = true;
            break;
        case "plot_filtered_jaccard_species":
            $("#thresholdGroupBrowser, #useShapesSpeciesGroupBrowser, #minOrthosGroupBrowser").show();
            hasParameters = true;
            break;
        case "not_implemented":
            hasParameters = false;
            console.log("Analysis not implemented yet");
            browserButton.prop("disabled", true);
            break;
        default:
            hasParameters = false;
    }
    if (hasParameters) {
        $("#browserParametersToggle").show();
    } else {
        $("#browserParametersToggle").hide();
    }
}

/**
 * Handle messages received from iframes.
 * @param {MessageEvent} event The message event.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function handleMessageEvent(event, baseUrl, sessionId) {
    if (event.data.action === "setImageList") {
        window.imageList = event.data.images;
    } else if (event.data.action === "openFancybox") {
        const imageName = event.data.imageName;
        const imageList = window.imageList || [imageName];
        let startIndex = imageList.indexOf(imageName);
        if (startIndex < 0) startIndex = 0;
        const items = imageList.map((name) =>
            ({ src: constructImageUrl(name, BASE_URL, OUTPUT_DIR, sessionId), type: "image" })
        );
        console.log("Received imageName:", imageName);
        console.log("Constructed URLs:", items);
        Fancybox.show(items, { startIndex: startIndex });
    } else if (event.data === "tableLoaded") {
        $("#infoDiv").fadeIn();
        $("#loadingBarInfo").hide();
    } else if (event.data.type === "tableData") {
        console.log("Table Data for further analysis:", event.data.data);
    } else if (
        event.data.type === "resizeIframe" ||
        event.data.type === "resizeBrowser"
    ) {
        const iframe = document.getElementById(event.data.target);
        if (iframe) {
            iframe.style.height = parseInt(event.data.height) + "px";
        }
    } else if (event.data.type === "resizeBrowserAnalysis") {
        const iframe = document.getElementById(event.data.target);
        if (iframe) {
            iframe.style.height = parseInt(event.data.height + 30) + "px";
        }
    }
}

/**
 * Send an analysis request and poll for status.
 * @param {string} baseUrl The base URL.
 * @param {Object} requestData The request data.
 * @param {string} type The analysis type ("browser" or "general").
 * @param {boolean} text Flag for text mode.
 */
export function sendAnalysisRequestUI(baseUrl, requestData, type, text) {
    const loadingBar = type === "browser"
        ? (text ? "#loadingText" : "#loadingBar3")
        : "#loadingBar4";
    const resultsDiv = type === "browser" ? "#browserAnalysisResults" : "#results";
    $(loadingBar).show();
    requestData.analysis_type = requestData.analysis_type || "";
    startAnalysis(baseUrl, requestData)
        .then((taskId) => {
            checkStatusUI(baseUrl, taskId, type, text);
        })
        .catch((error) => {
            console.error("Error starting analysis:", error);
            $(resultsDiv).html("<p>Error starting analysis: " + error + "</p>");
            $(loadingBar).hide();
        });
}

/**
 * Poll the analysis task status.
 * @param {string} baseUrl The base URL.
 * @param {string} taskId The task ID.
 * @param {string} type The analysis type.
 * @param {boolean} text Flag for text mode.
 */
export function checkStatusUI(baseUrl, taskId, type, text) {
    checkTaskStatus(baseUrl, taskId)
        .then((result) => {
            let resultsDiv, loadingBar, loadingText = "#loadingStatusText";
            if (type === "browser") {
                resultsDiv = "#browserAnalysisResults";
                loadingBar = text ? "#loadingText" : "#loadingBar3";
            } else {
                resultsDiv = "#results";
                loadingBar = "#loadingBar4";
            }
            if (result.state === "PENDING" || result.state === "STARTED") {
                setTimeout(() => checkStatusUI(baseUrl, taskId, type, text), 5000);
            } else if (result.state === "PROGRESS") {
                $(loadingText).text(result.status);
                setTimeout(() => checkStatusUI(baseUrl, taskId, type, text), 2000);
            } else if (result.state === "SUCCESS" && result.result) {
                console.log("Analysis result:", result.result);
                displayResultUI(result.result, type, text);
                $(loadingBar).hide();
            } else if (result.state === "FAILURE") {
                console.error("Task failed:", result.result.status);
                $(resultsDiv).html("<p>Error: " + result.result.message + "</p>");
                $(loadingBar).hide();
            } else {
                console.error("Unexpected task state:", result.state);
                $(resultsDiv).html("<p>Unexpected task state: " + result.state + "</p>");
                $(loadingBar).hide();
            }
        })
        .catch((error) => {
            console.error("Error checking task status:", error);
            const loadingBar = type === "browser" ? "#loadingBar3" : "#loadingBar4";
            $(loadingBar).hide();
        });
}

/**
 * Display analysis results.
 * @param {Object} result The analysis result.
 * @param {string} type The analysis type.
 * @param {boolean} text Flag for text mode.
 */
export function displayResultUI(result, type, text) {
    const timestamp = new Date().getTime();
    const newWindow = $("#displayInNewWindow").is(":checked");
    let resultsDiv, frameID, loadingBar, startButton;
    if (type === "browser") {
        resultsDiv = "#browserAnalysisResults";
        frameID = "browserResultsIframe";
        loadingBar = text ? "#loadingText" : "#loadingBar3";
        startButton = "#runBrowserAnalysis";
    } else {
        resultsDiv = "#results";
        frameID = "resultsIframe";
        loadingBar = "#loadingBar4";
        startButton = "#runAnalysis";
    }
    if (newWindow) {
        if (result.plot_url) {
            window.open(result.plot_url + "?t=" + timestamp, "_blank");
        } else if (result.message) {
            alert(result.message);
        } else {
            console.error("Invalid plot URL:", result.plot_url);
            alert("Error: Invalid plot URL");
        }
    } else {
        if (result.plot_url) {
            if (result.plot_url.endsWith(".html")) {
                console.log(
                    "Setting iFrame with URL:",
                    result.plot_url + "?t=" + timestamp
                );
                $(resultsDiv).html(
                    '<iframe id="' +
                    frameID +
                    '" scrolling="no" src="' +
                    result.plot_url +
                    "?t=" +
                    timestamp +
                    '" width="100%" height="600px" frameborder="0"></iframe>'
                );
            } else {
                console.log(
                    "Setting image with URL:",
                    result.plot_url + "?t=" + timestamp
                );
                $(resultsDiv).html(
                    '<img src="' +
                    result.plot_url +
                    "?t=" +
                    timestamp +
                    '" alt="Analysis Plot" width="100%" height="100%" />'
                );
            }
        } else if (result.message) {
            $(resultsDiv).html("<p>" + result.message + "</p>");
        } else {
            console.error("Invalid plot URL:", result.plot_url);
            $(resultsDiv).html("<p>Error: Invalid plot URL</p>");
        }
        $(resultsDiv).show();
    }

    if (sessionStorage.getItem("session_id")) {
        populateFileTypeDropdown(BASE_URL);
        updateFileListBasedOnExtension(BASE_URL);
    }

    $(loadingBar).hide();
    if (startButton) {
        $(startButton).prop("disabled", false);
    }
}

/**
 * Handle plant selection: show/hide sections and update options.
 */
export function selectPlants() {
    let plant = Array.from(
        document.querySelectorAll("#plant option:checked")
    ).map((option) => option.value);
    console.log("Plant:", plant);
    const singlePlantSelected = plant.length === 1;
    if (singlePlantSelected) {
        plant = plant[0];
    }
    console.log("Plant:", plant);
    console.log("selectPlants executed");
    $("#infoDiv, #settingsCard, #downstreamCard, #browserCard").show();
    adjustAnalysisOptions(singlePlantSelected);
    adjustBrowserOptions(singlePlantSelected, Array.isArray(plant) ? plant : [plant]);
    updateInputFields();
    updateBrowserInputFields();
}

/**
 * Adjust analysis options based on selection.
 * @param {boolean} single True if a single plant is selected.
 */
export function adjustAnalysisOptions(single) {
    const select = document.getElementById("analysisType");
    const filterType = single ? "single" : "multiple";
    select.innerHTML = "";
    window.originalOptions.forEach((option) => {
        if (
            option.value.endsWith(`_${filterType}`) ||
            option.value.endsWith("_common")
        ) {
            select.appendChild(option.cloneNode(true));
        }
    });
    if (select.options.length > 0) {
        select.value = select.options[0].value;
    }
}

/**
 * Adjust browser options based on selection.
 * @param {boolean} single True if a single plant is selected.
 * @param {Array} plants The selected plant(s).
 */
export function adjustBrowserOptions(single, plants) {
    const select = document.getElementById("browserAnalysisType");
    const filterType = single ? "single" : "multiple";
    select.innerHTML = "";
    window.originalBrowserOptions.forEach((option) => {
        const isVennOption = option.value.startsWith("plot_venn");
        const isFilterMatch =
            option.value.endsWith(`_${filterType}`) || option.value.endsWith("_common");
        if ((isFilterMatch && !isVennOption) || (isVennOption && plants.length === 3)) {
            select.appendChild(option.cloneNode(true));
        }
    });
    if (select.options.length > 0) {
        select.value = select.options[0].value;
    }
}

/**
 * Handle Run Browser Analysis button click.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function runBrowserAnalysisHandler(baseUrl, sessionId) {
    let plant = $("#plant").val();
    let analysisType = getBrowserAnalysisTypeWithoutSuffix(
        document.getElementById("browserAnalysisType").value
    );
    let analysisName = $("#browserAnalysisType option:selected").text();
    let nTopPercent = $("#nTopPercentBrowser").val();
    let nTop = $("#nTopBrowser").val();
    let threshold = $("#thresholdBrowser").val();
    let maxPval = $("#maxPvalBrowser").val();
    let maxNeighbors = $("#maxNeighborsBrowser").val();
    let minFe = $("#minFeBrowser").val();
    let minDepth = $("#minDepthBrowser").val();
    let useShapes = $("#useShapesBrowser").is(":checked");
    let useColors = $("#useColorsBrowser").is(":checked");
    let highlightList = $("#highlightListBrowser")
        .val()
        .split(/[\s,]+/)
        .map((item) => item.trim())
        .filter((item) => item.length > 0);
    let useShapesSpecies = $("#useShapesSpeciesBrowser").is(":checked");
    let interSpeciesOnly = $("#interSpeciesOnlyBrowser").is(":checked");
    let minOrthos = $("#minOrthosBrowser").val();
    let minClusterSize = $("#minClusterSizeBrowser").val();
    let prefix = $("#outputPrefix").val();
    let selectedPlotType = $('input[name="plotType"]:checked').val();
    let forceDetailedView = $("#forceDetailedViewBrowser").is(":checked");
    let text = false;
    if (prefix !== "" && !prefix.endsWith("_")) {
        prefix += "_";
    }
    if (plant.length === 1) {
        plant = plant[0];
    }
    const requestData = { plant, analysisName, session_id: sessionId, prefix };
    let isValid = true;
    let errorMessage = "";
    if (!plant || (Array.isArray(plant) && plant.length === 0)) {
        isValid = false;
        errorMessage += "Please select at least one plant.\n";
    }
    if (nTop === "" || !isNumeric(nTop) || parseInt(nTop) <= 0) {
        isValid = false;
        errorMessage += "Top n must be a positive number.\n";
    }
    if (
        nTopPercent === "" ||
        !isNumeric(nTopPercent) ||
        parseFloat(nTopPercent) < 0.1 ||
        parseFloat(nTopPercent) > 100
    ) {
        isValid = false;
        errorMessage += "Top percentage must be between 0.1 and 100.\n";
    }
    if (
        threshold === "" ||
        !isNumeric(threshold) ||
        parseFloat(threshold) < 0 ||
        parseFloat(threshold) > 1
    ) {
        isValid = false;
        errorMessage += "Threshold must be between 0 and 1.\n";
    }
    if (
        maxPval === "" ||
        !isNumeric(maxPval) ||
        parseFloat(maxPval) < 0 ||
        parseFloat(maxPval) > 1
    ) {
        isValid = false;
        errorMessage += "Max P-Value must be between 0 and 1.\n";
    }
    if (
        minFe === "" ||
        !isNumeric(minFe) ||
        parseFloat(minFe) < 0 ||
        parseFloat(minFe) > 1
    ) {
        isValid = false;
        errorMessage += "Min Fold Enrichment must be between 0 and 1.\n";
    }
    if (minDepth === "" || !isNumeric(minDepth) || parseInt(minDepth) < 1) {
        isValid = false;
        errorMessage += "Min Depth must be at least 1.\n";
    }
    if (maxNeighbors === "" || !isNumeric(maxNeighbors) || parseInt(maxNeighbors) < 0 || parseInt(maxNeighbors) > 50) {
        isValid = false;
        errorMessage += "Max Neighbors must be between 0 and 50.\n";
    }
    if (minClusterSize === "" || !isNumeric(minClusterSize) || parseInt(minClusterSize) < 2) {
        isValid = false;
        errorMessage += "Min Cluster Size must be at least 2.\n";
    }
    if (!isValid) {
        alert(errorMessage);
        return;
    }
    const iframe = document.getElementById("browserIframe");
    if (iframe) {
        iframe.contentWindow.postMessage({ type: "getTableParams" }, "*");
        window.addEventListener(
            "message",
            function handler(event) {
                if (event.data.type === "tableParams") {
                    const params = event.data.data;
                    fetch(`${baseUrl}/api/get_transcripts`, {
                        method: "POST",
                        headers: { "Content-Type": "application/json" },
                        body: JSON.stringify(params),
                    })
                        .then((response) => {
                            if (!response.ok)
                                throw new Error(`HTTP error! status: ${response.status}`);
                            return response.json();
                        })
                        .then((data) => {
                            const transcripts = data.transcripts;
                            requestData.transcripts = transcripts;
                            requestData.prefix = prefix;
                            if (analysisType === "plot_co_expression_network") {
                                requestData.selectedPlotType = selectedPlotType;
                                requestData.threshold = parseFloat(threshold);
                                requestData.highlightList = highlightList;
                                requestData.useShapesSpecies = useShapesSpecies;
                                requestData.useColors = useColors;
                                requestData.maxNeighbors = parseInt(maxNeighbors);
                                requestData.forceDetailedView = forceDetailedView;
                                requestData.minClusterSize = parseInt(minClusterSize);
                                text = true;
                            } else if (analysisType === "plot_go_terms") {
                                requestData.maxPval = parseFloat(maxPval);
                                requestData.minFe = parseFloat(minFe);
                                requestData.minDepth = parseInt(minDepth);
                                requestData.nTopPercent = parseFloat(nTopPercent);
                            } else if (analysisType === "plot_venn") {
                                console.log("Plotting Venn diagram");
                            } else if (analysisType === "plot_upset") {
                                console.log("Plotting UpSet plot");
                            } else if (analysisType === "plot_filtered_jaccard_modules") {
                                requestData.threshold = parseFloat(threshold);
                                requestData.useShapesSpecies = useShapesSpecies;
                                requestData.useColors = useColors;
                                requestData.interSpeciesOnly = interSpeciesOnly;
                                requestData.minOrthos = parseInt(minOrthos);
                            } else if (analysisType === "plot_filtered_jaccard_species") {
                                requestData.threshold = parseFloat(threshold);
                                requestData.useShapesSpecies = useShapesSpecies;
                                requestData.minOrthos = parseInt(minOrthos);
                            } else if (analysisType === "not_implemented") {
                                alert("Analysis not implemented yet.");
                                return;
                            }

                            requestData.analysis_type = analysisType;

                            sendAnalysisRequestUI(baseUrl, requestData, "browser", text);
                            $("#runBrowserAnalysis").prop("disabled", true);
                        })
                        .catch((error) => {
                            console.error("Error retrieving transcripts:", error);
                            $("#loadingBar3").hide();
                        });
                    window.removeEventListener("message", handler);
                } else {
                    $("#loadingBar3").hide();
                }
            },
            { once: true }
        );
    } else {
        $("#loadingBar3").hide();
    }
}

/**
 * Handle Run General Analysis button click.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function runGeneralAnalysisHandler(baseUrl, sessionId) {
    let plant = $("#plant").val();
    let analysisType = getAnalysisTypeWithoutSuffix(
        document.getElementById("analysisType").value
    );
    let analysisName = $("#analysisType option:selected").text();
    let nTop = $("#nTop").val();
    let nTopPercent = $("#nTopPercent").val();
    let threshold = $("#threshold").val();
    let maxPval = $("#maxPval").val();
    let minFe = $("#minFe").val();
    let minDepth = $("#minDepth").val();
    let positive = $("#positive").is(":checked");
    let useShapes = $("#useShapes").is(":checked");
    let useColors = $("#useColors").is(":checked");
    let interSpeciesOnly = $("#interSpeciesOnly").is(":checked");
    let prefix = $("#outputPrefix").val();
    let text = false;
    if (prefix !== "" && !prefix.endsWith("_")) {
        prefix += "_";
    }
    if (plant.length === 1) {
        plant = plant[0];
    }
    const requestData = { plant, analysisName, session_id: sessionId, prefix };
    let isValid = true;
    let errorMessage = "";
    if (!plant || (Array.isArray(plant) && plant.length === 0)) {
        isValid = false;
        errorMessage += "Please select at least one plant.\n";
    }
    if (nTop === "" || !isNumeric(nTop) || parseInt(nTop) <= 0) {
        isValid = false;
        errorMessage += "Top n must be a positive number.\n";
    }
    if (
        nTopPercent === "" ||
        !isNumeric(nTopPercent) ||
        parseFloat(nTopPercent) < 0.1 ||
        parseFloat(nTopPercent) > 100
    ) {
        isValid = false;
        errorMessage += "Top percentage must be between 0.1 and 100.\n";
    }
    if (
        threshold === "" ||
        !isNumeric(threshold) ||
        parseFloat(threshold) < 0 ||
        parseFloat(threshold) > 1
    ) {
        isValid = false;
        errorMessage += "Threshold must be between 0 and 1.\n";
    }
    if (
        maxPval === "" ||
        !isNumeric(maxPval) ||
        parseFloat(maxPval) < 0 ||
        parseFloat(maxPval) > 1
    ) {
        isValid = false;
        errorMessage += "Max P-Value must be between 0 and 1.\n";
    }
    if (
        minFe === "" ||
        !isNumeric(minFe) ||
        parseFloat(minFe) < 0 ||
        parseFloat(minFe) > 1
    ) {
        isValid = false;
        errorMessage += "Min Fold Enrichment must be between 0 and 1.\n";
    }
    if (minDepth === "" || !isNumeric(minDepth) || parseInt(minDepth) < 1) {
        isValid = false;
        errorMessage += "Min Depth must be at least 1.\n";
    }
    if (!isValid) {
        alert(errorMessage);
        return;
    }
    if (
        analysisType === "plot_correlation_network" ||
        analysisType === "plot_co_expression_network"
    ) {
        requestData.threshold = parseFloat(threshold);
        requestData.useShapes = useShapes;
        requestData.useColors = useColors;
    } else if (analysisType === "plot_species_correlation_network") {
        requestData.threshold = parseFloat(threshold);
        requestData.maxPval = parseFloat(maxPval);
        requestData.positive = positive;
        requestData.useShapes = useShapes;
        requestData.useColors = useColors;
        requestData.interSpeciesOnly = interSpeciesOnly;
    } else if (analysisType === "plot_module_goea") {
        requestData.maxPval = parseFloat(maxPval);
        requestData.minFe = parseFloat(minFe);
        requestData.minDepth = parseInt(minDepth);
        requestData.nTopPercent = parseFloat(nTopPercent);
    } else if (analysisType === "highest_expr_genes") {
        requestData.nTop = parseInt(nTop);
    } else if (analysisType === "plot_jaccard_tissues") {
        requestData.nTop = parseInt(nTop);
        requestData.threshold = parseFloat(threshold);
        requestData.useShapes = useShapes;
        requestData.interSpeciesOnly = interSpeciesOnly;
    } else if (analysisType === "plot_tissues_corr") {
        requestData.nTop = parseInt(nTop);
    } else if (analysisType === "plot_jaccard_modules") {
        requestData.threshold = parseFloat(threshold);
        requestData.useShapes = useShapes;
        requestData.useColors = useColors;
        requestData.interSpeciesOnly = interSpeciesOnly;
    }
    requestData.analysis_type = analysisType;
    sendAnalysisRequestUI(baseUrl, requestData, "general", text);
    $("#runAnalysis").prop("disabled", true);
}
