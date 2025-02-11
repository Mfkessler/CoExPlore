console.log(`Using BASE_URL: ${BASE_URL}`);
console.log(`Using OUTPUT_DIR: ${OUTPUT_DIR}`);

$(document).ready(function () {
    // Initialize the select elements
    const select = document.getElementById("analysisType");
    originalOptions = Array.from(select.options);
    originalBrowserOptions = Array.from(
        document.getElementById("browserAnalysisType").options
    );

    function sendHeartbeat() {
        const sessionId = sessionStorage.getItem("session_id");
        fetch(BASE_URL + "/api/heartbeat", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ session_id: sessionId }),
        })
            .then((response) => {
                return response.text();
            })
            .then((text) => {
                console.log("Heartbeat response:", text);
            })
            .catch((err) => {
                console.error("Error sending heartbeat:", err);
            });
    }

    if (!sessionStorage.getItem("session_id")) {
        sessionStorage.setItem("session_id", generateUUID());
    }
    const session_id = sessionStorage.getItem("session_id");

    console.log("Initial session ID:", session_id);
    sendHeartbeat();

    fetch(BASE_URL + "/api/init_session", {
        method: "POST",
        headers: {
            "Content-Type": "application/json",
            "Session-ID": session_id,
        },
    })
        .then((response) => response.json())
        .then((data) => {
            console.log("Success:", data);
        })
        .catch((error) => {
            console.error("Error:", error);
        });

    // Send heartbeat every 5 minutes
    setInterval(sendHeartbeat, 300000);

    function generateUUID() {
        // Public Domain/MIT
        var d = new Date().getTime(); //Timestamp
        var d2 = (performance && performance.now && performance.now() * 1000) || 0; //Time in microseconds since page-load or 0 if unsupported
        return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(
            /[xy]/g,
            function (c) {
                var r = Math.random() * 16; // random number between 0 and 16
                if (d > 0) {
                    //Use timestamp until depleted
                    r = (d + r) % 16 | 0;
                    d = Math.floor(d / 16);
                } else {
                    //Use microseconds since page-load if supported
                    r = (d2 + r) % 16 | 0;
                    d2 = Math.floor(d2 / 16);
                }
                return (c === "x" ? r : (r & 0x3) | 0x8).toString(16);
            }
        );
    }

    function populateFileTypeDropdown() {
        const sessionId = sessionStorage.getItem("session_id");

        fetch(BASE_URL + "/api/list_file_extensions", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ session_id: sessionId }),
        })
            .then((response) => response.json())
            .then((response) => {
                if (response.status === "success") {
                    const selectElement = document.getElementById("fileTypeSelect");
                    selectElement.innerHTML = "";
                    let allOption = document.createElement("option");

                    allOption.value = "all";
                    allOption.text = "All Results";
                    selectElement.appendChild(allOption);

                    const uniqueExtensions = new Set(response.extensions);

                    uniqueExtensions.forEach((extension) => {
                        let option = document.createElement("option");
                        option.value = extension;
                        option.text = extension;
                        selectElement.appendChild(option);
                    });

                    // if "html" is in the list, select it by default
                    if (uniqueExtensions.has(".html")) {
                        selectElement.value = ".html";
                        updateFileListBasedOnExtension();
                    } else {
                        selectElement.selectedIndex = 0;
                    }
                }
            })
            .catch((error) => {
                console.error("Error fetching file extensions:", error);
            });
    }

    function updateFileListBasedOnExtension() {
        const selectedExtension = document.getElementById("fileTypeSelect").value;
        const sessionId = sessionStorage.getItem("session_id");

        fetch(BASE_URL + "/api/list_filtered_files", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({
                session_id: sessionId,
                extension: selectedExtension,
            }),
        })
            .then((response) => response.json())
            .then((response) => {
                if (response.status === "success") {
                    var fileListHtml = '<ul class="list-group">';
                    response.files.forEach(function (file) {
                        const filePath = `${BASE_URL}/${OUTPUT_DIR}/${sessionId}/${file}`;
                        fileListHtml += `<li class="list-group-item"><a href="${filePath}" target="_blank">${file}</a></li>`;
                    });

                    fileListHtml += "</ul>";
                    document.getElementById("analysisList").innerHTML = fileListHtml;

                    $("#downloadResultsButton").prop("disabled", false);
                    $("#clearFilesButton").prop("disabled", false);
                } else {
                    document.getElementById("analysisList").innerHTML =
                        "<p>No files found.</p>";
                }
            })
            .catch((error) => {
                console.error("Error fetching filtered files:", error);
                document.getElementById("analysisList").innerHTML =
                    "<p>Error retrieving files.</p>";
            });
    }

    document
        .getElementById("fileTypeSelect")
        .addEventListener("change", updateFileListBasedOnExtension);

    $("#helpIcon").click(function () {
        $("#userGuideCard").toggle();
    });

    $("#helpIconBrowser").click(function () {
        $("#userGuideCardBrowser").toggle();
    });

    $("#helpIconDataset").click(function () {
        $("#userGuideCardDataset").toggle();
    });

    function selectPlants() {
        var plant = Array.from(
            document.querySelectorAll("#plant option:checked")
        ).map((option) => option.value);

        console.log("Plant:", plant);

        // Determine if a single plant is selected
        var singlePlantSelected = plant.length === 1;

        if (singlePlantSelected) {
            plant = plant[0];
        }

        console.log("Plant:", plant);
        console.log("selectPlants executed");

        // Show the four divs
        $("#infoDiv").show();
        $("#settingsCard").show();
        $("#downstreamCard").show();
        $("#browserCard").show();

        // Pass the correct boolean value to the functions
        adjustAnalysisOptions(singlePlantSelected);
        adjustBrowserOptions(singlePlantSelected, plant);
        updateInputFields();
        updateBrowserInputFields();
    }

    function adjustAnalysisOptions(single) {
        const select = document.getElementById("analysisType");
        const filterType = single ? "single" : "multiple";
        select.innerHTML = "";

        originalOptions.forEach((option) => {
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

    function adjustBrowserOptions(single, plants) {
        const select = document.getElementById("browserAnalysisType");
        const filterType = single ? "single" : "multiple";
        select.innerHTML = "";

        originalBrowserOptions.forEach((option) => {
            const isVennOption = option.value.startsWith("plot_venn");
            const isFilterMatch =
                option.value.endsWith(`_${filterType}`) ||
                option.value.endsWith("_common");

            if (
                (isFilterMatch && !isVennOption) ||
                (isVennOption && plants.length === 3)
            ) {
                select.appendChild(option.cloneNode(true));
            }
        });

        if (select.options.length > 0) {
            select.value = select.options[0].value;
        }
    }

    function updatePlantSelection() {
        var selectedPlant = $("#plant").val(); // Get the value of the plant selection
        var gallery = $("#results");

        // Ensure that selectedPlant is always an array
        if (typeof selectedPlant === "string") {
            selectedPlant = [selectedPlant];
        }

        console.log("Plant change event", selectedPlant); // Log the selection
        gallery.empty(); // Clear the gallery

        // If only one plant is selected, load the image data
        if (selectedPlant.length === 1) {
            var selectedCategory = $("#imageCategory").val(); // Get the category from the corresponding element
            loadImageData(selectedPlant[0], selectedCategory); // Load data for the selected plant
        }
    }

    $("#imageCategory").change(function () {
        var selectedCategory = $(this).val();
        var selectedPlant = $("#plant").val();
        if (selectedPlant.length == 1) {
            selectedPlant = selectedPlant[0];
        }
        loadImageData(selectedPlant, selectedCategory);
    });

    function loadImageData(plant, category) {
        fetch(BASE_URL + "/api/load_image_data", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ plant: plant, category: category }),
        })
            .then((response) => response.json())
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

    function displayImages(figures) {
        var gallery = $("#imageDiv");
        gallery.empty();
        figures.sort((a, b) => {
            if (!a.original || !b.original) {
                console.error("Undefined original property in:", a, b);
                return 0;
            }
            return a.original.localeCompare(b.original);
        });

        figures.forEach(function (imageData) {
            var thumb = $("<a>")
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
                        .css({
                            width: "100%",
                            height: "auto",
                            display: "block",
                        })
                );

            gallery.append(thumb);
        });

        adjustThumbnailSize();
        gallery.show();
    }

    function adjustThumbnailSize() {
        var gallery = $("#imageDiv");
        var thumbnails = $("#imageDiv a");
        var galleryWidth = gallery.width();
        var thumbnailWidth;

        if (window.innerWidth < 576) {
            thumbnailWidth = 55; // 1 image per row for extra small screens
        } else if (window.innerWidth < 768) {
            thumbnailWidth = 50; // 2 images per row for small screens
        } else if (window.innerWidth < 992) {
            thumbnailWidth = 33.33; // 3 images per row for medium screens
        } else {
            thumbnailWidth = 25; // 4 images per row for large and above screens
        }

        thumbnails.css("width", thumbnailWidth + "%"); // Set width in percentage
    }

    $(window).resize(function () {
        adjustThumbnailSize();
    });

    $(window).on("resize", adjustThumbnailSize);

    // Start the browser table
    function startBrowser() {
        var initialSpeciesFilter = $("#plant option:selected")
            .map(function () {
                return $(this).text();
            })
            .get();

        $("#loadingBarBrowser").show();
        $("#browserResults").hide();
        $("#startBrowser").prop("disabled", true);

        var speciesFilterParam = encodeURIComponent(
            JSON.stringify(initialSpeciesFilter)
        );
        $("#browserIframe").attr(
            "src",
            BASE_URL + "/api/browser?species_filter=" + speciesFilterParam
        );

        $("#browserIframe").on("load", function () {
            $("#loadingBarBrowser").hide();
            $("#browserResults").fadeIn();
            $("#browserAnalysisCard").fadeIn();
            $("#startBrowser").prop("disabled", false);
        });
    }

    // Start the info table
    function startInfo() {
        var selectedPlant = $("#plant").val();
        var initialSpeciesFilter = $("#plant option:selected")
            .map(function () {
                return $(this).text();
            })
            .get();

        $("#loadingBarInfo").show();
        $("#infoDiv").hide();

        var speciesFilterParam = encodeURIComponent(
            JSON.stringify(initialSpeciesFilter)
        );
        $("#infoIframe").attr(
            "src",
            BASE_URL + "/api/general_info?species_filter=" + speciesFilterParam
        );

        if (selectedPlant.length === 1) {
            $("#imageCategoryDiv").fadeIn(); // Show the image category selection
        }
    }

    window.addEventListener("message", function (event) {
        if (event.origin !== window.location.origin) {
            return;
        }
        if (event.data.action === "setImageList") {
            window.imageList = event.data.images;
        } else if (event.data.action === "openFancybox") {
            var imageName = event.data.imageName;
            var imageList = window.imageList || [imageName];
            var startIndex = imageList.indexOf(imageName);
            if (startIndex < 0) startIndex = 0;

            var items = imageList.map(function (name) {
                var url = constructImageUrl(name);
                return { src: url, type: "image" };
            });
            console.log("Received imageName:", event.data.imageName);
            console.log("Constructed URLs:", items);

            Fancybox.show(items, { startIndex: startIndex });
        }
        if (event.data === "tableLoaded") {
            $("#infoDiv").fadeIn();
            $("#loadingBarInfo").hide();
        }
        if (event.data.type === "tableData") {
            var tableData = event.data.data;
            console.log("Table Data for further analysis:", tableData);
        }
        if (event.data.type === "resizeIframe") {
            var iframe = document.getElementById(event.data.target);
            if (iframe) {
                iframe.style.height = parseInt(event.data.height) + "px";
            }
        }
        if (event.data.type === "resizeBrowser") {
            var iframe = document.getElementById(event.data.target);
            if (iframe) {
                iframe.style.height = parseInt(event.data.height) + "px";
            }
        }
        if (event.data.type === "resizeBrowserAnalysis") {
            var iframe = document.getElementById(event.data.target);
            if (iframe) {
                iframe.style.height = parseInt(event.data.height + 30) + "px";
            }
        }
    });

    function constructImageUrl(name) {
        if (BASE_URL === "") {
            return "/" + OUTPUT_DIR + "/" + session_id + "/" + name;
        }
        if (!name.startsWith(BASE_URL)) {
            return BASE_URL + "/" + OUTPUT_DIR + "/" + session_id + "/" + name;
        }
        return name;
    }

    $("#runBrowserAnalysis").click(function () {
        var plant = $("#plant").val();
        var analysisType = getBrowserAnalysisTypeWithoutSuffix();
        var analysisName = $("#browserAnalysisType option:selected").text();
        var nTopPercent = $("#nTopPercentBrowser").val();
        var nTop = $("#nTopBrowser").val();
        var threshold = $("#thresholdBrowser").val();
        var maxPval = $("#maxPvalBrowser").val();
        var minFe = $("#minFeBrowser").val();
        var minDepth = $("#minDepthBrowser").val();
        var useShapes = $("#useShapesBrowser").is(":checked");
        var useColors = $("#useColorsBrowser").is(":checked");
        var plotOnly = $("#plotOnlyBrowser").is(":checked");
        var highlightList = $("#highlightListBrowser")
            .val()
            .split(",")
            .map((item) => item.trim())
            .filter((item) => item.length > 0);
        var useShapesSpecies = $("#useShapesSpeciesBrowser").is(":checked");
        var interSpeciesOnly = $("#interSpeciesOnlyBrowser").is(":checked");
        var minOrthos = $("#minOrthosBrowser").val();
        var prefix = $("#outputPrefix").val();
        var selectedPlotType = $('input[name="plotType"]:checked').val();
        var text = false;

        if (prefix !== "" && !prefix.endsWith("_")) {
            prefix += "_";
        }

        if (plant.length == 1) {
            plant = plant[0];
        }

        var requestData = {
            plant: plant,
            analysisName: analysisName,
            session_id: session_id,
        };

        // Validate inputs
        var isValid = true;
        var errorMessage = "";

        if (!plant || plant.length === 0) {
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
            errorMessage += "Threshold must be between 0 und 1.\n";
        }

        if (
            maxPval === "" ||
            !isNumeric(maxPval) ||
            parseFloat(maxPval) < 0 ||
            parseFloat(maxPval) > 1
        ) {
            isValid = false;
            errorMessage += "Max P-Value must be between 0 und 1.\n";
        }

        if (
            minFe === "" ||
            !isNumeric(minFe) ||
            parseFloat(minFe) < 0 ||
            parseFloat(minFe) > 1
        ) {
            isValid = false;
            errorMessage += "Min Fold Enrichment must be between 0 und 1.\n";
        }

        if (minDepth === "" || !isNumeric(minDepth) || parseInt(minDepth) < 1) {
            isValid = false;
            errorMessage += "Min Depth must be at least 1.\n";
        }

        if (minOrthos === "" || !isNumeric(minOrthos) || parseInt(minOrthos) < 0) {
            isValid = false;
            errorMessage +=
                "Minimum Number of Orthogroups must be a positive number.\n";
        }

        if (!isValid) {
            alert(errorMessage);
            return;
        }

        var iframe = document.getElementById("browserIframe");
        if (iframe) {
            iframe.contentWindow.postMessage({ type: "getTableParams" }, "*");
            window.addEventListener(
                "message",
                async function (event) {
                    if (event.data.type === "tableParams") {
                        var params = event.data.data;

                        try {
                            const response = await fetch(BASE_URL + "/api/get_transcripts", {
                                method: "POST",
                                headers: {
                                    "Content-Type": "application/json",
                                },
                                body: JSON.stringify(params),
                            });

                            if (!response.ok) {
                                throw new Error(`HTTP error! status: ${response.status}`);
                            }

                            const data = await response.json();
                            var transcripts = data.transcripts;
                            requestData.transcripts = transcripts;
                            requestData.prefix = prefix;

                            // Add specific parameters based on the analysis type
                            if (analysisType === "plot_co_expression_network") {
                                requestData.selectedPlotType = selectedPlotType;
                                console.log("Selected plot type:", selectedPlotType);
                                requestData.threshold = parseFloat(threshold);
                                requestData.highlightList = highlightList;
                                requestData.useShapesSpecies = useShapesSpecies;
                                requestData.useColors = useColors;
                                requestData.plotOnly = plotOnly;
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

                            sendAnalysisRequest(analysisType, requestData, "browser", text);
                            $("#runBrowserAnalysis").prop("disabled", true);
                        } catch (error) {
                            console.error("Error retrieving transcripts:", error);
                            $("#loadingBar3").hide();
                        }
                    } else {
                        $("#loadingBar3").hide();
                    }
                },
                { once: true }
            );
        } else {
            $("#loadingBar3").hide();
        }
    });

    $("#runAnalysis").click(function () {
        var plant = $("#plant").val();
        var analysisType = getAnalysisTypeWithoutSuffix();
        var analysisName = $("#analysisType option:selected").text();
        var nTop = $("#nTop").val();
        var nTopPercent = $("#nTopPercent").val();
        var threshold = $("#threshold").val();
        var maxPval = $("#maxPval").val();
        var minFe = $("#minFe").val();
        var minDepth = $("#minDepth").val();
        var positive = $("#positive").is(":checked");
        var useShapes = $("#useShapes").is(":checked");
        var useColors = $("#useColors").is(":checked");
        var interSpeciesOnly = $("#interSpeciesOnly").is(":checked");
        var prefix = $("#outputPrefix").val();
        var text = false;

        if (prefix !== "" && !prefix.endsWith("_")) {
            prefix += "_";
        }

        if (plant.length == 1) {
            plant = plant[0];
        }

        var requestData = {
            plant: plant,
            analysisName: analysisName,
            session_id: session_id,
            prefix: prefix,
        };

        // Validate inputs
        var isValid = true;
        var errorMessage = "";

        if (!plant || plant.length === 0) {
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
            errorMessage += "Top percentage must be between 0.1 und 100.\n";
        }

        if (
            threshold === "" ||
            !isNumeric(threshold) ||
            parseFloat(threshold) < 0 ||
            parseFloat(threshold) > 1
        ) {
            isValid = false;
            errorMessage += "Threshold must be between 0 und 1.\n";
        }

        if (
            maxPval === "" ||
            !isNumeric(maxPval) ||
            parseFloat(maxPval) < 0 ||
            parseFloat(maxPval) > 1
        ) {
            isValid = false;
            errorMessage += "Max P-Value must be between 0 und 1.\n";
        }

        if (
            minFe === "" ||
            !isNumeric(minFe) ||
            parseFloat(minFe) < 0 ||
            parseFloat(minFe) > 1
        ) {
            isValid = false;
            errorMessage += "Min Fold Enrichment must be between 0 und 1.\n";
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

        sendAnalysisRequest(analysisType, requestData, "general", text);
        $("#runAnalysis").prop("disabled", true);
    });

    async function startAnalysis(requestData) {
        try {
            const response = await fetch(BASE_URL + "/api/start_analysis", {
                method: "POST",
                headers: {
                    "Content-Type": "application/json",
                },
                body: JSON.stringify(requestData),
            });

            if (!response.ok) {
                throw new Error(`HTTP error! Status: ${response.status}`);
            }

            const result = await response.json();
            return result.task_id;
        } catch (error) {
            console.error("Error starting analysis:", error);
            throw error;
        }
    }

    async function checkStatus(taskId, type, text) {
        try {
            const response = await fetch(`${BASE_URL}/api/check_status/${taskId}`);
            const result = await response.json();
            console.log("Task status:", result);

            let resultsDiv;
            let loadingBar;
            let loadingText = "#loadingStatusText";

            if (type === "browser") {
                resultsDiv = "#browserAnalysisResults";
                loadingBar = text ? "#loadingText" : "#loadingBar3";
            } else {
                resultsDiv = "#results";
                loadingBar = "#loadingBar4";
            }

            if (result.state === "PENDING" || result.state === "STARTED") {
                setTimeout(() => checkStatus(taskId, type, text), 5000);
            } else if (result.state === "PROGRESS") {
                $(loadingText).text(result.status);
                setTimeout(() => checkStatus(taskId, type, text), 2000);
            } else if (result.state === "SUCCESS" && result.result) {
                console.log("Analysis result:", result.result);
                displayResult(result.result, type, text);
                $(loadingBar).hide();
            } else if (result.state === "FAILURE") {
                console.error("Task failed:", result.result.status);
                $(resultsDiv).html("<p>Error: " + result.result.message + "</p>");
                $(loadingBar).hide();
            } else {
                console.error("Unexpected task state:", result.state);
                $(resultsDiv).html(
                    "<p>Unexpected task state: " + result.state + "</p>"
                );
                $(loadingBar).hide();
            }
        } catch (error) {
            console.error("Error checking task status:", error);
            let loadingBar = type === "browser" ? "#loadingBar3" : "#loadingBar4";
            $(loadingBar).hide();
        }
    }

    function displayResult(result, type, text) {
        console.log("Displaying result:", result);
        const timestamp = new Date().getTime();
        const newWindow = $("#displayInNewWindow").is(":checked");

        let resultsDiv;
        let loadingBar;
        let startButton;
        let frameID;

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

        populateFileTypeDropdown();
        updateFileListBasedOnExtension();

        $(loadingBar).hide();

        if (startButton) {
            $(startButton).prop("disabled", false);
        }
    }

    function sendAnalysisRequest(analysisType, requestData, type, text) {
        console.log("Sending request for analysis type:", analysisType);
        console.log("Request data:", requestData);

        let loadingBar;
        let resultsDiv;

        switch (type) {
            case "browser":
                loadingBar = text ? "#loadingText" : "#loadingBar3";
                resultsDiv = "#browserAnalysisResults";
                break;
            default:
                loadingBar = "#loadingBar4";
                resultsDiv = "#results";
                break;
        }

        $(loadingBar).show();
        requestData.analysis_type = analysisType;

        startAnalysis(requestData)
            .then((taskId) => {
                checkStatus(taskId, type, text);
            })
            .catch((error) => {
                console.error("Error:", error);
                $(resultsDiv).html("<p>Error starting analysis: " + error + "</p>");
                $(loadingBar).hide();
            });
    }

    function getAnalysisTypeWithoutSuffix() {
        var analysisType = $("#analysisType").val();
        var typeParts = analysisType.split("_");

        typeParts.pop();

        var analysisTypeWithoutSuffix = typeParts.join("_");

        return analysisTypeWithoutSuffix;
    }

    function getBrowserAnalysisTypeWithoutSuffix() {
        var analysisType = $("#browserAnalysisType").val();
        var typeParts = analysisType.split("_");

        typeParts.pop();

        var analysisTypeWithoutSuffix = typeParts.join("_");

        return analysisTypeWithoutSuffix;
    }

    function isNumeric(value) {
        console.log("Checking if value is numeric:", value);
        return !isNaN(value) && value.trim() !== "";
    }

    // Download results button click handler
    $("#downloadResultsButton").click(function () {
        const sessionId = sessionStorage.getItem("session_id");
        const selectedExtension = document.getElementById("fileTypeSelect").value;

        fetch(BASE_URL + "/api/download_results", {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({
                session_id: sessionId,
                extension: selectedExtension,
            }),
        })
            .then((response) => {
                if (!response.ok) throw new Error("Network response was not ok.");
                return response.blob();
            })
            .then((blob) => {
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement("a");
                a.style.display = "none";
                a.href = url;
                a.download = "analysis.zip";
                document.body.appendChild(a);
                a.click();
                window.URL.revokeObjectURL(url);
            })
            .catch((error) => {
                console.error("There was an issue downloading the files:", error);
            });
    });

    // Clean up session button click handler
    $("#clearFilesButton").click(function () {
        const sessionId = sessionStorage.getItem("session_id");
        fetch(`${BASE_URL}/api/cleanup_session_files`, {
            method: "POST",
            headers: {
                "Content-Type": "application/json",
            },
            body: JSON.stringify({ session_id: sessionId }),
        })
            .then((response) => {
                if (!response.ok) throw new Error("Network response was not ok.");
                return response.json();
            })
            .then((data) => {
                console.log("Cleanup response:", data);
                if (data.status === "success") {
                    $("#analysisList").empty();
                    $("#downloadResultsButton").prop("disabled", true);
                    $("#clearFilesButton").prop("disabled", true);
                } else {
                    alert("Error cleaning up session files: " + data.message);
                }
            })
            .catch((error) => {
                console.error(
                    "There was an issue cleaning up the session files:",
                    error
                );
            });
    });

    // Function to update input fields based on the selected analysis type
    function updateInputFields() {
        var shapeLabel = document.getElementById("useShapesLabel");
        var selectedAnalysis = getAnalysisTypeWithoutSuffix();
        var hasParameters = false;

        $(".general-analysis-description").hide();
        $("#generalHelpIcon").show();
        $("#generalHelpIcon")
            .off("click")
            .on("click", function () {
                $("#description_" + selectedAnalysis).toggle();
            });

        $("#nTopGroup").hide();
        $("#nTopPercentGroup").hide();
        $("#thresholdGroup").hide();
        $("#positiveGroup").hide();
        $("#maxPvalGroup").hide();
        $("#minFeGroup").hide();
        $("#minDepthGroup").hide();
        $("#useShapesGroup").hide();
        $("#useColorsGroup").hide();
        $("#interSpeciesOnlyGroup").hide();

        if (selectedAnalysis === "plot_correlation_network") {
            shapeLabel.innerHTML = "Use Shapes for Modules";
        } else {
            shapeLabel.innerHTML = "Use Shapes for Species";
        }

        // Determine which input fields to show
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
                $("#thresholdGroup").show();
                $("#useShapesGroup").show();
                $("#useColorsGroup").show();
                hasParameters = true;
                break;
            case "plot_species_correlation_network":
                $("#thresholdGroup").show();
                $("#maxPvalGroup").show();
                $("#positiveGroup").show();
                $("#useShapesGroup").show();
                $("#useColorsGroup").show();
                $("#interSpeciesOnlyGroup").show();
                hasParameters = true;
                break;
            case "plot_module_goea":
                $("#maxPvalGroup").show();
                $("#minFeGroup").show();
                $("#minDepthGroup").show();
                $("#nTopPercentGroup").show();
                hasParameters = true;
                break;
            case "plot_jaccard_tissues":
                $("#nTopGroup").show();
                $("#thresholdGroup").show();
                $("#useShapesGroup").show();
                $("#interSpeciesOnlyGroup").show();
                hasParameters = true;
                break;
            case "plot_tissues_corr":
                $("#nTopGroup").show();
                hasParameters = true;
                break;
            case "plot_jaccard_modules":
                $("#thresholdGroup").show();
                $("#useShapesGroup").show();
                $("#useColorsGroup").show();
                $("#interSpeciesOnlyGroup").show();
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

    // Analysis type change event handler
    $("#analysisType").change(function () {
        updateInputFields();
    });

    function updateBrowserInputFields() {
        var selectedAnalysis = getBrowserAnalysisTypeWithoutSuffix();
        var browserButton = $("#runBrowserAnalysis");
        var hasParameters = false;

        $(".analysis-description").hide();
        $("#browserHelpIcon").show();
        $("#browserHelpIcon")
            .off("click")
            .on("click", function () {
                $("#description_" + selectedAnalysis).toggle();
            });

        browserButton.prop("disabled", false);
        console.log("Selected analysis:", selectedAnalysis);

        $("#nTopPercentGroupBrowser").hide();
        $("#nTopGroupBrowser").hide();
        $("#thresholdGroupBrowser").hide();
        $("#maxPvalGroupBrowser").hide();
        $("#minFeGroupBrowser").hide();
        $("#minDepthGroupBrowser").hide();
        $("#useShapesGroupBrowser").hide();
        $("#useShapesSpeciesGroupBrowser").hide();
        $("#useColorsGroupBrowser").hide();
        $("#plotOnlyGroupBrowser").hide();
        $("#interSpeciesOnlyGroupBrowser").hide();
        $("#minOrthosGroupBrowser").hide();
        $("#highlightListGroupBrowser").hide();
        // $("#plotTypeGroupBrowser").hide();

        switch (selectedAnalysis) {
            case "plot_co_expression_network":
                // $("#plotTypeGroupBrowser").show();
                $("#thresholdGroupBrowser").show();
                $("#highlightListGroupBrowser").show();
                $("#useShapesSpeciesGroupBrowser").show();
                $("#useColorsGroupBrowser").show();
                $("#plotOnlyGroupBrowser").show();
                hasParameters = true;
                break;
            case "plot_go_terms":
                $("#maxPvalGroupBrowser").show();
                $("#minFeGroupBrowser").show();
                $("#minDepthGroupBrowser").show();
                $("#nTopPercentGroupBrowser").show();
                hasParameters = true;
                break;
            case "plot_filtered_jaccard_modules":
                $("#thresholdGroupBrowser").show();
                $("#useShapesSpeciesGroupBrowser").show();
                $("#useColorsGroupBrowser").show();
                $("#interSpeciesOnlyGroupBrowser").show();
                $("#minOrthosGroupBrowser").show();
                hasParameters = true;
                break;
            case "plot_filtered_jaccard_species":
                $("#thresholdGroupBrowser").show();
                $("#useShapesSpeciesGroupBrowser").show();
                $("#minOrthosGroupBrowser").show();
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

    // Browser analysis type change event handler
    $("#browserAnalysisType").change(function () {
        updateBrowserInputFields();
    });

    // Ensure the correct fields are displayed based on the initial selection when the page loads
    $("#analysisType").trigger("change");
    $("#browserAnalysisType").trigger("change");

    // Cleanup session on page unload
    window.addEventListener("beforeunload", function () {
        console.log("Sending cleanup request");
        const data = JSON.stringify({
            session_id: sessionStorage.getItem("session_id"),
        });
        const blob = new Blob([data], { type: "application/json" });
        navigator.sendBeacon(BASE_URL + "/api/cleanup", blob);
    });

    $("#plant").selectpicker({
        selectedTextFormat: "count > 1",
        countSelectedText: function (numSelected, numTotal) {
            if (numSelected === 1) {
                var selectedOption = this.$element.find("option:selected").text();
                return selectedOption;
            }
            return numSelected + " of " + numTotal + " selected";
        },
        selectAllText: "Select All",
        deselectAllText: "Deselect All",
        showTick: true,
    });

    $("#plant").on(
        "hidden.bs.select",
        function (e, clickedIndex, isSelected, previousValue) {
            $("#plant").prop("disabled", true);
            handleSelectionChange();
        }
    );

    // Initialize button text based on the collapse state for non-Navbar buttons
    $(".toggle-button")
        .not(".navbar-toggler")
        .each(function () {
            var target = $(this).data("bs-target");
            $(this)
                .find("span")
                .text($(target).hasClass("show") ? "Hide" : "Show");
        });

    // Update button text on collapse show/hide, excluding the Navbar toggle button
    $(".collapse")
        .on("show.bs.collapse", function () {
            $('button.toggle-button[data-bs-target="#' + this.id + '"]')
                .not(".navbar-toggler")
                .find("span")
                .text("Hide");
        })
        .on("hide.bs.collapse", function () {
            $('button.toggle-button[data-bs-target="#' + this.id + '"]')
                .not(".navbar-toggler")
                .find("span")
                .text("Show");
        });

    // Show the card section if a navbar link is clicked
    $(".navbar-nav .nav-link").click(function (event) {
        event.preventDefault(); // Prevent the default anchor behavior

        // Check the target ID from the href attribute of the clicked link
        var targetId = $(this).attr("href");

        // Find the corresponding card section and show it if hidden
        var card = $(targetId);

        var collapseSection = card.find(".collapse");

        if (!collapseSection.hasClass("show")) {
            collapseSection.collapse("show");
        }

        // Scroll to the target section smoothly
        $("html, body").animate(
            {
                scrollTop: $(targetId).offset().top,
            },
            500
        );
    });

    function handleSelectionChange() {
        var selectedOptions = $("#plant option:selected");
        if (selectedOptions.length > 0) {
            selectPlants();
            $("#browserAnalysisCard").hide();
            $("#browserAnalysisResults").hide();
            $("#browserResults").hide();
            $("#results").hide();
            startInfo();
            startBrowser();
        } else {
            $("#infoDiv").hide();
            $("#settingsCard").hide();
            $("#downstreamCard").hide();
            $("#browserCard").hide();
            $("#results").hide();
            $("#browserResults").hide();
            $("#browserAnalysisCard").hide();
            $("#browserAnalysisResults").hide();
            $("#imageCategoryDiv").hide();
        }

        if (selectedOptions.length > 1) {
            $("#imageCategoryDiv").hide();
            $("#plotTypePlotly").prop("disabled", true);
            $("#plotTypeCytoscape").prop("checked", true);
        } else {
            $("#plotTypePlotly").prop("disabled", false);
        }

        updatePlantSelection();
        $("#plant").prop("disabled", false);
    }

    handleSelectionChange();

    $("#selectPlantsCard").fadeIn();
    $("#infoDivCard").fadeIn();
    $("#settingsCard").fadeIn();
    $("#finishedCard").fadeIn();
    $("#downloadResultsButton").prop("disabled", true);
    $("#clearFilesButton").prop("disabled", true);
    $("#loadingSpinner").hide();
    $(".container").fadeIn();
});

function sendGoTermToIframe(goTerm) {
    var iframe = document.getElementById("browserIframe");
    if (iframe && iframe.contentWindow) {
        // Send the GO term to the iframe
        iframe.contentWindow.postMessage({ type: "goTerm", value: goTerm }, "*");
    }
}
