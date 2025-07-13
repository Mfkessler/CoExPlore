// events.js
import { downloadResults, clearSessionFiles } from "./api.js";
import {
    loadImageDataUI,
    runBrowserAnalysisHandler,
    runGeneralAnalysisHandler,
    adjustThumbnailSizeUI,
    handleMessageEvent,
    handleSelectionChange,
    updateInputFields,
    updateBrowserInputFields,
    updateFileListBasedOnExtension
} from "./ui.js";

/**
 * Bind all DOM events.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function bindEvents(baseUrl, sessionId) {
    // Help icon toggles
    const helpIconBrowser = document.getElementById("helpIconBrowser");
    if (helpIconBrowser) {
        helpIconBrowser.addEventListener("click", () => {
            $("#userGuideCardBrowser").toggle();
        });
    }

    const helpIconDataset = document.getElementById("helpIconDataset");
    if (helpIconDataset) {
        helpIconDataset.addEventListener("click", () => {
            $("#userGuideCardDataset").toggle();
        });
    }

    // File type dropdown change
    const fileTypeSelect = document.getElementById("fileTypeSelect");
    if (fileTypeSelect) {
      fileTypeSelect.addEventListener("change", () => {
        updateFileListBasedOnExtension(baseUrl);
      });
    }

    // Image category change event
    $("#imageCategory").change(function () {
        const selectedCategory = $(this).val();
        let selectedPlant = $("input[name='plant']:checked").map(function() {
            return this.value;
        }).get();
        if (selectedPlant.length === 1) {
            selectedPlant = selectedPlant[0];
        }
        loadImageDataUI(baseUrl, selectedPlant, selectedCategory);
    });

    // Run browser analysis
    $("#runBrowserAnalysis").click(function () {
        runBrowserAnalysisHandler(baseUrl, sessionId);
    });

    // Run general analysis
    $("#runAnalysis").click(function () {
        runGeneralAnalysisHandler(baseUrl, sessionId);
    });

    // Update input fields on analysisType change
    $("#analysisType").change(function () {
        updateInputFields();
    });

    // Update browser input fields on browserAnalysisType change
    $("#browserAnalysisType").change(function () {
        updateBrowserInputFields();
    });

    // Download results button
    $("#downloadResultsButton").click(function () {
        const selectedExtension =
            document.getElementById("fileTypeSelect").value;
        downloadResults(baseUrl, sessionId, selectedExtension);
    });

    // Clear files button
    $("#clearFilesButton").click(function () {
        clearSessionFiles(baseUrl, sessionId);
        $("#filesDiv").hide();
        $("#downloadResultsButton").prop("disabled", true);
        $("#clearFilesButton").prop("disabled", true);
        
    });

    // Toggle button text on collapse
    $(".collapse")
        .on("show.bs.collapse", function () {
            $(`button.toggle-button[data-bs-target="#${this.id}"]`)
                .not(".navbar-toggler")
                .find("span")
                .text("Hide");
        })
        .on("hide.bs.collapse", function () {
            $(`button.toggle-button[data-bs-target="#${this.id}"]`)
                .not(".navbar-toggler")
                .find("span")
                .text("Show");
        });

    // Navbar click: show card and scroll
    $(".navbar-nav .nav-link").click(function (event) {
        event.preventDefault();
        const targetId = $(this).attr("href");
        const card = $(targetId);
        const collapseSection = card.find(".collapse");
        if (!collapseSection.hasClass("show")) {
            collapseSection.collapse("show");
        }
        $("html, body").animate(
            { scrollTop: $(targetId).offset().top },
            500
        );
    });

    // Plant selection (Card Multi-Select)
    $("input[name='plant']").on("change", function () {
        $("#analysisForm").prop("disabled", true); // Optional, wenn du es brauchst
        handleSelectionChange(baseUrl, sessionId);
    });

    $("#selectAllSpecies").on("click", function() {
        $("input[name='plant']").prop("checked", true).trigger("change");
    });
    $("#deselectAllSpecies").on("click", function() {
        $("input[name='plant']").prop("checked", false).trigger("change");
    });

    // Window resize: adjust thumbnails
    $(window).on("resize", function () {
        adjustThumbnailSizeUI();
    });

    // Message events from iframes
    window.addEventListener("message", function (event) {
        if (event.origin !== window.location.origin) return;
        handleMessageEvent(event, baseUrl, sessionId);
    });

    // Cleanup on unload
    window.addEventListener("beforeunload", function () {
        const data = JSON.stringify({ session_id: sessionId });
        const blob = new Blob([data], { type: "application/json" });
        navigator.sendBeacon(`${baseUrl}/api/cleanup`, blob);
    });

    $(document).on('click', '#generalAiBtn', async function () {
        $('#generalAiTextResult').text('');
        $('#generalAiFallbackTable thead').empty();
        $('#generalAiFallbackTable tbody').empty();
        $('#generalAiResultContainer').hide();

        const question = $('#generalAiInput').val().trim();
        if (!question) {
            alert("Please enter a question.");
            return;
        }

        $('#generalAiBtn').prop('disabled', true);
        $('#generalAiSpinner').show();

        try {
            // 1. Get SQL query from AI (CyFISH) prompt
            const res = await fetch(BASE_URL + '/ai-search', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    question: question,
                    backend: "sql",
                    type: "general"
                })
            });
            const aiResult = await res.json();

            if (!aiResult.query) {
                // If only text is returned (no SQL, e.g. "Result not displayable as table"), show the text
                if (aiResult.text) {
                    $('#generalAiTextResult').text(aiResult.text);
                    $('#generalAiResultContainer').show();
                } else {
                    $('#generalAiTextResult').text('No results.');
                    $('#generalAiResultContainer').show();
                }
                return;
            }

            // 2. Execute query (on any table, no restriction)
            const dbRes = await fetch(BASE_URL + '/api/general-ai', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ ai_query: aiResult.query })
            });
            const dbResult = await dbRes.json();

            // If backend returns an error
            if (dbResult.error) {
                $('#generalAiTextResult').text('Error: ' + dbResult.error);
                $('#generalAiResultContainer').show();
                return;
            }

            // Display data
            if (dbResult.data && Array.isArray(dbResult.data) && dbResult.data.length > 0) {
                const keys = Object.keys(dbResult.data[0]);
                const headerRow = $('<tr>');
                keys.forEach(key => headerRow.append($('<th>').text(key)));
                $('#generalAiFallbackTable thead').empty().append(headerRow);
                $('#generalAiFallbackTable tbody').empty();
                dbResult.data.forEach(row => {
                    const rowEl = $('<tr>');
                    keys.forEach(key => {
                        const value = row[key] != null ? row[key] : '';
                        rowEl.append($('<td>').text(value));
                    });
                    $('#generalAiFallbackTable tbody').append(rowEl);
                });
                $('#generalAiFallbackTable').show();
                $('#generalAiResultContainer').show();
                // Optionally: also display/log the SQL query
                console.log("AI Query:", aiResult.query);
            } else {
                // If no table result, but maybe a text
                if (dbResult.text) {
                    $('#generalAiTextResult').text(dbResult.text);
                } else {
                    $('#generalAiTextResult').text('No results.');
                }
                $('#generalAiFallbackTable').hide();
                $('#generalAiResultContainer').show();
            }
        } catch (err) {
            console.error("General AI Search failed:", err);
            alert("AI Search failed: " + err.message);
        } finally {
            $('#generalAiBtn').prop('disabled', false);
            $('#generalAiBtnText').text('Run');
            $('#generalAiSpinner').hide();
        }
    });

}
