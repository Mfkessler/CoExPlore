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
        let selectedPlant = $("input[name='plant']:checked").map(function () {
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

    $("#selectAllSpecies").on("click", function () {
        $("input[name='plant']").prop("checked", true).trigger("change");
    });
    $("#deselectAllSpecies").on("click", function () {
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

    let generalAiDataTable = null;

    function showAiSqlQuery(query) {
        if (query) {
            let cleaned = query.replace(/^[\s\S]*?(?=SELECT)/i, '')
                .trim()
                .split('\n')
                .map(line => line.trim())
                .join('\n');
            $('#generalAiSqlQuery').text(cleaned);
            $('#generalAiSqlQueryCard').fadeIn(200);
            Prism.highlightElement(document.getElementById('generalAiSqlQuery'));
        } else {
            $('#generalAiSqlQueryCard').fadeOut(100, function() {
                $('#generalAiSqlQuery').text('');
            });
        }
    }

    $(document).on('click', '#copyAiSqlBtn', function () {
        const sql = $('#generalAiSqlQuery').text();
        if (sql) {
            navigator.clipboard.writeText(sql);
            $(this).tooltip({title: "Copied!", trigger: "manual"}).tooltip('show');
            setTimeout(() => $(this).tooltip('hide'), 800);
        }
    });

    $(document).on('click', '#generalAiBtn', async function () {
        // Hide result container and table before query (show only spinner)
        $('#generalAiResultContainer').hide();
        $('#generalAiTextResult').text('');
        $('#generalAiTableWrapper').hide(); // Always hide table at start

        const question = $('#generalAiInput').val().trim();
        if (!question) {
            alert("Please enter a question.");
            return;
        }

        $('#generalAiBtn').prop('disabled', true);
        $('#generalAiSpinner').show();

        try {
            // Step 1: Get SQL query from AI
            const res = await fetch(BASE_URL + '/ai-search', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ question: question, backend: "sql", mode: "general" })
            });
            const aiResult = await res.json();
            console.log("AI Result:", aiResult);

            // If no SQL query could be generated, show only error or info message (hide table!)
            if (!aiResult.query) {
                showAiSqlQuery(null);
                if (aiResult.error) {
                    $('#generalAiTextResult').text('Error: ' + aiResult.error);
                } else if (aiResult.text) {
                    $('#generalAiTextResult').text(aiResult.text);
                } else {
                    $('#generalAiTextResult').text('No results.');
                }
                // Destroy DataTable and hide table wrapper (no table visible)
                if ($.fn.DataTable.isDataTable('#generalAiResultTable')) {
                    $('#generalAiResultTable').DataTable().clear().draw();
                    $('#generalAiResultTable').DataTable().destroy();
                }
                $('#generalAiResultTable').html('<thead></thead><tbody></tbody>');
                $('#generalAiTableWrapper').hide(); // <--- Hide whole table including DataTables UI
                $('#generalAiResultContainer').show(); // Show section for error/info
                return;
            }

            showAiSqlQuery(aiResult.query);

            // Step 2: Execute the SQL query against the DB
            const dbRes = await fetch(BASE_URL + '/api/general-ai', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ ai_query: aiResult.query })
            });
            const dbResult = await dbRes.json();

            // If DB returns error, show only error message (hide table!)
            if (dbResult.error) {
                $('#generalAiTextResult').text('Error: ' + dbResult.error);
                if ($.fn.DataTable.isDataTable('#generalAiResultTable')) {
                    $('#generalAiResultTable').DataTable().clear().draw();
                    $('#generalAiResultTable').DataTable().destroy();
                }
                $('#generalAiResultTable').html('<thead></thead><tbody></tbody>');
                $('#generalAiTableWrapper').hide(); // <--- Hide table!
                showAiSqlQuery(null);
                $('#generalAiResultContainer').show();
                return;
            }

            // Success: Show result table if data is available
            if (dbResult.data && Array.isArray(dbResult.data) && dbResult.data.length > 0) {
                const data = dbResult.data;
                const keys = Object.keys(data[0]);
                let thead = '<tr>' + keys.map(k => `<th>${k}</th>`).join('') + '</tr>';

                if ($.fn.DataTable.isDataTable('#generalAiResultTable')) {
                    $('#generalAiResultTable').DataTable().clear().draw();
                    $('#generalAiResultTable').DataTable().destroy();
                }

                // Set new header/body
                $('#generalAiResultTable').html(`<thead>${thead}</thead><tbody></tbody>`);
                $('#generalAiTableWrapper').show(); // <--- Show table (incl. DataTables UI)
                $('#generalAiResultContainer').show();

                setTimeout(function() {
                    generalAiDataTable = $('#generalAiResultTable').DataTable({
                        data: data,
                        columns: keys.map(k => ({ data: k })),
                        responsive: true,
                        dom: 'Bfrtip',
                        pageLength: 10,
                        lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]],
                        buttons: [
                            { extend: 'colvis', text: 'Column Visibility' },
                            { extend: 'csvHtml5', text: 'Download CSV', filename: 'ai_general_result' }
                        ],
                        language: {
                            search: "Search:",
                            lengthMenu: "Show _MENU_ entries",
                            info: "Showing _START_ to _END_ of _TOTAL_ entries",
                            paginate: { previous: "Prev", next: "Next" }
                        }
                    });
                }, 0);

                $('#generalAiTextResult').text('');
            } else {
                // No results found: hide table and show only info
                if (dbResult.text) {
                    $('#generalAiTextResult').text(dbResult.text);
                } else {
                    $('#generalAiTextResult').text('No results.');
                }
                if ($.fn.DataTable.isDataTable('#generalAiResultTable')) {
                    $('#generalAiResultTable').DataTable().clear().draw();
                    $('#generalAiResultTable').DataTable().destroy();
                }
                $('#generalAiResultTable').html('<thead></thead><tbody></tbody>');
                $('#generalAiTableWrapper').hide(); // <--- Hide table!
            }

            // Always show result container for error/info/results
            $('#generalAiResultContainer').show();

        } catch (err) {
            // Fatal error: destroy table, show only error text
            console.error("General AI Search failed:", err);
            alert("AI Search failed: " + err.message);
            if ($.fn.DataTable.isDataTable('#generalAiResultTable')) {
                $('#generalAiResultTable').DataTable().clear().draw();
                $('#generalAiResultTable').DataTable().destroy();
            }
            $('#generalAiResultTable').html('<thead></thead><tbody></tbody>');
            $('#generalAiTableWrapper').hide(); // <--- Hide table!
            showAiSqlQuery(null);
            $('#generalAiResultContainer').show();
        } finally {
            $('#generalAiBtn').prop('disabled', false);
            $('#generalAiBtnText').text('Run');
            $('#generalAiSpinner').hide();
        }
    });

    $(document).on('click', '#generalAiRandomBtn', async function () {
        $('#generalAiRandomBtn').prop('disabled', true).text('...');
        try {
            const res = await fetch(BASE_URL + '/api/random-example', { method: 'GET' });
            const data = await res.json();
            if (data && data.question) {
                $('#generalAiInput').val(data.question);
            } else {
                alert("No example loaded!");
            }
        } catch (e) {
            alert("Failed to load random example.");
        } finally {
            $('#generalAiRandomBtn').prop('disabled', false).text('Random');
        }
    });

}
