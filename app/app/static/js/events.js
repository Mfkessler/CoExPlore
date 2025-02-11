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
    document.getElementById("helpIconBrowser").addEventListener("click", () => {
        $("#userGuideCardBrowser").toggle();
    });
    document.getElementById("helpIconDataset").addEventListener("click", () => {
        $("#userGuideCardDataset").toggle();
    });

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
        let selectedPlant = $("#plant").val();
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

    // Plant selection (selectpicker)
    $("#plant").on("hidden.bs.select", function () {
        $("#plant").prop("disabled", true);
        handleSelectionChange(baseUrl, sessionId);
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
}
