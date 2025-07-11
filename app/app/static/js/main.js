// main.js
import { generateUUID } from "./utils.js";
import { sendHeartbeat, initSession } from "./api.js";
import { bindEvents } from "./events.js";
import { handleSelectionChange } from "./ui.js";

console.log(`Using BASE_URL: ${BASE_URL}`);
console.log(`Using OUTPUT_DIR: ${OUTPUT_DIR}`);

$(document).ready(function () {
    // Store original options (for later filtering)
    window.originalOptions = Array.from(
        document.getElementById("analysisType").options
    );
    window.originalBrowserOptions = Array.from(
        document.getElementById("browserAnalysisType").options
    );

    // Initialize session
    if (!sessionStorage.getItem("session_id")) {
        sessionStorage.setItem("session_id", generateUUID());
    }
    const sessionId = sessionStorage.getItem("session_id");
    console.log("Initial session ID:", sessionId);
    sendHeartbeat(BASE_URL, sessionId);
    initSession(BASE_URL, sessionId);
    // Send heartbeat every 5 minutes
    setInterval(() => {
        sendHeartbeat(BASE_URL, sessionId);
    }, 300000);

    // Bind all events
    bindEvents(BASE_URL, sessionId);

    // Initialize toggle button text (excluding Navbar)
    $(".toggle-button")
        .not(".navbar-toggler")
        .each(function () {
            const target = $(this).data("bs-target");
            $(this)
                .find("span")
                .text($(target).hasClass("show") ? "Hide" : "Show");
        });

    handleSelectionChange(BASE_URL, sessionId);

    // Ensure the correct fields are displayed based on the initial selection when the page loads
    $("#analysisType").trigger("change");
    $("#browserAnalysisType").trigger("change");
    
    // Show relevant cards, disable buttons, and hide the spinner
    $("#selectPlantsCard, #infoDivCard, #settingsCard, #finishedCard").fadeIn();
    $("#downloadResultsButton, #clearFilesButton").prop("disabled", true);
    $("#loadingSpinner").hide();
    $(".container").fadeIn();
});

/**
 * Send a GO term to the browser iframe.
 * @param {string} goTerm The GO term to send.
 */
export function sendGoTermToIframe(goTerm) {
    const iframe = document.getElementById("browserIframe");
    if (iframe && iframe.contentWindow) {
        iframe.contentWindow.postMessage({ type: "goTerm", value: goTerm }, "*");
    }
}

window.sendGoTermToIframe = sendGoTermToIframe;

// Sidebar navigation and tab switching
document.querySelectorAll('.sidebar .nav-link').forEach(link => {
  link.addEventListener('click', function(e) {
    e.preventDefault();
    // Sidebar highlight
    document.querySelectorAll('.sidebar .nav-link').forEach(l => l.classList.remove('active'));
    this.classList.add('active');
    // Tabs content
    const tabId = this.getAttribute('data-tab');
    document.querySelectorAll('.tab-pane').forEach(pane => pane.classList.remove('show', 'active'));
    const tabPane = document.querySelector(tabId);
    if(tabPane) {
      tabPane.classList.add('show', 'active');
    }
    // Main tabs highlight
    document.querySelectorAll('.nav-tabs .nav-link').forEach(l => l.classList.remove('active'));
    const mainTab = document.querySelector('.nav-tabs .nav-link[data-bs-target="'+tabId+'"]');
    if(mainTab) mainTab.classList.add('active');
  });
});

// Synchronize sidebar links with main tabs
document.querySelectorAll('.nav-tabs .nav-link').forEach(tabLink => {
  tabLink.addEventListener('click', function() {
    const tabId = this.getAttribute('data-bs-target');
    document.querySelectorAll('.sidebar .nav-link').forEach(l => l.classList.remove('active'));
    const sidebarLink = document.querySelector('.sidebar .nav-link[data-tab="'+tabId+'"]');
    if(sidebarLink) sidebarLink.classList.add('active');
  });
});
