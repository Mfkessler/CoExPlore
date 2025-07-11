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

    $("input[name='analysisMode']").on("change", function() {
      if ($("#analysis-whole").is(":checked")) {
        $("#wholeAnalysisSection").show();
        $("#browserAnalysisSection").hide();
      } else if ($("#analysis-browser").is(":checked")) {
        $("#wholeAnalysisSection").hide();
        $("#browserAnalysisSection").show();
      }
    });

    if ($("#analysis-whole").is(":checked")) {
      $("#wholeAnalysisSection").show();
      $("#browserAnalysisSection").hide();
    } else if ($("#analysis-browser").is(":checked")) {
      $("#wholeAnalysisSection").hide();
      $("#browserAnalysisSection").show();
    }

    // Ensure the correct fields are displayed based on the initial selection when the page loads
    $("#analysisType").trigger("change");
    $("#browserAnalysisType").trigger("change");
    
    // Show relevant cards, disable buttons, and hide the spinner
    $("#selectPlantsCard, #infoDivCard, #settingsCard, #finishedCard", ).fadeIn();
    $("#downloadResultsButton, #clearFilesButton").prop("disabled", true);
    $("#loadingSpinner").hide();
    $(".container").show();
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

// Sidebar klickt einfach den richtigen Tab-Button an:
document.getElementById('sidebar-datasets').addEventListener('click', function() {
  document.getElementById('tabbtn-datasets').click();
});
document.getElementById('sidebar-info').addEventListener('click', function() {
  document.getElementById('tabbtn-info').click();
});
document.getElementById('sidebar-browser').addEventListener('click', function() {
  document.getElementById('tabbtn-browser').click();
});
document.getElementById('sidebar-analysis').addEventListener('click', function() {
  document.getElementById('tabbtn-analysis').click();
});
document.getElementById('sidebar-files').addEventListener('click', function() {
  document.getElementById('tabbtn-files').click();
});

// Sidebar-Highlight synchronisieren mit aktivem Tab:
document.querySelectorAll('.nav-tabs .nav-link').forEach(tabLink => {
  tabLink.addEventListener('shown.bs.tab', function() {
    const tabId = this.getAttribute('data-bs-target');
    document.querySelectorAll('.sidebar .sidebar-icon').forEach(sideIcon => {
      sideIcon.classList.remove('active');
    });
    const sideId = 'sidebar-' + tabId.replace('#tab-', '');
    const sideIcon = document.getElementById(sideId);
    if(sideIcon) sideIcon.classList.add('active');
  });
});
