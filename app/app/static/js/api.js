// api.js

/**
 * Send a heartbeat to the server.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function sendHeartbeat(baseUrl, sessionId) {
    fetch(`${baseUrl}/api/heartbeat`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: sessionId }),
    })
        .then((response) => response.text())
        .then((text) => console.log("Heartbeat response:", text))
        .catch((err) => console.error("Error sending heartbeat:", err));
}

/**
 * Initialize the session on the server.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function initSession(baseUrl, sessionId) {
    fetch(`${baseUrl}/api/init_session`, {
        method: "POST",
        headers: {
            "Content-Type": "application/json",
            "Session-ID": sessionId,
        },
    })
        .then((response) => response.json())
        .then((data) => console.log("Session initialized:", data))
        .catch((error) => console.error("Error initializing session:", error));
}

/**
 * Fetch available file extensions.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 * @returns {Promise<Object>} A promise with the JSON response.
 */
export function fetchFileExtensions(baseUrl, sessionId) {
    return fetch(`${baseUrl}/api/list_file_extensions`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: sessionId }),
    }).then((response) => response.json());
}

/**
 * Fetch filtered file list based on a file extension.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 * @param {string} extension The file extension.
 * @returns {Promise<Object>} A promise with the JSON response.
 */
export function fetchFilteredFiles(baseUrl, sessionId, extension) {
    return fetch(`${baseUrl}/api/list_filtered_files`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: sessionId, extension: extension }),
    }).then((response) => response.json());
}

/**
 * Load image data for a plant and category.
 * @param {string} baseUrl The base URL.
 * @param {string} plant The plant identifier.
 * @param {string} category The image category.
 * @returns {Promise<Object>} A promise with the JSON response.
 */
export function loadImageDataApi(baseUrl, plant, category) {
    return fetch(`${baseUrl}/api/load_image_data`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ plant: plant, category: category }),
    }).then((response) => response.json());
}

/**
 * Start an analysis on the server.
 * @param {string} baseUrl The base URL.
 * @param {Object} requestData The request data.
 * @returns {Promise<string>} A promise resolving to the task ID.
 */
export async function startAnalysis(baseUrl, requestData) {
    const response = await fetch(`${baseUrl}/api/start_analysis`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(requestData),
    });
    if (!response.ok) throw new Error(`HTTP error! Status: ${response.status}`);
    const result = await response.json();
    return result.task_id;
}

/**
 * Check the status of a task.
 * @param {string} baseUrl The base URL.
 * @param {string} taskId The task ID.
 * @returns {Promise<Object>} A promise with the JSON response.
 */
export function checkTaskStatus(baseUrl, taskId) {
    return fetch(`${baseUrl}/api/check_status/${taskId}`).then((response) =>
        response.json()
    );
}

/**
 * Download analysis results as a ZIP file.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 * @param {string} extension The file extension.
 */
export function downloadResults(baseUrl, sessionId, extension) {
    fetch(`${baseUrl}/api/download_results`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: sessionId, extension: extension }),
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
        .catch((error) => console.error("Error downloading results:", error));
}

/**
 * Clear session files on the server.
 * @param {string} baseUrl The base URL.
 * @param {string} sessionId The session ID.
 */
export function clearSessionFiles(baseUrl, sessionId) {
    fetch(`${baseUrl}/api/cleanup_session_files`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ session_id: sessionId }),
    })
        .then((response) => {
            if (!response.ok) throw new Error("Network response was not ok.");
            return response.json();
        })
        .then((data) => console.log("Cleanup response:", data))
        .catch((error) =>
            console.error("Error clearing session files:", error)
        );
}
