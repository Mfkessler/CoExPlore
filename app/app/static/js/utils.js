// utils.js

/**
 * Generate a new UUID.
 * @returns {string} A new UUID.
 */
export function generateUUID() {
    let d = new Date().getTime();
    let d2 = (performance && performance.now && performance.now() * 1000) || 0;
    return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g, (c) => {
        let r = Math.random() * 16;
        if (d > 0) {
            r = (d + r) % 16 | 0;
            d = Math.floor(d / 16);
        } else {
            r = (d2 + r) % 16 | 0;
            d2 = Math.floor(d2 / 16);
        }
        return (c === "x" ? r : (r & 0x3) | 0x8).toString(16);
    });
}

/**
 * Check if a string value is numeric.
 * @param {string} value The value to check.
 * @returns {boolean} True if numeric, false otherwise.
 */
export function isNumeric(value) {
    return !isNaN(value) && value.trim() !== "";
}

/**
 * Remove the suffix from an analysis type string.
 * @param {string} analysisType The analysis type string.
 * @returns {string} The analysis type without its suffix.
 */
export function getAnalysisTypeWithoutSuffix(analysisType) {
    const parts = analysisType.split("_");
    parts.pop();
    return parts.join("_");
}

/**
 * Remove the suffix from a browser analysis type string.
 * @param {string} analysisType The browser analysis type string.
 * @returns {string} The browser analysis type without its suffix.
 */
export function getBrowserAnalysisTypeWithoutSuffix(analysisType) {
    const parts = analysisType.split("_");
    parts.pop();
    return parts.join("_");
}

/**
 * Construct an image URL.
 * @param {string} name The image filename.
 * @param {string} baseUrl The base URL.
 * @param {string} outputDir The output directory.
 * @param {string} sessionId The session ID.
 * @returns {string} The complete image URL.
 */
export function constructImageUrl(name, baseUrl, outputDir, sessionId) {
    if (baseUrl === "") {
        return `/${outputDir}/${sessionId}/${name}`;
    }
    if (!name.startsWith(baseUrl)) {
        return `${baseUrl}/${outputDir}/${sessionId}/${name}`;
    }
    return name;
}
