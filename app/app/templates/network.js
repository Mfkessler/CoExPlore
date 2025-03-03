document.addEventListener('DOMContentLoaded', function() {
    // Initialize Cytoscape
    var filterEdges = {{ filter_edges|tojson }};  // Flag from Python
    var useBackgroundColor = {{ use_background_color|tojson }};  // Flag from Python
    var useShapes = {{ use_shapes|tojson }};
    var useClusterTooltip = {{ use_cluster_tooltip|tojson }};
    var nodeSize = {{node_size}};  // Initial node size
    var aggNodeSize = nodeSize * 4;  // Aggregated node size
    var edgeWidth = {{edge_width}};  // Initial edge width
    var highlightColor = '{{ highlight_color }}';  // Highlight color from Python
    var useEdgeTransparency = {{ use_edge_transparency|tojson }};  // Transparency flag
    var useEdgeColoring = {{ use_edge_coloring|tojson }};  // Edge coloring flag

    // Get the output path and custom filename from the template
    var customFilename = "{{ custom_filename }}";  // e.g., 'cyto_network'

    // Metadata
    const METADATA_MAPPING = {{ metadata_dict|tojson }};
    METADATA_MAPPING['gene'] = 'Transcript';
    const FORMATTED_KEYS = ['go_terms', 'ipr_id'];
    const HIDDEN_KEYS = ['description_long', 'modules_labels', 'module_count', 'unique_modules', 'unique_go_terms', 'transcript'];

    // Aggregated network
    let aggregated = true;

    let spinner = document.getElementById('loading-spinner');

    var cy = cytoscape({
        container: document.getElementById('cy'),
        elements: {{ network_data|tojson }},  // Data from Python
        style: [
            {
                selector: 'node',
                style: {
                    'background-color': function(node) {
                        if (node.data('aggregated')) {
                            return node.data('color');  // Aggregated node: use its own color
                        } else if (node.data('highlighted')) {
                            return highlightColor;
                        }
                        return useBackgroundColor ? node.data('module_colors') : 'black';
                    },
                    'width': function(node) {
                        return node.data('aggregated') ? aggNodeSize : nodeSize;
                    },
                    'height': function(node) {
                        return node.data('aggregated') ? aggNodeSize : nodeSize;
                    },
                    'text-valign': function(node) {
                        return node.data('aggregated') ? 'top' : 'center';
                    },
                    'text-halign': 'center',
                    'text-margin-y': function(node) {
                        return node.data('aggregated') ? -5 : 0;
                    },
                    'color': function(node) {
                        return node.data('aggregated') ? 'black' : 'white';
                    },
                    'font-size': function(node) {
                        return node.data('aggregated') ? 6 : 2;
                    },
                    'label': function(node) {
                        return node.data('aggregated') ? node.data('label') : '';
                    },
                    'border-color': 'black',
                    'border-width': 1
                }
            },
            {
                selector: 'node[is_neighbor]',
                style: {
                    'opacity': 0.40
                }
            },
            {
                selector: '.highlight',
                style: {
                    'background-color': 'magenta', 
                    'line-color': 'magenta',
                    'width': 2
                }
            }
        ],
        layout: {
            name: 'cose',
            fit: true,  // Ensures the entire graph is visible
            padding: 10,  // Adds padding around the layout
            nodeRepulsion: 400000,  // Increase this value to spread nodes further apart
            idealEdgeLength: 100,  // Sets an ideal edge length
            nodeDimensionsIncludeLabels: false,  // Doesn't include labels when calculating size
            animate: false  // Animates the layout
        }
    });

    function computeAggregatedFactor(count, minCount, maxCount) {
        const minScale = 1.33;
        const maxScale = 2.5;
        if (maxCount === minCount) {
            return (minScale + maxScale) / 2;
        }
        return minScale + ((count - minCount) / (maxCount - minCount)) * (maxScale - minScale);
    }
    

    /* SVG Shape Creation */
    function createShapeSVG(shape, color = 'black', size = 15) {
        const svgNS = "http://www.w3.org/2000/svg";
        const svg = document.createElementNS(svgNS, "svg");
        svg.setAttribute("width", size);
        svg.setAttribute("height", size);
        svg.setAttribute("viewBox", "0 0 20 20");
        let shapeElement;
        switch(shape) {
            case 'ellipse':
                shapeElement = document.createElementNS(svgNS, "circle");
                shapeElement.setAttribute("cx", "10");
                shapeElement.setAttribute("cy", "10");
                shapeElement.setAttribute("r", "9");
                break;
            case 'triangle':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "10,2 18,18 2,18");
                break;
            case 'rectangle':
                shapeElement = document.createElementNS(svgNS, "rect");
                shapeElement.setAttribute("x", "2");
                shapeElement.setAttribute("y", "2");
                shapeElement.setAttribute("width", "16");
                shapeElement.setAttribute("height", "16");
                break;
            case 'diamond':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "10,2 18,10 10,18 2,10");
                break;
            case 'pentagon':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "10,2 18,8 15,18 5,18 2,8");
                break;
            case 'hexagon':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "5,2 15,2 18,10 15,18 5,18 2,10");
                break;
            case 'barrel':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "4,2 16,2 18,6 18,14 16,18 4,18 2,14 2,6");
                break;
            case 'tag':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "2,2 14,2 18,10 14,18 2,18");
                break;
            case 'vee':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "2,2 10,8 18,2 10,18");
                break;
            case 'rhomboid':
                shapeElement = document.createElementNS(svgNS, "polygon");
                shapeElement.setAttribute("points", "4,18 16,18 12,2 0,2");
                break;
            default:
                shapeElement = document.createElementNS(svgNS, "rect");
                shapeElement.setAttribute("x", "2");
                shapeElement.setAttribute("y", "2");
                shapeElement.setAttribute("width", "16");
                shapeElement.setAttribute("height", "16");
                break;
        }
        shapeElement.setAttribute("fill", color);
        svg.appendChild(shapeElement);
        return svg;
    }    

    function assignDynamicNetworkIconToAggregatedNodes() {
        let aggregatedNodes = cy.nodes().filter(node => node.data('aggregated'));
        if (aggregatedNodes.empty()) return;
    
        // Determine the global minimum and maximum node_count (default value: 1)
        let counts = aggregatedNodes.map(node => node.data('node_count') || 1);
        let minCount = Math.min(...counts);
        let maxCount = Math.max(...counts);
    
        // For each aggregated node: Determine a scaled node count (between 2 and 10)
        aggregatedNodes.forEach(function(node) {
            let rawCount = node.data('node_count') || 1;
            let scaledCount;
            if (maxCount === minCount) {
                scaledCount = 6; // If all are the same, set a median value
            } else {
                scaledCount = Math.round(2 + ((rawCount - minCount) / (maxCount - minCount)) * (10 - 2));
            }
            // Base size: aggNodeSize (can be adjusted as needed)
            let baseSize = aggNodeSize;
            let newSize = baseSize; // Optionally, you can add additional scaling based on rawCount
    
            // Create the dynamic network icon with the scaled node count
            let svgElement = createDynamicNetworkIconSVG(node.data('color') || 'gray', newSize, scaledCount);
            
            // Serialize the SVG and create a Data URL from it
            let serializer = new XMLSerializer();
            let svgString = serializer.serializeToString(svgElement);
            let encodedData = encodeURIComponent(svgString);
            let dataUrl = 'data:image/svg+xml;charset=utf-8,' + encodedData;
            
            // Assign the new style to the node:
            node.style({
                'width': newSize,
                'height': newSize,
                'background-image': 'url(' + dataUrl + ')',
                'background-fit': 'cover',
                'background-opacity': 0.8,
                'background-color': 'lightgray',
                'border-width': 1,
                'border-color': 'darkgrey',
                'border-opacity': 0.8,
                'shape': 'ellipse',
                'label': '' // Remove the text if any
            });
        });

        updateEdgeStyles();
    }

    // Creates a dynamic network icon as SVG without a central circle.
    // Parameters:
    //   color: Fill color of the small nodes (taken from node data).
    //   size: Total size of the icon (in pixels).
    //   count: Number of small nodes (scaled value between 2 and 10).
    function createDynamicNetworkIconSVG(color, size, count) {
        const svgNS = "http://www.w3.org/2000/svg";
        // Create the SVG element.
        let svg = document.createElementNS(svgNS, "svg");
        svg.setAttribute("width", size);
        svg.setAttribute("height", size);
        svg.setAttribute("viewBox", "0 0 100 100");

        // Center and parameters for positioning
        const centerX = 50, centerY = 50;
        const placementRadius = 35; // Distance from the center
        const nodeRadius = 8;       // Radius of the small nodes

        // Calculate the positions of the small nodes (evenly distributed)
        let positions = [];
        for (let i = 0; i < count; i++) {
            let angle = (2 * Math.PI * i) / count;
            let cx = centerX + placementRadius * Math.cos(angle);
            let cy = centerY + placementRadius * Math.sin(angle);
            positions.push({ cx, cy });
        }

        // Draw lines from each small node to the center (for a network-like appearance)
        positions.forEach(pos => {
            let line = document.createElementNS(svgNS, "line");
            line.setAttribute("x1", pos.cx);
            line.setAttribute("y1", pos.cy);
            line.setAttribute("x2", centerX);
            line.setAttribute("y2", centerY);
            line.setAttribute("stroke", color);
            line.setAttribute("stroke-width", "2");
            svg.appendChild(line);
        });

        // Draw lines connecting adjacent nodes (close the circle)
        for (let i = 0; i < count; i++) {
            let start = positions[i];
            let end = positions[(i + 1) % count];
            let line = document.createElementNS(svgNS, "line");
            line.setAttribute("x1", start.cx);
            line.setAttribute("y1", start.cy);
            line.setAttribute("x2", end.cx);
            line.setAttribute("y2", end.cy);
            line.setAttribute("stroke", color);
            line.setAttribute("stroke-width", "2");
            svg.appendChild(line);
        }

        // Create the small nodes (circles) – with fill color and black border.
        positions.forEach(pos => {
            let circle = document.createElementNS(svgNS, "circle");
            circle.setAttribute("cx", pos.cx);
            circle.setAttribute("cy", pos.cy);
            circle.setAttribute("r", nodeRadius);
            circle.setAttribute("fill", "white");
            circle.setAttribute("stroke", color);
            circle.setAttribute("stroke-width", "1");
            svg.appendChild(circle);
        });

        return svg;
    }
    
    function debounce(func, wait) {
        let timeout;
        return function(...args) {
            clearTimeout(timeout);
            timeout = setTimeout(() => func.apply(this, args), wait);
        };
    }

    const metadataKeys = Object.keys(METADATA_MAPPING).filter(key => !HIDDEN_KEYS.includes(key));

    cy.ready(function() {
        /* Init edge visibility */
        const initialZoom = cy.zoom();
        cy.scratch('initialZoom', initialZoom);

        // Hide Legend
        document.getElementById('legend').style.display = 'none';
        
        const debouncedUpdateEdgesVisibility = debounce(updateEdgesVisibility, 100);
        document.getElementById('edge-visibility-range').addEventListener('input', debouncedUpdateEdgesVisibility);
        cy.on('zoom', debouncedUpdateEdgesVisibility);
        cy.on('pan', debouncedUpdateEdgesVisibility);

        /* Edge Filtering */
        if (filterEdges) {
            var thresholdDegree = 10; // Threshold for node degree
            var topEdges = 10;        // Number of edges to show for high-degree nodes

            // Calculate the degree for all nodes
            var degreeCount = {};
            cy.nodes().forEach(function(node) {
                degreeCount[node.id()] = node.degree();
            });

            // Create a map for edges between high-degree nodes
            var highDegreeEdges = {};

            cy.edges().forEach(function(edge) {
                var source = edge.source().id();
                var target = edge.target().id();

                // Check if both nodes exceed the degree threshold
                if (degreeCount[source] >= thresholdDegree && degreeCount[target] >= thresholdDegree) {
                    // Add this edge to the list of high-degree edges
                    if (!highDegreeEdges[source]) highDegreeEdges[source] = [];
                    if (!highDegreeEdges[target]) highDegreeEdges[target] = [];

                    highDegreeEdges[source].push(edge);
                    highDegreeEdges[target].push(edge);
                } else {
                    // Keep the edge as at least one node is below the threshold
                    edge.show();
                }
            });

            // Reduce the number of edges for high-degree nodes
            var hiddenEdges = new Set();
            Object.keys(highDegreeEdges).forEach(function(nodeId) {
                var edges = highDegreeEdges[nodeId];

                // Sort edges by weight
                var sortedEdges = edges.sort(function(a, b) {
                    return b.data('weight') - a.data('weight');
                });

                // Always show the strongest edge
                if (sortedEdges.length > 0) {
                    sortedEdges[0].show();
                }

                // Show only the top edges for this node (starting from index 1 as the strongest edge is already shown)
                for (var i = 1; i < sortedEdges.length && i < topEdges; i++) {
                    var edge = sortedEdges[i];
                    var edgeId = edge.id();

                    // Check if the edge has already been shown
                    if (!hiddenEdges.has(edgeId)) {
                        edge.show();
                        hiddenEdges.add(edgeId);  // Track shown edges to avoid duplicates
                    }
                }

                // Hide all remaining edges
                for (var i = topEdges; i < sortedEdges.length; i++) {
                    var edge = sortedEdges[i];
                    var edgeId = edge.id();

                    if (!hiddenEdges.has(edgeId)) {
                        edge.hide();
                        hiddenEdges.add(edgeId);
                    }
                }
            });
        }
    });

    function assignOrganismShapes() {
        if (useShapes && !aggregated) {
            const availableShapes = [
                'ellipse',
                'triangle',
                'rectangle',
                'diamond',
                'pentagon',
                'tag',
                'vee',
                'rhomboid',
                'hexagon',
                'barrel',
            ];

            let organismShapes = {};

            // Determine unique organisms
            let uniqueOrganisms = [...new Set(cy.nodes().map(node => node.data('species')))];
            // Assign shapes cyclically
            uniqueOrganisms.forEach((organism, index) => {
                organismShapes[organism] = availableShapes[index % availableShapes.length];
            });

            // Dynamically update node styles
            cy.batch(function () {
                cy.nodes().forEach(node => {
                    let organism = node.data('species');
                    let newShape = organismShapes[organism] || 'ellipse';
                    node.style('shape', newShape);
                });
            });

            // Dynamically build species legend with SVG symbols, a space, and italic species names
            let speciesLegendList = document.getElementById('species-legend-list');
            speciesLegendList.innerHTML = ''; // Clear the list
            Object.keys(organismShapes).forEach(species => {
                let shape = organismShapes[species];

                let li = document.createElement('li');
                li.style.marginBottom = '5px';

                // Create SVG icon for the respective shape
                let svgElement = createShapeSVG(shape);
                svgElement.classList.add('legend-shape');
                li.appendChild(svgElement);

                // Add species name
                let speciesSpan = document.createElement('span');
                speciesSpan.classList.add('legend-species');
                speciesSpan.textContent = species;
                li.appendChild(speciesSpan);

                speciesLegendList.appendChild(li);
            });
            // Show Legend
            document.getElementById('legend').style.display = 'block';
            // Show the species legend
            document.getElementById('species-legend').style.display = 'block';
        } else {
            // Hide the species legend if shapes are not used
            document.getElementById('species-legend').style.display = 'none';
        }
    }

    function clamp(val, min, max) {
        return Math.max(min, Math.min(max, val));
    }

    function lerp(a, b, t) {
        return a + (b - a) * t;
    }

    function updateEdgesVisibility() {
        if(!aggregated) {
            // Read the slider value (0 to 100)
            const sliderValue = parseInt(document.getElementById('edge-visibility-range').value, 10);
            
            // Parameters
            const exponentHigh = 4;    // For high-degree nodes – higher value = slower growth
            const thresholdLow = 15;   // Nodes with <= 10 edges are considered low-degree
            
            // Determine the current viewport (model coordinates)
            const extent = cy.extent();
            
            // Initialize sorted edges (once), sorted in descending order by weight
            if (!cy.scratch('sortedEdges')) {
                const sortedEdges = {};
                cy.nodes().forEach(node => {
                    const nodeEdges = node.connectedEdges();
                    const sorted = nodeEdges.sort((a, b) => (b.data('weight') || 0) - (a.data('weight') || 0));
                    sortedEdges[node.id()] = sorted.map(edge => edge.id());
                });
                cy.scratch('sortedEdges', sortedEdges);
            }
            const sortedEdges = cy.scratch('sortedEdges');
            
            // Hide all edges initially
            cy.edges().hide();
            
            function nodeInExtent(node) {
                const pos = node.position();
                return pos.x >= extent.x1 && pos.x <= extent.x2 &&
                    pos.y >= extent.y1 && pos.y <= extent.y2;
            }
            
            // For each node in the viewport: Calculate how many edges should be shown.
            cy.nodes().forEach(node => {
                if (nodeInExtent(node)) {
                    const edgeIds = sortedEdges[node.id()] || [];
                    const totalEdges = edgeIds.length;
                    let effectiveFraction;
                    
                    if (totalEdges <= thresholdLow) {
                        // For low-degree nodes: If slider ≥20%, show all edges; otherwise scale linearly.
                        effectiveFraction = sliderValue >= 20 ? 1 : sliderValue / 20;
                    } else {
                        // For high-degree nodes: Convex mapping for slower growth.
                        effectiveFraction = Math.pow(sliderValue / 100, exponentHigh);
                    }
                    
                    const allowedEdges = Math.ceil(effectiveFraction * totalEdges);
                    for (let i = 0; i < allowedEdges && i < totalEdges; i++) {
                        cy.getElementById(edgeIds[i]).show();
                    }
                }
            });
        }
    }

    // Calculate the minimum and maximum TOM value
    var minTOM = Math.min(...cy.edges().map(edge => edge.data('weight')));
    var maxTOM = Math.max(...cy.edges().map(edge => edge.data('weight')));

    /* Zoom Controls */

    // Store the initial zoom level and position
    var initialZoom = cy.zoom();
    var initialPan = cy.pan();
    var initialCenter = cy.center;

    // Zoom In button
    document.getElementById('zoom-in').addEventListener('click', function() {
        var currentZoom = cy.zoom();
        cy.zoom({
            level: currentZoom * 1.2,  // Zoom in by 20%
            renderedPosition: { x: cy.width() / 2, y: cy.height() / 2 }
        });
    });

    // Zoom Out button
    document.getElementById('zoom-out').addEventListener('click', function() {
        var currentZoom = cy.zoom();
        cy.zoom({
            level: currentZoom * 0.8,  // Zoom out by 20%
            renderedPosition: { x: cy.width() / 2, y: cy.height() / 2 }
        });
    });

    // Reset button
    document.getElementById('zoom-reset').addEventListener('click', function() {
        if (!aggregated) {
            var aggregatedFile = customFilename + "_aggregated.json";
            spinner.style.display = 'flex';
            fetch(aggregatedFile)
                .then(response => response.json())
                .then(aggregatedData => {
                    // Reload the network
                    cy.json({ elements: aggregatedData });
                    // Start layout
                    cy.layout({
                        name: 'cose',
                        fit: true,
                        padding: 10,
                        nodeRepulsion: 400000,
                        idealEdgeLength: 100,
                        animate: false,
                    }).run();
                    aggregated = true;
                    spinner.style.display = 'none';
                    // Apply dynamic aggregated node icons
                    assignDynamicNetworkIconToAggregatedNodes();
                    // Optionally adjust zoom and pan
                    cy.resize();
                    cy.fit();
                    // Hide legend
                    document.getElementById('legend').style.display = 'none';
                })
                .catch(error => {
                    console.error("Error loading aggregated network:", error);
                    spinner.style.display = 'none';
                });
        } 

        resetParameters();
    });    

    /* Toggle scroll zoom */

    // Variable to track the state of scroll zoom
    let scrollZoomActive = false; // Initial state off
    var rectangleZoomMode = false; // Rectangle zoom is off by default
    cy.userZoomingEnabled(false); // Ensure consistency

    document.getElementById('toggleZoom').addEventListener('click', function() {
        scrollZoomActive = !scrollZoomActive;  // Toggle scroll zoom

        if (scrollZoomActive) {
            this.classList.add('active'); // Update icon color
            cy.userZoomingEnabled(true);  // Enable scroll zoom

            // Disable rectangle zoom if active
            if (rectangleZoomMode) {
                document.getElementById('rectangle-zoom').classList.remove('active');
                rectangleZoomMode = false;
                cy.userPanningEnabled(true);
                cy.container().style.cursor = 'default';
            }
        } else {
            this.classList.remove('active'); // Update icon color
            cy.userZoomingEnabled(false);  // Disable scroll zoom
        }
    });

    // Event listener for the export button
    document.getElementById('export-json').addEventListener('click', function() {
        var cyJson = cy.json();  // Get the current Cytoscape network as JSON
        
        // Convert the JSON to a string
        var jsonString = JSON.stringify(cyJson, null, 2);

        // Create a blob from the JSON string
        var blob = new Blob([jsonString], { type: 'application/json' });

        // Create a link element to trigger the download
        var link = document.createElement('a');
        link.href = URL.createObjectURL(blob);
        link.download = 'cyto_network.json';  // Name of the exported file

        // Append the link, click it, and remove it
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    });

    /* Legend */

    // Toggle Legend
    document.getElementById('toggle-legend').addEventListener('click', function() {
        if(aggregated) {
            return;
        } else {
            // Show the legend if the detailed network is shown
            var legend = document.getElementById('legend');
            legend.classList.toggle('hidden');
            this.classList.toggle('active');
        }
      });

    /* Highlight Button */

    // Default: Highlighting is deactivated
    let highlightingActive = false;

    // Toggle button to activate/deactivate highlighting
    document.getElementById('toggle-highlighting').addEventListener('click', function() {
        highlightingActive = !highlightingActive;
        this.classList.toggle('active', highlightingActive);
    });

    function updateLegendDisplay() {
        var legend = document.getElementById('legend');
        var legendGradient = document.getElementById('legend-gradient');

        if (useEdgeColoring || useEdgeTransparency) {
            // Set the actual min and max TOM values in the legend
            document.getElementById('minTOM').textContent = minTOM.toFixed(2);
            document.getElementById('maxTOM').textContent = maxTOM.toFixed(2);

            // Update gradient for coloring or transparency
            if (useEdgeColoring) {
                // Set Viridis gradient if edge coloring is enabled
                legendGradient.style.background = 'linear-gradient(to right, #440154, #3b528b, #21918c, #5dc962, #fde725)';
            } else if (useEdgeTransparency) {
                // Set transparency gradient (gray-to-opacity) if only transparency is enabled
                legendGradient.style.background = 'linear-gradient(to right, rgba(169, 169, 169, 0.1), rgba(169, 169, 169, 1))';
            }
        } else {
            // Hide the legend if neither edge coloring nor transparency is enabled
            legend.style.display = 'none';
        }
    }

    function updateEdgeStyles() {
        cy.batch(function () {
            if (aggregated) {
                // For aggregated view: use fixed edge width and map edge opacity based on overlap weight
                // Filter aggregated edges (assumed to have the "overlap" attribute defined)
                let aggregatedEdges = cy.edges().filter(edge => edge.data('overlap') !== undefined);
                // Compute global min and max weight among these edges
                let weights = aggregatedEdges.map(edge => edge.data('weight'));
                let minWeight = Math.min(...weights);
                let maxWeight = Math.max(...weights);
                // Define a fixed edge width for all aggregated edges
                const fixedEdgeWidth = 1; // Adjust as needed
    
                aggregatedEdges.forEach(edge => {
                    let weight = edge.data('weight');
                    // Calculate opacity mapping:
                    // If all edges have the same weight, set opacity to 1.
                    // Otherwise, map minWeight to 0.25 and maxWeight to 1 linearly.
                    let opacity = (minWeight === maxWeight) 
                        ? 1 
                        : 0.25 + ((weight - minWeight) / (maxWeight - minWeight)) * (1 - 0.25);
                    edge.style({
                        'width': fixedEdgeWidth,
                        'line-color': 'black',
                        'opacity': opacity
                    });
                });
            } else {
                // Detail view: apply D3-based scaling for edge color and transparency
                cy.edges().forEach(function (edge) {
                    let styleUpdate = {};
                    if (useEdgeColoring) {
                        let colorScale = d3.scaleSequential()
                            .domain([minTOM, maxTOM])
                            .interpolator(d3.interpolateViridis);
                        styleUpdate['line-color'] = colorScale(edge.data('weight'));
                    }
                    if (useEdgeTransparency) {
                        let opacityScale = d3.scaleLinear()
                            .domain([minTOM, maxTOM])
                            .range([0.1, 1]);
                        styleUpdate['opacity'] = opacityScale(edge.data('weight'));
                    }
                    // Set default values if no scaling is applied
                    if (!styleUpdate['line-color']) {
                        styleUpdate['line-color'] = 'gray';
                    }
                    if (!styleUpdate['opacity']) {
                        styleUpdate['opacity'] = 1;
                    }
                    edge.style(styleUpdate);
                });
            }
        });
    }      

    function formatList(value) {
        if (typeof value !== 'string' || value.trim() === '') {
            return 'N/A'; 
        }
        if (!value.includes(',')) {
            return value; 
        }
        let valuesArray = value.split(',');
        return valuesArray[0] + ', +' + (valuesArray.length - 1);
    }    

    /* Tooltip */

    /* Tooltip on hover */
    cy.on('mouseover', 'node', function(evt) {
        var node = evt.target;
        var tooltip = document.getElementById('tooltip');
        var tooltipText = '';

        if (node.data('aggregated')) {
            // For aggregated nodes, show only cluster-specific metadata
            tooltipText += '<b>Cluster:</b> ' + node.data('label') + '<br>';
            tooltipText += '<b>Average Edge Weight:</b> ' + node.data('avg_weight').toFixed(2) + '<br>';
            tooltipText += '<b>Node Count:</b> ' + node.data('node_count') + '<br>';
        } else {
            // For detailed nodes, iterate through METADATA_MAPPING
            var dataDict = METADATA_MAPPING;
            for (var key in dataDict) {
                if (HIDDEN_KEYS.includes(key)) continue;
                var value = node.data(key);
                if (value !== undefined && value !== null) {
                    let formattedValue = FORMATTED_KEYS.includes(key) ? formatList(value) : value;
                    tooltipText += `<b>${dataDict[key]}:</b> ${formattedValue}<br>`;
                }
            }
            tooltipText += '<b>Degree:</b> ' + node.degree() + '<br>';
            if (useClusterTooltip) { 
                tooltipText += '<b>Cluster:</b> ' + node.data('cluster');
            }
        }

        tooltip.style.display = 'block';
        tooltip.innerHTML = tooltipText;
        var mouseX = evt.originalEvent.clientX;
        var mouseY = evt.originalEvent.clientY;
        tooltip.style.left = (mouseX + 15) + 'px';
        tooltip.style.top = (mouseY) + 'px';
    });

    cy.on('mouseout', 'node', function(evt) {
        document.getElementById('tooltip').style.display = 'none';
    });

    /* Zooming node and edge sizes */
    const nodeSizeRange = document.getElementById('node-size-range');
    const minNodeSize = 1; // Minimum node size
    const maxNodeSize = 30; // Maximum node size

    let currentNodeSize = nodeSize;  
    let currentEdgeWidth = edgeWidth;

    function adjustNodeAndEdgeSize() {
        var zoomLevel = cy.zoom();
    
        // Calculate scaled sizes
        var scaledNodeSize = currentNodeSize / zoomLevel;
        var scaledEdgeWidth = currentEdgeWidth / zoomLevel;
    
        // For non-aggregated nodes: Adjust size
        cy.nodes().forEach(function(node) {
            if (node.data('aggregated')) return;
            node.style('width', scaledNodeSize);
            node.style('height', scaledNodeSize);
        });
    
        // For edges: Do not overwrite aggregated edges
        cy.edges().forEach(function(edge) {
            // If it is an aggregated edge (assumed overlap is always defined) – skip
            if (aggregated && edge.data('overlap') !== undefined) return;
            edge.style('width', scaledEdgeWidth);
        });
    }

    // Jitter effect for all metadata
    let currentMetadataIndex = 0;
    let jitterIntervals = {};

    function stopJitter() {
        Object.values(jitterIntervals).forEach(interval => clearInterval(interval));
        jitterIntervals = {};
    }

    function jitterNodes(nodes) {
        stopJitter();
        nodes.forEach(node => {
            let nodeId = node.id();
            jitterIntervals[nodeId] = setInterval(() => {
                let currentPosition = node.position();
                let jitterX = (Math.random() - 0.5) * 5;
                let jitterY = (Math.random() - 0.5) * 5;

                node.position({
                    x: currentPosition.x + jitterX,
                    y: currentPosition.y + jitterY
                });

                setTimeout(() => node.position(currentPosition), 50);
            }, 100);
        });
    }

    function highlightNodesByMetadata(metadataKey, selectedNode = null) {
        stopJitter();
    
        // Retain the currently selected node if none is explicitly passed
        if (!selectedNode) {
            selectedNode = lastClickedNode ? cy.$(`node[id="${lastClickedNode}"]`) : null;
        }
        
        if (!selectedNode || selectedNode.empty()) return;  // Abort if no node exists
    
        let metadataValue = metadataKey === "degree" ? selectedNode.degree() : selectedNode.data(metadataKey);
        if (metadataValue === undefined) return;
    
        let nodesToJitter = cy.nodes().filter(node => 
            metadataKey === "degree" ? node.degree() === metadataValue : node.data(metadataKey) === metadataValue
        );
    
        jitterNodes(nodesToJitter);
    }    

    /* Fullscreen */
    const container = document.documentElement; // Entire document for fullscreen
    const fullscreenButton = document.getElementById("toggle-fullscreen");
    
    function toggleFullscreen() {
        if (!document.fullscreenElement) {
            // Enter fullscreen mode
            if (container.requestFullscreen) {
                container.requestFullscreen();
            } else if (container.mozRequestFullScreen) {
                container.mozRequestFullScreen();
            } else if (container.webkitRequestFullscreen) {
                container.webkitRequestFullscreen();
            } else if (container.msRequestFullscreen) {
                container.msRequestFullscreen();
            }
            fullscreenButton.classList.add("active"); // Activate button design
        } else {
            // Exit fullscreen mode
            if (document.exitFullscreen) {
                document.exitFullscreen();
            } else if (document.mozCancelFullScreen) {
                document.mozCancelFullScreen();
            } else if (document.webkitExitFullscreen) {
                document.webkitExitFullscreen();
            } else if (document.msExitFullscreen) {
                document.msExitFullscreen();
            }
            fullscreenButton.classList.remove("active"); // Deactivate button design
        }
    }
    
    // Event listener for the fullscreen button
    fullscreenButton.addEventListener("click", function() {
        toggleFullscreen();
        cy.resize();
        cy.fit();
    });

    function restoreNodeAndEdgeSize(initialNodeSize, initialEdgeWidth) {
        // Restore node and edge sizes
        cy.nodes().forEach(function(node) {
            if (initialNodeSize) {
                node.style('width', initialNodeSize);
                node.style('height', initialNodeSize);
            }
        });

        cy.edges().forEach(function(edge) {
            if (initialEdgeWidth) {
                edge.style('width', initialEdgeWidth);
            }
        });
    }
    
    // Event listener for fullscreen change -> Resize and fit the graph, adjust button status
    document.addEventListener("fullscreenchange", function() {
        cy.resize(); // Ensure the canvas is resized
        cy.fit(); // Fit the graph to the container
        adjustNodeAndEdgeSize(); // Apply the new size

        // If fullscreen is exited, reset the button
        if (!document.fullscreenElement) {
            cy.zoom(initialZoom);  // Restore the initial zoom level
            cy.pan(initialPan);    // Reset to the initial pan position
            cy.center();           // Ensure the graph is centered correctly
            fullscreenButton.classList.remove("active");
        } 
    });

    const toggleMetadataButton = document.getElementById("toggle-metadata");
    const buttonTooltip = document.getElementById("button-tooltip");

    let metadataActive = false; // Flag to track if metadata highlighting is active

    toggleMetadataButton.addEventListener('click', function() {
        if (!metadataActive) {
            // Activate metadata highlighting
            metadataActive = true;
            this.classList.add('active'); // Add the active class to the button
            highlightNodesByMetadata(metadataKeys[currentMetadataIndex]); // Show the first metadata
        } else {
            currentMetadataIndex = (currentMetadataIndex + 1) % metadataKeys.length;
            if (currentMetadataIndex === 0) {
                // Deactivate metadata highlighting after cycling through all metadata
                metadataActive = false;
                this.classList.remove('active'); // Remove the active class
                stopJitter(); // Stop any ongoing jitter effect
            } else {
                highlightNodesByMetadata(metadataKeys[currentMetadataIndex]);
            }
        }

        updateTooltipText();
    });    

    toggleMetadataButton.addEventListener("mouseover", function(event) {
        updateTooltipText();
    });

    function updateTooltipText() {
        if (metadataActive) { // Show metadata only if active
            buttonTooltip.innerHTML = "Current Metadata: <b>" + METADATA_MAPPING[metadataKeys[currentMetadataIndex]] + "</b>";
        } else {
            buttonTooltip.innerHTML = "Current Metadata: <b>None</b>"; // Show "None" when inactive
        }
        buttonTooltip.style.display = "block";
    
        let rect = toggleMetadataButton.getBoundingClientRect();
        buttonTooltip.style.left = (rect.left + window.scrollX + 20) + "px";
        buttonTooltip.style.top = (rect.top + window.scrollY - 10) + "px";
    }

    toggleMetadataButton.addEventListener("mouseout", function() {
        buttonTooltip.style.display = "none";
    });

    cy.on('zoom', function() {
        adjustNodeAndEdgeSize();
    });

    /* Adjust Node Size */
    nodeSizeRange.value = currentNodeSize;

    nodeSizeRange.addEventListener('input', function() {
        var oldNodeSize = currentNodeSize; // Store the old node size
    
        currentNodeSize = parseInt(this.value);  // Get node size from slider
    
        // Calculate the scaling factor based on the old and new node sizes
        var scaleFactor = currentNodeSize / oldNodeSize;
    
        // Scale the edge width by the same factor as the node size
        currentEdgeWidth *= scaleFactor; 
    
        adjustNodeAndEdgeSize();
    });
    
    /* Highlighting nodes */

    // Event listener to receive messages from the parent window (main page)
    window.addEventListener('message', function(event) {
        cy.batch(function() {
            if (event.data && event.data.tableSelectedNode !== undefined) {
                var selectedNode = event.data.tableSelectedNode;

                // Only reset table-highlight on affected nodes
                cy.nodes('.table-highlight').removeClass('table-highlight');

                if (selectedNode) {
                    // Highlight the node selected from the table in a different color
                    var node = cy.$("node[id='" + selectedNode + "']");
                    if (!node.hasClass('table-highlight')) {
                        node.addClass('table-highlight');
                    }
                }
            }

            if (event.data && event.data.visibleNodes) {
                var visibleNodes = event.data.visibleNodes;

                // Remove highlight only from previously highlighted nodes
                cy.nodes('.highlight').removeClass('highlight');

                // Highlight only the nodes that are visible based on the search in the DataTable
                visibleNodes.forEach(function(nodeID) {
                    var node = cy.$("node[id='" + nodeID + "']");
                    if (!node.hasClass('highlight')) {
                        node.addClass('highlight');
                    }
                });
            }
        });
    }, false);

    // Define CSS styles for different types of highlights
    cy.style().selector('.highlight')
        .style({
            'background-color': 'magenta', 
            'line-color': 'magenta',
            'width': 4
        })
        .update();

    cy.style().selector('.table-highlight')
        .style({
            'background-color': 'yellow', // Use yellow to differentiate
            'line-color': 'yellow',
            'width': 4
        })
        .update();

    // Track the state of clicked nodes
    let lastClickedNode = null; // Track the last clicked node

    cy.on('tap', 'node', function(evt) {
        var node = evt.target;
        var nodeId = node.id();
        lastClickedNode = nodeId;

        if (highlightingActive) {
            cy.batch(function() {
                if (node.hasClass('highlight')) {
                    // Highlight the neighborhood of the clicked node
                    node.neighborhood().addClass('highlight');
                } else {
                    // Remove highlight from all nodes and highlight the clicked node
                    cy.nodes('.highlight').removeClass('highlight');
                    node.addClass('highlight');
                }
            });

            // Collect highlighted nodes
            var highlightedNodes = cy.nodes('.highlight').map(function(node) {
                return node.data('id'); 
            });

            // Send the highlighted nodes back to the main page
            window.parent.postMessage({ highlightedNodes: highlightedNodes }, '*');
        }

        if (metadataActive) { // Only jitter if metadata highlighting is active
            highlightNodesByMetadata(metadataKeys[currentMetadataIndex], node);
        }

        // If the node is aggregated, load its detail network
        if (node.data('aggregated')) {
            // Extract cluster name from label (assumes "ClusterName (n nodes)")
            var label = node.data('label');
            var clusterName = label.split(" (")[0]; 
            var safeClusterName = clusterName.replace(/\s/g, "_");
            // Build file name WITHOUT outputPath prefix (files are in the same folder)
            var detailNetworkFile = customFilename + "_detail_" + safeClusterName + ".json";

            spinner.style.display = 'flex'; 
            fetch(detailNetworkFile)
                .then(function(response) { return response.json(); })
                .then(function(detailData) {
                    // Replace current network with detailed network and run layout
                    cy.json({ elements: detailData });
                    cy.layout({
                        name: 'cose',
                        fit: true,  // Ensures the entire graph is visible
                        padding: 10,  // Adds padding around the layout
                        nodeRepulsion: 400000,  // Increase this value to spread nodes further apart
                        idealEdgeLength: 100,  // Sets an ideal edge length
                        nodeDimensionsIncludeLabels: false,  // Doesn't include labels when calculating size
                        animate: false,  // Animates the layout
                        stop: () => {
                            spinner.style.display = 'none';
                        }
                    }).run();
                    aggregated = false;  // Detailed network is now shown
                    // Show Legend
                    document.getElementById('legend').style.display = 'block';
                    resetParameters();
                })
                .catch(function(error) {
                    console.error("Error loading detail network:", error);
                });
        }
    });

    // When clicking anywhere in the Cytoscape plot (but not on nodes), reset the table
    cy.on('tap', function(event) {
        if (event.target === cy) { // If we clicked on an empty space (not a node)
            cy.batch(function() {
                // Reset the highlights
                cy.nodes('.highlight').removeClass('highlight');
                lastClickedNode = null; // Clear the last clicked node tracking
            });

            stopJitter(); // Stop the jitter effect

            // Send an empty array to the main page to reset the DataTable
            window.parent.postMessage({ highlightedNodes: [] }, '*');
        }
    });

    /* Rectangle Zoom */

    // Variables for rectangle zoom
    var isMouseDown = false;
    var startX, startY;
    var rectangleOverlayElement;
    var containerOffset;

    // Disable default interactions during rectangle zoom
    if (rectangleZoomMode) {
        cy.userPanningEnabled(false);
        cy.userZoomingEnabled(false);
        cy.container().style.cursor = 'crosshair';
    }

    // Determine container offset for coordinate calculation
    var rect = cy.container().getBoundingClientRect();
    containerOffset = {
        left: rect.left + window.scrollX,
        top: rect.top + window.scrollY
    };

    // Event listener for the rectangle zoom button
    document.getElementById('rectangle-zoom').addEventListener('click', function() {
        if (rectangleZoomMode) {
            // Deactivate rectangle zoom mode
            rectangleZoomMode = false;
            this.classList.remove('active'); // Update icon color
            cy.userPanningEnabled(true);
            cy.userZoomingEnabled(scrollZoomActive); // Restore scroll zoom state
            cy.container().style.cursor = 'default';
        } else {
            // Activate rectangle zoom mode
            rectangleZoomMode = true;
            this.classList.add('active'); // Update icon color
            cy.userPanningEnabled(false);
            cy.userZoomingEnabled(false);
            cy.container().style.cursor = 'crosshair';
    
            var rect = cy.container().getBoundingClientRect();
            containerOffset = {
                left: rect.left + window.scrollX,
                top: rect.top + window.scrollY
            };

            // Deactivate scroll zoom if active
            if (scrollZoomActive) {
                document.getElementById('toggleZoom').classList.remove('active');
                cy.userZoomingEnabled(false);
                scrollZoomActive = false;
            }
        }
    });     

    // Mouse down event to start drawing
    cy.on('mousedown', function(event) {
        if (!rectangleZoomMode) return;
        if (event.target !== cy) return; // Only proceed if clicking on the background

        isMouseDown = true;

        startX = event.originalEvent.pageX - containerOffset.left;
        startY = event.originalEvent.pageY - containerOffset.top;

        // Create overlay element for the rectangle
        rectangleOverlayElement = document.createElement('div');
        rectangleOverlayElement.style.position = 'absolute';
        rectangleOverlayElement.style.border = '1px dashed #000';
        rectangleOverlayElement.style.backgroundColor = 'rgba(0, 0, 255, 0.1)';
        rectangleOverlayElement.style.pointerEvents = 'none';
        rectangleOverlayElement.style.zIndex = '9999';

        // Add to Cytoscape container
        cy.container().appendChild(rectangleOverlayElement);
    });

    // Mouse move event to update the rectangle dimensions
    cy.on('mousemove', function(event) {
        if (!rectangleZoomMode || !isMouseDown) return;

        var currentX = event.originalEvent.pageX - containerOffset.left;
        var currentY = event.originalEvent.pageY - containerOffset.top;

        var x = Math.min(startX, currentX);
        var y = Math.min(startY, currentY);
        var width = Math.abs(startX - currentX);
        var height = Math.abs(startY - currentY);

        rectangleOverlayElement.style.left = x + 'px';
        rectangleOverlayElement.style.top = y + 'px';
        rectangleOverlayElement.style.width = width + 'px';
        rectangleOverlayElement.style.height = height + 'px';
    });

    // Mouse up event to perform the zoom
    cy.on('mouseup', function(event) {
        if (!rectangleZoomMode || !isMouseDown) return;

        isMouseDown = false;

        var endX = event.originalEvent.pageX - containerOffset.left;
        var endY = event.originalEvent.pageY - containerOffset.top;

        // Remove the rectangle overlay
        if (rectangleOverlayElement && rectangleOverlayElement.parentNode) {
            rectangleOverlayElement.parentNode.removeChild(rectangleOverlayElement);
        }

        // Prevent zoom if rectangle is too small (e.g., a click instead of a drag)
        if (Math.abs(endX - startX) < 5 || Math.abs(endY - startY) < 5) {
            return;
        }

        // Calculate the bounding box in model coordinates
        var x1 = (Math.min(startX, endX) - cy.pan().x) / cy.zoom();
        var y1 = (Math.min(startY, endY) - cy.pan().y) / cy.zoom();
        var x2 = (Math.max(startX, endX) - cy.pan().x) / cy.zoom();
        var y2 = (Math.max(startY, endY) - cy.pan().y) / cy.zoom();

        var boundingBox = { x1: x1, y1: y1, x2: x2, y2: y2 };

        // Zoom into the bounding box
        cy.fit(boundingBox);
    });

    document.getElementById('export-png').addEventListener('click', function() {
        var exportScale = 2; // Higher resolution (2x)
        var cyImageData = cy.png({ full: false, scale: exportScale }); // Only visible area!
    
        var cyImage = new Image();
        cyImage.src = cyImageData;
    
        cyImage.onload = function() {
            createExportImage(cyImage);
        };
    });
    
    /* Export Image */
    function createExportImage(cyImage) {
        var exportCanvas = document.createElement('canvas');
        var ctx = exportCanvas.getContext('2d');
    
        // Define canvas size based on Cytoscape container
        var width = cy.container().clientWidth * 2;
        var height = cy.container().clientHeight * 2;
    
        exportCanvas.width = width;
        exportCanvas.height = height;
    
        // Set background color
        ctx.fillStyle = 'white';
        ctx.fillRect(0, 0, width, height);
    
        // Draw Cytoscape graph
        ctx.drawImage(cyImage, 0, 0, width, height);
    
        // Save as PNG
        var link = document.createElement('a');
        link.download = 'cytoscape_export.png';
        link.href = exportCanvas.toDataURL('image/png');
        link.click();
    }
    
    function resetParameters() {
        // Clear the sortedEdges cache so that edge filtering is recalculated
        cy.scratch('sortedEdges', null);
        
        // Reset global size variables to the initial values provided by Python
        currentNodeSize = nodeSize;
        currentEdgeWidth = edgeWidth;
        
        // Reset node sizes for all nodes:
        cy.nodes().forEach(function(node) {
            if (node.data('aggregated')) {
                // For aggregated nodes
                node.style('width', aggNodeSize);
                node.style('height', aggNodeSize);
            } else {
                // For detail nodes, use the default nodeSize
                node.style('width', nodeSize);
                node.style('height', nodeSize);
            }
        });
              
        // Update Species Legend
        assignOrganismShapes();

        // Re-run the edge visibility update (this will rebuild the sortedEdges cache)
        updateEdgesVisibility();
        
        // Calculate the minimum and maximum TOM value
        minTOM = Math.min(...cy.edges().map(edge => edge.data('weight')));
        maxTOM = Math.max(...cy.edges().map(edge => edge.data('weight')));

        updateEdgeStyles();
        updateLegendDisplay();

        adjustNodeAndEdgeSize();

        // Reset zoom and pan if not in fullscreen
        cy.resize();
        cy.fit();
    }

    resetParameters();
    assignDynamicNetworkIconToAggregatedNodes();

    // Wait 1 second than hide the spinner
    setTimeout(function() {
        spinner.style.display = 'none';
    }
    , 1000);

});