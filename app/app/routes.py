# routes.py

from flask import Blueprint, make_response, request, jsonify, send_from_directory, Response, current_app, render_template
import requests
import os
import uuid
import glob
import shutil
import tempfile
import json
import logging
import redis
import re
import scanpy as sc
import pandas as pd
import numpy as np
from PIL import Image
from sqlalchemy import create_engine, text
from celery import states
from datetime import datetime
from zipfile import ZipFile
from config import Config
from .task_mapping import analysis_tasks
from .factory import celery

# Initialize Blueprint
api_bp = Blueprint('api', __name__)

# Initialize Redis client
redis_client = redis.StrictRedis(host=os.getenv("REDIS_HOST", "redis"),
                                 port=int(os.getenv("REDIS_PORT", 6379)), 
                                 db=int(os.getenv("REDIS_DB", 0)), decode_responses=True)

# Configure logger
logger = logging.getLogger(__name__)

# Initialize database connection
engine = create_engine(f"postgresql+psycopg2://{os.getenv('DB_USER')}:{os.getenv('DB_PASSWORD')}@{os.getenv('DB_HOST')}:{os.getenv('DB_PORT')}/{os.getenv('DB_NAME')}")

"""
File serving route and index route
"""

@api_bp.route('/output/<env>/<path:filename>')
def serve_output_file(env, filename):
    parent_dir = os.path.dirname(current_app.root_path)
    output_dir = os.path.join(parent_dir, 'output', env)
    print(f"Serving from output_dir: {output_dir}")
    print(f"Filename: {filename}")
    return send_from_directory(output_dir, filename)

"""
Index route
"""

@api_bp.route('/')
def index():
    """
    Generates and renders the index HTML based on the mask type.
    """
    
    # Get mask_type from environment variable APP_NAME
    mask_type = os.environ.get("APP_NAME", "RanOmics")
    print(f"Generating index for mask type: {mask_type}")

    # Analysis options
    analysis_options = {
        "Common": {},
        "Single": {
            "plot_module_props": "Transcripts per Module",
            "plot_module_orthos": "Orthogroups per Module",
            "plot_module_hubs": "Top Hub Transcripts per Module",
            "plot_module_goea": "GOEA of each Module",
            "plot_correlation_network": "Module Eigengenes Correlation Network",
        },
        "Multiple": {
            "plot_jaccard_tissues": "Tissues Similarity Network",
            "plot_tissues_corr": "Tissues Similarity Heatmap",
            "plot_jaccard_modules": "Module Similarity Network",
            "plot_modules_corr": "Module Similarity Heatmap",
            "plot_species_correlation_network": "Module Membership Correlation Network",
            "plot_hog_level": "HOGs Level Distribution",
        },
    }

    # Browser options
    browser_options = {
        "Common": {
            "plot_co_expression_network": "Co-expression Network Analysis",
            "plot_filtered_jaccard_modules": "Module Similarity Network",
        },
        "Single": {
            "plot_go_terms": "Gene Ontology Enrichment Analysis",
        },
        "Multiple": {
            "plot_venn": "Venn Diagram",
            "plot_upset": "UpSet Plot",
            "plot_filtered_jaccard_species": "Species Similarity Network",
        },
    }

    # Flatten the analysis options
    flat_analysis_options = {}
    for category, options in analysis_options.items():
        for key, description in options.items():
            new_key = f"{key}_{category.lower()}"
            flat_analysis_options[new_key] = description

    # Flatten the browser options
    flat_browser_options = {}
    for category, options in browser_options.items():
        for key, description in options.items():
            new_key = f"{key}_{category.lower()}"
            flat_browser_options[new_key] = description

    # RanOmics species names
    species_names = [
        "Papaver somniferum",
        "Thalictrum thalictroides",
        "Aquilegia caerulea",
        "Staphisagria picta",
        "Capnoides sempervirens",
        "Hydrastis canadensis",
        "Eschscholzia californica",
        "Epimedium grandiflorum"
    ]

    # Function to generate codes and create the species list
    def get_species_list(names):
        return [("".join([word[0].upper() for word in name.split(" ")]), name) for name in names]

    # Load species information based on the mask type
    if mask_type == "RanOmics":
        species = get_species_list(species_names)
    elif mask_type == "Ceratopteris":
        species = [("CR", "Ceratopteris richardii")]
    elif mask_type == "Development":
        # Add "Ceratopteris richardii" to the species list
        species = get_species_list(species_names + ["Ceratopteris richardii"])
    elif mask_type == "Custom":
        # Load species names and abbreviations from the h5ad files
        species = []
        for file in glob.glob(os.path.join(Config.H5AD_DIR, "*.h5ad")):
            adata = sc.read_h5ad(file)
            species.append((adata.uns.get("name", "Unknown"), adata.uns.get("species", "Unknown")))
    else:
        raise ValueError(f"Unknown mask type: {mask_type}")

    # Sort species by name
    species = sorted(species, key=lambda x: x[1])

    # Render the template with the provided context
    return render_template(
        'template.html',
        species=species,
        analysis_options=flat_analysis_options,
        browser_options=flat_browser_options,
        BASE_URL=Config.BASE_URL,
        OUTPUT_DIR=Config.OUTPUT_DIR
    )

"""
Blueprints registration
"""

def register_routes(app):
    """
    Register all Blueprints with the Flask app.

    Parameters:
    - app (Flask): The Flask application instance.
    """

    app.register_blueprint(api_bp)

"""
Celery task management routes
"""

@api_bp.route('/api/check_status/<task_id>', methods=['GET'])
def check_status(task_id):
    """
    Endpoint to check the status of a Celery task.

    Parameters:
    - task_id (str): The ID of the Celery task.

    Returns:
    - JSON response with task state and result if available.
    """

    task = celery.AsyncResult(task_id)
    if task.state == states.PENDING:
        response = {
            'state': task.state,
            'status': 'Pending...'
        }
    elif task.state != states.FAILURE:
        response = {
            'state': task.state,
            'status': task.info.get('status', '')
        }
        if 'result' in task.info:
            response['result'] = task.info['result']
    else:
        response = {
            'state': task.state,
            'status': str(task.info),  # Exception info
        }

    return jsonify(response)

@api_bp.route('/api/start_analysis', methods=['POST'])
def start_analysis():
    """
    Endpoint to start various types of analyses based on the provided analysis_type.

    Expects JSON data with 'analysis_type' and other relevant parameters.

    Returns:
    - The Celery task ID.
    """

    data = request.json

    analysis_type = data.get('analysis_type')
    if not analysis_type:
        return jsonify({"status": "error", "message": "analysis_type is required"}), 400

    task = analysis_tasks.get(analysis_type)
    if not task:
        return jsonify({"status": "error", "message": "Invalid analysis type"}), 400

    celery_task = task.apply_async(args=[data])
    return jsonify({"task_id": celery_task.id}), 202

"""
Load image data
"""

@api_bp.route('/api/load_image_data', methods=['POST'])
def load_image_data():
    data = request.get_json()
    plant = data['plant']
    category = data['category']

    print(f"Loading image data for plant {plant}, category {category}")

    image_dir = os.path.join(current_app.static_folder, 'images', plant)
    thumbnail_dir = os.path.join(image_dir, 'thumbnails')
    os.makedirs(thumbnail_dir, exist_ok=True)

    if category == 'general':
        all_files = glob.glob(os.path.join(image_dir, '*'))
        figures = [
            f for f in all_files
            if os.path.isfile(f) and not os.path.basename(f).startswith(('module_barplot', 'module_heatmap'))
        ]
    else:
        figure_path = os.path.join(image_dir, f"{category}*")
        figures = [f for f in glob.glob(figure_path) if os.path.isfile(f)]

    if figures:
        figure_urls = []
        for fig in figures:
            original_filename = os.path.basename(fig)
            thumbnail_path = os.path.join(thumbnail_dir, original_filename)

            if not os.path.exists(thumbnail_path):
                with Image.open(fig) as img:
                    img.thumbnail((500, 500))
                    img.save(thumbnail_path, "PNG", quality=300)

            thumbnail_url = f"{Config.BASE_URL}/static/images/{plant}/thumbnails/{original_filename}"
            original_url = f"{Config.BASE_URL}/static/images/{plant}/{original_filename}"

            figure_urls.append({"thumbnail": thumbnail_url, "original": original_url})

        return jsonify({"status": "success", "figures": figure_urls})
    else:
        return jsonify({"status": "error", "message": "No images found"}), 404


"""
AI Search
"""

CYFISH_API_URL = os.environ.get("CYFISH_API_URL", "http://172.17.0.1:5004/api/transcripts")
logger.info(f"Using CyFISH API URL: {CYFISH_API_URL}")

@api_bp.route('/ai-search', methods=['POST'])
def ai_search():
    data = request.json
    res = requests.post(CYFISH_API_URL, json=data)

    return jsonify(res.json())

"""
Start table and plot tasks
"""

def start_task(task_name):
    """
    Starts a Celery task based on the given task_name.

    Parameters:
    - task_name (str): The key in the analysis_tasks dictionary.

    Returns:
    - JSON response with task_id or error message.
    """

    data = request.json
    task = analysis_tasks.get(task_name)

    if not task:
        return jsonify({"status": "error", "message": f"Invalid task type: {task_name}"}), 400

    if data:
        celery_task = task.apply_async(args=[data])
    else:
        celery_task = task.apply_async()
    
    celery_task = task.apply_async(args=[data])

    return jsonify({"task_id": celery_task.id}), 202

"""
Table routes helper functions
"""

def force_select(query: str) -> str:
    # Only apply if SELECT * is not already used
    if re.search(r'select\s+\*\s+from', query, re.IGNORECASE):
        return query

    # If aggregations or GROUP BY are present, do not modify
    if any(kw in query.lower() for kw in ["group by", "having", "count(", "sum(", "avg(", "min(", "max("]):
        return query

    # Replace SELECT ... FROM with SELECT * FROM
    select_pattern = re.compile(r'^select\s+.+?\s+from\s', re.IGNORECASE | re.DOTALL)
    return select_pattern.sub('SELECT * FROM ', query)


def parse_request_params(request) -> dict:
    """
    Parse request parameters for data tables.

    Parameters:
    - request (Request): The Flask request object.

    Returns:
    - dict: Parsed parameters including start, length, search value, columns, species filter,
            transcript filter (if provided), and order.
    """

    params = {}
    if request.is_json:
        values = request.get_json()
    else:
        values = request.values

    # Process start and length
    if values.get('get_all'):
        params['start'] = None  # No limit
        params['length'] = None  # No offset
    else:
        params['start'] = int(values.get('start', 0))
        params['length'] = int(values.get('length', 10))

    # Process global search value
    if request.is_json:
        params['search_value'] = values.get('search', {}).get('value', '')
    else:
        params['search_value'] = values.get('search[value]', '')

    # Process columns
    params['columns'] = []
    if request.is_json:
        columns = values.get('columns', [])
        for i, col in enumerate(columns):
            params['columns'].append({
                'data': col.get('data'),
                'name': col.get('data'),
                'searchable': col.get('searchable', False),
                'orderable': col.get('orderable', False),
                'search_value': col.get('search', {}).get('value', '')
            })
    else:
        i = 0
        while True:
            col_data = values.get(f'columns[{i}][data]')
            if col_data is None:
                break
            col = {
                'data': col_data,
                'name': col_data,
                'searchable': values.get(f'columns[{i}][searchable]') == 'true',
                'orderable': values.get(f'columns[{i}][orderable]') == 'true',
                'search_value': values.get(f'columns[{i}][search][value]', '')
            }
            params['columns'].append(col)
            i += 1

    # Process species filter
    species_filter_str = values.get('species_filter', '[]')
    try:
        params['species_filter'] = json.loads(species_filter_str)
    except json.JSONDecodeError:
        params['species_filter'] = []

    # Process transcript filter
    params['transcript_filter'] = values.get('transcript_filter', '')
    params['transcript_filter_column'] = values.get('transcript_filter_column', 'transcript')

    # Process order parameters
    params['order'] = []
    if request.is_json:
        orders = values.get('order', [])
        for order in orders:
            order_col_index = int(order.get('column'))
            order_dir = order.get('dir', 'asc')
            params['order'].append({
                'column': order_col_index,
                'dir': order_dir
            })
    else:
        i = 0
        while True:
            order_col_index = values.get(f'order[{i}][column]')
            if order_col_index is None:
                break
            order = {
                'column': int(order_col_index),
                'dir': values.get(f'order[{i}][dir]', 'asc')
            }
            params['order'].append(order)
            i += 1

    params['table_name'] = values.get('table_name', Config.BROWSER_DB)

    return params


def build_sql_query(params: dict, select_columns: str) -> tuple:
    """
    Build SQL query based on parsed parameters.

    Parameters:
    - params (dict): Parsed parameters including start, length, search value, columns, species filter,
                     transcript filter, and order.
    - select_columns (str): Columns to select in the SQL query.

    Returns:
    - tuple: SQL query string and search parameters.
    """
    
    where_conditions = []
    search_params = {}

    # Global search with OR for multiple tokens
    if params['search_value']:
        tokens = [token for token in re.split(r'[\s,]+', params['search_value'].strip()) if token]
        global_search_conditions = []
        for idx, col in enumerate(params['columns']):
            if col['searchable']:
                col_name = col['name']
                token_conditions = []
                for j, token in enumerate(tokens):
                    param_name = f"search_value_{idx}_{j}"
                    token_conditions.append(f"{col_name}::text ILIKE :{param_name}")
                    search_params[param_name] = f"%{token}%"
                global_search_conditions.append("(" + " OR ".join(token_conditions) + ")")
        if global_search_conditions:
            where_conditions.append('(' + ' OR '.join(global_search_conditions) + ')')

    # Column-specific search
    for idx, col in enumerate(params['columns']):
        if col['searchable'] and col['search_value']:
            col_name = col['name']
            param_name = f"col_search_value_{idx}"
            where_conditions.append(f"{col_name}::text ILIKE :{param_name}")
            search_params[param_name] = f"%{col['search_value']}%"

    # Add species filter if provided
    if params['species_filter']:
        species_placeholders = ', '.join([f":species_{i}" for i in range(len(params['species_filter']))])
        where_conditions.append(f"species IN ({species_placeholders})")
        for i, species in enumerate(params['species_filter']):
            search_params[f"species_{i}"] = species

    # Extended filter: transcript_filter (general filter_column) with list logic
    transcript_filter = params.get('transcript_filter', '').strip()
    if transcript_filter:
        # Use filter column, default: 'transcript'
        filter_column = params.get('transcript_filter_column', 'transcript')
        # Validation: only use allowed columns
        allowed_columns = [col['name'] for col in params['columns']]
        if filter_column not in allowed_columns:
            filter_column = 'transcript'

        tokens = [token for token in re.split(r'[\s,]+', transcript_filter) if token]

        if tokens:
            # Pass search tokens as an array (in lowercase)
            search_params['transcript_tokens'] = [token.lower() for token in tokens]
            # If commas are present, use array overlap (&&) - it is sufficient if at least one token matches.
            condition = (
                f"((POSITION(',' IN {filter_column}) > 0 AND "
                f"string_to_array(LOWER({filter_column}), ',') && :transcript_tokens) OR "
                f"(POSITION(',' IN {filter_column}) = 0 AND LOWER({filter_column}) = ANY(:transcript_tokens)))"
            )
            where_conditions.append(condition)

    where_clause = 'WHERE ' + ' AND '.join(where_conditions) if where_conditions else ''

    # Create ORDER BY clause
    order_params = []
    for order in params['order']:
        col = params['columns'][order['column']]
        if col['orderable']:
            col_name = col['name']
            order_dir = 'ASC' if order['dir'].lower() == 'asc' else 'DESC'
            order_params.append(f"{col_name} {order_dir}")

    order_clause = 'ORDER BY ' + ', '.join(order_params) if order_params else ''

    # Add LIMIT and OFFSET only if 'length' and 'start' are set
    limit_offset = ''
    if params.get('length') is not None and params.get('start') is not None:
        limit_offset = "LIMIT :length OFFSET :start"
        search_params['length'] = params['length']
        search_params['start'] = params['start']

    query = f"SELECT {select_columns} FROM {params['table_name']} {where_clause} {order_clause} {limit_offset}"

    return query, search_params


def execute_sql_query(query: str, params: dict) -> pd.DataFrame:
    """
    Execute SQL query and return results as a DataFrame.

    Parameters:
    - query (str): SQL query string.
    - params (dict): Search parameters for the SQL query.

    Returns:
    - pd.DataFrame: Query results as a DataFrame.
    """

    with engine.connect() as conn:
        df = pd.read_sql(text(query), conn, params=params)

    return df


"""
Table Routes
"""

@api_bp.route('/api/data', methods=['GET', 'POST'])
def data():
    try:
        params = parse_request_params(request)
        table = params['table_name']
        select_columns = ', '.join([col['name'] for col in params['columns']])

        ai_query = request.json.get("ai_query")
        order_clause = ""
        limit_offset_clause = ""

        if ai_query:
            # Security: Only allow SELECT statements
            if not ai_query.lower().strip().startswith("select"):
                return jsonify({"error": "Only SELECT statements allowed."}), 400

            base_query = ai_query.strip().rstrip(";")
            base_query = force_select(base_query)

            # Optionally extract ORDER BY from request
            order_params = []
            for order in params.get("order", []):
                col = params["columns"][order["column"]]
                if col.get("orderable", True):
                    col_name = col["name"]
                    dir_sql = "ASC" if order["dir"].lower() == "asc" else "DESC"
                    order_params.append(f"{col_name} {dir_sql}")
            if order_params:
                order_clause = "ORDER BY " + ", ".join(order_params)

            # LIMIT + OFFSET
            limit = int(params.get("length", 10))
            offset = int(params.get("start", 0))

            # Check if LIMIT is already in the base query
            if re.search(r'\blimit\b', base_query, re.IGNORECASE):
                limit_offset_clause = ""
            else:
                limit_offset_clause = f"LIMIT {limit} OFFSET {offset}"

            query = f"{base_query} {order_clause} {limit_offset_clause}"

            query = f"{base_query} {order_clause} {limit_offset_clause}"
            query_params = {}

            logger.info("AI query: %s", query)
            logger.info("AI query params: %s", query_params)

        else:
            # Default case
            query, query_params = build_sql_query(params, select_columns)

        logger.info("Received transcript_filter: %s", params.get('transcript_filter'))
        logger.info("Received transcript_filter_column: %s", params.get('transcript_filter_column'))
        logger.info("Columns in request: %s", [col['name'] for col in params['columns']])

        # --- total_records + total_filtered_records ---
        if ai_query:
            # AI queries: count via subquery
            count_query = f"SELECT COUNT(*) FROM ({base_query}) AS subquery"
            with engine.connect() as conn:
                total_filtered_records = pd.read_sql(text(count_query), conn).iloc[0, 0]
                total_records = total_filtered_records
        else:
            with engine.connect() as conn:
                total_records = pd.read_sql(text(f"SELECT COUNT(*) FROM {table}"), conn).iloc[0, 0]
            count_query = f"SELECT COUNT(*) FROM {table} {query.partition(f'FROM {table}')[2].partition('ORDER BY')[0]}"
            with engine.connect() as conn:
                total_filtered_records = pd.read_sql(text(count_query), conn, params=query_params).iloc[0, 0]

        # --- Fetch data ---
        df = execute_sql_query(query, query_params)
        data = df.replace({np.nan: None}).to_dict(orient='records')

        return jsonify({
            'draw': int(request.values.get('draw', 1)),
            'recordsTotal': int(total_records),
            'recordsFiltered': int(total_filtered_records),
            'data': data
        })

    except Exception as e:
        print(f"Error: {e}")
        return jsonify({"error": "Server Error"}), 500


@api_bp.route('/api/get_transcripts', methods=['POST'])
def get_transcripts():
    try:
        params = parse_request_params(request)
        select_columns = 'transcript, species'
        params['start'] = None
        params['length'] = None

        ai_query = request.json.get("ai_query")
        if ai_query:
            ai_query = force_select(ai_query.strip().rstrip(";"))
            query = ai_query
            query_params = {}
        else:
            query, query_params = build_sql_query(params, select_columns)

        df = execute_sql_query(query, query_params)
        transcripts_dict = df.groupby('species')['transcript'].apply(list).to_dict()

        return jsonify({'transcripts': transcripts_dict})
    
    except Exception as e:
        print(f"Error: {e}")
        return jsonify({'error': 'Server Error'}), 500


@api_bp.route('/api/get_column_entries', methods=['POST'])
def get_column_entries():
    try:
        params = parse_request_params(request)
        selected_column = request.json.get('selected_column')
        unique = request.json.get('unique', False)
        params['start'] = None
        params['length'] = None

        ai_query = request.json.get("ai_query")
        if ai_query:
            ai_query = force_select(ai_query.strip().rstrip(";"))
            query = ai_query
            query_params = {}
        else:
            query, query_params = build_sql_query(params, selected_column)

        df = execute_sql_query(query, query_params)
        entries = df[selected_column].dropna().astype(str)
        entries = entries[entries.str.strip() != '']

        if unique:
            entries = entries.drop_duplicates()

        return jsonify({'entries': entries.tolist()})

    except Exception as e:
        print(f"Error: {e}")
        return jsonify({'error': 'Server Error'}), 500


@api_bp.route('/api/export_data', methods=['POST'])
def export_data():
    try:
        params = parse_request_params(request)
        select_columns = ', '.join([col['name'] for col in params['columns']])
        params['start'] = None
        params['length'] = None

        ai_query = request.json.get("ai_query")
        if ai_query:
            ai_query = force_select(ai_query.strip().rstrip(";"))
            query = ai_query
            query_params = {}
        else:
            query, query_params = build_sql_query(params, select_columns)

        df = execute_sql_query(query, query_params)
        csv_data = df.to_csv(index=False, sep='\t')
        response = make_response(csv_data)
        response.headers['Content-Disposition'] = 'attachment; filename=export.csv'
        response.headers['Content-Type'] = 'text/csv'

        return response

    except Exception as e:
        print(f"Error: {e}")
        return jsonify({'error': 'Server Error'}), 500


@api_bp.route('/api/browser')
def browser():
    context = {
        "BASE_URL": Config.BASE_URL,
        "BROWSER_DB": Config.BROWSER_DB,
        "METADATA_DICT": Config.METADATA_DICT
    }

    return render_template('browser_template.html', **context)


@api_bp.route('/api/general_info')
def general_info():
    context = {
        "BASE_URL": Config.BASE_URL,
        "INFO_DB": Config.INFO_DB
    }
    return render_template('general_info_template.html', **context)

"""
Browser Routes
"""

@api_bp.route('/api/plot_co_expression_network', methods=['POST'])
def plot_co_expression_network():
    return start_task('plot_co_expression_network')


@api_bp.route('/api/plot_go_terms', methods=['POST'])
def plot_go_terms():
    return start_task('plot_go_terms')


@api_bp.route('/api/plot_venn', methods=['POST'])
def plot_venn():
    return start_task('plot_venn')


@api_bp.route('/api/plot_upset', methods=['POST'])
def plot_upset():
    return start_task('plot_upset')


@api_bp.route('/api/plot_filtered_jaccard_modules', methods=['POST'])
def plot_filtered_jaccard_modules():
    return start_task('plot_filtered_jaccard_modules')


@api_bp.route('/api/plot_filtered_jaccard_species', methods=['POST'])
def plot_filtered_jaccard_species():
    return start_task('plot_filtered_jaccard_species')


"""
Single Plant Routes
"""

@api_bp.route('/api/plot_module_goea', methods=['POST'])
def plot_module_goea():
    return start_task('plot_module_goea')


@api_bp.route('/api/plot_correlation_network', methods=['POST'])
def plot_correlation_network():
    return start_task('plot_correlation_network')


@api_bp.route('/api/plot_module_props', methods=['POST'])
def plot_module_props():
    return start_task('plot_module_props')


@api_bp.route('/api/plot_module_orthos', methods=['POST'])
def plot_module_orthos():
    return start_task('plot_module_orthos')


@api_bp.route('/api/plot_module_hubs', methods=['POST'])
def plot_module_hubs():
    return start_task('plot_module_hubs')


"""
Multiple Plant Routes
"""

@api_bp.route('/api/plot_jaccard_tissues', methods=['POST'])
def plot_jaccard_tissues():
    return start_task('plot_jaccard_tissues')


@api_bp.route('/api/plot_tissues_corr', methods=['POST'])
def plot_tissues_corr():
    return start_task('plot_tissues_corr')


@api_bp.route('/api/plot_jaccard_modules', methods=['POST'])
def plot_jaccard_modules():
    return start_task('plot_jaccard_modules')


@api_bp.route('/api/plot_modules_corr', methods=['POST'])
def plot_modules_corr():
    return start_task('plot_modules_corr')


@api_bp.route('/api/plot_hog_level', methods=['POST'])
def plot_hog_level():
    return start_task('plot_hog_level')


@api_bp.route('/api/plot_species_correlation_network', methods=['POST'])
def plot_species_correlation_network():
    return start_task('plot_species_correlation_network')

"""
Session Management
"""

@api_bp.route('/api/init_session', methods=['POST'])
def init_session():
    session_id = request.headers.get('Session-ID')

    if not session_id:
        session_id = str(uuid.uuid4()) 

    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)
    os.makedirs(user_dir, exist_ok=True)

    logger.info(f"User directory: {user_dir}, Session ID: {session_id}")

    return jsonify({"status": "success", "message": "Session and user directory created successfully.", "session_id": session_id})

@api_bp.route('/api/heartbeat', methods=['POST'])
def heartbeat():
    data = request.get_json()
    session_id = data.get('session_id')

    if session_id:
        current_time = datetime.now().timestamp()
        redis_client.hset('last_heartbeat', session_id, current_time)
        
        print(f"Received heartbeat from session: {session_id}")
        return '', 204
    else:
        print("Heartbeat received with no session ID")
        return jsonify({"status": "error", "message": "Session ID missing"}), 400

# Cleanup session directory before unloading
@api_bp.route('/api/cleanup', methods=['POST'])
def cleanup():
    data = request.get_json(force=True)  # force=True to ignore content type
    session_id = data.get('session_id')
    
    if session_id:
        # Check if the session ID exists in Redis
        if redis_client.hexists('last_heartbeat', session_id):
            user_dir = os.path.join(Config.OUTPUT_DIR, session_id)

            if os.path.exists(user_dir):
                shutil.rmtree(user_dir)
                print(user_dir, "removed")

            # Remove the session from Redis
            redis_client.hdel('last_heartbeat', session_id)
            return jsonify({"status": "success", "message": "Cleanup successful"})
        else:
            return jsonify({"status": "error", "message": "Session not found"}), 404
    else:
        return jsonify({"status": "error", "message": "Session ID missing"}), 400

# Clear session files button
@api_bp.route('/api/cleanup_session_files', methods=['POST'])
def cleanup_session_files():
    data = request.get_json()
    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)

    if os.path.exists(user_dir):
        for file in os.listdir(user_dir):
            file_path = os.path.join(user_dir, file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                print(e)
        return jsonify({"status": "success", "message": "Session files cleaned up"})
    return jsonify({"status": "error", "message": "Session ID missing"}), 400

"""
List and download files
"""

# List all files in a session directory
@api_bp.route('/api/list_session_files', methods=['POST'])
def list_session_files():
    data = request.get_json()
    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)

    if not os.path.exists(user_dir):
        return jsonify({"status": "error", "message": "No files found"}), 404

    files = [f for f in os.listdir(user_dir) if os.path.isfile(os.path.join(user_dir, f))]
    filtered_files = [f for f in files if not (f.startswith('browser') or f.startswith('general'))]

    return jsonify({"status": "success", "files": filtered_files})

# Download results as a ZIP file
@api_bp.route('/api/download_results', methods=['POST'])
def download_results():
    data = request.get_json()
    session_id = data.get('session_id')
    extension = data.get('extension', 'all')  # Standardmäßig alle Dateien
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)

    if not os.path.exists(user_dir):
        return jsonify({"status": "error", "message": "No files to download"}), 404

    temp_zip = tempfile.NamedTemporaryFile(mode='w+b', suffix='.zip', delete=False)
    
    with ZipFile(temp_zip.name, 'w') as zipf:
        for root, dirs, files in os.walk(user_dir):
            for file in files:
                # Datei darf nicht mit 'browser' oder 'general' anfangen
                if not (file.startswith('browser') or file.startswith('general')):
                    # Dateien nach der gewählten Erweiterung filtern
                    if extension == 'all' or file.endswith(extension):
                        zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), user_dir))

    temp_zip.close()

    def generate():
        with open(temp_zip.name, 'rb') as f:
            yield from f  # Inhalt der Datei streamen
        os.remove(temp_zip.name)  # Temporäre Datei nach dem Senden löschen

    return Response(generate(), mimetype='application/zip',
                    headers={"Content-Disposition": "attachment; filename=analysis.zip"})

@api_bp.route('/api/list_file_extensions', methods=['POST'])
def list_file_extensions():
    data = request.get_json()
    session_id = data.get('session_id')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)

    if not os.path.exists(user_dir):
        return jsonify({"status": "error", "message": "No files found"}), 404

    files = [f for f in os.listdir(user_dir) if os.path.isfile(os.path.join(user_dir, f))]
    extensions = set([os.path.splitext(f)[1] for f in files if not (f.startswith('browser') or f.startswith('general'))])

    return jsonify({"status": "success", "extensions": list(extensions)})

@api_bp.route('/api/list_filtered_files', methods=['POST'])
def list_filtered_files():
    data = request.get_json()
    session_id = data.get('session_id')
    extension = data.get('extension', 'all')
    user_dir = os.path.join(Config.OUTPUT_DIR, session_id)

    if not os.path.exists(user_dir):
        return jsonify({"status": "error", "message": "No files found"}), 404

    files = [f for f in os.listdir(user_dir) if os.path.isfile(os.path.join(user_dir, f))]
    filtered_files = [f for f in files if not (f.startswith('browser') or f.startswith('general'))]

    if extension != 'all':
        filtered_files = [f for f in filtered_files if f.endswith(extension)]

    return jsonify({"status": "success", "files": filtered_files})