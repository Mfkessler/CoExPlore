# factory.py

from flask import Flask
from celery import Celery
from celery.signals import worker_ready
from config import Config
from flask_cors import CORS
import os
import logging

def create_app():
    """
    Create and configure the Flask application.
    """

    app_dir = os.path.dirname(os.path.abspath(__file__))
    static_dir = os.path.join(app_dir, 'static')
    templates_dir = os.path.join(app_dir, 'templates')

    app = Flask(
        __name__,
        static_folder=static_dir,
        template_folder=templates_dir
    )

    env_template_dir = os.path.join(templates_dir, Config.TEMPLATE_OUTPUT_DIR)
    app.jinja_loader.searchpath.insert(0, env_template_dir)

    app.secret_key = "my_key_dev"
    app.config['SESSION_COOKIE_NAME'] = 'flask_session_dev'
    app.config['PERMANENT_SESSION_LIFETIME'] = 86400

    CORS(app)

    # Set Flask and Celery configurations
    app.config.update(
        broker_url=Config.CELERY_BROKER_URL,
        result_backend=Config.CELERY_RESULT_BACKEND,
        broker_connection_retry_on_startup=True,
        task_default_queue=Config.CELERY_TASK_DEFAULT_QUEUE,
        beat_schedule={
            'cleanup-sessions-every-60-minutes': {
                'task': 'cleanup_old_sessions_task',
                'schedule': 3600.0,
                'options': {'queue': Config.CELERY_TASK_DEFAULT_QUEUE}
            },
        }
    )

    # Initialize logger
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Ensure output directory exists
    if not os.path.exists(Config.OUTPUT_DIR):
        try:
            os.makedirs(Config.OUTPUT_DIR)
        except Exception as e:
            logger.error(f"Error creating output directory: {e}")

    return app

def make_celery(app):
    """
    Create a new Celery instance for the given Flask app.

    Parameters:
    - app (Flask): The Flask application instance.

    Returns:
    - Celery: The Celery instance.
    """

    celery = Celery(
        app.import_name,
        broker=app.config['broker_url'],
        backend=app.config['result_backend']
    )
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask

    return celery

app = create_app()
celery = make_celery(app)

# Initialize AnnData cache
@worker_ready.connect
def at_worker_start(sender, **kwargs):
    from .tasks import init_adata_cache_task
    try:
        init_adata_cache_task.apply_async()
        logging.info("init_adata_cache_task was successfully started at worker start.")
    except Exception as e:
        logging.error(f"Error starting init_adata_cache_task: {e}")