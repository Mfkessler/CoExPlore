import os
import json

class Config:
    ENV_NAME = os.getenv("ENV_NAME", "default")
    APP_NAME = os.getenv("APP_NAME", "DefaultApp")
    APP_VERSION = os.getenv("APP_VERSION", "0.1.0")
    BASE_URL = os.getenv("BASE_URL", "")

    CELERY_TASK_DEFAULT_QUEUE = os.getenv("CELERY_TASK_DEFAULT_QUEUE", f"{ENV_NAME}_queue")
    BROWSER_DB = os.getenv("BROWSER_DB", "default_browser_db")
    INFO_DB = os.getenv("INFO_DB", "default_info_db")
    TEMPLATE_OUTPUT_DIR = os.getenv("TEMPLATE_OUTPUT_DIR", ENV_NAME)

    OUTPUT_DIR = os.getenv("OUTPUT_DIR", f"output/{ENV_NAME}")
    H5AD_DIR = os.getenv("H5AD_DIR", "/data/h5ad")
    DATA_DIR = os.getenv("DATA_DIR", "/data")
    STATIC_DIR = os.getenv("STATIC_DIR", "static")

    SECRET_KEY = os.getenv("SECRET_KEY", "default_secret_key")
    SESSION_COOKIE_NAME = os.getenv("SESSION_COOKIE_NAME", "flask_session")
    PERMANENT_SESSION_LIFETIME = int(os.getenv("PERMANENT_SESSION_LIFETIME", 86400))

    REDIS_HOST = os.getenv("REDIS_HOST", "localhost")
    REDIS_PORT = int(os.getenv("REDIS_PORT", 6379))
    REDIS_DB = int(os.getenv("REDIS_DB", 0))

    CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL")
    CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND")
    BROKER_CONNECTION_RETRY_ON_STARTUP = os.getenv("BROKER_CONNECTION_RETRY_ON_STARTUP", "True") == "True"

    TASK_DEFAULT_QUEUE = CELERY_TASK_DEFAULT_QUEUE
    BEAT_SCHEDULE = {
        os.getenv("CLEANUP_TASK_NAME", "cleanup_task"): {
            "task": os.getenv("CLEANUP_TASK_NAME", "cleanup_task"),
            "schedule": float(os.getenv("CLEANUP_SCHEDULE", 3600.0)),
            "options": {"queue": CELERY_TASK_DEFAULT_QUEUE},
        },
    }

    METADATA_DICT = json.loads(os.getenv("METADATA_DICT", "{}"))