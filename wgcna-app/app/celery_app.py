# celery_app.py

from celery import Celery
from config import Config

def make_celery():
    """
    Create a Celery instance for the application

    Returns:
    - Celery: The Celery instance.
    """

    celery = Celery(Config.APP_NAME, broker=Config.CELERY_BROKER_URL, backend=Config.CELERY_RESULT_BACKEND)
    celery.conf.update({
        'task_default_queue': Config.CELERY_TASK_DEFAULT_QUEUE
    })

    return celery

celery = make_celery()
