# extensions.py

import redis
from config import Config

# Initialize Redis client
redis_client = redis.StrictRedis(host=Config.REDIS_HOST, port=Config.REDIS_PORT, db=Config.REDIS_DB)
